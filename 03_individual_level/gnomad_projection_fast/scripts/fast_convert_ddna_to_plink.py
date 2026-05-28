#!/usr/bin/env python3
"""
Vectorized, parallel DDNA → PLINK bed/bim/fam converter.

Byte-identical drop-in for the slow pipeline's two-step
  convert_ddna_to_plink.py (tped, Python row loop)
  + plink2 --tped/--tfam --make-bed (tped → bed)
  + Python biallelic filter.

This script:

  1. Discovers sample dirs and reads DDNA TXTs in parallel via
     multiprocessing.Pool. pd.read_csv is used for the parse only — all
     filtering is numpy bytes ops, ~10x faster than pandas .str methods.
  2. Builds the SNP reference from the first sample's autosomal,
     ACGT-only, gs>=min_gs SNPs (deduped by (chrom, pos)).
  3. Looks up each subsequent sample's genotypes against that reference.
  4. Per variant, sets A1 = minor allele in the cohort, A2 = major.
     Mono-allelic variants kept as A1='.', A2=observed (matches plink2
     output from tped with `len(observed) <= 2` pre-filter).
     Tie-break (counts equal): A2 = first non-missing sample's gt[0],
     A1 = the other observed allele. Empirically verified against
     plink2 --make-bed at N=10: 100% match on the ~95k tied variants.
  5. Pads the bed file's unused trailing bit-slots with 0b00 (matching
     plink2's behaviour, undocumented but consistent).

Usage:
    fast_convert_ddna_to_plink.py <data_dir> <out_prefix> [--min-gs F]
                                  [--workers N]
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

for _parent in Path(__file__).resolve().parents:
    if (_parent / "tools" / "genotype_normalizer.py").exists():
        sys.path.insert(0, str(_parent))
        break
sys.path.append("/opt/biovault")
from tools import genotype_normalizer as _genoio  # noqa: E402


AUTOSOMES_S2 = np.array([str(c).encode() for c in range(1, 23)], dtype="S2")
BASES_S1 = np.frombuffer(b"ACGT", dtype="S1")
BASE_CODE = {b"A": 1, b"C": 2, b"G": 3, b"T": 4}


def read_ddna_fast(path: str, min_gs: float):
    """Return (rsid, chrom_bytes, pos, gt_bytes) numpy arrays for one DDNA TXT.

    Filtered to: autosomes 1..22, gs >= min_gs, len(gt)==2 with both bases in
    ACGT. Bytes dtype (S1/S2) is used throughout for cheap allele ops.
    """
    if _genoio.sniff_format(path) == "illumina":
        # GSGT synthetic carries no GenCall score; treat every call as
        # confident (gs = 1.0) so the gs >= min_gs filter is a no-op.
        gdf = _genoio.read_pipeline_genotypes(path)
        rsid = gdf["rsid"].to_numpy(dtype=object)
        chrom = gdf["chrom"].to_numpy(dtype=object).astype("S2")
        pos = gdf["pos"].astype("int64").to_numpy()
        gt = gdf["gt"].to_numpy(dtype=object).astype("S2")
        gs = np.ones(len(rsid), dtype=np.float32)
    else:
        df = pd.read_csv(
            path, sep="\t", comment="#", header=None,
            names=["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"],
            usecols=["rsid", "chrom", "pos", "gt", "gs"],
            dtype={"rsid": str, "chrom": str, "pos": str, "gt": str, "gs": str},
            engine="c", na_filter=False,
        )
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
        df["gs"] = pd.to_numeric(df["gs"], errors="coerce")
        rsid = df["rsid"].to_numpy(dtype=object)
        chrom = df["chrom"].to_numpy(dtype=object).astype("S2")
        pos = df["pos"].to_numpy(dtype=np.float64)
        gt = df["gt"].to_numpy(dtype=object).astype("S2")
        gs = df["gs"].to_numpy(dtype=np.float32)

    # Drop rows where rsid == "rsid" (any embedded header row).
    mask = rsid != "rsid"
    mask &= np.isfinite(pos)
    mask &= np.isfinite(gs)
    mask &= gs >= np.float32(min_gs)
    mask &= np.isin(chrom, AUTOSOMES_S2)

    # gt length 2 (raw "S2" zero-pads short strings, so check char[1] != \0).
    gt_view = gt.view("S1").reshape(-1, 2)
    mask &= gt_view[:, 1] != b""
    mask &= np.isin(gt_view[:, 0], BASES_S1)
    mask &= np.isin(gt_view[:, 1], BASES_S1)

    keep = np.flatnonzero(mask)
    return rsid[keep], chrom[keep], pos[keep].astype(np.int64), gt[keep]


# Globals populated in worker init / main process; readable across forks.
_REF_RSID_TO_IDX: dict | None = None
_MIN_GS: float = 0.15


def _worker_init(rsid_to_idx, min_gs):
    global _REF_RSID_TO_IDX, _MIN_GS
    _REF_RSID_TO_IDX = rsid_to_idx
    _MIN_GS = min_gs


def _worker_process(args):
    i, path = args
    rsid, _chrom, _pos, gt = read_ddna_fast(path, _MIN_GS)
    # Map rsid -> reference index; dropna -> kept rows.
    # Use Series.map for vectorized dict lookup.
    idx_series = pd.Series(rsid).map(_REF_RSID_TO_IDX)
    valid = idx_series.notna().to_numpy()
    if not valid.any():
        return i, np.empty(0, dtype=np.int64), np.empty(0, dtype=np.uint8), np.empty(0, dtype=np.uint8)
    ix = idx_series[valid].to_numpy(dtype=np.int64)
    gt_view = gt[valid].view("S1").reshape(-1, 2)
    # Encode allele bytes to 1..4. Use cumulative mask trick — faster than dict.
    a1_codes = np.zeros(gt_view.shape[0], dtype=np.uint8)
    a2_codes = np.zeros(gt_view.shape[0], dtype=np.uint8)
    for code, base in enumerate(BASES_S1, start=1):
        m1 = gt_view[:, 0] == base
        m2 = gt_view[:, 1] == base
        a1_codes[m1] = code
        a2_codes[m2] = code
    # If same rsid appears multiple times in a sample, keep first (matches
    # the slow pipeline, which uses a dict-based one-shot fill in row order).
    # pd.Series.map preserves order so this is naturally first-seen; but if
    # there were duplicate rsids in the sample DDNA, our final write of A1[i,
    # ix] = a1v would let last-occurrence win. Use np.unique on ix to dedup.
    if np.unique(ix).size != ix.size:
        # Keep first occurrence.
        _, first_idx = np.unique(ix, return_index=True)
        first_idx.sort()
        ix = ix[first_idx]
        a1_codes = a1_codes[first_idx]
        a2_codes = a2_codes[first_idx]
    return i, ix, a1_codes, a2_codes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("--min-gs", type=float, default=0.15)
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 4)))
    args = ap.parse_args()

    sample_dirs = sorted(
        d for d in glob.glob(os.path.join(args.data_dir, "*"))
        if os.path.isdir(d)
    )
    if not sample_dirs:
        sys.exit(f"ERROR: no sample dirs in {args.data_dir}")
    n_samples = len(sample_dirs)
    sample_ids = [os.path.basename(d) for d in sample_dirs]
    print(f"Found {n_samples} samples")

    def txt_for(d: str) -> str:
        txts = glob.glob(os.path.join(d, "*.txt"))
        if not txts:
            raise FileNotFoundError(f"no .txt in {d}")
        return txts[0]

    # --- Build SNP reference from first sample ----------------------------
    t0 = time.time()
    first_path = txt_for(sample_dirs[0])
    print(f"Building SNP reference from: {os.path.basename(first_path)}")
    rsid0, chrom0, pos0, gt0 = read_ddna_fast(first_path, args.min_gs)

    # Dedup by (chrom, pos), keep first.
    n_before_dedup = rsid0.size
    # Encode chrom/pos to a single uint64 sort/group key (chrom in 1..22 fits 5 bits).
    chrom_int = np.array([int(c) for c in chrom0.astype(str)], dtype=np.int32)
    sort_key = chrom_int.astype(np.int64) * (10**11) + pos0
    order = np.argsort(sort_key, kind="stable")
    rsid0 = rsid0[order]
    chrom0 = chrom0[order]
    pos0 = pos0[order]
    gt0 = gt0[order]
    chrom_int = chrom_int[order]
    sort_key = sort_key[order]
    # First occurrence of each (chrom, pos)
    _, first_idx = np.unique(sort_key, return_index=True)
    first_idx.sort()
    rsid0 = rsid0[first_idx]
    chrom0 = chrom0[first_idx]
    pos0 = pos0[first_idx]
    gt0 = gt0[first_idx]

    ref_chrom = chrom0
    ref_rsid = rsid0
    ref_pos = pos0
    n_snps = ref_rsid.size
    rsid_to_idx = {str(r): i for i, r in enumerate(ref_rsid.tolist())}
    print(f"Reference SNP set: {n_snps:,} autosomal SNPs ({time.time()-t0:.2f}s)")

    # --- Read all samples in parallel -------------------------------------
    A1i = np.zeros((n_samples, n_snps), dtype=np.uint8)
    A2i = np.zeros((n_samples, n_snps), dtype=np.uint8)

    t0 = time.time()
    work = list(enumerate([txt_for(d) for d in sample_dirs]))
    if args.workers > 1 and n_samples > 1:
        ctx = mp.get_context("fork")
        with ctx.Pool(min(args.workers, n_samples),
                      initializer=_worker_init,
                      initargs=(rsid_to_idx, args.min_gs)) as pool:
            for i, ix, a1c, a2c in pool.imap_unordered(_worker_process, work, chunksize=1):
                A1i[i, ix] = a1c
                A2i[i, ix] = a2c
    else:
        _worker_init(rsid_to_idx, args.min_gs)
        for w in work:
            i, ix, a1c, a2c = _worker_process(w)
            A1i[i, ix] = a1c
            A2i[i, ix] = a2c
    print(f"Read {n_samples} samples in {time.time()-t0:.2f}s")

    # --- Per-variant allele counts, biallelic filter -----------------------
    t0 = time.time()
    n_var = n_snps
    counts = np.zeros((n_var, 4), dtype=np.int32)
    both_present = (A1i != 0) & (A2i != 0)
    for code in range(1, 5):
        counts[:, code - 1] = (((A1i == code) & both_present).sum(axis=0)
                               + ((A2i == code) & both_present).sum(axis=0))
    present_mask = counts > 0
    n_distinct = present_mask.sum(axis=1)
    keep_variant = (n_distinct >= 1) & (n_distinct <= 2)
    print(f"Allele counts in {time.time()-t0:.2f}s; "
          f"mono+biallelic kept: {int(keep_variant.sum()):,}/{n_var:,} "
          f"(mono: {int((n_distinct==1).sum()):,}, bi: {int((n_distinct==2).sum()):,})")

    # --- A1/A2 assignment --------------------------------------------------
    t0 = time.time()
    keep_idx = np.where(keep_variant)[0]
    k = len(keep_idx)
    cnt_k = counts[keep_idx]
    pres_k = present_mask[keep_idx]
    nd_k = n_distinct[keep_idx]
    order = np.argsort(~pres_k, axis=1, kind="stable")
    allele_lo_code = (order[:, 0] + 1).astype(np.uint8)
    allele_hi_code = (order[:, 1] + 1).astype(np.uint8)
    mono = nd_k == 1
    if mono.any():
        allele_hi_code = np.where(mono, 0, allele_hi_code).astype(np.uint8)
    cnt_lo = cnt_k[np.arange(k), allele_lo_code - 1]
    safe_hi_idx = np.where(allele_hi_code > 0, allele_hi_code - 1, 0)
    cnt_hi = np.where(allele_hi_code > 0, cnt_k[np.arange(k), safe_hi_idx], 0)
    tied = (cnt_lo == cnt_hi) & (~mono)

    a1_code = np.where(cnt_lo < cnt_hi, allele_lo_code, allele_hi_code).astype(np.uint8)
    a2_code = np.where(cnt_lo < cnt_hi, allele_hi_code, allele_lo_code).astype(np.uint8)

    if mono.any():
        a1_code = np.where(mono, 0, a1_code).astype(np.uint8)
        a2_code = np.where(mono, allele_lo_code, a2_code).astype(np.uint8)

    if tied.any():
        A1i_k = A1i[:, keep_idx]
        nonmiss = A1i_k != 0
        first_sample = nonmiss.argmax(axis=0)
        anchor = A1i_k[first_sample, np.arange(k)]
        anchor_is_lo = anchor == allele_lo_code
        tied_a2 = np.where(anchor_is_lo, allele_lo_code, allele_hi_code)
        tied_a1 = np.where(anchor_is_lo, allele_hi_code, allele_lo_code)
        a1_code = np.where(tied, tied_a1, a1_code).astype(np.uint8)
        a2_code = np.where(tied, tied_a2, a2_code).astype(np.uint8)
    print(f"Per-variant A1/A2 in {time.time()-t0:.2f}s "
          f"(tied: {int(tied.sum()):,}, mono: {int(mono.sum()):,})")

    # --- Genotype codes (PLINK bed encoding) -------------------------------
    t0 = time.time()
    A1i_k = A1i[:, keep_idx]
    A2i_k = A2i[:, keep_idx]
    a1c = a1_code[None, :]
    a2c = a2_code[None, :]
    missing = (A1i_k == 0) | (A2i_k == 0)
    n_a1 = ((A1i_k == a1c).astype(np.uint8) + (A2i_k == a1c).astype(np.uint8))
    code = np.where(n_a1 == 2, 0b00,
            np.where(n_a1 == 1, 0b10,
            np.where(n_a1 == 0, 0b11, 0b01))).astype(np.uint8)
    code[missing] = 0b01

    out_of_set = (~missing) & (
        ((A1i_k != a1c) & (A1i_k != a2c)) |
        ((A2i_k != a1c) & (A2i_k != a2c))
    )
    if out_of_set.any():
        code[out_of_set] = 0b01
    print(f"Genotype codes in {time.time()-t0:.2f}s; matrix shape {code.shape}")

    # --- Write .fam ---------------------------------------------------------
    t0 = time.time()
    with open(args.out_prefix + ".fam", "w") as f:
        for sid in sample_ids:
            f.write(f"0\t{sid}\t0\t0\t0\t-9\n")

    # --- Write .bim ---------------------------------------------------------
    base_chars = np.array([".", "A", "C", "G", "T"])
    a1_chars = base_chars[a1_code]
    a2_chars = base_chars[a2_code]
    chrom_strs = ref_chrom[keep_idx].astype(str)
    rsid_strs = ref_rsid[keep_idx].astype(str)
    pos_arr = ref_pos[keep_idx]
    with open(args.out_prefix + ".bim", "w") as f:
        # Build all lines, write once. ~6x faster than per-row f.write.
        lines = (f"{c}\t{r}\t0\t{int(p)}\t{a1}\t{a2}"
                 for c, r, p, a1, a2 in zip(chrom_strs, rsid_strs, pos_arr,
                                            a1_chars, a2_chars))
        f.write("\n".join(lines))
        f.write("\n")

    # --- Write .bed ---------------------------------------------------------
    code_vmaj = code.T
    n_samp = code_vmaj.shape[1]
    bytes_per_var = (n_samp + 3) // 4
    pad = bytes_per_var * 4 - n_samp
    if pad:
        # plink2 pads unused slots with 0b00 (verified empirically).
        pad_arr = np.zeros((code_vmaj.shape[0], pad), dtype=np.uint8)
        code_padded = np.concatenate([code_vmaj, pad_arr], axis=1)
    else:
        code_padded = code_vmaj
    reshaped = code_padded.reshape(code_padded.shape[0], bytes_per_var, 4)
    bed_bytes = (
        reshaped[:, :, 0]
        | (reshaped[:, :, 1] << 2)
        | (reshaped[:, :, 2] << 4)
        | (reshaped[:, :, 3] << 6)
    ).astype(np.uint8)
    with open(args.out_prefix + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01")
        f.write(bed_bytes.tobytes())
    print(f"Wrote bed/bim/fam in {time.time()-t0:.2f}s")
    print(f"  {args.out_prefix}.fam ({n_samples} samples)")
    print(f"  {args.out_prefix}.bim ({k} SNPs)")
    print(f"  {args.out_prefix}.bed")


if __name__ == "__main__":
    main()

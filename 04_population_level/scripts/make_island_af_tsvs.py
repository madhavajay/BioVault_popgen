#!/usr/bin/env python3
"""
Compute per-island allele frequencies from DDNA TXT genotype files and write
one TSV per island in the format that 04_population_level/fst_islands/
scripts/01_load_merge.py consumes.

Inputs
------
--mapping   TSV with columns: participant_id, island
            (produced by 01_mock_data_generation/scripts/make_island_mapping.py)
--data-dir  Dir containing <pid>/<file>.txt for each participant.
--out-dir   Where to write the per-island AF TSVs. Files inside the dir are
            named allele_freq_<filename_label>.tsv with case matching the
            existing ISLAND_FILES map in 01_load_merge.py.

Output columns (one TSV per island)
-----------------------------------
locus_key   rsid   allele_freq   allele_number

locus_key   = "<chrom>-<pos>-<ref>-<alt>", with ref/alt set to the two non-zero
              alleles observed across the whole cohort, sorted alphabetically.
              Multi-allelic loci (>2 distinct alleles observed cohort-wide)
              are dropped.
allele_freq = alt allele count / total called alleles in this island.
allele_number = total called alleles in this island (≈ 2 × n_genotyped).

Notes
-----
* Autosomes only (chr1..22).
* GS (gencall) below --min-gs is treated as missing.
* No structure exists in biosynth's --alt-frequency 0.5 mock, so all islands
  will land near AF=0.5 ± sqrt(0.25/n_alleles). That's a useful correctness
  check; pipelines downstream will report FST ≈ 0.
"""

from __future__ import annotations

import argparse
import logging
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# Island label (in mapping) → filename suffix (matches 01_load_merge.py).
ISLAND_FILENAMES = {
    "BVI":       "BVI",
    "TT":        "TT",
    "Bahamas":   "bahamas",
    "Barbados":  "barbados",
    "Bermuda":   "bermuda",
    "StLucia":   "stlucia",
}

BASES = ("A", "C", "G", "T")
_BASE_COL = {b: i for i, b in enumerate(BASES)}


def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")


def read_genotypes(txt_path: Path, min_gs: float) -> pd.DataFrame:
    """Read one DDNA TXT, return [rsid, chrom, pos, a1, a2] for valid SNP rows."""
    df = pd.read_csv(
        txt_path,
        sep="\t",
        header=None,
        comment="#",
        names=["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"],
        usecols=["rsid", "chrom", "pos", "gt", "gs"],
        dtype={"rsid": str, "chrom": str, "pos": np.int32, "gt": str, "gs": float},
        engine="c",
    )
    df = df[df["rsid"] != "rsid"]                                    # drop header row if any
    df = df[df["chrom"].isin([str(c) for c in range(1, 23)])]        # autosomes only
    df = df[df["gs"] >= min_gs]                                      # gencall filter
    df = df[df["gt"].str.len() == 2]                                 # 2-char genotype
    df["a1"] = df["gt"].str[0]
    df["a2"] = df["gt"].str[1]
    df = df[df["a1"].isin(BASES) & df["a2"].isin(BASES)]
    return df


def main():
    setup_logging()
    log = logging.getLogger(__name__)

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mapping", type=Path, required=True)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-gs", type=float, default=0.15)
    args = ap.parse_args()

    mapping = pd.read_csv(args.mapping, sep="\t", dtype={"participant_id": str})
    by_island: dict[str, list[str]] = defaultdict(list)
    for _, row in mapping.iterrows():
        by_island[row["island"]].append(row["participant_id"])
    islands = sorted(by_island.keys())
    log.info(f"Mapping: {len(mapping)} participants across {len(islands)} islands")
    for island in islands:
        log.info(f"  {island}: {len(by_island[island])} participants")

    # SNP universe from the first available participant.
    first_pid = mapping["participant_id"].iloc[0]
    first_txts = list((args.data_dir / first_pid).glob("*.txt"))
    if not first_txts:
        raise SystemExit(f"No .txt in {args.data_dir / first_pid}")
    log.info(f"Building SNP universe from {first_txts[0].name} ...")
    universe = read_genotypes(first_txts[0], min_gs=0.0)
    universe = universe[["rsid", "chrom", "pos"]].drop_duplicates(subset="rsid").reset_index(drop=True)
    universe["snp_idx"] = np.arange(len(universe), dtype=np.int32)
    n_snps = len(universe)
    rsid_to_idx = dict(zip(universe["rsid"], universe["snp_idx"]))
    log.info(f"SNP universe: {n_snps:,}")

    counts = {island: np.zeros((n_snps, 4), dtype=np.uint32) for island in islands}

    t0 = time.time()
    total_files = sum(len(v) for v in by_island.values())
    done = 0
    for island in islands:
        for pid in by_island[island]:
            txts = list((args.data_dir / pid).glob("*.txt"))
            if not txts:
                log.warning(f"  skip {pid}: no .txt found")
                continue
            df = read_genotypes(txts[0], min_gs=args.min_gs)
            idx_series = df["rsid"].map(rsid_to_idx)
            valid = idx_series.notna()
            idx = idx_series[valid].to_numpy(dtype=np.int32)
            a1 = df.loc[valid, "a1"].map(_BASE_COL).to_numpy(dtype=np.int32)
            a2 = df.loc[valid, "a2"].map(_BASE_COL).to_numpy(dtype=np.int32)
            np.add.at(counts[island], (idx, a1), 1)
            np.add.at(counts[island], (idx, a2), 1)
            done += 1
            if done % 50 == 0 or done == total_files:
                rate = done / max(time.time() - t0, 1e-6)
                eta = (total_files - done) / max(rate, 1e-6)
                log.info(f"  ingested {done}/{total_files} files "
                         f"({rate:.1f}/s, ETA {eta:.0f}s)")

    # Cohort-wide alleles → biallelic mask + ref/alt assignment (lex order).
    total = sum(counts.values())
    nonzero = total > 0
    n_alleles_per_snp = nonzero.sum(axis=1)
    biallelic_mask = n_alleles_per_snp == 2
    biallelic_idx = np.where(biallelic_mask)[0]
    log.info(f"Cohort biallelic loci: {biallelic_idx.size:,} / {n_snps:,}")

    ref_col = np.full(n_snps, -1, dtype=np.int8)
    alt_col = np.full(n_snps, -1, dtype=np.int8)
    for i in biallelic_idx:
        cols = np.where(nonzero[i])[0]
        ref_col[i] = cols[0]
        alt_col[i] = cols[1]

    chrom_arr = universe["chrom"].to_numpy()
    pos_arr = universe["pos"].to_numpy()
    rsid_arr = universe["rsid"].to_numpy()
    bases_arr = np.array(BASES)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for island in islands:
        c = counts[island]
        ref_count = c[biallelic_idx, ref_col[biallelic_idx]]
        alt_count = c[biallelic_idx, alt_col[biallelic_idx]]
        total_calls = ref_count + alt_count
        keep = total_calls > 0
        sub_idx = biallelic_idx[keep]
        ref_count = ref_count[keep]
        alt_count = alt_count[keep]
        total_calls = total_calls[keep]

        ref_b = bases_arr[ref_col[sub_idx]]
        alt_b = bases_arr[alt_col[sub_idx]]
        ch = chrom_arr[sub_idx]
        po = pos_arr[sub_idx]
        rs = rsid_arr[sub_idx]

        locus_key = [f"{ch[i]}-{po[i]}-{ref_b[i]}-{alt_b[i]}" for i in range(len(sub_idx))]
        af = alt_count.astype(np.float64) / total_calls.astype(np.float64)

        out_df = pd.DataFrame({
            "locus_key": locus_key,
            "rsid": rs,
            "allele_freq": af,
            "allele_number": total_calls.astype(np.int64),
        })
        label = ISLAND_FILENAMES.get(island, island)
        out_path = args.out_dir / f"allele_freq_{label}.tsv"
        out_df.to_csv(out_path, sep="\t", index=False)
        log.info(f"  {island} -> {out_path}  ({len(out_df):,} rows, mean AF {af.mean():.3f})")

    log.info(f"Done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()

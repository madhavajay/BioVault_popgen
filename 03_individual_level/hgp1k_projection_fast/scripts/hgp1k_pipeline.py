#!/usr/bin/env python3
"""Joint PCA of BioVault samples and filtered 1000 Genomes high-coverage VCFs.

This pipeline deliberately does not project onto precomputed gnomAD loadings.
It builds a shared SNP matrix from:

  * 1000 Genomes high-coverage VCFs filtered to tools/locus_map.tsv positions.
  * BioVault/DDNA/Illumina genotype text files normalized by tools.genotype_normalizer.

It then computes allele frequencies for both groups, imputes missing genotypes
to the combined mean dosage, standardizes by the combined allele frequency, and
fits a PCA on the combined 1KGP + study matrix.

Outputs in out_dir:
  pca_scores.tsv              all 1KGP and study samples
  study_pca_projection.tsv    study-only legacy shape: s<TAB>scores
  allele_freqs.tsv            per-variant 1KGP and study allele frequencies
  qc_report.txt               run summary
  pca_projection.png          PC1/PC2 scatter
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

for _parent in Path(__file__).resolve().parents:
    if (_parent / "tools" / "genotype_normalizer.py").exists():
        sys.path.insert(0, str(_parent))
        break
sys.path.append("/opt/biovault")
from tools import genotype_normalizer as genoio  # noqa: E402


BASES = {"A", "C", "G", "T"}
AUTOSOMES = [str(i) for i in range(1, 23)]
MISSING = np.uint8(255)
REGION_CMAPS = {
    "AFR": "Oranges",
    "AMR": "Blues",
    "CSA": "YlGn",
    "EAS": "Greens",
    "EUR": "Reds",
    "MID": "Purples",
    "OCE": "PuRd",
}
REGION_ORDER = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]
REGION_ALIASES = {
    "SAS": "CSA",
    "AMR_HGDP": "AMR",
}

_WORKER_KEY_TO_IDX: dict[tuple[str, int], int] = {}
_WORKER_REF: np.ndarray | None = None
_WORKER_ALT: np.ndarray | None = None
_WORKER_MIN_GS: float = 0.15


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str

    @property
    def key(self) -> tuple[str, int]:
        return self.chrom, self.pos


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")


def normalize_chrom(value: str) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return "MT" if text == "M" else text.upper()


def chromosome_list(value: str) -> list[str]:
    if value == "all":
        return AUTOSOMES
    out: list[str] = []
    for raw in value.split(","):
        chrom = normalize_chrom(raw)
        if chrom not in AUTOSOMES:
            raise ValueError(f"unsupported chromosome for PCA: {raw!r}")
        out.append(chrom)
    return out


def vcf_name_for_chr(chrom: str) -> str:
    if chrom == "X":
        return "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    return f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def find_vcfs(vcf_dir: Path, chromosomes: list[str]) -> list[Path]:
    vcfs: list[Path] = []
    for chrom in chromosomes:
        candidate = vcf_dir / vcf_name_for_chr(chrom)
        if not candidate.exists():
            logging.warning("missing filtered 1KGP VCF for chr%s: %s", chrom, candidate)
            continue
        if (candidate.parent / f"{candidate.name}.aria2").exists():
            logging.warning("skipping active/incomplete aria2 download: %s", candidate)
            continue
        if not (candidate.parent / f"{candidate.name}.tbi").exists():
            logging.warning("skipping unindexed VCF: %s", candidate)
            continue
        vcfs.append(candidate)
    if not vcfs:
        raise FileNotFoundError(f"no usable filtered VCFs found in {vcf_dir}")
    return vcfs


def run_text(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)


def vcf_samples(vcf: Path) -> list[str]:
    text = run_text(["bcftools", "query", "-l", str(vcf)])
    return [line.strip() for line in text.splitlines() if line.strip()]


def parse_gt_dosage(gt: str) -> int:
    gt = gt.split(":", 1)[0]
    if "." in gt:
        return int(MISSING)
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return int(MISSING)
    dosage = 0
    for allele in alleles:
        if allele == "0":
            continue
        if allele == "1":
            dosage += 1
            continue
        return int(MISSING)
    return dosage


def load_locus_rsids(locus_map: Path) -> dict[tuple[str, int], str]:
    df = pd.read_csv(locus_map, sep="\t", dtype={"chrom": str, "pos": int, "rsid": str})
    df["chrom"] = df["chrom"].map(normalize_chrom)
    out: dict[tuple[str, int], str] = {}
    conflicts: set[tuple[str, int]] = set()
    for chrom, pos, rsid in df[["chrom", "pos", "rsid"]].itertuples(index=False, name=None):
        key = (chrom, int(pos))
        if key in conflicts:
            continue
        existing = out.get(key)
        if existing is not None and existing != rsid:
            out.pop(key, None)
            conflicts.add(key)
            continue
        out[key] = str(rsid)
    return out


def stream_vcf_matrix(
    vcfs: list[Path],
    locus_rsids: dict[tuple[str, int], str],
    max_variants: int | None = None,
) -> tuple[list[str], list[Variant], np.ndarray]:
    log = logging.getLogger(__name__)
    samples = vcf_samples(vcfs[0])
    rows: list[np.ndarray] = []
    variants: list[Variant] = []
    seen_sample_tuple = tuple(samples)

    for vcf in vcfs:
        current_samples = vcf_samples(vcf)
        if tuple(current_samples) != seen_sample_tuple:
            raise ValueError(f"{vcf}: sample list/order differs from first VCF")

        log.info("reading 1KGP VCF: %s", vcf)
        cmd = [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n",
            str(vcf),
        ]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )
        assert proc.stdout is not None
        n_read = 0
        n_kept = 0
        stopped_early = False
        for line in proc.stdout:
            n_read += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 + len(samples):
                continue
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
                continue
            if "," in alt:
                continue
            rsid = locus_rsids.get((chrom, pos), fields[2])
            dosage = np.fromiter(
                (parse_gt_dosage(gt) for gt in fields[5:]),
                dtype=np.uint8,
                count=len(samples),
            )
            rows.append(dosage)
            variants.append(Variant(chrom=chrom, pos=pos, rsid=rsid, ref=ref, alt=alt))
            n_kept += 1
            if max_variants is not None and len(variants) >= max_variants:
                stopped_early = True
                proc.terminate()
                break
        if stopped_early:
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        rc = proc.wait()
        if rc != 0 and not stopped_early:
            raise RuntimeError(f"bcftools query failed for {vcf}:\n{stderr}")
        if stderr.strip():
            log.warning("bcftools warning for %s: %s", vcf.name, stderr.strip().splitlines()[-1])
        log.info("  read %s VCF rows, kept %s biallelic SNP rows", f"{n_read:,}", f"{n_kept:,}")
        if max_variants is not None and len(variants) >= max_variants:
            break

    if not variants:
        raise ValueError("no biallelic SNP variants loaded from 1KGP VCFs")
    matrix = np.stack(rows, axis=1)  # samples x variants, uint8, 255 = missing
    return samples, variants, matrix


def _decode_array_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def load_reference_matrix_npz(path: Path) -> tuple[list[str], list[Variant], np.ndarray]:
    arc = np.load(path, allow_pickle=True)
    required = {"dosage", "samples", "chrom", "pos", "rsid", "ref", "alt"}
    missing = required - set(arc.files)
    if missing:
        raise ValueError(f"{path}: missing npz keys: {sorted(missing)}")
    dosage = arc["dosage"].astype(np.uint8, copy=True)
    if "missing_mask_packed" in arc.files and "missing_mask_shape" in arc.files:
        shape = tuple(int(x) for x in arc["missing_mask_shape"].tolist())
        if shape != dosage.shape:
            raise ValueError(f"{path}: missing mask shape {shape} does not match dosage shape {dosage.shape}")
        missing = np.unpackbits(arc["missing_mask_packed"], count=dosage.size).reshape(dosage.shape).astype(bool)
        if missing.any():
            dosage[missing] = MISSING
    samples = [_decode_array_value(value) for value in arc["samples"].tolist()]
    chrom = [_decode_array_value(value) for value in arc["chrom"].tolist()]
    pos = arc["pos"].astype(np.int64)
    rsid = [_decode_array_value(value) for value in arc["rsid"].tolist()]
    ref = [_decode_array_value(value) for value in arc["ref"].tolist()]
    alt = [_decode_array_value(value) for value in arc["alt"].tolist()]
    variants = [
        Variant(
            chrom=normalize_chrom(c),
            pos=int(p),
            rsid=r,
            ref=rf,
            alt=al,
        )
        for c, p, r, rf, al in zip(chrom, pos, rsid, ref, alt)
    ]
    if dosage.shape != (len(samples), len(variants)):
        raise ValueError(
            f"{path}: dosage shape {dosage.shape} does not match "
            f"{len(samples)} samples x {len(variants)} variants"
        )
    return samples, variants, dosage


def drop_duplicate_loci(variants: list[Variant], matrix: np.ndarray) -> tuple[list[Variant], np.ndarray, list[Variant]]:
    counts: dict[tuple[str, int], int] = {}
    for variant in variants:
        counts[variant.key] = counts.get(variant.key, 0) + 1
    keep = np.array([counts[v.key] == 1 for v in variants], dtype=bool)
    kept = [v for v, ok in zip(variants, keep) if ok]
    dropped = [v for v, ok in zip(variants, keep) if not ok]
    return kept, matrix[:, keep], dropped


def discover_study_files(data_dir: Path) -> list[tuple[str, Path]]:
    sample_dirs = sorted(d for d in data_dir.iterdir() if d.is_dir())
    out: list[tuple[str, Path]] = []
    for sample_dir in sample_dirs:
        txts = sorted(sample_dir.glob("*.txt")) + sorted(sample_dir.glob("*.txt.gz"))
        if txts:
            out.append((sample_dir.name, txts[0]))
        else:
            logging.warning("skip %s: no .txt or .txt.gz genotype file", sample_dir.name)
    if not out:
        raise FileNotFoundError(f"no participant .txt or .txt.gz files found under {data_dir}")
    return out


def _init_study_worker(
    variants: list[Variant],
    min_gs: float,
) -> None:
    global _WORKER_KEY_TO_IDX, _WORKER_REF, _WORKER_ALT, _WORKER_MIN_GS
    _WORKER_KEY_TO_IDX = {variant.key: idx for idx, variant in enumerate(variants)}
    _WORKER_REF = np.array([v.ref for v in variants], dtype=object)
    _WORKER_ALT = np.array([v.alt for v in variants], dtype=object)
    _WORKER_MIN_GS = min_gs


def _parse_study_sample(task: tuple[int, str, str]) -> tuple[int, str, list[tuple[int, int]], list[dict[str, str]]]:
    row_idx, sample_id, path_text = task
    path = Path(path_text)
    errors: list[dict[str, str]] = []
    hits: list[tuple[int, int]] = []
    try:
        # PCA matching below uses chrom/pos/ref/alt, so resolving probe IDs
        # through the default locus map only adds repeated parse overhead.
        df = genoio.read_pipeline_genotypes(path, locus_map={})
    except Exception as exc:
        errors.append({
            "participant_id": sample_id,
            "file": str(path),
            "code": "PARSE_FAILED",
            "message": str(exc),
        })
        return row_idx, sample_id, hits, errors
    if df.empty:
        return row_idx, sample_id, hits, errors

    df["chrom"] = df["chrom"].map(normalize_chrom)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["gs"] = pd.to_numeric(df["gs"], errors="coerce").fillna(1.0)
    df = df[df["pos"].notna() & df["chrom"].isin(AUTOSOMES)]
    df = df[df["gs"] >= _WORKER_MIN_GS]
    df = df[df["gt"].astype(str).str.fullmatch(r"[ACGT]{2}", na=False)]

    if _WORKER_REF is None or _WORKER_ALT is None:
        raise RuntimeError("study worker was not initialized")

    seen: set[int] = set()
    for chrom, pos_raw, gt in df[["chrom", "pos", "gt"]].itertuples(index=False, name=None):
        key = (chrom, int(pos_raw))
        idx = _WORKER_KEY_TO_IDX.get(key)
        if idx is None or idx in seen:
            continue
        seen.add(idx)
        a1, a2 = str(gt)[0], str(gt)[1]
        ref = _WORKER_REF[idx]
        alt = _WORKER_ALT[idx]
        if (a1 not in {ref, alt}) or (a2 not in {ref, alt}):
            continue
        hits.append((idx, int((a1 == alt) + (a2 == alt))))
    return row_idx, sample_id, hits, errors


def load_study_matrix(
    sample_files: list[tuple[str, Path]],
    variants: list[Variant],
    min_gs: float,
    workers: int,
) -> tuple[list[str], np.ndarray, list[dict[str, str]]]:
    key_to_idx = {variant.key: idx for idx, variant in enumerate(variants)}
    ref = np.array([v.ref for v in variants], dtype=object)
    alt = np.array([v.alt for v in variants], dtype=object)
    matrix = np.full((len(sample_files), len(variants)), MISSING, dtype=np.uint8)
    sample_ids: list[str] = []
    errors: list[dict[str, str]] = []

    if workers > 1 and len(sample_files) > 1:
        sample_ids = [sample_id for sample_id, _path in sample_files]
        tasks = [
            (idx, sample_id, str(path))
            for idx, (sample_id, path) in enumerate(sample_files)
        ]
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_study_worker,
            initargs=(variants, min_gs),
        ) as pool:
            futures = [pool.submit(_parse_study_sample, task) for task in tasks]
            completed = 0
            for fut in as_completed(futures):
                row_idx, _sample_id, hits, sample_errors = fut.result()
                errors.extend(sample_errors)
                if hits:
                    cols = np.fromiter((idx for idx, _dosage in hits), dtype=np.int64, count=len(hits))
                    vals = np.fromiter((dosage for _idx, dosage in hits), dtype=np.uint8, count=len(hits))
                    matrix[row_idx, cols] = vals
                completed += 1
                if completed % 25 == 0 or completed == len(sample_files):
                    logging.info("parsed study samples: %s/%s", completed, len(sample_files))
        return sample_ids, matrix, errors

    for row_idx, (sample_id, path) in enumerate(sample_files):
        sample_ids.append(sample_id)
        try:
            # PCA matching below uses chrom/pos/ref/alt, so resolving probe IDs
            # through the default locus map only adds repeated parse overhead.
            df = genoio.read_pipeline_genotypes(path, locus_map={})
        except Exception as exc:
            errors.append({
                "participant_id": sample_id,
                "file": str(path),
                "code": "PARSE_FAILED",
                "message": str(exc),
            })
            continue
        if df.empty:
            continue
        df["chrom"] = df["chrom"].map(normalize_chrom)
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
        df["gs"] = pd.to_numeric(df["gs"], errors="coerce").fillna(1.0)
        df = df[df["pos"].notna() & df["chrom"].isin(AUTOSOMES)]
        df = df[df["gs"] >= min_gs]
        df = df[df["gt"].astype(str).str.fullmatch(r"[ACGT]{2}", na=False)]

        seen: set[int] = set()
        for chrom, pos_raw, gt in df[["chrom", "pos", "gt"]].itertuples(index=False, name=None):
            key = (chrom, int(pos_raw))
            idx = key_to_idx.get(key)
            if idx is None or idx in seen:
                continue
            seen.add(idx)
            a1, a2 = str(gt)[0], str(gt)[1]
            if (a1 not in {ref[idx], alt[idx]}) or (a2 not in {ref[idx], alt[idx]}):
                continue
            matrix[row_idx, idx] = np.uint8((a1 == alt[idx]) + (a2 == alt[idx]))
    return sample_ids, matrix, errors


def allele_stats(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    called = matrix != MISSING
    called_samples = called.sum(axis=0).astype(np.int64)
    alt_sum = np.where(called, matrix, 0).sum(axis=0).astype(np.float64)
    allele_n = called_samples * 2
    af = np.divide(alt_sum, allele_n, out=np.full_like(alt_sum, np.nan), where=allele_n > 0)
    missing_rate = 1.0 - (called_samples / max(matrix.shape[0], 1))
    return af, allele_n.astype(np.int64), missing_rate


def _parse_bvs_locus_key(locus_key: str) -> tuple[str, int, str, str] | None:
    text = str(locus_key)
    parts = text.split(":") if ":" in text else text.split("-")
    if len(parts) != 4:
        return None
    chrom, pos_raw, ref, alt = parts
    try:
        pos = int(pos_raw)
    except ValueError:
        return None
    return normalize_chrom(chrom), pos, ref.upper(), alt.upper()


def load_bvs_study_allele_freqs(
    path: Path,
    variants: list[Variant],
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"locus_key", "allele_number", "allele_freq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing bvs allele frequency column(s): {sorted(missing)}")

    key_to_idx = {
        (variant.chrom, variant.pos, variant.ref.upper(), variant.alt.upper()): idx
        for idx, variant in enumerate(variants)
    }
    af = np.full(len(variants), np.nan, dtype=np.float64)
    allele_n = np.zeros(len(variants), dtype=np.int64)
    stats = {"rows": 0, "matched": 0, "unmatched": 0, "duplicates": 0}

    seen: set[int] = set()
    for row in df.itertuples(index=False):
        stats["rows"] += 1
        row_map = row._asdict()
        parsed = _parse_bvs_locus_key(row_map["locus_key"])
        if parsed is None:
            stats["unmatched"] += 1
            continue
        idx = key_to_idx.get(parsed)
        if idx is None:
            stats["unmatched"] += 1
            continue
        if idx in seen:
            stats["duplicates"] += 1
            continue
        seen.add(idx)
        try:
            allele_n[idx] = int(row_map["allele_number"])
            af[idx] = float(row_map["allele_freq"])
        except (TypeError, ValueError):
            stats["unmatched"] += 1
            continue
        stats["matched"] += 1
    return af, allele_n, stats


def load_bvs_study_dosage(
    dosage_npy: Path,
    samples_tsv: Path,
    variants: list[Variant],
) -> tuple[list[str], np.ndarray]:
    samples_df = pd.read_csv(samples_tsv, sep="\t", dtype=str)
    if "sample_id" not in samples_df.columns:
        raise ValueError(f"{samples_tsv}: missing sample_id column")
    sample_ids = samples_df["sample_id"].astype(str).tolist()
    matrix = np.load(dosage_npy, mmap_mode=None)
    if matrix.dtype != np.uint8:
        raise ValueError(f"{dosage_npy}: expected uint8 matrix, got {matrix.dtype}")
    expected_shape = (len(sample_ids), len(variants))
    if matrix.shape != expected_shape:
        raise ValueError(f"{dosage_npy}: shape {matrix.shape} does not match expected {expected_shape}")
    return sample_ids, np.asarray(matrix, dtype=np.uint8)


def write_errors(path: Path, errors: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["participant_id", "file", "severity", "code", "message"])
        for row in errors:
            writer.writerow([
                row.get("participant_id", ""),
                row.get("file", ""),
                "ERROR",
                row.get("code", ""),
                row.get("message", ""),
            ])


def write_allele_freqs(
    path: Path,
    variants: list[Variant],
    ref_af: np.ndarray,
    ref_an: np.ndarray,
    study_af: np.ndarray,
    study_an: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chrom", "pos", "rsid", "ref", "alt",
            "hgp1k_af", "hgp1k_allele_number",
            "study_af", "study_allele_number",
        ])
        for i, variant in enumerate(variants):
            writer.writerow([
                variant.chrom,
                variant.pos,
                variant.rsid,
                variant.ref,
                variant.alt,
                "" if np.isnan(ref_af[i]) else f"{ref_af[i]:.8g}",
                int(ref_an[i]),
                "" if np.isnan(study_af[i]) else f"{study_af[i]:.8g}",
                int(study_an[i]),
            ])


def _float_or_blank(values: np.ndarray, idx: int, precision: str = ".8g") -> str:
    value = float(values[idx])
    if not np.isfinite(value):
        return ""
    return format(value, precision)


def write_pca_variant_files(
    used_path: Path,
    dropped_path: Path,
    variants: list[Variant],
    duplicate_dropped: list[Variant],
    keep: np.ndarray,
    af_ok: np.ndarray,
    ref_missing_ok: np.ndarray,
    combined_af: np.ndarray,
    combined_an: np.ndarray,
    ref_missing: np.ndarray,
    study_missing: np.ndarray,
    min_af: float,
    max_ref_missing: float,
) -> None:
    fieldnames = [
        "chrom",
        "pos",
        "rsid",
        "ref",
        "alt",
        "stage",
        "reason",
        "combined_af",
        "combined_allele_number",
        "ref_missing_rate",
        "study_missing_rate",
    ]

    with used_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for idx in np.flatnonzero(keep):
            variant = variants[int(idx)]
            writer.writerow({
                "chrom": variant.chrom,
                "pos": variant.pos,
                "rsid": variant.rsid,
                "ref": variant.ref,
                "alt": variant.alt,
                "stage": "pca",
                "reason": "used",
                "combined_af": _float_or_blank(combined_af, int(idx)),
                "combined_allele_number": int(combined_an[int(idx)]),
                "ref_missing_rate": _float_or_blank(ref_missing, int(idx), ".6f"),
                "study_missing_rate": _float_or_blank(study_missing, int(idx), ".6f"),
            })

    with dropped_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for variant in duplicate_dropped:
            writer.writerow({
                "chrom": variant.chrom,
                "pos": variant.pos,
                "rsid": variant.rsid,
                "ref": variant.ref,
                "alt": variant.alt,
                "stage": "duplicate_position",
                "reason": "duplicate_chrom_pos",
                "combined_af": "",
                "combined_allele_number": "",
                "ref_missing_rate": "",
                "study_missing_rate": "",
            })
        for idx in np.flatnonzero(~keep):
            i = int(idx)
            variant = variants[i]
            reasons: list[str] = []
            if not bool(af_ok[i]):
                if not np.isfinite(combined_af[i]):
                    reasons.append("combined_af_missing")
                elif combined_af[i] < min_af:
                    reasons.append("combined_af_below_min")
                elif combined_af[i] > (1.0 - min_af):
                    reasons.append("combined_af_above_max")
                else:
                    reasons.append("combined_af_failed")
            if not bool(ref_missing_ok[i]):
                reasons.append("ref_missing_above_max")
            writer.writerow({
                "chrom": variant.chrom,
                "pos": variant.pos,
                "rsid": variant.rsid,
                "ref": variant.ref,
                "alt": variant.alt,
                "stage": "pca_filter",
                "reason": ";".join(reasons) if reasons else "unknown",
                "combined_af": _float_or_blank(combined_af, i),
                "combined_allele_number": int(combined_an[i]),
                "ref_missing_rate": _float_or_blank(ref_missing, i, ".6f"),
                "study_missing_rate": _float_or_blank(study_missing, i, ".6f"),
            })


def fill_standardized_block(
    out: np.ndarray,
    start_row: int,
    matrix: np.ndarray,
    keep_idx: np.ndarray,
    mean: np.ndarray,
    sd: np.ndarray,
    chunk_size: int,
) -> float:
    total_ss = 0.0
    n_cols = keep_idx.size
    for start in range(0, n_cols, chunk_size):
        end = min(start + chunk_size, n_cols)
        idx = keep_idx[start:end]
        block = matrix[:, idx].astype(np.float32, copy=True)
        miss = block == float(MISSING)
        block -= mean[start:end][None, :].astype(np.float32)
        block /= sd[start:end][None, :].astype(np.float32)
        block[miss] = 0.0
        out[start_row:start_row + matrix.shape[0], start:end] = block
        total_ss += float(np.square(block, dtype=np.float64).sum())
    return total_ss


def load_metadata(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".csv":
        sep = ","
    else:
        sep = r"\s+"
    df = pd.read_csv(path, sep=sep, dtype=str)
    cols = {c.lower(): c for c in df.columns}
    sample_col = (
        cols.get("sample_id")
        or cols.get("sample")
        or cols.get("sampleid")
        or cols.get("iid")
    )
    if not sample_col:
        raise ValueError(f"{path}: expected sample_id/sample/SampleID/IID column")
    df = df.rename(columns={sample_col: "sample_id"})
    rename = {}
    for canonical, options in {
        "population": ("population", "pop"),
        "superpopulation": ("superpopulation", "super_pop", "genetic_region"),
    }.items():
        for option in options:
            col = cols.get(option)
            if col and col != canonical:
                rename[col] = canonical
                break
    if rename:
        df = df.rename(columns=rename)
    if "superpopulation" in df.columns:
        df["region"] = df["superpopulation"].map(lambda x: REGION_ALIASES.get(str(x).strip(), str(x).strip()))
    elif "region" not in df.columns:
        df["region"] = ""
    if "population" not in df.columns:
        df["population"] = ""
    return df


def subpopulation_colors(ref: pd.DataFrame):
    from matplotlib import cm

    colors = {}
    for region in REGION_ORDER:
        pops = sorted(ref.loc[ref["region"] == region, "population"].dropna().unique())
        pops = [pop for pop in pops if str(pop).strip()]
        if not pops:
            continue
        cmap = cm.get_cmap(REGION_CMAPS.get(region, "Greys"))
        for i, pop in enumerate(pops):
            shade = 0.35 + 0.60 * (i / max(1, len(pops) - 1))
            colors[(region, pop)] = cmap(shade)
    return colors


def plot_scores(scores_df: pd.DataFrame, out_path: Path, metadata: pd.DataFrame | None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    df = scores_df.copy()
    if metadata is not None:
        df = df.merge(metadata, on="sample_id", how="left")

    fig, ax = plt.subplots(figsize=(12, 9))
    ref = df[df["source"] == "1KGP"]
    study = df[df["source"] == "BioVault"]

    handles = []
    if metadata is not None and {"region", "population"}.issubset(ref.columns):
        colors = subpopulation_colors(ref)
        for region in REGION_ORDER:
            region_ref = ref[ref["region"] == region]
            if region_ref.empty:
                continue
            for pop in sorted(region_ref["population"].dropna().unique()):
                pop_ref = region_ref[region_ref["population"] == pop]
                if pop_ref.empty:
                    continue
                color = colors.get((region, pop), "#8a8f98")
                ax.scatter(
                    pop_ref["PC1"],
                    pop_ref["PC2"],
                    s=13,
                    alpha=0.65,
                    color=color,
                    edgecolors="none",
                )
                handles.append(Line2D(
                    [0], [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=color,
                    markeredgecolor="none",
                    markersize=5,
                    label=f"[{region}] {pop} (n={len(pop_ref)})",
                ))
    else:
        ax.scatter(ref["PC1"], ref["PC2"], c="#8a8f98", s=9, alpha=0.35, label="1KGP")
        handles.append(Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor="#8a8f98",
            markeredgecolor="none", markersize=5, label=f"1KGP (n={len(ref)})",
        ))

    ax.scatter(
        study["PC1"], study["PC2"],
        c="#111111", marker="*", s=180,
        edgecolors="#f2c94c", linewidths=0.8,
    )
    handles.append(Line2D(
        [0], [0],
        marker="*",
        linestyle="",
        markerfacecolor="#111111",
        markeredgecolor="#f2c94c",
        markeredgewidth=0.8,
        markersize=13,
        label=f"BioVault (n={len(study)})",
    ))
    for _, row in study.iterrows():
        ax.annotate(
            str(row["sample_id"]),
            (row["PC1"], row["PC2"]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.axhline(0, color="#b8bcc4", linewidth=0.6)
    ax.axvline(0, color="#b8bcc4", linewidth=0.6)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("BioVault samples and 1000 Genomes high-coverage PCA")
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
        ncol=2 if len(handles) > 30 else 1,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_closest_population_labels(
    scores_df: pd.DataFrame,
    metadata: pd.DataFrame | None,
    out_path: Path,
    summary_path: Path,
    n_pcs: int = 2,
    top_n: int = 5,
    tie_ratio: float = 1.15,
) -> None:
    if metadata is None:
        with out_path.open("w", encoding="utf-8") as handle:
            handle.write("sample_id\trank\tregion\tpopulation\tlabel\tn_ref\tdistance\n")
        with summary_path.open("w", encoding="utf-8") as handle:
            handle.write("sample_id\tclosest_labels\n")
        return

    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1) if f"PC{i}" in scores_df.columns]
    if len(pc_cols) < 2:
        raise ValueError("closest population labels need at least PC1 and PC2")

    df = scores_df.merge(metadata, on="sample_id", how="left")
    ref = df[
        (df["source"] == "1KGP")
        & df["region"].fillna("").astype(str).str.strip().ne("")
        & df["population"].fillna("").astype(str).str.strip().ne("")
    ].copy()
    study = df[df["source"] != "1KGP"].copy()

    if ref.empty:
        raise ValueError("metadata did not match any 1KGP sample ids for closest-label output")

    centroids = (
        ref.groupby(["region", "population"], dropna=False)[pc_cols]
        .mean()
        .reset_index()
    )
    counts = ref.groupby(["region", "population"], dropna=False).size().rename("n_ref").reset_index()
    centroids = centroids.merge(counts, on=["region", "population"], how="left")

    rows: list[dict[str, object]] = []
    summary_rows: list[tuple[str, str]] = []
    centroid_values = centroids[pc_cols].to_numpy(dtype=np.float64)
    for _, sample in study.iterrows():
        sample_id = str(sample["sample_id"])
        point = sample[pc_cols].to_numpy(dtype=np.float64)
        distances = np.sqrt(((centroid_values - point[None, :]) ** 2).sum(axis=1))
        order = np.argsort(distances)
        best = float(distances[order[0]])
        close_cutoff = best * tie_ratio if best > 0 else np.nextafter(0.0, 1.0)
        close_labels: list[str] = []
        for rank, idx in enumerate(order[:top_n], start=1):
            row = centroids.iloc[int(idx)]
            region = str(row["region"])
            pop = str(row["population"])
            label = f"[{region}] {pop}"
            distance = float(distances[idx])
            rows.append({
                "sample_id": sample_id,
                "rank": rank,
                "region": region,
                "population": pop,
                "label": label,
                "n_ref": int(row["n_ref"]),
                "distance": distance,
            })
            if distance <= close_cutoff:
                close_labels.append(label)
        summary_rows.append((sample_id, ";".join(close_labels)))

    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False, float_format="%.8g")
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("sample_id\tclosest_labels\n")
        for sample_id, labels in summary_rows:
            handle.write(f"{sample_id}\t{labels}\n")


def main() -> None:
    setup_logging()
    log = logging.getLogger(__name__)

    repo_root = next(
        (p for p in Path(__file__).resolve().parents if (p / "tools" / "locus_map.tsv").exists()),
        Path.cwd(),
    )
    default_vcf_dir = repo_root / "data" / "1kgp_high_coverage" / "filtered"
    default_metadata = repo_root / "data" / "1kgp_high_coverage" / "20130606_g1k_3202_samples_ped_population.txt"

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--vcf-dir", type=Path, default=Path(os.environ.get("HGP1K_VCF_DIR", default_vcf_dir)))
    ap.add_argument("--matrix-npz", type=Path, default=Path(os.environ["HGP1K_MATRIX_NPZ"]) if os.environ.get("HGP1K_MATRIX_NPZ") else None)
    ap.add_argument("--locus-map", type=Path, default=Path(os.environ.get("LOCUS_MAP", repo_root / "tools" / "locus_map.tsv")))
    ap.add_argument("--chromosomes", default=os.environ.get("CHROMOSOMES", "all"))
    ap.add_argument(
        "--metadata-tsv",
        type=Path,
        default=Path(os.environ["HGP1K_METADATA_TSV"]) if os.environ.get("HGP1K_METADATA_TSV") else (
            default_metadata if default_metadata.exists() else None
        ),
    )
    ap.add_argument("--n-components", type=int, default=10)
    ap.add_argument("--min-gs", type=float, default=0.15)
    ap.add_argument("--min-af", type=float, default=0.01)
    ap.add_argument("--max-ref-missing", type=float, default=0.05)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--study-workers", type=int, default=max(1, (os.cpu_count() or 1)))
    ap.add_argument("--max-variants", type=int, default=None, help="debug/testing cap after loading VCF SNP rows")
    ap.add_argument("--study-af-tsv", type=Path, default=None, help="Optional bvs fast-allele-freq TSV for study AF output")
    ap.add_argument("--study-dosage-npy", type=Path, default=None, help="Optional bvs target-panel-aligned uint8 study dosage .npy")
    ap.add_argument("--study-samples-tsv", type=Path, default=None, help="Sample metadata TSV for --study-dosage-npy")
    ap.add_argument("--random-state", type=int, default=7)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    chromosomes = chromosome_list(args.chromosomes)
    vcfs: list[Path] = []
    if args.matrix_npz:
        log.info("loading compact 1KGP matrix: %s", args.matrix_npz)
        log.info("compact 1KGP matrix filename: %s", args.matrix_npz.name)
        ref_sample_ids, variants, ref_matrix = load_reference_matrix_npz(args.matrix_npz)
        if args.max_variants is not None:
            variants = variants[:args.max_variants]
            ref_matrix = ref_matrix[:, :args.max_variants]
        log.info("using compact matrix instead of VCF dir")
    else:
        vcfs = find_vcfs(args.vcf_dir, chromosomes)
        log.info("using chromosomes: %s", ",".join(chromosomes))
        log.info("using %d filtered 1KGP VCF(s)", len(vcfs))
        locus_rsids = load_locus_rsids(args.locus_map)
        ref_sample_ids, variants, ref_matrix = stream_vcf_matrix(vcfs, locus_rsids, args.max_variants)
    loaded_variant_count = len(variants)
    variants, ref_matrix, duplicate_dropped = drop_duplicate_loci(variants, ref_matrix)
    n_dup_loci = len(duplicate_dropped)
    log.info("1KGP matrix after duplicate-locus drop: %s samples x %s variants", ref_matrix.shape[0], ref_matrix.shape[1])

    sample_files = discover_study_files(args.data_dir)
    study_workers = max(1, min(args.study_workers, len(sample_files)))
    if args.study_dosage_npy:
        if not args.study_samples_tsv:
            raise ValueError("--study-samples-tsv is required with --study-dosage-npy")
        log.info("loading bvs study dosage matrix: %s", args.study_dosage_npy)
        study_sample_ids, study_matrix = load_bvs_study_dosage(args.study_dosage_npy, args.study_samples_tsv, variants)
        errors = []
    else:
        log.info("study parse workers: %s", study_workers)
        study_sample_ids, study_matrix, errors = load_study_matrix(sample_files, variants, args.min_gs, study_workers)
    write_errors(args.out_dir / "errors.tsv", errors)
    log.info("study matrix: %s samples x %s variants", study_matrix.shape[0], study_matrix.shape[1])

    ref_af, ref_an, ref_missing = allele_stats(ref_matrix)
    study_af_from_matrix, study_an_from_matrix, study_missing = allele_stats(study_matrix)
    if args.study_af_tsv:
        study_af, study_an, bvs_af_stats = load_bvs_study_allele_freqs(args.study_af_tsv, variants)
        log.info(
            "bvs study AF: rows=%s matched=%s unmatched=%s duplicates=%s",
            f"{bvs_af_stats['rows']:,}",
            f"{bvs_af_stats['matched']:,}",
            f"{bvs_af_stats['unmatched']:,}",
            f"{bvs_af_stats['duplicates']:,}",
        )
    else:
        study_af, study_an = study_af_from_matrix, study_an_from_matrix
    write_allele_freqs(args.out_dir / "allele_freqs.tsv", variants, ref_af, ref_an, study_af, study_an)

    combined = np.concatenate([ref_matrix, study_matrix], axis=0)
    combined_af, combined_an, _combined_missing = allele_stats(combined)
    af_ok = np.isfinite(combined_af) & (combined_af >= args.min_af) & (combined_af <= (1.0 - args.min_af))
    ref_missing_ok = ref_missing <= args.max_ref_missing
    keep = af_ok & ref_missing_ok
    keep_idx = np.flatnonzero(keep)
    if keep_idx.size < 2:
        raise SystemExit("ERROR: fewer than two variants passed PCA filters")
    write_pca_variant_files(
        args.out_dir / "pca_variants_used.tsv",
        args.out_dir / "pca_variants_dropped.tsv",
        variants,
        duplicate_dropped,
        keep,
        af_ok,
        ref_missing_ok,
        combined_af,
        combined_an,
        ref_missing,
        study_missing,
        args.min_af,
        args.max_ref_missing,
    )

    n_components = min(args.n_components, combined.shape[0] - 1, keep_idx.size)
    if n_components < 2:
        raise SystemExit("ERROR: need at least two samples and two variants for PC1/PC2")

    mean = (2.0 * combined_af[keep_idx]).astype(np.float32)
    sd = np.sqrt(2.0 * combined_af[keep_idx] * (1.0 - combined_af[keep_idx])).astype(np.float32)
    sd = np.where(sd > 0, sd, 1.0).astype(np.float32)

    matrix_path = args.out_dir / "combined_standardized.float32.mmap"
    X = np.memmap(matrix_path, dtype=np.float32, mode="w+", shape=(combined.shape[0], keep_idx.size))
    ss_ref = fill_standardized_block(X, 0, ref_matrix, keep_idx, mean, sd, args.chunk_size)
    ss_study = fill_standardized_block(X, ref_matrix.shape[0], study_matrix, keep_idx, mean, sd, args.chunk_size)
    X.flush()

    from sklearn.utils.extmath import randomized_svd
    log.info("running randomized SVD PCA: %s samples x %s variants", X.shape[0], X.shape[1])
    U, S, _Vt = randomized_svd(
        X,
        n_components=n_components,
        n_iter=5,
        random_state=args.random_state,
    )
    scores = U * S[None, :]
    explained_var = (S ** 2) / max(X.shape[0] - 1, 1)
    total_var = (ss_ref + ss_study) / max(X.shape[0] - 1, 1)
    explained_ratio = explained_var / total_var if total_var > 0 else np.full_like(explained_var, np.nan)

    all_sample_ids = ref_sample_ids + study_sample_ids
    sources = ["1KGP"] * len(ref_sample_ids) + ["BioVault"] * len(study_sample_ids)
    score_cols = [f"PC{i}" for i in range(1, n_components + 1)]
    scores_df = pd.DataFrame(scores, columns=score_cols)
    scores_df.insert(0, "source", sources)
    scores_df.insert(0, "sample_id", all_sample_ids)
    scores_df.to_csv(args.out_dir / "pca_scores.tsv", sep="\t", index=False)

    with (args.out_dir / "study_pca_projection.tsv").open("w", encoding="utf-8") as handle:
        handle.write("s\tscores\n")
        for sid, row in zip(study_sample_ids, scores[len(ref_sample_ids):]):
            handle.write(f"{sid}\t[{','.join(repr(float(v)) for v in row)}]\n")

    metadata = load_metadata(args.metadata_tsv)
    plot_scores(scores_df, args.out_dir / "pca_projection.png", metadata)
    write_closest_population_labels(
        scores_df,
        metadata,
        args.out_dir / "closest_population_labels.tsv",
        args.out_dir / "closest_population_summary.tsv",
    )

    with (args.out_dir / "qc_report.txt").open("w", encoding="utf-8") as handle:
        handle.write("=== hgp1k_projection_fast QC ===\n")
        handle.write(f"VCF dir: {args.vcf_dir}\n")
        handle.write(f"Matrix NPZ: {args.matrix_npz or ''}\n")
        handle.write(f"Matrix NPZ filename: {args.matrix_npz.name if args.matrix_npz else ''}\n")
        handle.write(f"Matrix NPZ resolved: {args.matrix_npz.resolve() if args.matrix_npz else ''}\n")
        handle.write(f"BVS study AF TSV: {args.study_af_tsv or ''}\n")
        handle.write(f"BVS study dosage NPY: {args.study_dosage_npy or ''}\n")
        handle.write(f"BVS study samples TSV: {args.study_samples_tsv or ''}\n")
        handle.write(f"Locus map: {args.locus_map}\n")
        handle.write(f"Chromosomes requested: {','.join(chromosomes)}\n")
        handle.write(f"VCFs used: {len(vcfs)}\n")
        handle.write(f"1KGP samples: {len(ref_sample_ids):,}\n")
        handle.write(f"Study samples: {len(study_sample_ids):,}\n")
        handle.write(f"Study parse workers: {0 if args.study_dosage_npy else study_workers}\n")
        handle.write(f"Loaded biallelic SNP variants: {loaded_variant_count:,}\n")
        handle.write(f"Duplicate-position variants dropped: {n_dup_loci:,}\n")
        handle.write(f"Variants after duplicate-position drop: {len(variants):,}\n")
        handle.write(f"PCA variants after AF/ref-missing filters: {keep_idx.size:,}\n")
        handle.write(f"PCA variants dropped by AF/ref-missing filters: {int((~keep).sum()):,}\n")
        handle.write("PCA variant audit files: pca_variants_used.tsv, pca_variants_dropped.tsv\n")
        handle.write(f"min_af={args.min_af}, max_ref_missing={args.max_ref_missing}, min_gs={args.min_gs}\n")
        handle.write(f"Mean 1KGP missing rate across kept variants: {float(ref_missing[keep_idx].mean()):.6f}\n")
        handle.write(f"Mean study missing rate across kept variants: {float(study_missing[keep_idx].mean()):.6f}\n")
        if args.study_af_tsv:
            handle.write(
                "BVS study AF rows/matched/unmatched/duplicates: "
                f"{bvs_af_stats['rows']}/{bvs_af_stats['matched']}/"
                f"{bvs_af_stats['unmatched']}/{bvs_af_stats['duplicates']}\n"
            )
        handle.write("Explained variance ratio:\n")
        for i, value in enumerate(explained_ratio, start=1):
            handle.write(f"  PC{i}: {float(value):.8g}\n")
        handle.write(f"Elapsed seconds: {time.time() - t_start:.1f}\n")

    try:
        matrix_path.unlink()
    except OSError:
        pass

    log.info("done: outputs in %s", args.out_dir)


if __name__ == "__main__":
    main()

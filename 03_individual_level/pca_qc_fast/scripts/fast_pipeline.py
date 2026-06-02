#!/usr/bin/env python3
"""
Fast pca_qc implementation for large cohorts.

The original pca_qc path was intentionally simple: merge every genotype into
SNP x sample TSVs, encode that TSV into another dense TSV, then run PCA. That is
fine for the documented 10-sample smoke test, but it is the wrong storage model
for 1000+ samples. This implementation keeps the same published outputs while
using compact on-disk arrays internally:

  * pass 1 builds the full cohort SNP universe from all readable files;
  * pass 2 fills uint8 allele memmaps, so every valid file contributes without
    materialising a pandas string matrix;
  * QC and LD pruning run in chunks over an int8 dosage memmap;
  * PCA is computed from a chunked sample Gram matrix instead of a dense
    samples x SNPs float matrix;
  * PLINK BED/BIM/FAM files are written directly for external inspection.

Set BV_WRITE_MATRICES=1 to also emit the legacy giant TSV matrices.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import sys
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

for _parent in Path(__file__).resolve().parents:
    if (_parent / "tools" / "genotype_normalizer.py").exists():
        sys.path.insert(0, str(_parent))
        break
sys.path.append("/opt/biovault")
from tools import genotype_normalizer as _genoio  # noqa: E402


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("BIOVAULT_DATA_DIR", BASE_DIR.parent)).resolve()
MERGED_DIR = BASE_DIR / "data" / "merged"
PLINK_DIR = BASE_DIR / "data" / "plink"
PCA_DIR = BASE_DIR / "data" / "pca"
PLOTS_DIR = BASE_DIR / "plots"
LOG_DIR = BASE_DIR / "logs"
QC_DIR = BASE_DIR / "data" / "qc"
WORK_DIR = BASE_DIR / "data" / "work"
ERRORS_TSV = LOG_DIR / "errors.tsv"
WARNINGS_TSV = LOG_DIR / "warnings.tsv"
FILTERED_SNPS_TSV = QC_DIR / "filtered_snps.tsv"

N_PCS = int(os.environ.get("BV_N_PCS", "20"))
GENO = float(os.environ.get("BV_GENO", "0.05"))
MIND = float(os.environ.get("BV_MIND", "0.10"))
MAF = float(os.environ.get("BV_MAF", "0.01"))
HWE_P = float(os.environ.get("BV_HWE_P", "1e-4"))
LD_WINDOW = int(os.environ.get("BV_LD_WINDOW", "50"))
LD_STEP = int(os.environ.get("BV_LD_STEP", "5"))
LD_R2 = float(os.environ.get("BV_LD_R2", "0.2"))
CHUNK_VARIANTS = int(os.environ.get("BV_CHUNK_VARIANTS", "4096"))
WRITE_MATRICES = os.environ.get("BV_WRITE_MATRICES", "0").strip().lower() in {"1", "true", "yes"}
PLINK_BACKEND = os.environ.get("BV_PCA_BACKEND", "auto").strip().lower()
PARSE_MODE = os.environ.get("BV_PARSE_MODE", "process").strip().lower()
DEFAULT_WORKERS = os.cpu_count() or 1
WORKERS = max(1, int(os.environ.get("BV_WORKERS", str(DEFAULT_WORKERS))))

BASES = np.array([b"A", b"C", b"G", b"T"], dtype="S1")
BASE_LABELS = np.array([".", "A", "C", "G", "T"], dtype=object)
BASE_TO_CODE = {b"A": 1, b"C": 2, b"G": 3, b"T": 4}
MISSING = np.int8(-1)
FILTERED_SNPS_COLUMNS = [
    "variant_index",
    "rsid",
    "chromosome",
    "position",
    "filter",
    "n_homref",
    "n_het",
    "n_homalt",
    "n_missing",
    "frac_homref",
    "frac_het",
    "frac_homalt",
    "call_rate",
    "maf",
    "hwe_p",
]


def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s  %(message)s",
        handlers=[
            logging.FileHandler(LOG_DIR / "fast_pipeline.log", mode="w"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger("pca_qc_fast")


log = setup_logging()
os.environ.setdefault("BIOVAULT_WARNINGS_TSV", str(WARNINGS_TSV))
os.environ.setdefault("BIOVAULT_FAST_NORMALIZE", "1")


def timed(label: str):
    class Timer:
        def __enter__(self):
            self.start = time.perf_counter()
            log.info("%s ...", label)
            return self

        def __exit__(self, exc_type, exc, tb):
            if exc_type is None:
                log.info("%s done in %.2fs", label, time.perf_counter() - self.start)

    return Timer()


def write_errors(rows: list[dict[str, str]]) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with ERRORS_TSV.open("w", encoding="utf-8") as handle:
        handle.write("participant_id\tfile\tseverity\tcode\tmessage\n")
        for row in rows:
            msg = row["message"].replace("\t", " ").replace("\n", " ")
            handle.write(f"{row['participant_id']}\t{row['file']}\tERROR\t{row['code']}\t{msg}\n")


def append_error(participant_id: str, file: str, code: str, message: str) -> None:
    if not ERRORS_TSV.exists():
        write_errors([])
    with ERRORS_TSV.open("a", encoding="utf-8") as handle:
        msg = message.replace("\t", " ").replace("\n", " ")
        handle.write(f"{participant_id}\t{file}\tERROR\t{code}\t{msg}\n")


def write_empty_pca_outputs(reason: str) -> None:
    PCA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    (PCA_DIR / "pca.eigenvec").write_text("", encoding="utf-8")
    (PCA_DIR / "pca.eigenval").write_text("", encoding="utf-8")
    append_error("COHORT", str(DATA_DIR), "INSUFFICIENT_USABLE_DATA", reason)
    log.error(reason)


def ensure_filtered_snps_file() -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    if not FILTERED_SNPS_TSV.exists():
        pd.DataFrame(columns=FILTERED_SNPS_COLUMNS).to_csv(FILTERED_SNPS_TSV, sep="\t", index=False)


def write_filtered_snps(df: pd.DataFrame | None = None) -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    if df is None or df.empty:
        pd.DataFrame(columns=FILTERED_SNPS_COLUMNS).to_csv(FILTERED_SNPS_TSV, sep="\t", index=False)
        return
    for col in FILTERED_SNPS_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df[FILTERED_SNPS_COLUMNS].to_csv(FILTERED_SNPS_TSV, sep="\t", index=False)


def discover_samples() -> list[tuple[str, Path]]:
    samples: list[tuple[str, Path]] = []
    for sample_dir in sorted(DATA_DIR.iterdir()):
        if not sample_dir.is_dir():
            continue
        txt_files = sorted(sample_dir.glob("*.txt"))
        if txt_files:
            samples.append((sample_dir.name, txt_files[0]))
    if not samples:
        raise RuntimeError(f"No sample directories with .txt files found under {DATA_DIR}")
    return samples


def read_sample(sample_id: str, path: Path) -> pd.DataFrame:
    with path.open("rb") as handle:
        if b"\x00" in handle.read(65536):
            raise ValueError("file contains NUL bytes; refusing binary-like genotype input")
    g = _genoio.read_pipeline_genotypes(path)
    if g.empty:
        return pd.DataFrame(columns=["rsid", "chromosome", "position", "genotype"])
    df = g[["rsid", "chrom", "pos", "gt"]].rename(
        columns={"chrom": "chromosome", "pos": "position", "gt": "genotype"}
    )
    df = df.dropna(subset=["rsid", "genotype"]).copy()
    df["rsid"] = df["rsid"].astype(str)
    df["chromosome"] = df["chromosome"].astype(str)
    df["position"] = pd.to_numeric(df["position"], errors="coerce").fillna(0).astype(np.int64)
    df["genotype"] = df["genotype"].astype(str).str.upper().str.strip()
    df = df.drop_duplicates("rsid", keep="first")
    return df


def _cache_path_for(cache_dir: Path, idx: int, sample_id: str) -> Path:
    safe = "".join(c if c.isalnum() or c in {"_", "-", "."} else "_" for c in sample_id)
    return cache_dir / f"{idx:06d}_{safe}.npz"


def parse_and_cache_sample(task: tuple[int, str, str, str]) -> dict[str, str | int | None]:
    idx, sample_id, path_s, cache_path_s = task
    path = Path(path_s)
    cache_path = Path(cache_path_s)
    try:
        df = read_sample(sample_id, path)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            cache_path,
            rsid=df["rsid"].to_numpy(dtype="S32"),
            chrom=df["chromosome"].to_numpy(dtype="S2"),
            pos=df["position"].to_numpy(dtype=np.int64),
            gt=df["genotype"].to_numpy(dtype="S2"),
        )
        return {
            "idx": idx,
            "sample_id": sample_id,
            "file": str(path),
            "cache": str(cache_path),
            "rows": int(len(df)),
            "error": None,
        }
    except Exception as exc:
        return {
            "idx": idx,
            "sample_id": sample_id,
            "file": str(path),
            "cache": str(cache_path),
            "rows": 0,
            "error": str(exc),
        }


def parse_samples_to_cache(samples: list[tuple[str, Path]]) -> tuple[list[dict[str, str | int]], list[dict[str, str]]]:
    cache_dir = WORK_DIR / "parsed_samples"
    cache_dir.mkdir(parents=True, exist_ok=True)
    for old in cache_dir.glob("*.npz"):
        old.unlink()

    workers = max(1, min(len(samples), WORKERS))
    parse_mode = PARSE_MODE if PARSE_MODE in {"process", "thread", "serial"} else "process"
    if parse_mode != PARSE_MODE:
        log.warning("Unknown BV_PARSE_MODE=%r; using process", PARSE_MODE)

    tasks = [
        (idx, sample_id, str(path), str(_cache_path_for(cache_dir, idx, sample_id)))
        for idx, (sample_id, path) in enumerate(samples)
    ]
    results: list[dict[str, str | int]] = []
    errors: list[dict[str, str]] = []
    start = time.perf_counter()

    with timed(f"Parsing {len(samples)} samples to cache with {workers} {parse_mode} workers"):
        completed = 0
        if workers > 1 and parse_mode != "serial":
            executor_cls = ProcessPoolExecutor if parse_mode == "process" else ThreadPoolExecutor
            try:
                pool_ctx = executor_cls(max_workers=workers)
            except PermissionError as exc:
                if parse_mode != "process":
                    raise
                log.warning("Process workers unavailable (%s); falling back to thread workers", exc)
                parse_mode = "thread"
                pool_ctx = ThreadPoolExecutor(max_workers=workers)
            with pool_ctx as pool:
                futures = {pool.submit(parse_and_cache_sample, task): task for task in tasks}
                for future in as_completed(futures):
                    row = future.result()
                    completed += 1
                    if row["error"]:
                        log.error("Skipping %s: %s", row["sample_id"], row["error"])
                        errors.append({
                            "participant_id": str(row["sample_id"]),
                            "file": str(row["file"]),
                            "code": "PARSE_FAILED",
                            "message": str(row["error"]),
                        })
                    else:
                        results.append(row)  # type: ignore[arg-type]
                    log_progress("Parsed", completed, len(samples), start)
        else:
            for task in tasks:
                row = parse_and_cache_sample(task)
                completed += 1
                if row["error"]:
                    log.error("Skipping %s: %s", row["sample_id"], row["error"])
                    errors.append({
                        "participant_id": str(row["sample_id"]),
                        "file": str(row["file"]),
                        "code": "PARSE_FAILED",
                        "message": str(row["error"]),
                    })
                else:
                    results.append(row)  # type: ignore[arg-type]
                log_progress("Parsed", completed, len(samples), start)

    results.sort(key=lambda row: int(row["idx"]))
    write_errors(errors)
    if not results:
        raise RuntimeError("No genotype files could be parsed; see errors.tsv")
    log.info("Parsed cache rows: %d samples, %d genotype rows", len(results), sum(int(r["rows"]) for r in results))
    return results, errors


def log_progress(label: str, done: int, total: int, start: float) -> None:
    if done % 50 != 0 and done != total:
        return
    elapsed = max(time.perf_counter() - start, 1e-6)
    rate = done / elapsed
    eta = (total - done) / max(rate, 1e-6)
    log.info("%s %d/%d (%.1f files/s, ETA %.0fs)", label, done, total, rate, eta)


def genotype_codes(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gt = values.astype("S2", copy=False)
    chars = gt.view("S1").reshape(-1, 2)
    a1 = np.zeros(chars.shape[0], dtype=np.uint8)
    a2 = np.zeros(chars.shape[0], dtype=np.uint8)
    for code, base in enumerate(BASES, start=1):
        a1[chars[:, 0] == base] = code
        a2[chars[:, 1] == base] = code
    called = (a1 != 0) & (a2 != 0)
    return a1, a2, called


def build_snp_universe(samples: list[tuple[str, Path]]) -> tuple[list[str], pd.DataFrame, list[dict[str, str | int]], list[dict[str, str]]]:
    snp_to_idx: dict[str, int] = {}
    records: list[tuple[str, str, int]] = []
    parsed, errors = parse_samples_to_cache(samples)
    start = time.perf_counter()
    with timed(f"Pass 1: building full SNP universe from {len(parsed)} cached samples"):
        for i, row in enumerate(parsed, 1):
            data = np.load(str(row["cache"]), allow_pickle=False)
            rsids = data["rsid"].astype(str)
            chroms = data["chrom"].astype(str)
            for rsid, chrom, pos in zip(rsids, chroms, data["pos"]):
                if rsid not in snp_to_idx:
                    snp_to_idx[rsid] = len(records)
                    records.append((rsid, str(chrom), int(pos)))
            log_progress("Scanned", i, len(parsed), start)

    if not records:
        raise RuntimeError("Parsed files contained no usable genotype rows")

    snp_info = pd.DataFrame(records, columns=["rsid", "chromosome", "position"])
    MERGED_DIR.mkdir(parents=True, exist_ok=True)
    snp_info.to_csv(MERGED_DIR / "snp_info.tsv", sep="\t", index=False)
    log.info("Full SNP universe: %d SNPs x %d samples", len(records), len(parsed))
    return [r[0] for r in records], snp_info, parsed, errors


def fill_allele_memmaps(
    samples: list[dict[str, str | int]],
    snp_ids: list[str],
) -> tuple[np.memmap, np.memmap, list[str]]:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    for name in ("allele1.uint8.mmap", "allele2.uint8.mmap", "dosage.int8.mmap"):
        path = WORK_DIR / name
        if path.exists():
            path.unlink()
    n_samples = len(samples)
    n_snps = len(snp_ids)
    a1_mm = np.memmap(WORK_DIR / "allele1.uint8.mmap", dtype=np.uint8, mode="w+", shape=(n_samples, n_snps))
    a2_mm = np.memmap(WORK_DIR / "allele2.uint8.mmap", dtype=np.uint8, mode="w+", shape=(n_samples, n_snps))
    a1_mm[:] = 0
    a2_mm[:] = 0
    snp_to_idx = {rsid: idx for idx, rsid in enumerate(snp_ids)}
    sample_ids = [str(row["sample_id"]) for row in samples]

    start = time.perf_counter()
    with timed(f"Pass 2: filling compact allele memmaps for {n_samples} samples"):
        for i, row in enumerate(samples):
            data = np.load(str(row["cache"]), allow_pickle=False)
            idx = pd.Series(data["rsid"].astype(str), copy=False).map(snp_to_idx)
            valid = idx.notna().to_numpy()
            if valid.any():
                ix = idx[valid].to_numpy(dtype=np.int64)
                c1, c2, _called = genotype_codes(data["gt"][valid])
                a1_mm[i, ix] = c1
                a2_mm[i, ix] = c2
            log_progress("Loaded", i + 1, n_samples, start)
    a1_mm.flush()
    a2_mm.flush()
    return a1_mm, a2_mm, sample_ids


def allele_counts(a1_mm: np.memmap, a2_mm: np.memmap) -> np.ndarray:
    n_samples, n_snps = a1_mm.shape
    counts = np.zeros((n_snps, 4), dtype=np.int32)
    with timed("Counting alleles across full matrix"):
        for start in range(0, n_snps, CHUNK_VARIANTS):
            end = min(start + CHUNK_VARIANTS, n_snps)
            a1 = np.asarray(a1_mm[:, start:end])
            a2 = np.asarray(a2_mm[:, start:end])
            called = (a1 != 0) & (a2 != 0)
            for code in range(1, 5):
                counts[start:end, code - 1] = (
                    ((a1 == code) & called).sum(axis=0)
                    + ((a2 == code) & called).sum(axis=0)
                )
    return counts


def choose_biallelic_variants(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    present = counts > 0
    n_distinct = present.sum(axis=1)
    keep = (n_distinct >= 1) & (n_distinct <= 2)
    keep_idx = np.flatnonzero(keep)
    if keep_idx.size == 0:
        return keep_idx, np.empty(0, dtype=np.uint8), np.empty(0, dtype=np.uint8)

    kept_counts = counts[keep_idx]
    order = np.argsort(-kept_counts, axis=1, kind="stable")
    major = (order[:, 0] + 1).astype(np.uint8)
    minor = (order[:, 1] + 1).astype(np.uint8)
    mono = (kept_counts[np.arange(keep_idx.size), order[:, 1]] == 0)
    minor[mono] = major[mono]
    log.info(
        "Biallelic/monomorphic variants retained for dosage/PLINK: %d/%d (multiallelic or empty dropped: %d)",
        keep_idx.size,
        counts.shape[0],
        counts.shape[0] - keep_idx.size,
    )
    return keep_idx, major, minor


def build_dosage_memmap(
    a1_mm: np.memmap,
    a2_mm: np.memmap,
    keep_idx: np.ndarray,
    major: np.ndarray,
    minor: np.ndarray,
) -> np.memmap:
    n_samples = a1_mm.shape[0]
    dosage = np.memmap(WORK_DIR / "dosage.int8.mmap", dtype=np.int8, mode="w+", shape=(keep_idx.size, n_samples))
    with timed("Encoding compact int8 dosage memmap"):
        for out_start in range(0, keep_idx.size, CHUNK_VARIANTS):
            out_end = min(out_start + CHUNK_VARIANTS, keep_idx.size)
            src = keep_idx[out_start:out_end]
            a1 = np.asarray(a1_mm[:, src]).T
            a2 = np.asarray(a2_mm[:, src]).T
            maj = major[out_start:out_end, None]
            minr = minor[out_start:out_end, None]
            called = (a1 != 0) & (a2 != 0)
            valid = called & (
                ((a1 == maj) | (a1 == minr))
                & ((a2 == maj) | (a2 == minr))
            )
            d = ((a1 == minr).astype(np.int8) + (a2 == minr).astype(np.int8))
            d[~valid] = MISSING
            dosage[out_start:out_end, :] = d
    dosage.flush()
    return dosage


def _chrom_code(value: str) -> str:
    chrom = str(value).replace("chr", "").replace("CHR", "")
    return (
        chrom.replace("XY", "25")
        .replace("MT", "26")
        .replace("M", "26")
        .replace("X", "23")
        .replace("Y", "24")
    )


def write_plink_binary(
    dosage: np.memmap,
    snp_info: pd.DataFrame,
    keep_idx: np.ndarray,
    major: np.ndarray,
    minor: np.ndarray,
    sample_ids: list[str],
) -> None:
    PLINK_DIR.mkdir(parents=True, exist_ok=True)
    prefix = PLINK_DIR / "genotypes"
    with timed("Writing PLINK BED/BIM/FAM directly"):
        with (prefix.with_suffix(".fam")).open("w", encoding="utf-8") as handle:
            for sid in sample_ids:
                handle.write(f"{sid}\t{sid}\t0\t0\t0\t-9\n")

        kept_info = snp_info.iloc[keep_idx].reset_index(drop=True)
        with (prefix.with_suffix(".bim")).open("w", encoding="utf-8") as handle:
            for row, maj, minr in zip(kept_info.itertuples(index=False), major, minor):
                handle.write(
                    f"{_chrom_code(row.chromosome)}\t{row.rsid}\t0\t{int(row.position)}\t"
                    f"{BASE_LABELS[minr]}\t{BASE_LABELS[maj]}\n"
                )

        n_var, n_samples = dosage.shape
        bytes_per_var = (n_samples + 3) // 4
        with (prefix.with_suffix(".bed")).open("wb") as handle:
            handle.write(b"\x6c\x1b\x01")
            for start in range(0, n_var, CHUNK_VARIANTS):
                end = min(start + CHUNK_VARIANTS, n_var)
                d = np.asarray(dosage[start:end, :], dtype=np.int8)
                code = np.where(
                    d == MISSING,
                    0b01,
                    np.where(d == 2, 0b00, np.where(d == 1, 0b10, 0b11)),
                ).astype(np.uint8)
                pad = bytes_per_var * 4 - n_samples
                if pad:
                    code = np.concatenate([code, np.zeros((code.shape[0], pad), dtype=np.uint8)], axis=1)
                reshaped = code.reshape(code.shape[0], bytes_per_var, 4)
                bed_bytes = (
                    reshaped[:, :, 0]
                    | (reshaped[:, :, 1] << 2)
                    | (reshaped[:, :, 2] << 4)
                    | (reshaped[:, :, 3] << 6)
                ).astype(np.uint8)
                handle.write(bed_bytes.tobytes())

        # Legacy scripts expect PED/MAP names to exist; keep MAP cheap and make
        # PED an explicit placeholder instead of generating a multi-GB text file.
        with (PLINK_DIR / "genotypes.map").open("w", encoding="utf-8") as handle:
            for row in kept_info.itertuples(index=False):
                handle.write(f"{_chrom_code(row.chromosome)}\t{row.rsid}\t0\t{int(row.position)}\n")
        (PLINK_DIR / "genotypes.ped").write_text(
            "PED text output is intentionally not generated by pca_qc_fast at large scale; "
            "use genotypes.bed/.bim/.fam.\n",
            encoding="utf-8",
        )


def read_bim_snp_info(prefix: Path) -> pd.DataFrame:
    bim = prefix.with_suffix(".bim")
    if not bim.exists():
        return pd.DataFrame(columns=["rsid", "chromosome", "position"])
    df = pd.read_csv(
        bim,
        sep=r"\s+",
        header=None,
        names=["chromosome", "rsid", "cm", "position", "a1", "a2"],
        dtype={"chromosome": str, "rsid": str},
    )
    return df[["rsid", "chromosome", "position"]].copy()


def write_plink_filtered_snps(input_prefix: Path, qc_prefix: Path) -> None:
    input_info = read_bim_snp_info(input_prefix)
    if input_info.empty:
        write_filtered_snps(None)
        return
    if qc_prefix.with_suffix(".bim").exists():
        kept_rsids = set(read_bim_snp_info(qc_prefix)["rsid"].astype(str))
    else:
        kept_rsids = set()
    fail = ~input_info["rsid"].astype(str).isin(kept_rsids)
    if not fail.any():
        write_filtered_snps(None)
        return
    rows = input_info.loc[fail].copy()
    rows.insert(0, "variant_index", rows.index.astype(np.int64))
    rows["filter"] = "plink_qc"
    write_filtered_snps(rows.reset_index(drop=True))


def write_legacy_matrices(
    dosage: np.memmap,
    snp_info: pd.DataFrame,
    keep_idx: np.ndarray,
    sample_ids: list[str],
) -> None:
    if not WRITE_MATRICES:
        log.info("Skipping legacy genotype_matrix_*.tsv; set BV_WRITE_MATRICES=1 to emit them")
        return
    with timed("Writing legacy numeric matrix TSV"):
        path = MERGED_DIR / "genotype_matrix_numeric.tsv"
        with path.open("w", encoding="utf-8") as handle:
            handle.write("rsid\t" + "\t".join(sample_ids) + "\n")
            kept_rsids = snp_info.iloc[keep_idx]["rsid"].to_numpy()
            for start in range(0, dosage.shape[0], CHUNK_VARIANTS):
                end = min(start + CHUNK_VARIANTS, dosage.shape[0])
                block = np.asarray(dosage[start:end, :])
                for rsid, row in zip(kept_rsids[start:end], block):
                    vals = ["NA" if v == MISSING else str(int(v)) for v in row]
                    handle.write(str(rsid) + "\t" + "\t".join(vals) + "\n")
    # Raw string matrix cannot be reconstructed exactly after allele coding
    # without another pass through every genotype file; leave a clear marker.
    (MERGED_DIR / "genotype_matrix_raw.tsv").write_text(
        "Raw matrix TSV skipped by compact pca_qc_fast backend. "
        "Set BV_WRITE_MATRICES=1 for numeric matrix; use PLINK BED/BIM/FAM for genotype data.\n",
        encoding="utf-8",
    )


def class_counts(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return (
        (block == 0).sum(axis=1),
        (block == 1).sum(axis=1),
        (block == 2).sum(axis=1),
        (block == MISSING).sum(axis=1),
    )


def qc_filter(dosage: np.memmap, snp_info: pd.DataFrame | None = None) -> tuple[np.ndarray, np.ndarray]:
    n_var, n_samples = dosage.shape
    with timed("Applying chunked sample missingness filter"):
        sample_missing = np.zeros(n_samples, dtype=np.int64)
        for start in range(0, n_var, CHUNK_VARIANTS):
            end = min(start + CHUNK_VARIANTS, n_var)
            sample_missing += (np.asarray(dosage[start:end, :]) == MISSING).sum(axis=0)
        sample_keep = (sample_missing / max(n_var, 1)) <= MIND
    log.info("After --mind %.3g: %d/%d samples retained", MIND, int(sample_keep.sum()), n_samples)
    if sample_keep.sum() == 0:
        return np.empty(0, dtype=np.int64), sample_keep

    kept_samples = np.flatnonzero(sample_keep)
    keep_variants: list[np.ndarray] = []
    filtered_rows: list[pd.DataFrame] = []
    with timed("Applying chunked call-rate, MAF, and HWE filters"):
        for start in range(0, n_var, CHUNK_VARIANTS):
            end = min(start + CHUNK_VARIANTS, n_var)
            block = np.asarray(dosage[start:end, :])[:, kept_samples]
            n_homref, n_het, n_homalt, n_missing = class_counts(block)
            called = n_homref + n_het + n_homalt
            call_rate = np.divide(called, block.shape[1], out=np.zeros_like(called, dtype=np.float64), where=block.shape[1] > 0)
            alt_freq = np.divide(
                (2 * n_homalt + n_het),
                2 * called,
                out=np.zeros_like(called, dtype=np.float64),
                where=called > 0,
            )
            maf = np.minimum(alt_freq, 1.0 - alt_freq)

            p_ref = 1.0 - alt_freq
            exp_homref = called * p_ref**2
            exp_het = called * 2 * p_ref * alt_freq
            exp_homalt = called * alt_freq**2
            chi2 = np.zeros_like(alt_freq)
            for obs, exp in ((n_homref, exp_homref), (n_het, exp_het), (n_homalt, exp_homalt)):
                chi2 += np.divide((obs - exp) ** 2, exp, out=np.zeros_like(exp, dtype=np.float64), where=exp > 0)
            hwe_p = stats.chi2.sf(chi2, df=1)
            zero_expected = (exp_homref < 1e-6) | (exp_het < 1e-6) | (exp_homalt < 1e-6)

            pass_call = (1.0 - call_rate) <= GENO
            pass_maf = maf >= MAF
            pass_hwe = (hwe_p >= HWE_P) | zero_expected
            keep = pass_call & pass_maf & pass_hwe
            keep_variants.append(np.arange(start, end, dtype=np.int64)[keep])

            fail = ~keep
            if fail.any():
                reason = np.where(~pass_call, "call_rate", np.where(~pass_maf, "maf", "hwe"))
                idx = np.arange(start, end, dtype=np.int64)[fail]
                denom = np.where(called[fail] > 0, called[fail], 1)
                rows = pd.DataFrame({
                    "variant_index": idx,
                    "filter": reason[fail],
                    "n_homref": n_homref[fail],
                    "n_het": n_het[fail],
                    "n_homalt": n_homalt[fail],
                    "n_missing": n_missing[fail],
                    "frac_homref": np.round(n_homref[fail] / denom, 4),
                    "frac_het": np.round(n_het[fail] / denom, 4),
                    "frac_homalt": np.round(n_homalt[fail] / denom, 4),
                    "call_rate": np.round(call_rate[fail], 6),
                    "maf": np.round(maf[fail], 6),
                    "hwe_p": np.round(hwe_p[fail], 6),
                })
                if snp_info is not None and not snp_info.empty:
                    info = snp_info.iloc[idx].reset_index(drop=True)
                    rows.insert(1, "rsid", info["rsid"].astype(str).to_numpy())
                    rows.insert(2, "chromosome", info["chromosome"].astype(str).to_numpy())
                    rows.insert(3, "position", info["position"].to_numpy())
                filtered_rows.append(rows)

    variant_keep = np.concatenate(keep_variants) if keep_variants else np.empty(0, dtype=np.int64)
    write_filtered_snps(pd.concat(filtered_rows, ignore_index=True) if filtered_rows else None)
    log.info("After QC filters: %d/%d SNPs retained", variant_keep.size, n_var)
    return variant_keep, sample_keep


def standardize_row(row: np.ndarray) -> np.ndarray:
    x = row.astype(np.float32, copy=True)
    called = x != MISSING
    if called.any():
        mean = float(x[called].mean())
        x[~called] = mean
    else:
        x[:] = 0.0
        return x
    x -= float(x.mean())
    std = float(x.std())
    if std > 0:
        x /= std
    return x


def ld_prune_streaming(dosage: np.memmap, qc_idx: np.ndarray, sample_keep: np.ndarray) -> np.ndarray:
    if qc_idx.size == 0:
        return qc_idx
    kept_samples = np.flatnonzero(sample_keep)
    selected: list[int] = []
    recent: deque[tuple[int, np.ndarray]] = deque()
    with timed("LD pruning with streaming window"):
        for pos, variant_idx in enumerate(qc_idx):
            while recent and pos - recent[0][0] > LD_WINDOW:
                recent.popleft()
            x = standardize_row(np.asarray(dosage[variant_idx, kept_samples]))
            drop = False
            for _prev_pos, prev in recent:
                r = float(np.dot(prev, x) / max(len(x), 1))
                if r * r > LD_R2:
                    drop = True
                    break
            if not drop:
                selected.append(int(variant_idx))
                recent.append((pos, x))
    out = np.asarray(selected, dtype=np.int64)
    log.info("After LD pruning: %d SNPs retained from %d", out.size, qc_idx.size)
    return out


def run_chunked_pca(dosage: np.memmap, pruned_idx: np.ndarray, sample_keep: np.ndarray, sample_ids: list[str]) -> None:
    PCA_DIR.mkdir(parents=True, exist_ok=True)
    kept_samples = np.flatnonzero(sample_keep)
    kept_ids = [sample_ids[i] for i in kept_samples]
    if pruned_idx.size == 0:
        write_empty_pca_outputs("No SNPs remained after QC/LD pruning; PCA cannot run.")
        return
    if len(kept_ids) < 2:
        write_empty_pca_outputs(f"Only {len(kept_ids)} usable sample(s) remained after QC; PCA requires at least 2.")
        return

    n_samples = len(kept_ids)
    gram = np.zeros((n_samples, n_samples), dtype=np.float64)
    with timed("Building chunked sample covariance for PCA"):
        for start in range(0, pruned_idx.size, CHUNK_VARIANTS):
            end = min(start + CHUNK_VARIANTS, pruned_idx.size)
            block = np.asarray(dosage[pruned_idx[start:end], :])[:, kept_samples].astype(np.float32)
            missing = block == MISSING
            called = ~missing
            sums = np.where(called, block, 0).sum(axis=1)
            counts = called.sum(axis=1)
            means = np.divide(sums, counts, out=np.zeros_like(sums, dtype=np.float32), where=counts > 0)
            block[missing] = np.take(means, np.where(missing)[0])
            block -= block.mean(axis=1, keepdims=True)
            x = block.T
            gram += x @ x.T

    with timed("Eigen-decomposing sample covariance"):
        eigvals, eigvecs = np.linalg.eigh(gram)
        order = np.argsort(eigvals)[::-1]
        eigvals = np.maximum(eigvals[order], 0)
        eigvecs = eigvecs[:, order]
        n_pcs = min(N_PCS, n_samples - 1, pruned_idx.size)
        eigvals = eigvals[:n_pcs]
        eigvecs = eigvecs[:, :n_pcs]
        scores = eigvecs * np.sqrt(eigvals)[None, :]
        explained_variance = eigvals / max(n_samples - 1, 1)

    with timed("Writing PCA outputs"):
        with (PCA_DIR / "pca.eigenvec").open("w", encoding="utf-8") as handle:
            for sid, row in zip(kept_ids, scores):
                handle.write(f"{sid} {sid} " + " ".join(f"{v:.6f}" for v in row) + "\n")
        with (PCA_DIR / "pca.eigenval").open("w", encoding="utf-8") as handle:
            for val in explained_variance:
                handle.write(f"{val:.6f}\n")

    total = float(explained_variance.sum())
    if total > 0:
        for i, val in enumerate((explained_variance / total * 100)[:5], 1):
            log.info("PC%s: %.2f%% variance explained among emitted PCs", i, val)


def run_plink_backend_if_requested() -> bool:
    if PLINK_BACKEND == "python":
        return False
    if PLINK_BACKEND not in {"auto", "plink", "plink2"}:
        log.warning("Unknown BV_PCA_BACKEND=%r; using chunked Python backend", PLINK_BACKEND)
        return False
    plink2 = shutil.which("plink2")
    if not plink2:
        if PLINK_BACKEND == "auto":
            log.info("plink2 not found; using chunked Python QC/PCA backend")
            return False
        raise RuntimeError("BV_PCA_BACKEND=plink requested, but plink2 was not found on PATH")
    prefix = PLINK_DIR / "genotypes"
    qc = PLINK_DIR / "qc_pass"
    pruned = PLINK_DIR / "pruned"
    out = PCA_DIR / "pca"
    PCA_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    commands = [
        [plink2, "--bfile", str(prefix), "--geno", str(GENO), "--mind", str(MIND), "--maf", str(MAF), "--hwe", str(HWE_P), "--make-bed", "--out", str(qc), "--allow-no-sex"],
        [plink2, "--bfile", str(qc), "--indep-pairwise", str(LD_WINDOW), str(LD_STEP), str(LD_R2), "--out", str(pruned), "--allow-no-sex"],
        [plink2, "--bfile", str(qc), "--extract", str(pruned) + ".prune.in", "--pca", str(N_PCS), "--out", str(out), "--allow-no-sex"],
    ]
    labels = ["plink_qc", "plink_prune", "plink_pca"]
    for label, cmd in zip(labels, commands):
        with timed(f"Running {label}"):
            with (LOG_DIR / f"{label}.log").open("w", encoding="utf-8") as log_handle:
                try:
                    subprocess.run(cmd, check=True, stdout=log_handle, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError:
                    if PLINK_BACKEND == "auto":
                        log.warning("%s failed under auto backend; falling back to chunked Python QC/PCA", label)
                        return False
                    raise
    write_plink_filtered_snps(prefix, qc)
    return True


def load_eigenvec() -> pd.DataFrame:
    df = pd.read_csv(PCA_DIR / "pca.eigenvec", sep=r"\s+", comment="#", header=None)
    if df.empty:
        return df
    n_pcs = df.shape[1] - 2
    df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]
    df["sample_id"] = df["IID"].astype(str)
    return df


def load_eigenval() -> list[float]:
    vals = pd.read_csv(PCA_DIR / "pca.eigenval", header=None)[0].tolist()
    total = sum(vals)
    return [v / total * 100 for v in vals] if total > 0 else vals


def scatter_pca(df: pd.DataFrame, pc_x: str, pc_y: str, var_exp: list[float], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    palette = plt.get_cmap("tab10", max(len(df), 1))
    colors = [palette(i) for i in range(len(df))]
    ax.scatter(df[pc_x], df[pc_y], c=colors, s=45, edgecolors="k", linewidths=0.3, alpha=0.9)
    if len(df) <= 50:
        for _i, row in df.iterrows():
            ax.annotate(row["sample_id"], (row[pc_x], row[pc_y]), textcoords="offset points", xytext=(6, 4), fontsize=7)
    x_num = int(pc_x.replace("PC", ""))
    y_num = int(pc_y.replace("PC", ""))
    x_label = f"{pc_x} ({var_exp[x_num - 1]:.1f}% var)" if len(var_exp) >= x_num else pc_x
    y_label = f"{pc_y} ({var_exp[y_num - 1]:.1f}% var)" if len(var_exp) >= y_num else pc_y
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title("Ancestry PCA", fontsize=13, fontweight="bold")
    ax.axhline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved plot -> %s", os.path.basename(str(out_path)))


def plot_pca() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    eigenvec = PCA_DIR / "pca.eigenvec"
    if not eigenvec.exists() or eigenvec.stat().st_size == 0:
        log.warning("PCA plot skipped: no PCA coordinates were produced")
        return
    with timed("Plotting PCA"):
        df = load_eigenvec()
        if df.empty or "PC2" not in df.columns:
            log.warning("PCA plot skipped: fewer than two PCs were produced")
            return
        var_exp = load_eigenval() if (PCA_DIR / "pca.eigenval").exists() else []
        scatter_pca(df, "PC1", "PC2", var_exp, PLOTS_DIR / "pca_pc1_pc2.png")
        if "PC4" in df.columns:
            scatter_pca(df, "PC3", "PC4", var_exp, PLOTS_DIR / "pca_pc3_pc4.png")


def read_bed_to_dosage(prefix: Path) -> tuple[np.memmap, list[str]]:
    """Read a PLINK .bed/.bim/.fam into the (n_var x n_samples) int8 dosage memmap
    the Python QC/PCA backend expects. 2-bit codes: 00->2, 10->1, 11->0, 01->missing
    (matches build_dosage_memmap's encoding, dosage = count of A1/minor)."""
    sample_ids = [
        line.split()[1]
        for line in prefix.with_suffix(".fam").read_text().splitlines()
        if line.strip()
    ]
    n_samples = len(sample_ids)
    n_var = sum(1 for line in prefix.with_suffix(".bim").read_text().splitlines() if line.strip())
    raw = np.fromfile(prefix.with_suffix(".bed"), dtype=np.uint8)
    if raw.size < 3 or raw[0] != 0x6C or raw[1] != 0x1B or raw[2] != 0x01:
        raise ValueError(f"{prefix}.bed is not a valid variant-major PLINK bed")
    body = raw[3:]
    bytes_per_var = (n_samples + 3) // 4
    body = body.reshape(n_var, bytes_per_var)
    # unpack 2-bit codes -> (n_var, bytes_per_var*4)
    two_bit = np.zeros((n_var, bytes_per_var * 4), dtype=np.uint8)
    for k in range(4):
        two_bit[:, k::4] = (body >> (2 * k)) & 0b11
    two_bit = two_bit[:, :n_samples]
    dosage = np.memmap(WORK_DIR / "dosage.int8.mmap", dtype=np.int8, mode="w+", shape=(n_var, n_samples))
    # 00->2, 10->1, 11->0, 01->missing(-1)
    dosage[:] = np.select(
        [two_bit == 0b00, two_bit == 0b10, two_bit == 0b11],
        [np.int8(2), np.int8(1), np.int8(0)],
        default=MISSING,
    ).astype(np.int8)
    dosage.flush()
    return dosage, sample_ids


def run_from_prebuilt_bed(prefix: Path) -> None:
    """Skip parse/universe/memmap/bed-build: use a prebuilt PLINK bed (e.g. from
    `bvs cohort-bed`) and run the exact same QC/PCA backend + plots. Keeps outputs
    byte-identical to the full pipeline because they run on the identical bed."""
    log.info("Using prebuilt PLINK bed: %s.{bed,bim,fam}", prefix)
    PLINK_DIR.mkdir(parents=True, exist_ok=True)
    dest = PLINK_DIR / "genotypes"
    for ext in ("bed", "bim", "fam"):
        shutil.copy(str(prefix.with_suffix(f".{ext}")), str(dest.with_suffix(f".{ext}")))
    if not run_plink_backend_if_requested():
        dosage, sample_ids = read_bed_to_dosage(dest)
        snp_info = read_bim_snp_info(dest)
        qc_idx, sample_keep = qc_filter(dosage, snp_info)
        if qc_idx.size == 0:
            write_empty_pca_outputs("No SNPs remained after call-rate/MAF/HWE filtering; PCA cannot run.")
            return
        pruned_idx = ld_prune_streaming(dosage, qc_idx, sample_keep)
        run_chunked_pca(dosage, pruned_idx, sample_keep, sample_ids)
    plot_pca()


def main() -> None:
    total = time.perf_counter()
    for path in (MERGED_DIR, PLINK_DIR, PCA_DIR, PLOTS_DIR, LOG_DIR, QC_DIR, WORK_DIR):
        path.mkdir(parents=True, exist_ok=True)
    ensure_filtered_snps_file()

    prebuilt = os.environ.get("BV_PREBUILT_BED", "").strip()
    if prebuilt:
        run_from_prebuilt_bed(Path(prebuilt))
        log.info("Fast pipeline (prebuilt bed) complete in %.2fs", time.perf_counter() - total)
        return

    samples = discover_samples()
    log.info("Discovered %d samples", len(samples))

    snp_ids, snp_info, usable_samples, _errors = build_snp_universe(samples)
    a1_mm, a2_mm, sample_ids = fill_allele_memmaps(usable_samples, snp_ids)
    counts = allele_counts(a1_mm, a2_mm)
    keep_idx, major, minor = choose_biallelic_variants(counts)
    if keep_idx.size == 0:
        write_empty_pca_outputs("No mono/biallelic ACGT SNPs were available for PCA.")
        return

    dosage = build_dosage_memmap(a1_mm, a2_mm, keep_idx, major, minor)
    write_plink_binary(dosage, snp_info, keep_idx, major, minor, sample_ids)
    write_legacy_matrices(dosage, snp_info, keep_idx, sample_ids)

    if not run_plink_backend_if_requested():
        dosage_snp_info = snp_info.iloc[keep_idx].reset_index(drop=True)
        qc_idx, sample_keep = qc_filter(dosage, dosage_snp_info)
        if qc_idx.size == 0:
            write_empty_pca_outputs("No SNPs remained after call-rate/MAF/HWE filtering; PCA cannot run.")
            return
        pruned_idx = ld_prune_streaming(dosage, qc_idx, sample_keep)
        run_chunked_pca(dosage, pruned_idx, sample_keep, sample_ids)
    plot_pca()

    log.info("Fast pipeline complete in %.2fs", time.perf_counter() - total)


if __name__ == "__main__":
    main()

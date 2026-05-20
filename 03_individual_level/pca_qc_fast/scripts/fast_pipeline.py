#!/usr/bin/env python3
"""
Fast pca_qc implementation.

Outputs intentionally mirror ../pca_qc:
  data/merged/genotype_matrix_raw.tsv
  data/merged/genotype_matrix_numeric.tsv
  data/merged/snp_info.tsv
  data/plink/genotypes.ped
  data/plink/genotypes.map
  data/pca/pca.eigenvec
  data/pca/pca.eigenval
  plots/pca_pc1_pc2.png
  plots/pca_pc3_pc4.png
"""

from __future__ import annotations

import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import genoio as _genoio  # noqa: E402  (synced fork, see genoio.py header)


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("BIOVAULT_DATA_DIR", BASE_DIR.parent)).resolve()
MERGED_DIR = BASE_DIR / "data" / "merged"
PLINK_DIR = BASE_DIR / "data" / "plink"
PCA_DIR = BASE_DIR / "data" / "pca"
PLOTS_DIR = BASE_DIR / "plots"
LOG_DIR = BASE_DIR / "logs"
QC_DIR = BASE_DIR / "data" / "qc"

COLS = ["rsid", "chromosome", "position", "genotype", "gs", "baf", "lrr"]
VALID_BASES = np.array([b"A", b"C", b"G", b"T"], dtype="S1")
BASE_LABELS = np.array(["A", "C", "G", "T"], dtype=object)

N_PCS = 20
GENO = 0.05
MIND = 0.10
MAF = 0.01
HWE_P = 1e-4
LD_WINDOW = 50
LD_R2 = 0.2


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


def timed(label: str):
    class Timer:
        def __enter__(self):
            self.start = time.perf_counter()
            log.info("%s ...", label)
            return self

        def __exit__(self, exc_type, exc, tb):
            elapsed = time.perf_counter() - self.start
            if exc_type is None:
                log.info("%s done in %.2fs", label, elapsed)

    return Timer()


def discover_samples() -> list[tuple[str, Path]]:
    samples: list[tuple[str, Path]] = []
    for sample_dir in sorted(DATA_DIR.iterdir()):
        if not sample_dir.is_dir() or not sample_dir.name.isdigit():
            continue
        txt_files = sorted(sample_dir.glob("*.txt"))
        if txt_files:
            samples.append((sample_dir.name, txt_files[0]))
    if not samples:
        raise RuntimeError(f"No numeric sample directories with .txt files found under {DATA_DIR}")
    return samples


def read_sample(task: tuple[str, Path]) -> tuple[str, pd.DataFrame]:
    sample_id, path = task
    if _genoio.sniff_format(path) == "illumina":
        g = _genoio.read_genotypes(path)
        df = g[["rsid", "chrom", "pos", "gt"]].rename(columns={
            "chrom": "chromosome", "pos": "position", "gt": "genotype"})
        df["position"] = df["position"].astype(np.int64)
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            names=COLS,
            usecols=["rsid", "chromosome", "position", "genotype"],
            dtype={"rsid": str, "chromosome": str, "position": np.int64, "genotype": str},
        )
    df = df.dropna(subset=["rsid", "genotype"])
    df["genotype"] = df["genotype"].str.upper().str.strip()
    return sample_id, df


def load_and_merge(samples: list[tuple[str, Path]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    workers = min(len(samples), os.cpu_count() or 1)
    with timed(f"Reading {len(samples)} samples with {workers} workers"):
        if workers > 1:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                loaded = list(pool.map(read_sample, samples))
        else:
            loaded = [read_sample(s) for s in samples]

    loaded.sort(key=lambda item: item[0])
    snp_info = (
        loaded[0][1][["rsid", "chromosome", "position"]]
        .drop_duplicates("rsid")
        .reset_index(drop=True)
    )

    with timed("Merging genotype matrix"):
        series = [
            df.drop_duplicates("rsid").set_index("rsid")["genotype"].rename(sample_id)
            for sample_id, df in loaded
        ]
        matrix = pd.concat(series, axis=1, join="outer")
        matrix.index.name = "rsid"

    MERGED_DIR.mkdir(parents=True, exist_ok=True)
    with timed("Writing merged outputs"):
        matrix.to_csv(MERGED_DIR / "genotype_matrix_raw.tsv", sep="\t")
        snp_info.to_csv(MERGED_DIR / "snp_info.tsv", sep="\t", index=False)

    log.info("Loaded matrix: %s SNPs x %s samples", matrix.shape[0], matrix.shape[1])
    return matrix, snp_info


def genotype_chars(matrix: pd.DataFrame) -> np.ndarray:
    values = matrix.fillna("--").to_numpy(dtype="S2", copy=True)
    values = np.ascontiguousarray(values)
    return values.view("S1").reshape(values.shape[0], values.shape[1], 2)


def infer_major_minor(chars: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with timed("Inferring major/minor alleles"):
        counts = np.empty((chars.shape[0], 4), dtype=np.int16)
        first_pos = np.empty((chars.shape[0], 4), dtype=np.int16)
        flat = chars.reshape(chars.shape[0], -1)
        positions = np.arange(flat.shape[1], dtype=np.int16)
        for i, base in enumerate(VALID_BASES):
            matches = flat == base
            counts[:, i] = matches.sum(axis=1)
            first = np.where(matches, positions, np.iinfo(np.int16).max).min(axis=1)
            first_pos[:, i] = first

        # Match Counter.most_common() tie behavior from the original pipeline:
        # highest count first, then first allele observed while scanning samples.
        order = np.lexsort((first_pos, -counts), axis=1)
        major_idx = order[:, 0]
        minor_idx = order[:, 1]
        no_minor = counts[np.arange(counts.shape[0]), minor_idx] == 0
        minor_idx[no_minor] = major_idx[no_minor]
        has_observed = counts.sum(axis=1) > 0

    return major_idx, minor_idx, has_observed


def encode_numeric(
    matrix: pd.DataFrame,
    snp_info: pd.DataFrame,
    chars: np.ndarray,
    major_idx: np.ndarray,
    minor_idx: np.ndarray,
    has_observed: np.ndarray,
) -> pd.DataFrame:
    n_snps, n_samples, _ = chars.shape
    major = VALID_BASES[major_idx]
    minor = VALID_BASES[minor_idx]

    with timed("Encoding numeric dosage"):
        valid = (
            ((chars[:, :, 0] == major[:, None]) | (chars[:, :, 0] == minor[:, None]))
            & ((chars[:, :, 1] == major[:, None]) | (chars[:, :, 1] == minor[:, None]))
            & has_observed[:, None]
        )
        dosage = (chars == minor[:, None, None]).sum(axis=2).astype(np.float32)
        dosage[~valid] = np.nan
        numeric = pd.DataFrame(dosage, index=matrix.index, columns=matrix.columns)
        numeric.index.name = "rsid"

    with timed("Writing numeric matrix"):
        numeric.to_csv(MERGED_DIR / "genotype_matrix_numeric.tsv", sep="\t")

    write_plink_files(matrix, snp_info, chars, major, minor, has_observed)
    return numeric


def write_plink_files(
    matrix: pd.DataFrame,
    snp_info: pd.DataFrame,
    chars: np.ndarray,
    major: np.ndarray,
    minor: np.ndarray,
    has_observed: np.ndarray,
) -> None:
    PLINK_DIR.mkdir(parents=True, exist_ok=True)
    valid = (
        ((chars[:, :, 0] == major[:, None]) | (chars[:, :, 0] == minor[:, None]))
        & ((chars[:, :, 1] == major[:, None]) | (chars[:, :, 1] == minor[:, None]))
        & has_observed[:, None]
    )

    with timed("Writing PLINK PED"):
        ped_path = PLINK_DIR / "genotypes.ped"
        with ped_path.open("w") as f:
            for sample_idx, sid in enumerate(matrix.columns):
                sample_chars = chars[:, sample_idx, :].astype("U1")
                sample_valid = valid[:, sample_idx]
                a1 = np.where(sample_valid, sample_chars[:, 0], "0")
                a2 = np.where(sample_valid, sample_chars[:, 1], "0")
                alleles = np.empty(a1.size * 2, dtype=object)
                alleles[0::2] = a1
                alleles[1::2] = a2
                f.write(f"{sid} {sid} 0 0 0 -9 {' '.join(alleles)}\n")

    with timed("Writing PLINK MAP"):
        info = snp_info.drop_duplicates("rsid").set_index("rsid").reindex(matrix.index)
        chrom = info["chromosome"].fillna("0").astype(str)
        chrom = (
            chrom.str.replace("XY", "25", regex=False)
            .str.replace("X", "23", regex=False)
            .str.replace("Y", "24", regex=False)
            .str.replace("MT", "26", regex=False)
        )
        pos = info["position"].fillna(0).astype("int64")
        out = pd.DataFrame(
            {
                "chromosome": chrom.to_numpy(),
                "rsid": matrix.index.to_numpy(),
                "genetic_distance": 0,
                "position": pos.to_numpy(),
            }
        )
        out.to_csv(PLINK_DIR / "genotypes.map", sep="\t", header=False, index=False)


def _class_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Per-SNP genotype-class counts from the numeric matrix
    (0 = hom REF/major, 1 = het, 2 = hom ALT/minor, NaN = missing)."""
    v = df.to_numpy(dtype=np.float64, copy=False)
    n_homref = np.nansum(v == 0, axis=1).astype(int)
    n_het = np.nansum(v == 1, axis=1).astype(int)
    n_homalt = np.nansum(v == 2, axis=1).astype(int)
    n_miss = np.isnan(v).sum(axis=1).astype(int)
    called = (n_homref + n_het + n_homalt)
    denom = np.where(called > 0, called, 1)
    return pd.DataFrame({
        "n_homref": n_homref, "n_het": n_het, "n_homalt": n_homalt,
        "n_missing": n_miss,
        "frac_homref": np.round(n_homref / denom, 4),
        "frac_het": np.round(n_het / denom, 4),
        "frac_homalt": np.round(n_homalt / denom, 4),
    }, index=df.index)


def _write_filtered(records: pd.DataFrame, mode: str) -> None:
    """Append/initialise QC_DIR/filtered_snps.tsv — one row per dropped SNP
    with the filter that removed it and its genotype-class breakdown.
    Aggregate per-variant stats only (no per-individual data)."""
    QC_DIR.mkdir(parents=True, exist_ok=True)
    out = QC_DIR / "filtered_snps.tsv"
    header = mode == "w"
    records.to_csv(out, sep="\t", mode=mode, header=header, index=False)


def filter_qc(numeric: pd.DataFrame) -> pd.DataFrame:
    cc = _class_counts(numeric)  # per-SNP class counts on the full matrix
    dropped: list[pd.DataFrame] = []

    def _record(idx, reason, **extra):
        if len(idx) == 0:
            return
        r = cc.loc[idx].reset_index()  # 'rsid' + class counts
        r.insert(1, "filter", reason)
        for k, val in extra.items():
            r[k] = np.round(np.asarray(val), 6)
        dropped.append(r)

    with timed("Applying call-rate and MAF filters"):
        mat = numeric
        snp_missing = mat.isna().mean(axis=1)
        ind_missing = mat.isna().mean(axis=0)

        geno_fail = snp_missing.index[snp_missing > GENO]
        _record(geno_fail, "call_rate",
                 call_rate=(1.0 - snp_missing.loc[geno_fail]).to_numpy())
        n_ind_dropped = int((ind_missing > MIND).sum())

        mat = mat.loc[snp_missing <= GENO, ind_missing <= MIND]

        values = mat.to_numpy(dtype=np.float64, copy=False)
        alt_freq = np.nanmean(values, axis=1) / 2.0
        maf_vals = pd.Series(np.minimum(alt_freq, 1.0 - alt_freq),
                             index=mat.index)
        maf_fail = maf_vals.index[maf_vals < MAF]
        _record(maf_fail, "maf", maf=maf_vals.loc[maf_fail].to_numpy())
        mat = mat.loc[maf_vals >= MAF]
        log.info("After call-rate/MAF filters: %s SNPs x %s samples "
                 "(call_rate dropped %d, maf dropped %d, %d samples dropped "
                 "by --mind)", mat.shape[0], mat.shape[1],
                 len(geno_fail), len(maf_fail), n_ind_dropped)

    with timed("Applying vectorized HWE filter"):
        values = mat.to_numpy(dtype=np.float64, copy=False)
        n_hom_ref = np.nansum(values == 0, axis=1)
        n_het = np.nansum(values == 1, axis=1)
        n_hom_alt = np.nansum(values == 2, axis=1)
        n = (n_hom_ref + n_het + n_hom_alt).astype(np.float64)
        n_hom_ref = n_hom_ref.astype(np.float64)
        n_het = n_het.astype(np.float64)
        n_hom_alt = n_hom_alt.astype(np.float64)

        p_alt = np.divide(2 * n_hom_alt + n_het, 2 * n, out=np.zeros_like(n), where=n > 0)
        p_ref = 1.0 - p_alt
        exp_hom_ref = n * p_ref**2
        exp_het = n * 2 * p_ref * p_alt
        exp_hom_alt = n * p_alt**2

        chi2 = np.zeros_like(n, dtype=np.float64)
        for obs, exp in (
            (n_hom_ref, exp_hom_ref),
            (n_het, exp_het),
            (n_hom_alt, exp_hom_alt),
        ):
            chi2 += np.divide((obs - exp) ** 2, exp, out=np.zeros_like(exp, dtype=np.float64), where=exp > 0)

        p_values = stats.chi2.sf(chi2, df=1)
        zero_expected = (exp_hom_ref < 1e-6) | (exp_het < 1e-6) | (exp_hom_alt < 1e-6)
        keep = (p_values >= HWE_P) | zero_expected
        hwe_fail = mat.index[~keep]
        _record(hwe_fail, "hwe",
                hwe_p=p_values[~keep])
        mat = mat.loc[keep]
        log.info("After HWE filter: %s SNPs x %s samples (hwe dropped %d)",
                 mat.shape[0], mat.shape[1], len(hwe_fail))

    if dropped:
        allcols = ["rsid", "filter", "n_homref", "n_het", "n_homalt",
                   "n_missing", "frac_homref", "frac_het", "frac_homalt",
                   "call_rate", "maf", "hwe_p"]
        out = pd.concat(dropped, ignore_index=True)
        for c in allcols:
            if c not in out.columns:
                out[c] = ""
        _write_filtered(out[allcols], "w")
        summary = out["filter"].value_counts().to_dict()
        log.info("QC drops by filter: %s -> %s", summary,
                 os.path.relpath(QC_DIR / "filtered_snps.tsv", BASE_DIR))

    return mat


def ld_prune_fast(mat: pd.DataFrame) -> pd.DataFrame:
    with timed("LD pruning"):
        imputer = SimpleImputer(strategy="mean")
        x = imputer.fit_transform(mat.to_numpy(dtype=np.float64, copy=False).T).T
        x -= x.mean(axis=1, keepdims=True)
        std = x.std(axis=1, keepdims=True)
        std[std == 0] = 1.0
        x /= std

        n_snps, n_samples = x.shape
        selected = np.ones(n_snps, dtype=bool)
        for i in range(n_snps):
            if not selected[i]:
                continue
            end = min(i + LD_WINDOW + 1, n_snps)
            candidates = np.arange(i + 1, end)
            candidates = candidates[selected[candidates]]
            if candidates.size == 0:
                continue
            r = (x[candidates] @ x[i]) / n_samples
            selected[candidates[r * r > LD_R2]] = False

        pruned = mat.iloc[selected]
        ld_dropped = mat.index[~selected]
        if len(ld_dropped) > 0:
            r = _class_counts(mat.loc[ld_dropped]).reset_index()
            r.insert(1, "filter", "ld_prune")
            for c in ("call_rate", "maf", "hwe_p"):
                r[c] = ""
            _write_filtered(r[["rsid", "filter", "n_homref", "n_het",
                               "n_homalt", "n_missing", "frac_homref",
                               "frac_het", "frac_homalt", "call_rate",
                               "maf", "hwe_p"]], "a")
        log.info("After LD pruning: %s SNPs retained from %s "
                 "(ld_prune dropped %d, appended to filtered_snps.tsv)",
                 pruned.shape[0], n_snps, len(ld_dropped))
        return pruned


def run_pca(mat: pd.DataFrame) -> None:
    PCA_DIR.mkdir(parents=True, exist_ok=True)
    with timed("Imputing and running PCA"):
        imputer = SimpleImputer(strategy="mean")
        x = imputer.fit_transform(mat.to_numpy(dtype=np.float64, copy=False).T)
        n_pcs = min(N_PCS, x.shape[0] - 1, x.shape[1])
        pca = PCA(n_components=n_pcs)
        scores = pca.fit_transform(x)

    with timed("Writing PCA outputs"):
        with (PCA_DIR / "pca.eigenvec").open("w") as f:
            for i, sid in enumerate(mat.columns):
                pc_vals = " ".join(f"{v:.6f}" for v in scores[i])
                f.write(f"{sid} {sid} {pc_vals}\n")

        with (PCA_DIR / "pca.eigenval").open("w") as f:
            for val in pca.explained_variance_:
                f.write(f"{val:.6f}\n")

    var_exp = pca.explained_variance_ratio_ * 100
    for i, val in enumerate(var_exp[:5], 1):
        log.info("PC%s: %.2f%% variance explained", i, val)


def load_eigenvec() -> pd.DataFrame:
    df = pd.read_csv(PCA_DIR / "pca.eigenvec", sep=r"\s+", header=None)
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
    palette = plt.get_cmap("tab10", len(df))
    colors = [palette(i) for i in range(len(df))]

    ax.scatter(df[pc_x], df[pc_y], c=colors, s=80, edgecolors="k", linewidths=0.5, alpha=0.9)
    for i, row in df.iterrows():
        ax.annotate(row["sample_id"], (row[pc_x], row[pc_y]), textcoords="offset points", xytext=(6, 4), fontsize=7)

    x_num = int(pc_x.replace("PC", ""))
    y_num = int(pc_y.replace("PC", ""))
    x_label = f"{pc_x} ({var_exp[x_num - 1]:.1f}% var)" if var_exp else pc_x
    y_label = f"{pc_y} ({var_exp[y_num - 1]:.1f}% var)" if var_exp else pc_y

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[i], markersize=9, label=df["sample_id"].iloc[i])
        for i in range(len(df))
    ]
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title("Ancestry PCA - Mock Synthetic Data", fontsize=13, fontweight="bold")
    ax.legend(handles=handles, title="Sample", fontsize=8, loc="best", framealpha=0.7)
    ax.axhline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved plot -> %s", os.path.basename(str(out_path)))


def plot_pca() -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    with timed("Plotting PCA"):
        df = load_eigenvec()
        var_exp = load_eigenval() if (PCA_DIR / "pca.eigenval").exists() else []
        scatter_pca(df, "PC1", "PC2", var_exp, PLOTS_DIR / "pca_pc1_pc2.png")
        pc_cols = [c for c in df.columns if c.startswith("PC")]
        if len(pc_cols) >= 4:
            scatter_pca(df, "PC3", "PC4", var_exp, PLOTS_DIR / "pca_pc3_pc4.png")


def main() -> None:
    total = time.perf_counter()
    for path in (MERGED_DIR, PLINK_DIR, PCA_DIR, PLOTS_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)

    samples = discover_samples()
    log.info("Discovered %d samples", len(samples))

    matrix, snp_info = load_and_merge(samples)
    chars = genotype_chars(matrix)
    major_idx, minor_idx, has_observed = infer_major_minor(chars)
    numeric = encode_numeric(matrix, snp_info, chars, major_idx, minor_idx, has_observed)
    qc = filter_qc(numeric)
    pruned = ld_prune_fast(qc)
    run_pca(pruned)
    plot_pca()

    log.info("Fast pipeline complete in %.2fs", time.perf_counter() - total)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Project per-participant DDNA genotypes onto the gnomAD HGDP+1kGP PCA space
using a pre-exported loadings matrix. No PLINK, no Hail, no Spark — just
pandas + numpy.

Math (same as hl.experimental.pc_project):
    PC[i, k] = Σ_j (G[i, j] − 2·pca_af[j]) · loadings[j, k]
with missing G values contributing 0 (i.e. centered to the population mean).

Match key: (chrom, pos). Variants whose observed alleles don't match either
ref or alt at the loadings locus are treated as missing (no strand-flip).

Outputs (in <out_dir>):
    study_pca_projection.tsv     s\\tscores       (scores = "[v1,v2,...]")
    pca_projection.png            PC1 vs PC2 scatter
    qc_report.txt                 small summary

Usage:
    fast_pipeline.py <data_dir> <out_dir>
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


AUTOSOMES = [str(c) for c in range(1, 23)]


def setup_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")


def load_loadings(path):
    arc = np.load(path)
    chrom = arc["chrom"].astype(str)
    pos = arc["pos"].astype(np.int64)
    ref = arc["ref"].astype(str)
    alt = arc["alt"].astype(str)
    loadings = arc["loadings"].astype(np.float32)
    pca_af = arc["pca_af"].astype(np.float32)
    key = np.char.add(np.char.add(chrom, ":"), pos.astype(str))
    key_to_idx = dict(zip(key.tolist(), np.arange(len(key), dtype=np.int64).tolist()))
    return ref, alt, loadings, pca_af, key_to_idx


def compute_dosage(txt_path, key_to_idx, ref_arr, alt_arr, n_var, min_gs):
    df = pd.read_csv(
        txt_path, sep="\t", header=None, comment="#",
        names=["rsid", "chrom", "pos", "gt", "gs", "baf", "lrr"],
        usecols=["chrom", "pos", "gt", "gs"],
        dtype={"chrom": str, "pos": np.int64, "gt": str, "gs": float},
        engine="c",
    )
    df = df[df["rsid"] != "rsid"] if "rsid" in df.columns else df
    df = df[df["chrom"].isin(AUTOSOMES)]
    df = df[df["gs"] >= min_gs]
    df = df[df["gt"].str.len() == 2]

    keys = df["chrom"].to_numpy(dtype=object) + ":" + df["pos"].astype(str).to_numpy(dtype=object)
    idx_series = pd.Series(keys).map(key_to_idx)
    valid_key = idx_series.notna()
    idx = idx_series[valid_key].to_numpy(dtype=np.int64)

    if idx.size == 0:
        return np.full(n_var, np.nan, dtype=np.float32)

    gt = df.loc[valid_key.values, "gt"].to_numpy(dtype=object)
    a1 = np.array([g[0] for g in gt])
    a2 = np.array([g[1] for g in gt])

    ref = ref_arr[idx]
    alt = alt_arr[idx]

    m_a1_alt = (a1 == alt)
    m_a2_alt = (a2 == alt)
    m_a1_ref = (a1 == ref)
    m_a2_ref = (a2 == ref)
    valid_alleles = (m_a1_alt | m_a1_ref) & (m_a2_alt | m_a2_ref)

    dosage = np.full(n_var, np.nan, dtype=np.float32)
    if valid_alleles.any():
        dosage[idx[valid_alleles]] = (
            m_a1_alt[valid_alleles].astype(np.float32)
            + m_a2_alt[valid_alleles].astype(np.float32)
        )
    return dosage


def process_sample(args):
    pid, txt_path, key_to_idx, ref_arr, alt_arr, n_var, min_gs = args
    dosage = compute_dosage(txt_path, key_to_idx, ref_arr, alt_arr, n_var, min_gs)
    n_obs = int(np.isfinite(dosage).sum())
    return pid, dosage, n_obs


def main():
    setup_logging()
    log = logging.getLogger(__name__)

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--loadings", type=Path,
                    default=Path(os.environ.get(
                        "LOADINGS_NPZ",
                        "/opt/biovault/reference/pca_loadings/loadings.npz")))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--min-gs", type=float, default=0.15)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Loading {args.loadings}")
    ref_arr, alt_arr, loadings, pca_af, key_to_idx = load_loadings(args.loadings)
    n_var, n_pcs = loadings.shape
    log.info(f"Loadings: {n_var:,} variants × {n_pcs} PCs")

    sample_dirs = sorted(
        d for d in args.data_dir.iterdir()
        if d.is_dir() and d.name.isdigit()
    )
    log.info(f"Samples: {len(sample_dirs)}")

    work_items = []
    for d in sample_dirs:
        txts = list(d.glob("*.txt"))
        if not txts:
            log.warning(f"  skip {d.name}: no .txt")
            continue
        work_items.append((
            d.name, txts[0], key_to_idx, ref_arr, alt_arr, n_var, args.min_gs,
        ))

    t0 = time.time()
    results = {}
    obs_counts = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(process_sample, it) for it in work_items]
        for f in as_completed(futures):
            pid, dosage, n_obs = f.result()
            results[pid] = dosage
            obs_counts[pid] = n_obs
    ingest_s = time.time() - t0
    log.info(f"Ingested {len(results)} samples in {ingest_s:.1f}s")

    pids = sorted(results)
    G = np.stack([results[p] for p in pids])  # (n_samples, n_var)
    log.info(f"Dosage matrix: {G.shape}")

    # Normalize and project. Matches hl.experimental.pc_project:
    #   gt_norm = (G - 2*af) / sqrt(n_variants * 2 * af * (1 - af))
    #   score   = Σ_j loadings[j,k] * gt_norm[i,j]
    # Missing G stays 0 (no contribution).
    t0 = time.time()
    af = pca_af.astype(np.float64)
    af_valid = (af > 0.0) & (af < 1.0)
    n_var_total = int(n_var)  # total loadings rows, like Hail's loadings_ht.count()
    sd = np.sqrt(n_var_total * 2.0 * af * (1.0 - af))
    sd_safe = np.where(sd > 0, sd, 1.0)
    loadings_scaled = loadings.astype(np.float64) / sd_safe[:, None]
    loadings_scaled[~af_valid] = 0.0

    deviation = G.astype(np.float64) - 2.0 * af[None, :]
    np.copyto(deviation, 0.0, where=~np.isfinite(deviation))
    deviation[:, ~af_valid] = 0.0
    scores = (deviation @ loadings_scaled).astype(np.float64)
    proj_s = time.time() - t0
    log.info(f"Projection in {proj_s:.3f}s")

    out_tsv = args.out_dir / "study_pca_projection.tsv"
    with out_tsv.open("w") as f:
        f.write("s\tscores\n")
        for pid, row in zip(pids, scores):
            f.write(f"{pid}\t[{','.join(repr(float(v)) for v in row)}]\n")
    log.info(f"Wrote {out_tsv}")

    qc = args.out_dir / "qc_report.txt"
    qc.write_text(
        "=== gnomad_projection_fast QC ===\n"
        f"Loadings variants: {n_var:,}\n"
        f"Samples projected: {len(pids)}\n"
        f"Mean observed loadings/sample: {np.mean(list(obs_counts.values())):.0f}\n"
        f"Ingest time: {ingest_s:.1f}s\n"
        f"Projection time: {proj_s:.3f}s\n"
    )
    log.info(f"Wrote {qc}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(scores[:, 0], scores[:, 1], c="black", marker="*", s=200,
                   edgecolors="gold", linewidth=0.8)
        for i, pid in enumerate(pids):
            ax.annotate(str(pid), (scores[i, 0], scores[i, 1]),
                        fontsize=8, alpha=0.7, xytext=(4, 4),
                        textcoords="offset points")
        ax.set_xlabel("PC1 (gnomAD HGDP+1kGP space)")
        ax.set_ylabel("PC2 (gnomAD HGDP+1kGP space)")
        ax.grid(True, alpha=0.2)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)
        ax.set_title("Study samples projected onto gnomAD reference PCs (fast)")
        plot_path = args.out_dir / "pca_projection.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200, bbox_inches="tight")
        plt.close()
        log.info(f"Plot -> {plot_path}")
    except Exception as e:
        log.warning(f"plot skipped: {e}")


if __name__ == "__main__":
    main()

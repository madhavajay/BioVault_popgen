#!/usr/bin/env python3
"""Find ADMIXTURE K and label components from self-reported ancestry anchors.

This stage is intentionally autosome-only. It prepares a QC/LD-pruned PLINK
dataset, runs ADMIXTURE across K, and uses participants who reported exactly
one ancestry label as anchors to map anonymous components to labels such as
AFR/EUR/SAS. A later stage can reuse the selected K and component map for
autosome-vs-X sex-biased admixture.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_LABEL_ALIASES = {
    "african": "AFR",
    "africa": "AFR",
    "afr": "AFR",
    "black": "AFR",
    "european": "EUR",
    "europe": "EUR",
    "eur": "EUR",
    "white": "EUR",
    "caucasian": "EUR",
    "indian": "SAS",
    "south asian": "SAS",
    "south-asian": "SAS",
    "sas": "SAS",
    "asian indian": "SAS",
}


def run(cmd: list[str], log_path: Path | None = None, cwd: Path | None = None) -> str:
    text = " ".join(cmd)
    print(f"[find_k] {text}", flush=True)
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if log_path:
        log_path.write_text(proc.stdout, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {text}\n{proc.stdout}")
    return proc.stdout


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"ERROR: required tool not found on PATH: {name}")
    return path


def split_ancestry(raw: str) -> list[str]:
    text = str(raw or "").strip()
    if not text:
        return []
    parts = re.split(r"\s*(?:[;,/|+]|\band\b|\&)\s*", text, flags=re.IGNORECASE)
    return [part.strip() for part in parts if part.strip()]


def canonical_label(raw: str, aliases: dict[str, str]) -> str:
    key = re.sub(r"\s+", " ", str(raw or "").strip().lower())
    return aliases.get(key, str(raw or "").strip())


def load_ancestry_map(path: Path, aliases: dict[str, str]) -> pd.DataFrame:
    sep = "," if path.suffix == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str).fillna("")
    cols = {c.lower(): c for c in df.columns}
    sample_col = (
        cols.get("participant_id")
        or cols.get("sample_id")
        or cols.get("sample")
        or cols.get("iid")
    )
    ancestry_col = (
        cols.get("self_reported_ancestry")
        or cols.get("reported_ancestry")
        or cols.get("ancestry")
        or cols.get("ethnicity")
    )
    if not sample_col or not ancestry_col:
        raise ValueError(
            f"{path}: expected participant_id/sample_id and ancestry columns, got {list(df.columns)}"
        )
    rows = []
    for row in df[[sample_col, ancestry_col]].itertuples(index=False, name=None):
        sample_id, raw = row
        labels = [canonical_label(part, aliases) for part in split_ancestry(raw)]
        labels = [label for label in labels if label]
        rows.append({
            "sample_id": str(sample_id).strip(),
            "reported_ancestry": str(raw).strip(),
            "normalized_ancestries": ";".join(labels),
            "n_ancestries": len(set(labels)),
            "anchor_label": labels[0] if len(set(labels)) == 1 else "",
        })
    return pd.DataFrame(rows)


def read_fam_samples(prefix: Path) -> list[str]:
    fam = prefix.with_suffix(".fam")
    samples = []
    with fam.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                samples.append(fields[1])
    return samples


def parse_cv_error(log_text: str) -> float:
    matches = re.findall(r"CV error[^:]*:\s*([0-9.eE+-]+)", log_text)
    if not matches:
        matches = re.findall(r"CV error.*?([0-9]+(?:\.[0-9]+)?)", log_text)
    return float(matches[-1]) if matches else float("nan")


def prepare_plink(
    input_prefix: Path,
    out_prefix: Path,
    threads: int,
    geno: float,
    mind: float,
    maf: float,
    hwe: float,
    ld_window: int,
    ld_step: int,
    ld_r2: float,
) -> Path:
    plink2 = require_tool("plink2")
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    qc = out_prefix.parent / "autosomes_qc"
    pruned = out_prefix.parent / "autosomes_pruned"
    run([
        plink2,
        "--bfile", str(input_prefix),
        "--chr", "1-22",
        "--rm-dup", "exclude-all",
        "--geno", str(geno),
        "--mind", str(mind),
        "--maf", str(maf),
        "--hwe", str(hwe),
        "--make-bed",
        "--out", str(qc),
        "--allow-no-sex",
        "--threads", str(threads),
    ], out_prefix.parent / "plink_qc.log")
    run([
        plink2,
        "--bfile", str(qc),
        "--indep-pairwise", str(ld_window), str(ld_step), str(ld_r2),
        "--out", str(pruned),
        "--allow-no-sex",
        "--threads", str(threads),
    ], out_prefix.parent / "plink_prune.log")
    run([
        plink2,
        "--bfile", str(qc),
        "--extract", str(pruned) + ".prune.in",
        "--make-bed",
        "--out", str(out_prefix),
        "--allow-no-sex",
        "--threads", str(threads),
    ], out_prefix.parent / "plink_extract_pruned.log")
    return out_prefix


def run_admixture(
    prefix: Path,
    out_dir: Path,
    k_values: list[int],
    threads: int,
    reps: int,
) -> pd.DataFrame:
    admixture = require_tool("admixture")
    rows = []
    for k in k_values:
        for rep in range(1, reps + 1):
            seed = 1729 + (k * 1000) + rep
            log_path = out_dir / f"admixture_K{k}_rep{rep}.log"
            t0 = time.perf_counter()
            text = run([
                admixture,
                "--cv",
                f"-s{seed}",
                str(prefix.with_suffix(".bed")),
                str(k),
                f"-j{threads}",
            ], log_path, cwd=out_dir)
            elapsed = time.perf_counter() - t0
            q_src = out_dir / f"{prefix.name}.{k}.Q"
            p_src = out_dir / f"{prefix.name}.{k}.P"
            q_dst = out_dir / f"admixture_K{k}_rep{rep}.Q"
            p_dst = out_dir / f"admixture_K{k}_rep{rep}.P"
            if q_src.exists():
                q_src.replace(q_dst)
            if p_src.exists():
                p_src.replace(p_dst)
            rows.append({
                "K": k,
                "rep": rep,
                "seed": seed,
                "cv_error": parse_cv_error(text),
                "elapsed_seconds": elapsed,
                "q_file": q_dst.name,
                "p_file": p_dst.name,
                "log_file": log_path.name,
            })
    return pd.DataFrame(rows)


def label_components(
    fam_samples: list[str],
    ancestry_df: pd.DataFrame,
    q_path: Path,
    k: int,
    min_anchor_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    q = np.loadtxt(q_path, dtype=float)
    if q.ndim == 1:
        q = q.reshape(1, -1)
    q_df = pd.DataFrame(q, columns=[f"component_{i + 1}" for i in range(k)])
    q_df.insert(0, "sample_id", fam_samples)

    anchors = ancestry_df[ancestry_df["anchor_label"] != ""].copy()
    anchors = anchors.merge(q_df, on="sample_id", how="inner")
    means = []
    labels = []
    for label, group in anchors.groupby("anchor_label", sort=True):
        if len(group) < min_anchor_n:
            continue
        values = group[[f"component_{i + 1}" for i in range(k)]].mean(axis=0)
        best_component = int(np.argmax(values.to_numpy())) + 1
        for i, value in enumerate(values, start=1):
            means.append({
                "K": k,
                "anchor_label": label,
                "component": i,
                "mean_q": float(value),
                "n_anchor": int(len(group)),
            })
        labels.append({
            "K": k,
            "anchor_label": label,
            "component": best_component,
            "mean_q": float(values.iloc[best_component - 1]),
            "n_anchor": int(len(group)),
        })
    return pd.DataFrame(means), pd.DataFrame(labels)


def write_labeled_q(fam_samples: list[str], q_path: Path, component_labels: pd.DataFrame, out_path: Path) -> None:
    q = np.loadtxt(q_path, dtype=float)
    if q.ndim == 1:
        q = q.reshape(1, -1)
    label_by_component = {
        int(row.component): str(row.anchor_label)
        for row in component_labels.itertuples(index=False)
    }
    cols = [label_by_component.get(i + 1, f"component_{i + 1}") for i in range(q.shape[1])]
    used = {}
    unique_cols = []
    for col in cols:
        used[col] = used.get(col, 0) + 1
        unique_cols.append(col if used[col] == 1 else f"{col}_{used[col]}")
    out = pd.DataFrame(q, columns=unique_cols)
    out.insert(0, "sample_id", fam_samples)
    out.to_csv(out_path, sep="\t", index=False, float_format="%.8g")


def summarize_cv(cv: pd.DataFrame) -> pd.DataFrame:
    if cv.empty:
        return pd.DataFrame(columns=[
            "K", "n_reps", "mean_cv_error", "sd_cv_error", "best_cv_error",
            "best_rep", "mean_elapsed_seconds", "total_elapsed_seconds",
        ])
    rows = []
    for k, group in cv.groupby("K", sort=True):
        best = group.sort_values("cv_error", na_position="last").iloc[0]
        rows.append({
            "K": int(k),
            "n_reps": int(len(group)),
            "mean_cv_error": float(group["cv_error"].mean()),
            "sd_cv_error": float(group["cv_error"].std(ddof=0)) if len(group) > 1 else 0.0,
            "best_cv_error": float(best["cv_error"]),
            "best_rep": int(best["rep"]),
            "mean_elapsed_seconds": float(group["elapsed_seconds"].mean()),
            "total_elapsed_seconds": float(group["elapsed_seconds"].sum()),
        })
    return pd.DataFrame(rows)


def plot_cv_summary(summary: pd.DataFrame, out_path: Path) -> None:
    if summary.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.errorbar(
        summary["K"],
        summary["mean_cv_error"],
        yerr=summary["sd_cv_error"],
        color="#1f4e79",
        marker="o",
        linewidth=1.5,
        capsize=3,
    )
    ax1.set_xlabel("K")
    ax1.set_ylabel("CV error")
    ax1.grid(True, alpha=0.25)

    ax2 = ax1.twinx()
    ax2.bar(
        summary["K"],
        summary["mean_elapsed_seconds"],
        color="#d9a441",
        alpha=0.35,
        width=0.55,
    )
    ax2.set_ylabel("Mean runtime per replicate (seconds)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed-prefix", type=Path, required=True)
    ap.add_argument("--ancestry-map", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--k-min", type=int, default=int(os.environ.get("BV_ADMIXTURE_K_MIN", "2")))
    ap.add_argument("--k-max", type=int, default=int(os.environ.get("BV_ADMIXTURE_K_MAX", "8")))
    ap.add_argument("--threads", type=int, default=int(os.environ.get("BV_THREADS", "8")))
    ap.add_argument("--reps", type=int, default=int(os.environ.get("BV_ADMIXTURE_REPS", "1")))
    ap.add_argument("--min-anchor-n", type=int, default=int(os.environ.get("BV_MIN_ANCHOR_N", "5")))
    ap.add_argument("--skip-plink-prep", action="store_true", help="Use --bed-prefix directly as ADMIXTURE input")
    ap.add_argument("--geno", type=float, default=0.05)
    ap.add_argument("--mind", type=float, default=0.10)
    ap.add_argument("--maf", type=float, default=0.01)
    ap.add_argument("--hwe", type=float, default=1e-6)
    ap.add_argument("--ld-window", type=int, default=200)
    ap.add_argument("--ld-step", type=int, default=50)
    ap.add_argument("--ld-r2", type=float, default=0.2)
    args = ap.parse_args()

    if args.k_min < 2 or args.k_max < args.k_min:
        raise SystemExit("ERROR: expected 2 <= k_min <= k_max")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    aliases = DEFAULT_LABEL_ALIASES.copy()
    ancestry = load_ancestry_map(args.ancestry_map, aliases)
    ancestry.to_csv(args.out_dir / "reported_ancestry_normalized.tsv", sep="\t", index=False)
    ancestry[ancestry["anchor_label"] != ""].to_csv(
        args.out_dir / "ancestry_anchor_samples.tsv", sep="\t", index=False
    )

    if args.skip_plink_prep:
        pruned_prefix = args.bed_prefix
    else:
        pruned_prefix = prepare_plink(
            args.bed_prefix,
            args.out_dir / "admixture_pruned",
            args.threads,
            args.geno,
            args.mind,
            args.maf,
            args.hwe,
            args.ld_window,
            args.ld_step,
            args.ld_r2,
        )
    fam_samples = read_fam_samples(pruned_prefix)
    cv = run_admixture(pruned_prefix, args.out_dir, list(range(args.k_min, args.k_max + 1)), args.threads, args.reps)
    cv.to_csv(args.out_dir / "admixture_cv_errors.tsv", sep="\t", index=False, float_format="%.8g")
    summary = summarize_cv(cv)
    summary.to_csv(args.out_dir / "admixture_k_summary.tsv", sep="\t", index=False, float_format="%.8g")
    plot_cv_summary(summary, args.out_dir / "admixture_k_summary.png")

    all_means = []
    all_labels = []
    for row in cv.itertuples(index=False):
        q_path = args.out_dir / row.q_file
        if not q_path.exists():
            continue
        means, labels = label_components(fam_samples, ancestry, q_path, int(row.K), args.min_anchor_n)
        all_means.append(means)
        all_labels.append(labels)
        if not labels.empty:
            write_labeled_q(fam_samples, q_path, labels, args.out_dir / f"admixture_K{int(row.K)}_rep{int(row.rep)}_labeled_Q.tsv")

    means_df = pd.concat(all_means, ignore_index=True) if all_means else pd.DataFrame()
    labels_df = pd.concat(all_labels, ignore_index=True) if all_labels else pd.DataFrame()
    means_df.to_csv(args.out_dir / "component_anchor_means.tsv", sep="\t", index=False, float_format="%.8g")
    labels_df.to_csv(args.out_dir / "component_labels.tsv", sep="\t", index=False, float_format="%.8g")

    best_k = int(summary.sort_values("mean_cv_error", na_position="last").iloc[0]["K"]) if not summary.empty else -1
    with (args.out_dir / "selected_k.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["selected_k", "criterion", "note"])
        writer.writerow([best_k, "minimum_mean_cv_error", "Review component_labels.tsv for anchor label collisions before final use."])

    with (args.out_dir / "find_k_report.txt").open("w", encoding="utf-8") as handle:
        handle.write("=== sex_biased_admixture_find_k ===\n")
        handle.write(f"Input BED prefix: {args.bed_prefix}\n")
        handle.write(f"Ancestry map: {args.ancestry_map}\n")
        handle.write(f"K range: {args.k_min}-{args.k_max}\n")
        handle.write(f"Replicates per K: {args.reps}\n")
        handle.write(f"Skipped PLINK prep: {args.skip_plink_prep}\n")
        handle.write(f"Samples after pruning: {len(fam_samples)}\n")
        handle.write(f"Single-ancestry anchors: {(ancestry['anchor_label'] != '').sum()}\n")
        handle.write(f"Selected K: {best_k}\n")


if __name__ == "__main__":
    main()

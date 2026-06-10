#!/usr/bin/env python3
"""Regenerate sex-biased ADMIXTURE reports from labelled Q tables.

This is a recovery tool for flow 09 when ADMIXTURE completed but the final
step-8 report/figure writing did not finish. It needs only the labelled
autosome and X Q TSVs.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pandas as pd


META_COLUMNS = {"sample_id", "superpopulation", "sex", "group"}


def infer_k(*paths: Path) -> int:
    for path in paths:
        if not path:
            continue
        match = re.search(r"_K(\d+)_", path.name)
        if match:
            return int(match.group(1))
    return 5


def read_labeled_q(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{name} labelled Q file not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype={"sample_id": str, "sex": str, "group": str})
    missing = {"sample_id", "sex", "group"} - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    ancestry = ancestry_columns(df)
    if not ancestry:
        raise ValueError(f"{path} has no ancestry/component columns")
    for col in ancestry:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def ancestry_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if col not in META_COLUMNS]


def compare_x_vs_auto(auto_df: pd.DataFrame, x_df: pd.DataFrame, k: int,
                      labels: list[str], group: str = "study") -> pd.DataFrame:
    rows = []
    auto_group = auto_df[auto_df["group"] == group]
    x_group = x_df[x_df["group"] == group]
    for lab in labels:
        a = float(auto_group[lab].mean()) if lab in auto_group.columns else np.nan
        xall = float(x_group[lab].mean()) if lab in x_group.columns else np.nan
        xf = float(x_group.loc[x_group["sex"] == "2", lab].mean()) if lab in x_group.columns else np.nan
        xm = float(x_group.loc[x_group["sex"] == "1", lab].mean()) if lab in x_group.columns else np.nan
        rows.append({
            "K": k,
            "ancestry": lab,
            "mean_auto": round(a, 4),
            "mean_x": round(xall, 4),
            "delta_x_minus_auto": round(xall - a, 4),
            "mean_x_female": round(xf, 4),
            "mean_x_male_haploid": round(xm, 4),
        })
    return pd.DataFrame(rows)


def build_per_sample(auto_df: pd.DataFrame, x_df: pd.DataFrame, k: int,
                     labels: list[str], group: str = "study") -> pd.DataFrame:
    a = auto_df[auto_df["group"] == group].copy()
    x = x_df[x_df["group"] == group].copy()
    acols = ["sample_id", "sex"] + [label for label in labels if label in a.columns]
    xcols = ["sample_id"] + [label for label in labels if label in x.columns]
    merged = a[acols].merge(x[xcols], on="sample_id", suffixes=("_auto", "_x"))
    for label in labels:
        ca, cx = f"{label}_auto", f"{label}_x"
        if ca in merged.columns and cx in merged.columns:
            merged[f"{label}_delta"] = (merged[cx] - merged[ca]).round(6)
    merged.insert(1, "K", k)
    merged["sex_label"] = merged["sex"].map({"1": "M", "2": "F"}).fillna("?")
    return merged


def plot_x_vs_auto(comp: pd.DataFrame, out_png: Path) -> None:
    if comp.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ks = sorted(comp["K"].unique())
    fig, axes = plt.subplots(1, len(ks), figsize=(4.2 * len(ks), 4.0), squeeze=False)
    for ax, k in zip(axes[0], ks):
        sub = comp[comp["K"] == k]
        x = np.arange(len(sub))
        width = 0.38
        ax.bar(x - width / 2, sub["mean_auto"], width, label="Autosomes", color="#aaaaaa")
        ax.bar(x + width / 2, sub["mean_x"], width, label="X chromosome", color="#B2182B")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["ancestry"], rotation=30, ha="right", fontsize=8)
        ax.set_title(f"K={k}", fontsize=10)
        ax.set_ylabel("Mean ancestry proportion (study)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0][0].legend(fontsize=8, frameon=False)
    fig.suptitle("Sex-biased admixture: study ancestry on autosomes vs X", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_figure(per_sample: pd.DataFrame, cohort: pd.DataFrame, k: int,
                labels: list[str], out_png: Path, focal_ancestry: str | None = None) -> str | None:
    if per_sample.empty or not labels:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    col_f, col_m = "#B2182B", "#2166AC"
    palette = ["#1A9641", "#4575B4", "#F46D43", "#984EA3", "#FF7F00"]
    sex_color = {"F": col_f, "M": col_m, "?": "#888888"}

    csub = cohort[cohort["K"] == k]
    fallback_focal = (csub.reindex(csub["delta_x_minus_auto"].abs().sort_values(ascending=False).index)
                      ["ancestry"].iloc[0]) if not csub.empty else labels[0]
    focal = focal_ancestry if (
        focal_ancestry
        and f"{focal_ancestry}_auto" in per_sample.columns
        and f"{focal_ancestry}_x" in per_sample.columns
    ) else fallback_focal
    fa, fx, fd = f"{focal}_auto", f"{focal}_x", f"{focal}_delta"

    ps = per_sample.sort_values([fa]).reset_index(drop=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(f"HGP1K sex-biased ADMIXTURE - K={k} (focal ancestry: {focal})", fontsize=12)

    ax = axes[0][0]
    lim = [0, 1]
    ax.plot(lim, lim, ls="--", lw=0.8, color="#aaa")
    for _, row in ps.iterrows():
        ax.scatter(row.get(fa, np.nan), row.get(fx, np.nan), s=45,
                   color=sex_color.get(row["sex_label"], "#888"),
                   edgecolors="k", linewidths=0.4)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(f"{focal} - autosomes")
    ax.set_ylabel(f"{focal} - X")
    ax.set_title("Autosomal vs X ancestry (per individual)")
    ax.legend(handles=[mpatches.Patch(color=col_f, label="Female"),
                       mpatches.Patch(color=col_m, label="Male")], fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[0][1]
    ax.axvline(0, ls="--", lw=0.8, color="#aaa")
    for i, (_, row) in enumerate(ps.iterrows()):
        delta = row.get(fd, 0.0)
        color = sex_color.get(row["sex_label"], "#888")
        ax.hlines(i, 0, delta, color=color, lw=1.2)
        ax.scatter(delta, i, s=35, color=color, edgecolors="k", linewidths=0.4)
    ax.set_yticks(range(len(ps)))
    ax.set_yticklabels(ps["sample_id"], fontsize=6)
    ax.set_xlabel(f"X - autosomal {focal} (delta proportion)")
    ax.set_title("Sex-biased signal per individual")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1][0]
    psx = per_sample.sort_values(["sex_label", fa]).reset_index(drop=True)
    n = len(psx)
    xs = np.arange(n)
    bar_width = 0.38
    for offset, suffix in [(-bar_width / 2, "_auto"), (bar_width / 2, "_x")]:
        bottom = np.zeros(n)
        for idx, label in enumerate(labels):
            col = f"{label}{suffix}"
            if col not in psx.columns:
                continue
            vals = psx[col].fillna(0).to_numpy()
            ax.bar(xs + offset, vals, bar_width, bottom=bottom,
                   color=palette[idx % len(palette)], edgecolor="white", linewidth=0.2,
                   label=label if suffix == "_auto" else "_nolegend_")
            bottom += vals
    ax.set_xticks(xs)
    ax.set_xticklabels(psx["sample_id"], rotation=60, ha="right", fontsize=5)
    ax.set_ylabel("Component proportion")
    ax.set_ylim(0, 1.05)
    ax.set_title("Ancestry: autosomes (A) vs X - per individual")
    ax.legend(fontsize=7, ncol=len(labels), frameon=False, loc="lower center")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1][1]
    cx = np.arange(len(csub))
    bar_width = 0.38
    ax.bar(cx - bar_width / 2, csub["mean_auto"], bar_width, color="#aaaaaa", label="Autosomes")
    ax.bar(cx + bar_width / 2, csub["mean_x"], bar_width, color=col_f, label="X chromosome")
    ax.set_xticks(cx)
    ax.set_xticklabels(csub["ancestry"], fontsize=8)
    ax.set_ylabel("Mean proportion (study)")
    ax.set_title("Cohort mean: autosomes vs X")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return focal


def resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path | None]:
    if args.input_dir:
        input_dir = args.input_dir.expanduser()
        auto = input_dir / "admixture_auto_K5_labeled_Q.tsv"
        x = input_dir / "admixture_x_K5_labeled_Q.tsv"
        combined = input_dir / "admixture_combined_K5_labeled_Q.tsv"
        return auto, x, combined if combined.exists() else None
    if not args.auto or not args.x:
        raise SystemExit("Either --input-dir or both --auto and --x are required")
    return args.auto.expanduser(), args.x.expanduser(), args.combined.expanduser() if args.combined else None


def write_outputs(auto_path: Path, x_path: Path, combined_path: Path | None,
                  out_dir: Path, k: int | None, with_reference: bool) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    k = k or infer_k(auto_path, x_path, combined_path)

    auto_df = read_labeled_q(auto_path, "autosome")
    x_df = read_labeled_q(x_path, "X")
    labels = sorted(set(ancestry_columns(auto_df)) & set(ancestry_columns(x_df)))
    if not labels:
        raise ValueError("No shared labelled ancestry columns between auto and X files")

    written: list[Path] = []
    sex_bias = compare_x_vs_auto(auto_df, x_df, k, labels, group="study")
    sex_bias_path = out_dir / f"sex_bias_x_vs_auto.tsv"
    sex_bias.to_csv(sex_bias_path, sep="\t", index=False)
    written.append(sex_bias_path)

    per_sample = build_per_sample(auto_df, x_df, k, labels, group="study")
    per_sample_path = out_dir / f"sex_bias_per_sample_K{k}.tsv"
    per_sample.to_csv(per_sample_path, sep="\t", index=False, float_format="%.6g")
    written.append(per_sample_path)

    plot_path = out_dir / "sex_bias_x_vs_auto.png"
    plot_x_vs_auto(sex_bias, plot_path)
    written.extend([plot_path, plot_path.with_suffix(".pdf")])

    figure_path = out_dir / f"figure_sex_biased_admixture_K{k}.png"
    focal = plot_figure(per_sample, sex_bias, k, labels, figure_path)
    written.extend([figure_path, figure_path.with_suffix(".pdf")])

    if with_reference and "reference" in set(auto_df["group"]) and "reference" in set(x_df["group"]):
        ref_bias = compare_x_vs_auto(auto_df, x_df, k, labels, group="reference")
        ref_bias_path = out_dir / "sex_bias_x_vs_auto_reference.tsv"
        ref_bias.to_csv(ref_bias_path, sep="\t", index=False)
        written.append(ref_bias_path)
        ref_per_sample = build_per_sample(auto_df, x_df, k, labels, group="reference")
        ref_per_sample_path = out_dir / f"sex_bias_per_sample_K{k}_reference.tsv"
        ref_per_sample.to_csv(ref_per_sample_path, sep="\t", index=False, float_format="%.6g")
        written.append(ref_per_sample_path)
        ref_figure_path = out_dir / f"figure_sex_biased_admixture_K{k}_reference.png"
        plot_figure(ref_per_sample, ref_bias, k, labels, ref_figure_path, focal_ancestry=focal)
        written.extend([ref_figure_path, ref_figure_path.with_suffix(".pdf")])

    for src in [auto_path, x_path, combined_path]:
        if src and src.exists():
            dst = out_dir / src.name
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)
            written.append(dst)

    missing = [path for path in written if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError("Some expected outputs were not written: " + ", ".join(map(str, missing)))
    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate flow-09 sex-biased ADMIXTURE TSVs and figures from labelled Q files."
    )
    parser.add_argument("--input-dir", type=Path,
                        help="Directory containing admixture_auto_K5_labeled_Q.tsv and admixture_x_K5_labeled_Q.tsv")
    parser.add_argument("--auto", type=Path, help="Path to admixture_auto_K*_labeled_Q.tsv")
    parser.add_argument("--x", type=Path, help="Path to admixture_x_K*_labeled_Q.tsv")
    parser.add_argument("--combined", type=Path,
                        help="Optional path to admixture_combined_K*_labeled_Q.tsv; copied through for provenance")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory to write regenerated outputs")
    parser.add_argument("--k", type=int, default=None, help="K value; defaults to parsing filenames, then 5")
    parser.add_argument("--no-reference", action="store_true",
                        help="Skip regenerating reference negative-control TSV/figure outputs")
    args = parser.parse_args()

    auto_path, x_path, combined_path = resolve_inputs(args)
    written = write_outputs(
        auto_path=auto_path,
        x_path=x_path,
        combined_path=combined_path,
        out_dir=args.out_dir.expanduser(),
        k=args.k,
        with_reference=not args.no_reference,
    )

    print("Regenerated outputs:")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()

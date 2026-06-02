#!/usr/bin/env python3
"""
Plot PGx call accumulation across cohort facets.

Input is the flow-level pgx_participant_results.tsv produced from PharmCAT
report TSVs. PharmCAT reports do not carry population allele frequencies, so
"rare/non-reference" here is an operational PGx burden signal: gene calls whose
diplotype is not an obvious reference/wildtype call and whose row is not a
plain no-call.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


FIG_DPI = 300
MIN_TOTAL_CALLS = 1


def is_reference_like(value: str) -> bool:
    text = (value or "").strip().lower()
    if not text:
        return False
    if text in {"no call", "none", "na", "n/a"}:
        return True
    compact = re.sub(r"\s+", " ", text)
    reference_patterns = [
        r"^reference/reference$",
        r"^\*1/\*1$",
        r"^b \(reference\)/b \(reference\)$",
        r"^rs\d+ reference \([acgt]\)/rs\d+ reference \([acgt]\)$",
    ]
    return any(re.match(pattern, compact) for pattern in reference_patterns)


def load_participant_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"participant_id", "country", "sex", "gene", "source_diplotype"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(missing)}")
    for col in ["country", "sex", "gene", "source_diplotype", "phenotype"]:
        df[col] = df[col].astype(str).str.strip()
    df = df[df["gene"] != ""].copy()
    df["is_non_reference_call"] = ~df["source_diplotype"].map(is_reference_like)
    df.loc[df["source_diplotype"].str.lower().eq("no call"), "is_non_reference_call"] = False
    return df


def make_gene_country_tables(df: pd.DataFrame, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = (
        df.groupby(["country", "gene"], dropna=False)
        .agg(
            total_calls=("participant_id", "nunique"),
            non_reference_calls=("is_non_reference_call", "sum"),
        )
        .reset_index()
    )
    grouped = grouped[grouped["total_calls"] >= MIN_TOTAL_CALLS].copy()
    grouped["non_reference_rate"] = grouped["non_reference_calls"] / grouped["total_calls"]
    grouped.to_csv(out_dir / "pgx_gene_country_burden.tsv", sep="\t", index=False, float_format="%.6f")

    by_sex = (
        df.groupby(["country", "sex", "gene"], dropna=False)
        .agg(
            total_calls=("participant_id", "nunique"),
            non_reference_calls=("is_non_reference_call", "sum"),
        )
        .reset_index()
    )
    by_sex = by_sex[by_sex["total_calls"] >= MIN_TOTAL_CALLS].copy()
    by_sex["non_reference_rate"] = by_sex["non_reference_calls"] / by_sex["total_calls"]
    by_sex.to_csv(out_dir / "pgx_gene_country_sex_burden.tsv", sep="\t", index=False, float_format="%.6f")
    return grouped, by_sex


def plot_country_gene_heatmap(grouped: pd.DataFrame, out_dir: Path) -> None:
    matrix = grouped.pivot(index="gene", columns="country", values="non_reference_rate").fillna(0.0)
    counts = grouped.pivot(index="gene", columns="country", values="non_reference_calls").fillna(0).astype(int)
    if matrix.empty:
        return

    gene_order = counts.sum(axis=1).sort_values(ascending=False).index.tolist()
    matrix = matrix.loc[gene_order]
    counts = counts.loc[gene_order]

    fig_h = max(5.5, len(matrix) * 0.34)
    fig_w = max(7.5, len(matrix.columns) * 0.85)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="mako",
        vmin=0,
        vmax=max(1.0, float(matrix.max().max())),
        annot=counts,
        fmt="d",
        linewidths=0.4,
        linecolor="#f2f2f2",
        cbar_kws={"label": "Non-reference PGx call rate"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("Gene")
    ax.set_title("PGx non-reference call accumulation by country", fontsize=11, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    fig.tight_layout()
    for suffix, kwargs in {
        "png": {"dpi": FIG_DPI, "bbox_inches": "tight"},
        "pdf": {"bbox_inches": "tight"},
    }.items():
        fig.savefig(out_dir / f"pgx_non_reference_country_heatmap.{suffix}", **kwargs)
    plt.close(fig)


def plot_country_sex_heatmap(by_sex: pd.DataFrame, out_dir: Path) -> None:
    data = by_sex.copy()
    data["country_sex"] = data["country"] + " | " + data["sex"]
    matrix = data.pivot(index="gene", columns="country_sex", values="non_reference_rate").fillna(0.0)
    counts = data.pivot(index="gene", columns="country_sex", values="non_reference_calls").fillna(0).astype(int)
    if matrix.empty:
        return

    gene_order = counts.sum(axis=1).sort_values(ascending=False).index.tolist()
    matrix = matrix.loc[gene_order]
    counts = counts.loc[gene_order]

    fig_h = max(5.5, len(matrix) * 0.34)
    fig_w = max(8.5, len(matrix.columns) * 0.72)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="rocket_r",
        vmin=0,
        vmax=max(1.0, float(matrix.max().max())),
        annot=counts,
        fmt="d",
        linewidths=0.35,
        linecolor="#f2f2f2",
        cbar_kws={"label": "Non-reference PGx call rate"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("Gene")
    ax.set_title("PGx non-reference call accumulation by country and sex", fontsize=11, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
    fig.tight_layout()
    for suffix, kwargs in {
        "png": {"dpi": FIG_DPI, "bbox_inches": "tight"},
        "pdf": {"bbox_inches": "tight"},
    }.items():
        fig.savefig(out_dir / f"pgx_non_reference_country_sex_heatmap.{suffix}", **kwargs)
    plt.close(fig)


def plot_top_gene_bar(grouped: pd.DataFrame, out_dir: Path, top_n: int = 15) -> None:
    totals = (
        grouped.groupby("gene", as_index=False)["non_reference_calls"]
        .sum()
        .sort_values("non_reference_calls", ascending=False)
        .head(top_n)
    )
    if totals.empty:
        return
    fig, ax = plt.subplots(figsize=(8, max(4.2, len(totals) * 0.32)))
    ax.barh(totals["gene"], totals["non_reference_calls"], color="#3B7A78")
    ax.invert_yaxis()
    ax.set_xlabel("Non-reference PGx calls")
    ax.set_ylabel("")
    ax.set_title("Top PGx genes by non-reference call burden", fontsize=11, fontweight="bold")
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    for suffix, kwargs in {
        "png": {"dpi": FIG_DPI, "bbox_inches": "tight"},
        "pdf": {"bbox_inches": "tight"},
    }.items():
        fig.savefig(out_dir / f"pgx_top_gene_burden.{suffix}", **kwargs)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--participant-results", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_participant_results(args.participant_results)
    grouped, by_sex = make_gene_country_tables(df, args.out_dir)
    plot_country_gene_heatmap(grouped, args.out_dir)
    plot_country_sex_heatmap(by_sex, args.out_dir)
    plot_top_gene_bar(grouped, args.out_dir)


if __name__ == "__main__":
    main()

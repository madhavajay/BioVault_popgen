"""
AIMs step 06 (flow fork) - AIMs panels + clustermaps + combined PCA.

Forked from
04_population_level/aims_differential_snps/scripts/06_AIMs_dendrogram.py.
Only change: population set + display labels + per-population colours are
derived from the country facet (popset) instead of hardcoded dicts. AIMs are
still picked purely from gnomAD reference contrasts (AFR/NFE, AFR/SAS), so the
panel selection is unaffected by which countries the cohort contains.

Env:
  BV_WORK_DIR     AIMs working tree (reads data/master_af_table.tsv)
  BV_RAW_DIR      per-country AF dir (population resolution / fail-loud)
  BV_POPULATIONS  comma-separated normalized country labels

Outputs:
  data/aims/aims_AFR_NFE.tsv  aims_AFR_SAS.tsv  aims_combined.tsv
  plots/aims_AFR_NFE_clustermap.{png,pdf}
  plots/aims_AFR_SAS_clustermap.{png,pdf}
  plots/aims_combined_pca.{png,pdf}
"""

import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from popset import display, require_columns, resolve_populations

BASE = Path(os.environ.get("BV_WORK_DIR", Path(__file__).resolve().parents[1]))
RAW_DIR = Path(os.environ.get("BV_RAW_DIR", BASE.parent / "raw_allele_freq_country"))
MASTER = BASE / "data" / "master_af_table.tsv"
AIM_DIR = BASE / "data" / "aims"
PLOTS_DIR = BASE / "plots"

GNOMAD_REFS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]
TOP_N = 200
MAF_FLOOR = 0.05

POPULATIONS = resolve_populations(RAW_DIR)

_pop_cmap = plt.cm.get_cmap("tab10", max(len(POPULATIONS), 1))
DISPLAY = {p: display(p) for p in POPULATIONS}
DISPLAY.update(
    {
        "gnomAD_global": "gnomAD global",
        "gnomAD_AFR": "gnomAD AFR",
        "gnomAD_NFE": "gnomAD NFE",
        "gnomAD_SAS": "gnomAD SAS",
    }
)
COLORS = {p: matplotlib.colors.to_hex(_pop_cmap(i)) for i, p in enumerate(POPULATIONS)}
COLORS.update(
    {
        "gnomAD_AFR": "#000000",
        "gnomAD_NFE": "#7F7F7F",
        "gnomAD_SAS": "#999999",
        "gnomAD_global": "#444444",
    }
)


def maf(x: pd.Series) -> pd.Series:
    return x.clip(0, 1).combine(1 - x, min)


def pick_aims(df: pd.DataFrame, ref_a: str, ref_b: str, top_n: int) -> pd.DataFrame:
    common = (maf(df[ref_a]) >= MAF_FLOOR) & (maf(df[ref_b]) >= MAF_FLOOR)
    sub = df.loc[common].copy()
    sub["delta_abs"] = (sub[ref_a] - sub[ref_b]).abs()
    sub["pooled_af"] = (sub[ref_a] + sub[ref_b]) / 2
    sub = (
        sub.sort_values(["delta_abs", "pooled_af"], ascending=[False, True])
        .head(top_n)
        .reset_index(drop=True)
    )
    sub["rank"] = sub.index + 1
    return sub


def write_panel(panel: pd.DataFrame, name: str) -> Path:
    cols = ["locus_key", "rsid", "rank", "delta_abs", "pooled_af"] + POPULATIONS + GNOMAD_REFS
    f = AIM_DIR / f"{name}.tsv"
    panel[cols].to_csv(f, sep="\t", index=False, float_format="%.6f")
    return f


def plot_clustermap(panel: pd.DataFrame, refs_in_view, label_a, label_b, fname_stem) -> None:
    cols = POPULATIONS + refs_in_view
    M = panel.set_index("locus_key")[cols].copy()
    M.columns = [DISPLAY[c] for c in cols]
    col_colors = [COLORS[c] for c in cols]

    g = sns.clustermap(
        M,
        cmap="RdBu_r",
        center=0.5,
        vmin=0,
        vmax=1,
        figsize=(6.0, 9.0),
        method="average",
        metric="euclidean",
        col_cluster=True,
        row_cluster=True,
        col_colors=[col_colors],
        cbar_kws={"label": "Allele frequency"},
        xticklabels=True,
        yticklabels=False,
        dendrogram_ratio=(0.18, 0.12),
        colors_ratio=0.025,
        linewidths=0,
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel(
        f"AIMs (top {TOP_N} {label_a}<->{label_b}, |dAF|-ranked)", fontsize=8
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    g.fig.suptitle(
        f"AIMs from gnomAD {label_a} vs {label_b} - applied to cohort populations",
        fontsize=9,
        y=1.0,
    )
    pdf = PLOTS_DIR / f"{fname_stem}.pdf"
    png = PLOTS_DIR / f"{fname_stem}.png"
    g.savefig(pdf, bbox_inches="tight")
    g.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(g.fig)
    print(f"[06] wrote {pdf}")
    print(f"[06] wrote {png}")


def pca_on_aims(combined: pd.DataFrame) -> None:
    cohorts = POPULATIONS + ["gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS", "gnomAD_global"]
    X = combined[cohorts].to_numpy().T
    X = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    PC = U * S
    var_explained = (S**2) / np.sum(S**2) * 100

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for i, c in enumerate(cohorts):
        is_ref = c.startswith("gnomAD")
        ax.scatter(
            PC[i, 0],
            PC[i, 1],
            s=110 if not is_ref else 80,
            color=COLORS.get(c, "#444"),
            edgecolor="black",
            linewidth=0.6,
            marker="o" if not is_ref else "s",
            label=DISPLAY[c].replace("\n", " "),
            zorder=3,
        )
        ax.annotate(
            DISPLAY[c].replace("\n", " "),
            (PC[i, 0], PC[i, 1]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=6.5,
            zorder=4,
        )
    ax.axhline(0, color="grey", lw=0.4, ls="--", zorder=1)
    ax.axvline(0, color="grey", lw=0.4, ls="--", zorder=1)
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    ax.set_title(
        f"PCA of cohorts on combined AIMs\n"
        f"(top {TOP_N} AFR<->NFE U top {TOP_N} AFR<->SAS, deduplicated)",
        fontsize=8,
    )
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    pdf = PLOTS_DIR / "aims_combined_pca.pdf"
    png = PLOTS_DIR / "aims_combined_pca.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[06] wrote {pdf}")
    print(f"[06] wrote {png}")


def main() -> None:
    AIM_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(MASTER, sep="\t")
    require_columns(master.columns, POPULATIONS, f"master_af_table {MASTER}")
    print(f"[06] master rows: {len(master):,}")

    panel_ne = pick_aims(master, "gnomAD_AFR", "gnomAD_NFE", TOP_N)
    panel_sa = pick_aims(master, "gnomAD_AFR", "gnomAD_SAS", TOP_N)
    write_panel(panel_ne, "aims_AFR_NFE")
    write_panel(panel_sa, "aims_AFR_SAS")
    print(f"[06] AFR/NFE panel: {len(panel_ne)}; AFR/SAS panel: {len(panel_sa)}")

    plot_clustermap(
        panel_ne,
        refs_in_view=["gnomAD_AFR", "gnomAD_NFE", "gnomAD_global"],
        label_a="AFR",
        label_b="NFE",
        fname_stem="aims_AFR_NFE_clustermap",
    )
    plot_clustermap(
        panel_sa,
        refs_in_view=["gnomAD_AFR", "gnomAD_SAS", "gnomAD_global"],
        label_a="AFR",
        label_b="SAS",
        fname_stem="aims_AFR_SAS_clustermap",
    )

    combined = pd.concat([panel_ne, panel_sa]).drop_duplicates(subset="locus_key")
    combined.to_csv(
        AIM_DIR / "aims_combined.tsv", sep="\t", index=False, float_format="%.6f"
    )
    print(f"[06] combined AIMs (dedup): {len(combined)}")
    pca_on_aims(combined)


if __name__ == "__main__":
    main()

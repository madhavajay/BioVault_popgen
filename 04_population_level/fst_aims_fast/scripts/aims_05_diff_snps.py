"""
AIMs step 05 (flow fork) - per-population differential SNPs vs gnomAD.

Forked from
04_population_level/aims_differential_snps/scripts/05_differential_snps_per_island.py.
Only change: the population set is the country facet (popset) instead of a
hardcoded ISLANDS list, and DISPLAY labels are derived from popset.display.

Env:
  BV_WORK_DIR     AIMs working tree (reads data/master_af_table.tsv)
  BV_RAW_DIR      per-country AF dir (population resolution / fail-loud)
  BV_POPULATIONS  comma-separated normalized country labels

Outputs:
  data/differential_snps/<pop>_vs_<ref>_<enriched|depleted>.tsv
  data/differential_snps/all_outliers_long.tsv
  plots/diff_snps_heatmap_vs_AFR.{png,pdf}
  plots/diff_snps_heatmap_vs_global.{png,pdf}
"""

import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import leaves_list, linkage

from popset import display, require_columns, resolve_populations

BASE = Path(os.environ.get("BV_WORK_DIR", Path(__file__).resolve().parents[1]))
RAW_DIR = Path(os.environ.get("BV_RAW_DIR", BASE.parent / "raw_allele_freq_country"))
MASTER = BASE / "data" / "master_af_table.tsv"
DIFF_DIR = BASE / "data" / "differential_snps"
PLOTS_DIR = BASE / "plots"

GNOMAD_REFS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]
TOP_N = 25
MIN_REF_AF = 0.01

POPULATIONS = resolve_populations(RAW_DIR)
DISPLAY = {p: display(p) for p in POPULATIONS}
DISPLAY.update(
    {
        "gnomAD_global": "gnomAD\nglobal",
        "gnomAD_AFR": "gnomAD\nAFR",
        "gnomAD_NFE": "gnomAD\nNFE",
        "gnomAD_SAS": "gnomAD\nSAS",
    }
)


def collect_outliers(df: pd.DataFrame, reference: str) -> pd.DataFrame:
    rows = []
    af_ref = df[reference]
    common = af_ref.clip(0, 1).combine(1 - af_ref, min) >= MIN_REF_AF
    sub = df.loc[common].copy()
    for pop in POPULATIONS:
        d = sub[pop] - af_ref.loc[sub.index]
        s = sub.assign(delta=d, island=pop, reference=reference)
        top = s.nlargest(TOP_N, "delta").assign(
            direction="enriched", rank=range(1, TOP_N + 1)
        )
        bot = s.nsmallest(TOP_N, "delta").assign(
            direction="depleted", rank=range(1, TOP_N + 1)
        )
        rows.append(top)
        rows.append(bot)
    return pd.concat(rows, ignore_index=True)


def write_per_pop_tables(out: pd.DataFrame, ref_short: str) -> None:
    cols = ["locus_key", "rsid", "delta", "rank"] + POPULATIONS + GNOMAD_REFS
    for pop in POPULATIONS:
        for direc in ("enriched", "depleted"):
            sub = (
                out.query("island == @pop and direction == @direc")
                .sort_values("rank")[cols]
            )
            f = DIFF_DIR / f"{pop}_vs_{ref_short}_{direc}.tsv"
            sub.to_csv(f, sep="\t", index=False, float_format="%.6f")


def diverging_cmap():
    return LinearSegmentedColormap.from_list(
        "div", ["#2166AC", "#F7F7F7", "#B2182B"]
    )


def plot_heatmap(out: pd.DataFrame, master: pd.DataFrame, reference: str, ref_short: str) -> None:
    keys = (
        out.sort_values(["direction", "island", "rank"])["locus_key"]
        .drop_duplicates()
        .tolist()
    )
    M = master.set_index("locus_key").loc[keys, POPULATIONS + GNOMAD_REFS]

    delta = M[POPULATIONS].sub(M[reference], axis=0).to_numpy()
    if len(delta) > 2:
        Z = linkage(delta, method="average", metric="euclidean")
        order = leaves_list(Z)
        M = M.iloc[order]
        keys = [keys[i] for i in order]

    rsid_lookup = master.set_index("locus_key")["rsid"].to_dict()
    row_labels = [rsid_lookup.get(k) or k for k in keys]

    n_rows = len(M)
    fig, ax = plt.subplots(figsize=(7.5, max(6, n_rows * 0.13)))
    cmap = diverging_cmap()
    norm = TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=1.0)
    im = ax.imshow(M.to_numpy(), aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(
        [DISPLAY.get(c, c) for c in M.columns], rotation=0, fontsize=6.5
    )
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(row_labels, fontsize=4.5)
    ax.tick_params(axis="x", which="both", length=0)
    ax.tick_params(axis="y", which="both", length=0)

    ref_idx = M.columns.get_loc(reference)
    ax.axvline(ref_idx - 0.5, color="black", lw=0.4)
    ax.axvline(ref_idx + 0.5, color="black", lw=0.4)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Allele frequency", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.set_title(
        f"Per-population top differential SNPs vs {DISPLAY[reference].replace(chr(10), ' ')}\n"
        f"(top {TOP_N} enriched + {TOP_N} depleted per population; "
        f"AF_ref >= {MIN_REF_AF*100:.0f}% common)",
        fontsize=8,
    )
    fig.tight_layout()
    pdf = PLOTS_DIR / f"diff_snps_heatmap_vs_{ref_short}.pdf"
    png = PLOTS_DIR / f"diff_snps_heatmap_vs_{ref_short}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[05] wrote {pdf}")
    print(f"[05] wrote {png}")


def main() -> None:
    DIFF_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(MASTER, sep="\t")
    require_columns(master.columns, POPULATIONS, f"master_af_table {MASTER}")
    print(f"[05] master table rows: {len(master):,}")

    long_rows = []
    for ref, short in [("gnomAD_AFR", "AFR"), ("gnomAD_global", "global")]:
        print(f"[05] -- differential SNPs vs {ref} --")
        out = collect_outliers(master, ref)
        write_per_pop_tables(out, short)
        long_rows.append(out)
        plot_heatmap(out, master, ref, short)

    long = pd.concat(long_rows, ignore_index=True)
    long_path = DIFF_DIR / "all_outliers_long.tsv"
    long.to_csv(long_path, sep="\t", index=False, float_format="%.6f")
    print(f"[05] wrote {long_path}")


if __name__ == "__main__":
    main()

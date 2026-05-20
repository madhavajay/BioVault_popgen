"""
sex_biased_admixture.py  —  v3 (facet-driven)
==============================================
Sex-biased admixture analysis over a cohort of GSA genotype files
(DDNA or Illumina GSGT — read via the shared genoio adapter).

What is computed
----------------
• Per-sample, per-chromosome SNP counts, heterozygosity rate, BAF
  mean/variance, and LRR mean  (read directly from the genotype files)
• Ancestry-like components via NMF(k=3) on subsampled BAF values,
  autosomes and X chromosome treated separately
• Per-sex contrast of the above; the key panel (d) compares
  autosomal vs X heterozygosity by sex.

Sex assignment
--------------
Sex is the **`sex` participant facet**, passed through exactly like the
`country`/`island` facet: `assign_sex()` reads a `sex_mapping.tsv`
(`participant_id<TAB>sex`, values M/F/Male/Female/1/2) located via
$BIOVAULT_SEX_MAPPING or next to the genotype data. If no mapping is
found it falls back to a deterministic positional split (legacy mock
behaviour) and logs a loud warning — that fallback is NOT real sex.

Synthetic positive control
--------------------------
Plain `--alt-frequency 0.5` data carries no sex signal, so panel d is
flat by construction. To get a verifiable known-truth signal the
synthetic generator injects an X-chromosome block that mimics male
X-hemizygosity: Male → homozygous (low X heterozygosity), Female →
heterozygous (normal diploid X het). That reproduces the canonical
sex-biased signature panel d is designed to detect (male X-het ≪
female X-het, autosomal het ≈ equal). With no injection / real data the
flat result is still the correct output.

Outputs
-------
  results/sex_bias_results.tsv
  plots/figure4_sex_biased_admixture.pdf
  plots/figure4_sex_biased_admixture.png
  logs/sex_biased_admixture.log
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
from sklearn.decomposition import NMF

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genoio as _genoio  # noqa: E402  (synced fork, see genoio.py header)

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
BASE        = SCRIPT_DIR.parent                          # sex_biased_admixture/
REPO_ROOT   = BASE.parents[1]                            # BioVault_popgen/
DATA_DIR    = REPO_ROOT / "01_mock_data_generation" / "output"
PLOTS_DIR   = BASE / "plots"
RESULTS_DIR = BASE / "results"
LOGS_DIR    = BASE / "logs"

for _d in [PLOTS_DIR, RESULTS_DIR, LOGS_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "sex_biased_admixture.log"),
    ],
)
log = logging.getLogger(__name__)

# ── Plot style ─────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "axes.labelsize":    7,
    "axes.titlesize":    7.5,
    "xtick.labelsize":   6.5,
    "ytick.labelsize":   6.5,
    "legend.fontsize":   6.5,
    "axes.linewidth":    0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
})

COL_FEMALE = "#B2182B"
COL_MALE   = "#2166AC"
COL_C1     = "#1A9641"   # Component 1
COL_C2     = "#4575B4"   # Component 2
COL_C3     = "#F46D43"   # Component 3

AUTOSOME_CHROMS = [str(c) for c in range(1, 23)]
X_CHROM         = "X"
NMF_SEED        = 42
SUBSAMPLE_STEP  = 14   # autosomes: ~686 K SNPs → ~49 K after subsampling


# ── Step 1: Load raw genotype files ───────────────────────────────────────────
def read_txt(path: Path) -> pd.DataFrame:
    """Read one GSA .txt file; return DataFrame with chrom/baf/lrr/is_het.

    Illumina GSGT synthetic carries no BAF/LRR, so those come back NaN and
    the BAF-derived NMF/stat columns are uninformative for Illumina samples
    (genotype-based het_rate is still valid)."""
    if _genoio.sniff_format(path) == "illumina":
        g = _genoio.read_genotypes(path)
        df = g.rename(columns={"gt": "genotype"})[
            ["rsid", "chrom", "pos", "genotype", "gs", "baf", "lrr"]]
    else:
        df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            header=None,
            names=["rsid", "chrom", "pos", "genotype", "gs", "baf", "lrr"],
            dtype={"chrom": str, "genotype": str},
        )
    df["baf"] = pd.to_numeric(df["baf"], errors="coerce")
    df["lrr"] = pd.to_numeric(df["lrr"], errors="coerce")
    # Heterozygous = two different, non-missing alleles
    df["is_het"] = df["genotype"].apply(
        lambda g: len(g) == 2 and g[0] != g[1] and "0" not in g
    )
    return df


def load_all_samples(data_dir: Path) -> dict:
    """Return {sample_id: DataFrame} for all numeric-named subdirectories."""
    dirs = sorted(d for d in data_dir.iterdir() if d.is_dir() and d.name.isdigit())
    samples = {}
    for d in dirs:
        txt = list(d.glob("*.txt"))
        if not txt:
            log.warning(f"  No .txt file in {d}; skipping")
            continue
        log.info(f"  Reading {d.name} …")
        samples[d.name] = read_txt(txt[0])
    return samples


# ── Step 2: Sex assignment ─────────────────────────────────────────────────────
_SEX_NORM = {
    "m": "M", "male": "M", "1": "M",
    "f": "F", "female": "F", "2": "F",
}


def _load_sex_mapping() -> dict:
    """Read participant_id → {M,F} from a sex facet mapping file, mirroring
    how the island/country facet is consumed elsewhere.

    Lookup order: $BIOVAULT_SEX_MAPPING, then sex_mapping.tsv next to the
    genotype data (DATA_DIR) or one level up. Returns {} if none found.
    """
    import os

    candidates = []
    env = os.environ.get("BIOVAULT_SEX_MAPPING")
    if env:
        candidates.append(Path(env))
    candidates += [
        DATA_DIR / "sex_mapping.tsv",
        DATA_DIR.parent / "sex_mapping.tsv",
    ]
    for path in candidates:
        if not path.is_file():
            continue
        mapping = {}
        with path.open() as f:
            for i, line in enumerate(f):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                pid, raw = parts[0].strip(), parts[1].strip()
                if i == 0 and raw.lower() in ("sex", "gender"):
                    continue  # header row
                sx = _SEX_NORM.get(raw.lower())
                if sx:
                    mapping[pid] = sx
        if mapping:
            log.info(f"Sex facet: loaded {len(mapping)} labels from {path}")
            return mapping
    return {}


def assign_sex(sample_ids: list) -> dict:
    """Sex per participant from the sex facet mapping file. Falls back to a
    deterministic positional split (first half → F) only if no mapping is
    available, so legacy mock runs still work — emits a clear warning."""
    mapping = _load_sex_mapping()
    if mapping:
        missing = [s for s in sample_ids if s not in mapping]
        if missing:
            log.warning(
                f"Sex facet: {len(missing)} sample(s) absent from mapping; "
                f"defaulting them to 'M' (e.g. {missing[:3]})")
        return {sid: mapping.get(sid, "M") for sid in sample_ids}

    log.warning(
        "Sex facet: no sex_mapping.tsv found — falling back to positional "
        "split (first half sorted = Female). This is NOT real sex; provide "
        "a sex_mapping.tsv (participant_id<TAB>sex) for valid results.")
    sids = sorted(sample_ids)
    n = len(sids)
    return {sid: ("F" if i < n // 2 else "M") for i, sid in enumerate(sids)}


# ── Step 3: Per-chromosome statistics ─────────────────────────────────────────
def chrom_stats_for_sample(df: pd.DataFrame, chroms: list) -> dict:
    """Return per-chromosome dict of n_snps/het_rate/baf_mean/baf_var/lrr_mean."""
    result = {}
    for ch in chroms:
        sub = df[df["chrom"] == ch]
        n   = len(sub)
        if n == 0:
            result[ch] = dict(n=0, het_rate=np.nan, baf_mean=np.nan,
                              baf_var=np.nan, lrr_mean=np.nan)
        else:
            result[ch] = dict(
                n        = n,
                het_rate = float(sub["is_het"].mean()),
                baf_mean = float(sub["baf"].mean()),
                baf_var  = float(sub["baf"].var()),
                lrr_mean = float(sub["lrr"].mean()),
            )
    return result


# ── Step 4: NMF ancestry components ───────────────────────────────────────────
def nmf_components(samples: dict, chroms: list,
                   subsample_step: int = 1, seed: int = NMF_SEED) -> pd.DataFrame:
    """
    Build a BAF matrix (n_samples × n_SNPs), run NMF(k=3), normalise
    each row to sum to 1, and return a DataFrame with columns c1/c2/c3.

    Note: with --alt-frequency 0.5 mock data every row is ~[0.33, 0.33, 0.33].
    """
    sample_ids = sorted(samples)
    ref        = samples[sample_ids[0]]

    # Build SNP index (same positions for all samples — same array)
    mask    = ref["chrom"].isin(chroms)
    snp_idx = np.where(mask)[0][::subsample_step]
    log.info(f"    NMF: {len(snp_idx):,} SNPs from chroms {chroms[:3]}{'…' if len(chroms)>3 else ''}")

    # Matrix: rows = samples, cols = SNPs
    mat = np.empty((len(sample_ids), len(snp_idx)), dtype=np.float32)
    for i, sid in enumerate(sample_ids):
        baf         = samples[sid]["baf"].values[snp_idx]
        mat[i]      = np.where(np.isnan(baf), 0.5, baf)

    model = NMF(n_components=3, random_state=seed, max_iter=1000)
    W     = model.fit_transform(mat)          # (n_samples, 3)
    row_s = W.sum(axis=1, keepdims=True)
    row_s = np.where(row_s == 0, 1.0, row_s)
    W_norm = W / row_s

    log.info(f"    NMF reconstruction error: {model.reconstruction_err_:.4f}")
    return pd.DataFrame(W_norm, index=sample_ids, columns=["c1", "c2", "c3"])


# ── Step 5: Assemble results table ────────────────────────────────────────────
def build_results_table(samples: dict) -> pd.DataFrame:
    sample_ids = sorted(samples)
    sex_map    = assign_sex(sample_ids)

    log.info("Running NMF on autosomes …")
    auto_nmf = nmf_components(samples, AUTOSOME_CHROMS,
                              subsample_step=SUBSAMPLE_STEP)

    log.info("Running NMF on X chromosome …")
    x_nmf    = nmf_components(samples, [X_CHROM], subsample_step=1)

    rows = []
    for sid in sample_ids:
        df  = samples[sid]
        ast = chrom_stats_for_sample(df, AUTOSOME_CHROMS)
        xst = chrom_stats_for_sample(df, [X_CHROM]).get(X_CHROM, {})

        auto_snps = sum(v["n"]        for v in ast.values())
        auto_het  = float(np.nanmean([v["het_rate"] for v in ast.values()]))
        auto_lrr  = float(np.nanmean([v["lrr_mean"] for v in ast.values()]))

        rows.append({
            "sample_id":       sid,
            "sex":             sex_map[sid],
            # Raw statistics from mock data
            "auto_snps":       auto_snps,
            "x_snps":          xst.get("n", 0),
            "auto_het":        round(auto_het, 4),
            "x_het":           round(xst.get("het_rate", np.nan), 4),
            "auto_lrr_mean":   round(auto_lrr, 4),
            "x_lrr_mean":      round(xst.get("lrr_mean", np.nan), 4),
            # NMF ancestry components
            "auto_c1":         round(auto_nmf.loc[sid, "c1"], 4),
            "auto_c2":         round(auto_nmf.loc[sid, "c2"], 4),
            "auto_c3":         round(auto_nmf.loc[sid, "c3"], 4),
            "x_c1":            round(x_nmf.loc[sid, "c1"], 4),
            "x_c2":            round(x_nmf.loc[sid, "c2"], 4),
            "x_c3":            round(x_nmf.loc[sid, "c3"], 4),
            "x_minus_auto_c1": round(
                x_nmf.loc[sid, "c1"] - auto_nmf.loc[sid, "c1"], 4),
        })

    return pd.DataFrame(rows)


# ── Step 6: Figure ─────────────────────────────────────────────────────────────
def plot_figure(df: pd.DataFrame):
    """
    4-panel figure mirroring the sex-biased admixture layout.

    a  Scatter: auto_c1 vs x_c1, by sex
    b  Lollipop: x_c1 − auto_c1 per individual
    c  Grouped stacked bars: NMF components on autosomes vs X
    d  Heterozygosity rate: autosomes vs X per individual
       (in real data males show lower X het; mock data shows no difference)
    """
    sex_color = {"F": COL_FEMALE, "M": COL_MALE}
    colors    = [COL_C1, COL_C2, COL_C3]
    comp_lbl  = ["Component 1", "Component 2", "Component 3"]

    females = df[df["sex"] == "F"].sort_values("auto_c1")
    males   = df[df["sex"] == "M"].sort_values("auto_c1")
    ordered = pd.concat([females, males]).reset_index(drop=True)
    n_f     = len(females)
    n_m     = len(males)

    fig = plt.figure(figsize=(9.0, 7.2))
    gs  = gridspec.GridSpec(
        2, 2, figure=fig,
        hspace=0.48, wspace=0.38,
        left=0.08, right=0.97, top=0.93, bottom=0.10,
    )

    def _panel_label(ax, letter, x=-0.13, y=1.10):
        ax.text(x, y, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="left")

    # ── Panel a: Scatter ───────────────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, 0])

    all_c1 = pd.concat([df["auto_c1"], df["x_c1"]]).dropna()
    margin = 0.05
    lo     = max(0.0, all_c1.min() - margin)
    hi     = min(1.0, all_c1.max() + margin)
    lims   = [lo, hi]

    ax_a.plot(lims, lims, color="#aaaaaa", lw=0.9, ls="--", zorder=1,
              label="_nolegend_")
    ax_a.fill_between(lims, lims, [hi, hi],
                      color=COL_FEMALE, alpha=0.06, zorder=0)
    ax_a.fill_between(lims, [lo, lo], lims,
                      color=COL_MALE,   alpha=0.06, zorder=0)

    for _, row in df.iterrows():
        ax_a.scatter(row.auto_c1, row.x_c1, s=55,
                     color=sex_color[row.sex], edgecolors="k",
                     linewidths=0.55, zorder=4)
        ax_a.annotate(row.sample_id, (row.auto_c1, row.x_c1),
                      textcoords="offset points", xytext=(5, 3),
                      fontsize=5.0, color="#333333")

    ax_a.set_xlim(lims); ax_a.set_ylim(lims)
    ax_a.set_xlabel("NMF Component 1 — autosomes",   fontsize=7)
    ax_a.set_ylabel("NMF Component 1 — X chromosome", fontsize=7)
    ax_a.set_title("Autosomal vs X-linked ancestry component",
                   fontsize=7.5, pad=3, fontweight="normal")
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.text(0.97, 0.97, "X enriched", transform=ax_a.transAxes,
              fontsize=5.5, color=COL_FEMALE, ha="right", va="top",
              style="italic")
    ax_a.text(0.03, 0.03, "Auto enriched", transform=ax_a.transAxes,
              fontsize=5.5, color=COL_MALE, ha="left", va="bottom",
              style="italic")
    ax_a.legend(
        handles=[
            mpatches.Patch(color=COL_FEMALE, label=f"Female (n={n_f})"),
            mpatches.Patch(color=COL_MALE,   label=f"Male (n={n_m})"),
            Line2D([0], [0], color="#aaaaaa", ls="--", lw=0.9,
                   label="X = Auto (null)"),
        ],
        fontsize=5.8, framealpha=0.7, edgecolor="#cccccc", loc="upper left",
    )
    _panel_label(ax_a, "a")

    # ── Panel b: Lollipop ─────────────────────────────────────────────────────
    ax_b  = fig.add_subplot(gs[0, 1])
    y_pos = np.arange(len(ordered))

    ax_b.axvline(0, color="#aaaaaa", lw=0.8, ls="--", zorder=1)
    for i, row in ordered.iterrows():
        col  = sex_color[row.sex]
        diff = row.x_minus_auto_c1
        ax_b.hlines(i, 0, diff, color=col, lw=1.2, zorder=2)
        ax_b.scatter(diff, i, s=45, color=col,
                     edgecolors="k", linewidths=0.5, zorder=3)

    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(
        [f"{r.sample_id}  {'♀' if r.sex=='F' else '♂'}"
         for _, r in ordered.iterrows()],
        fontsize=5.8,
    )
    for ytl, (_, row) in zip(ax_b.get_yticklabels(), ordered.iterrows()):
        ytl.set_color(sex_color[row.sex])

    ax_b.set_xlabel("X − Autosomal Component 1 (Δ proportion)", fontsize=7)
    ax_b.set_title("Sex-biased signal: X vs autosomal ancestry component",
                   fontsize=7.5, pad=3, fontweight="normal")
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
    _panel_label(ax_b, "b")

    # ── Panel c: Stacked bars ─────────────────────────────────────────────────
    ax_c     = fig.add_subplot(gs[1, 0])
    n_ind    = len(ordered)
    x_coords = np.arange(n_ind)
    bar_w    = 0.35
    gap      = 0.10

    for col_idx, (c1_col, c2_col, c3_col, region) in enumerate([
        ("auto_c1", "auto_c2", "auto_c3", "Autosomes"),
        ("x_c1",   "x_c2",   "x_c3",   "X chromosome"),
    ]):
        xpos   = x_coords + (col_idx - 0.5) * (bar_w + gap)
        bottom = np.zeros(n_ind)
        for col, color, lbl in zip([c1_col, c2_col, c3_col], colors, comp_lbl):
            vals = ordered[col].values
            ax_c.bar(xpos, vals, bar_w, bottom=bottom, color=color,
                     edgecolor="white", linewidth=0.3,
                     label=lbl if col_idx == 0 else "_nolegend_")
            bottom += vals
        letter = "A" if col_idx == 0 else "X"
        for xp in xpos:
            ax_c.text(xp, 1.015, letter, ha="center", va="bottom",
                      fontsize=4.5, color="#555555", fontstyle="italic")

    ax_c.axvline(n_f - 0.5, color="#555555", lw=0.7, ls=":", zorder=5)
    ax_c.text(n_f / 2 - 0.5, 1.07, "Females", ha="center", fontsize=6,
              color=COL_FEMALE, transform=ax_c.get_xaxis_transform())
    ax_c.text(n_f + n_m / 2 - 0.5, 1.07, "Males", ha="center", fontsize=6,
              color=COL_MALE, transform=ax_c.get_xaxis_transform())

    ax_c.set_xticks(x_coords)
    ax_c.set_xticklabels(ordered["sample_id"].values,
                          rotation=35, ha="right", fontsize=5.0)
    ax_c.set_ylim(0, 1.14)
    ax_c.set_ylabel("NMF component proportion", fontsize=7)
    ax_c.set_title("NMF ancestry decomposition: autosomes (A) vs X (X)",
                   fontsize=7.5, pad=3, fontweight="normal")
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.legend(
        handles=[mpatches.Patch(color=c, label=l)
                 for c, l in zip(colors, comp_lbl)],
        fontsize=5.8, framealpha=0.7, edgecolor="#cccccc",
        loc="lower right", ncol=3,
    )
    _panel_label(ax_c, "c")

    # ── Panel d: Heterozygosity auto vs X ─────────────────────────────────────
    ax_d     = fig.add_subplot(gs[1, 1])
    x_coords = np.arange(n_ind)
    bar_w    = 0.35

    ax_d.bar(x_coords - bar_w / 2, ordered["auto_het"].values, bar_w,
             color="#aaaaaa", edgecolor="white", linewidth=0.3,
             label="Autosomes")
    for xi, (_, row) in zip(x_coords + bar_w / 2, ordered.iterrows()):
        ax_d.bar(xi, row["x_het"], bar_w,
                 color=sex_color[row.sex], edgecolor="white", linewidth=0.3)

    ax_d.axvline(n_f - 0.5, color="#555555", lw=0.7, ls=":", zorder=5)
    ax_d.text(n_f / 2 - 0.5, 1.03, "Females", ha="center", fontsize=6,
              color=COL_FEMALE, transform=ax_d.get_xaxis_transform())
    ax_d.text(n_f + n_m / 2 - 0.5, 1.03, "Males", ha="center", fontsize=6,
              color=COL_MALE, transform=ax_d.get_xaxis_transform())

    ax_d.set_xticks(x_coords)
    ax_d.set_xticklabels(ordered["sample_id"].values,
                          rotation=35, ha="right", fontsize=5.0)
    ax_d.set_ylabel("Heterozygosity rate", fontsize=7)
    # Data-driven annotation: the panel must describe what THIS cohort
    # actually shows (biased data carries a real X signal; plain
    # --alt-frequency data does not), not a hardcoded assumption.
    _mean_auto = float(ordered["auto_het"].mean())
    _xf = float(ordered["x_het"].iloc[:n_f].mean()) if n_f else float("nan")
    _xm = float(ordered["x_het"].iloc[n_f:].mean()) if n_m else float("nan")
    ax_d.set_title(
        "Heterozygosity: autosomes vs X chromosome\n"
        f"(mean X het — female {_xf:.2f}, male {_xm:.2f}; "
        f"autosomal {_mean_auto:.2f})",
        fontsize=7.0, pad=3, fontweight="normal",
    )
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.legend(
        handles=[
            mpatches.Patch(color="#aaaaaa",  label="Autosomes"),
            mpatches.Patch(color=COL_FEMALE, label="X – Female"),
            mpatches.Patch(color=COL_MALE,   label="X – Male"),
        ],
        fontsize=5.8, framealpha=0.7, edgecolor="#cccccc", loc="upper right",
    )
    _dx = _xf - _xm
    if np.isfinite(_dx) and abs(_dx) >= 0.05:
        _note = (
            f"Female X het {'>' if _dx > 0 else '<'} male X het by "
            f"{abs(_dx):.2f} — sex-biased X signal present"
        )
    else:
        _note = "No sex difference in X heterozygosity (expected for unbiased data)"
    ax_d.text(
        0.5, -0.22, _note,
        transform=ax_d.transAxes, ha="center", fontsize=5.5,
        color="#888888", style="italic",
    )
    _panel_label(ax_d, "d")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_pdf = PLOTS_DIR / "figure4_sex_biased_admixture.pdf"
    out_png = PLOTS_DIR / "figure4_sex_biased_admixture.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"Saved PDF → {out_pdf}")
    log.info(f"Saved PNG → {out_png}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    log.info("=== Sex-biased admixture v2 — mock-data-driven ===")
    log.info(f"Data dir : {DATA_DIR}")

    if not DATA_DIR.exists():
        log.error(f"Data directory not found: {DATA_DIR}")
        sys.exit(1)

    log.info("Loading genotype files …")
    samples = load_all_samples(DATA_DIR)
    if not samples:
        log.error("No samples found.")
        sys.exit(1)
    log.info(f"Loaded {len(samples)} samples: {sorted(samples)}")

    sex_map = assign_sex(sorted(samples))
    log.info("Sex assignment (sorted IDs, first 5 = F, last 5 = M):")
    for sid, sex in sex_map.items():
        log.info(f"  {sid} → {sex}")

    log.info("Building results table …")
    df = build_results_table(samples)

    out_tsv = RESULTS_DIR / "sex_bias_results.tsv"
    df.to_csv(out_tsv, sep="\t", index=False, float_format="%.4f")
    log.info(f"Results → {out_tsv}")
    log.info("\n" + df[[
        "sample_id", "sex", "auto_het", "x_het",
        "auto_c1", "x_c1", "x_minus_auto_c1",
    ]].to_string(index=False))

    log.info("Generating figure …")
    plot_figure(df)
    log.info("Done.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""HGP1K-anchored sex-biased ADMIXTURE.

Joint ADMIXTURE of the study cohort + a frozen 1000 Genomes reference
(AFR/EUR/SAS, baked as PLINK BED) at K = 3,4,5, run on:
  * the combined genome   (the headline ancestry result)
  * autosomes only        }  compared, per labelled component, to detect
  * X chromosome only      }  sex-biased admixture (component over-represented
                              on X vs autosomes -> female-biased contribution).

Components are anonymous out of ADMIXTURE, so each is LABELLED by which
reference superpopulation loads highest on it. Labelling the autosome and X
runs against the SAME reference makes their components directly comparable.

X handling (see flow README / Kathy Liu): males are hemizygous, so the study
`sex` facet is materialised and applied (`plink2 --update-sex`) before the X
run; PLINK then encodes male non-PAR X as haploid. All samples are kept; a
per-sex breakdown of the X estimates is emitted.

Pipeline (all validated against real 1KGP + bvs output):
  study_raw (bvs cohort-bed, mixed rsID/chr:pos IDs, monomorphic rows)
    -> plink1.9 round-trip + fix monomorphic ALT==REF -> A1='0'
    -> plink2 re-key chr:pos, --update-sex, --rm-dup
    -> split autosomes / X
    -> per compartment: intersect with reference (chr:pos), drop strand-ambiguous
       (A/T,C/G), plink1.9 --bmerge with .missnp retry
    -> QC (geno/mind/maf/hwe) + LD prune (indep-pairwise) per compartment
    -> combined = concat(auto_pruned, x_pruned)
    -> ADMIXTURE K on combined / auto / x
    -> label components by reference superpop; X-vs-auto comparison; plots.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --- config / defaults ------------------------------------------------------
GENO = float(os.environ.get("BV_GENO", "0.05"))
MIND = float(os.environ.get("BV_MIND", "0.10"))
MAF = float(os.environ.get("BV_MAF", "0.01"))
HWE = float(os.environ.get("BV_HWE", "1e-6"))
LD_WINDOW = int(os.environ.get("BV_LD_WINDOW", "200"))
LD_STEP = int(os.environ.get("BV_LD_STEP", "50"))
LD_R2 = float(os.environ.get("BV_LD_R2", "0.2"))
SEED = int(os.environ.get("BV_ADMIXTURE_SEED", "1729"))
BUILD = os.environ.get("BV_BUILD", "hg38")

SEX_TO_PLINK = {"m": "1", "male": "1", "1": "1", "f": "2", "female": "2", "2": "2"}


def run(cmd, cwd=None, log=None, check=True):
    cmd = [str(c) for c in cmd]
    print(f"[hgp1k] $ {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    if log:
        Path(log).write_text(proc.stdout, encoding="utf-8")
    if check and proc.returncode != 0:
        sys.stderr.write(proc.stdout + "\n")
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc


def tool(name):
    p = shutil.which(name)
    if not p:
        raise SystemExit(f"ERROR: required tool not on PATH: {name}")
    return p


def n_variants(prefix: Path) -> int:
    bim = Path(f"{prefix}.bim")
    return sum(1 for _ in bim.open()) if bim.exists() else 0


def read_fam_ids(prefix: Path) -> list[str]:
    return [ln.split()[1] for ln in Path(f"{prefix}.fam").open() if ln.strip()]


# --- study prep -------------------------------------------------------------
def load_sex_mapping(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path or not Path(path).is_file():
        return out
    for i, line in enumerate(Path(path).open()):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        pid, raw = parts[0].strip(), parts[1].strip()
        if i == 0 and raw.lower() in ("sex", "gender"):
            continue
        sx = SEX_TO_PLINK.get(raw.lower())
        if sx:
            out[pid] = sx
    return out


def write_sex_file(sample_ids, sex_map, path: Path):
    with path.open("w") as h:
        h.write("#IID\tSEX\n")
        for sid in sample_ids:
            h.write(f"{sid}\t{sex_map.get(sid, '0')}\n")


def prep_study(study_bed: Path, sex_file: Path, work: Path, plink, plink2, threads):
    """bvs study_raw -> chr:pos-keyed, sex-coded study_chrpos."""
    clean = work / "study_clean"
    run([plink, "--bfile", study_bed, "--make-bed", "--allow-extra-chr",
         "--out", clean], log=work / "study_clean.log")
    # fix monomorphic ALT==REF (col5==col6) -> A1='0' so plink2 will load it
    bim = Path(f"{clean}.bim")
    lines = []
    for ln in bim.open():
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 6 and f[4] == f[5]:
            f[4] = "0"
        lines.append("\t".join(f))
    bim.write_text("\n".join(lines) + "\n")

    chrpos = work / "study_chrpos"
    run([plink2, "--bfile", clean, "--set-all-var-ids", "@:#",
         "--rm-dup", "exclude-all", "--update-sex", sex_file,
         "--make-bed", "--out", chrpos, "--threads", threads, "--memory", "6000"],
        log=work / "study_chrpos.log")
    return chrpos


def split_compartment(study_chrpos: Path, sex_file: Path, region: str, work: Path, plink2, threads):
    out = work / f"study_{region}"
    args = [plink2, "--bfile", study_chrpos, "--make-bed", "--out", out,
            "--threads", threads, "--memory", "6000"]
    if region == "auto":
        args += ["--chr", "1-22"]
    else:
        args += ["--chr", "X", "--update-sex", sex_file, "--split-par", BUILD]
    run(args, log=work / f"study_{region}.log", check=False)
    return out if Path(f"{out}.bed").exists() else None


# --- merge study + reference (validated recipe) -----------------------------
def merge_with_reference(study_pref: Path, ref_pref: Path, region: str,
                         work: Path, plink, plink2, threads):
    sids = set(ln.split()[1] for ln in Path(f"{study_pref}.bim").open())
    rids = set(ln.split()[1] for ln in Path(f"{ref_pref}.bim").open())
    common = sorted(sids & rids)
    cfile = work / f"common_{region}.txt"
    cfile.write_text("\n".join(common) + "\n")
    print(f"[hgp1k] {region}: study {len(sids)} ∩ ref {len(rids)} = {len(common)} common", flush=True)
    if len(common) < 100:
        raise RuntimeError(f"{region}: only {len(common)} common variants — check build/IDs")

    s_c = work / f"s_{region}_common"
    r_c = work / f"r_{region}_common"
    run([plink2, "--bfile", study_pref, "--extract", cfile, "--make-bed", "--out", s_c,
         "--threads", threads, "--memory", "6000"], log=work / f"s_{region}_common.log")
    run([plink2, "--bfile", ref_pref, "--extract", cfile, "--make-bed", "--out", r_c,
         "--threads", threads, "--memory", "6000"], log=work / f"r_{region}_common.log")

    # drop strand-ambiguous (A/T, C/G) from study, apply to both
    ambi = work / f"ambi_{region}.txt"
    amb = []
    for ln in Path(f"{s_c}.bim").open():
        f = ln.split()
        a, b = f[4], f[5]
        if {a, b} in ({"A", "T"}, {"C", "G"}):
            amb.append(f[1])
    ambi.write_text("\n".join(amb) + "\n")
    s_na, r_na = s_c, r_c
    if amb:
        s_na = work / f"s_{region}_na"
        r_na = work / f"r_{region}_na"
        run([plink2, "--bfile", s_c, "--exclude", ambi, "--make-bed", "--out", s_na,
             "--threads", threads, "--memory", "6000"], log=work / f"s_{region}_na.log")
        run([plink2, "--bfile", r_c, "--exclude", ambi, "--make-bed", "--out", r_na,
             "--threads", threads, "--memory", "6000"], log=work / f"r_{region}_na.log")

    merged = work / f"merged_{region}"
    p = run([plink, "--bfile", s_na, "--bmerge", r_na, "--make-bed", "--allow-no-sex",
             "--out", merged, "--threads", threads], log=work / f"merge_{region}.log", check=False)
    missnp = work / f"merged_{region}-merge.missnp"
    if p.returncode != 0:
        if not missnp.exists():
            sys.stderr.write(Path(work / f"merge_{region}.log").read_text())
            raise RuntimeError(f"{region}: merge failed without a .missnp file")
        s2 = work / f"s_{region}_na2"
        r2 = work / f"r_{region}_na2"
        run([plink, "--bfile", s_na, "--exclude", missnp, "--make-bed", "--allow-no-sex",
             "--out", s2, "--threads", threads], log=work / f"s_{region}_na2.log")
        run([plink, "--bfile", r_na, "--exclude", missnp, "--make-bed", "--allow-no-sex",
             "--out", r2, "--threads", threads], log=work / f"r_{region}_na2.log")
        run([plink, "--bfile", s2, "--bmerge", r2, "--make-bed", "--allow-no-sex",
             "--out", merged, "--threads", threads], log=work / f"merge_{region}2.log")
    return merged


def qc_and_prune(merged: Path, region: str, work: Path, plink2, threads):
    qc = work / f"qc_{region}"
    run([plink2, "--bfile", merged, "--geno", GENO, "--mind", MIND, "--maf", MAF,
         "--hwe", HWE, "--make-bed", "--out", qc, "--allow-no-sex",
         "--threads", threads, "--memory", "6000"], log=work / f"qc_{region}.log")
    prune = work / f"prune_{region}"
    run([plink2, "--bfile", qc, "--indep-pairwise", LD_WINDOW, LD_STEP, LD_R2,
         "--out", prune, "--allow-no-sex", "--threads", threads, "--memory", "6000"],
        log=work / f"prune_{region}.log")
    pruned = work / f"pruned_{region}"
    run([plink2, "--bfile", qc, "--extract", f"{prune}.prune.in", "--make-bed",
         "--out", pruned, "--allow-no-sex", "--threads", threads, "--memory", "6000"],
        log=work / f"pruned_{region}.log")
    return pruned


def concat_compartments(prefixes, out: Path, work: Path, plink2, threads):
    mlist = work / "combine_list.txt"
    mlist.write_text("\n".join(str(p) for p in prefixes[1:]) + "\n")
    run([plink2, "--bfile", prefixes[0], "--pmerge-list", mlist, "bfile",
         "--make-bed", "--out", out, "--threads", threads, "--memory", "8000"],
        log=work / "combine.log")
    return out


# --- ADMIXTURE --------------------------------------------------------------
def admixture_ready(pruned: Path, work: Path, tag: str) -> Path:
    """Copy bed/fam, rewrite .bim chromosome col to '1' (ADMIXTURE ignores
    chrom/pos but refuses non-integer codes like 'X')."""
    dst = work / f"adx_{tag}"
    shutil.copy(f"{pruned}.bed", f"{dst}.bed")
    shutil.copy(f"{pruned}.fam", f"{dst}.fam")
    out = []
    for ln in Path(f"{pruned}.bim").open():
        f = ln.rstrip("\n").split("\t")
        f[0] = "1"
        out.append("\t".join(f))
    Path(f"{dst}.bim").write_text("\n".join(out) + "\n")
    return dst


def run_admixture(pruned: Path, tag: str, k: int, work: Path, out_dir: Path,
                  admixture, threads) -> tuple[Path, list[str]]:
    adx = admixture_ready(pruned, work, f"{tag}_K{k}")
    seed = SEED + k
    run([admixture, f"-s{seed}", f"{adx}.bed", k, f"-j{threads}"],
        cwd=work, log=out_dir / f"admixture_{tag}_K{k}.log")
    q_src = work / f"{adx.name}.{k}.Q"
    p_src = work / f"{adx.name}.{k}.P"
    q_dst = out_dir / f"admixture_{tag}_K{k}.Q"
    p_dst = out_dir / f"admixture_{tag}_K{k}.P"
    shutil.copy(q_src, q_dst)
    if p_src.exists():
        shutil.copy(p_src, p_dst)
    return q_dst, read_fam_ids(pruned)


# --- labelling + comparison -------------------------------------------------
def label_and_frame(q_path: Path, sample_ids, k, labels_df, sex_df) -> tuple[pd.DataFrame, dict]:
    q = np.loadtxt(q_path, dtype=float)
    if q.ndim == 1:
        q = q.reshape(1, -1)
    comp_cols = [f"comp{i+1}" for i in range(k)]
    df = pd.DataFrame(q, columns=comp_cols)
    df.insert(0, "sample_id", sample_ids)
    df = df.merge(labels_df, on="sample_id", how="left")  # adds superpopulation (ref only)
    df = df.merge(sex_df, on="sample_id", how="left")      # adds sex
    df["group"] = np.where(df["superpopulation"].notna(), "reference", "study")

    # component -> ancestry label: superpop with highest mean Q on that component
    ref = df[df["group"] == "reference"]
    comp_label = {}
    if not ref.empty:
        for c in comp_cols:
            means = ref.groupby("superpopulation")[c].mean()
            comp_label[c] = means.idxmax() if not means.empty else c
    else:
        comp_label = {c: c for c in comp_cols}
    # dedup duplicate labels (comp -> AFR, AFR_2, ...)
    seen: dict[str, int] = {}
    rename = {}
    for c in comp_cols:
        lab = comp_label[c]
        seen[lab] = seen.get(lab, 0) + 1
        rename[c] = lab if seen[lab] == 1 else f"{lab}_{seen[lab]}"
    df = df.rename(columns=rename)
    return df, rename


def compare_x_vs_auto(auto_df, x_df, k, auto_map, x_map, group="study") -> pd.DataFrame:
    """Per labelled ancestry, mean proportion on autosomes vs X for one group
    (group='study' is the result; group='reference' is the unadmixed-founder
    negative control — its deltas should sit at ~0)."""
    rows = []
    auto_labels = list(auto_map.values())
    x_labels = list(x_map.values())
    auto_study = auto_df[auto_df["group"] == group]
    x_study = x_df[x_df["group"] == group]
    for lab in sorted(set(auto_labels) & set(x_labels)):
        a = float(auto_study[lab].mean()) if lab in auto_study else np.nan
        xall = float(x_study[lab].mean()) if lab in x_study else np.nan
        xf = float(x_study.loc[x_study["sex"] == "2", lab].mean()) if lab in x_study else np.nan
        xm = float(x_study.loc[x_study["sex"] == "1", lab].mean()) if lab in x_study else np.nan
        rows.append({
            "K": k, "ancestry": lab,
            "mean_auto": round(a, 4),
            "mean_x": round(xall, 4),
            "delta_x_minus_auto": round(xall - a, 4),
            "mean_x_female": round(xf, 4),
            "mean_x_male_haploid": round(xm, 4),
        })
    return pd.DataFrame(rows)


def plot_x_vs_auto(comp: pd.DataFrame, out_png: Path):
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
        w = 0.38
        ax.bar(x - w / 2, sub["mean_auto"], w, label="Autosomes", color="#aaaaaa")
        ax.bar(x + w / 2, sub["mean_x"], w, label="X chromosome", color="#B2182B")
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


def build_per_sample(auto_df, x_df, k, auto_map, x_map, group="study") -> tuple[pd.DataFrame, list]:
    """Per individual (of the given group): autosomal vs X proportion (+ delta)
    for each labelled ancestry — the row-level data behind the per-individual
    panels. group='reference' gives the negative-control figure."""
    labels = sorted(set(auto_map.values()) & set(x_map.values()))
    a = auto_df[auto_df["group"] == group].copy()
    x = x_df[x_df["group"] == group].copy()
    acols = ["sample_id", "sex"] + [l for l in labels if l in a.columns]
    xcols = ["sample_id"] + [l for l in labels if l in x.columns]
    m = a[acols].merge(x[xcols], on="sample_id", suffixes=("_auto", "_x"))
    for l in labels:
        ca, cx = f"{l}_auto", f"{l}_x"
        if ca in m.columns and cx in m.columns:
            m[f"{l}_delta"] = (m[cx] - m[ca]).round(6)
    m.insert(1, "K", k)
    m["sex_label"] = m["sex"].map({"1": "M", "2": "F"}).fillna("?")
    return m, labels


def plot_figure(per_sample: pd.DataFrame, cohort: pd.DataFrame, k: int,
                labels: list, out_png: Path, focal_ancestry: str | None = None):
    """4-panel figure (mirrors the old figure4): per-individual auto-vs-X
    scatter + lollipop for the most sex-biased ancestry, stacked component
    bars by sex, and the cohort mean bar. PNG + PDF."""
    if per_sample.empty or not labels:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    COL_F, COL_M = "#B2182B", "#2166AC"
    palette = ["#1A9641", "#4575B4", "#F46D43", "#984EA3", "#FF7F00"]
    sex_color = {"F": COL_F, "M": COL_M, "?": "#888888"}

    # Default focal ancestry = largest |mean delta| this K. Reference-control
    # plots can pass the study focal ancestry so the control directly tests the
    # same signal instead of highlighting a tiny unrelated residual.
    csub = cohort[cohort["K"] == k]
    fallback_focal = (csub.reindex(csub["delta_x_minus_auto"].abs().sort_values(ascending=False).index)
                      ["ancestry"].iloc[0]) if not csub.empty else labels[0]
    focal = focal_ancestry if (
        focal_ancestry
        and f"{focal_ancestry}_auto" in per_sample.columns
        and f"{focal_ancestry}_x" in per_sample.columns
    ) else fallback_focal
    fa, fx, fd = f"{focal}_auto", f"{focal}_x", f"{focal}_delta"

    ps = per_sample.sort_values([f"{focal}_auto"]).reset_index(drop=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(f"HGP1K sex-biased admixture — K={k} (focal ancestry: {focal})", fontsize=12)

    # a) scatter auto vs X for the focal ancestry, by sex
    ax = axes[0][0]
    lim = [0, 1]
    ax.plot(lim, lim, ls="--", lw=0.8, color="#aaa")
    for _, r in ps.iterrows():
        ax.scatter(r.get(fa, np.nan), r.get(fx, np.nan), s=45,
                   color=sex_color.get(r["sex_label"], "#888"), edgecolors="k", linewidths=0.4)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(f"{focal} — autosomes"); ax.set_ylabel(f"{focal} — X")
    ax.set_title("Autosomal vs X ancestry (per individual)")
    ax.legend(handles=[mpatches.Patch(color=COL_F, label="Female"),
                       mpatches.Patch(color=COL_M, label="Male")], fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    # b) lollipop of (X - auto) for the focal ancestry, per individual
    ax = axes[0][1]
    ax.axvline(0, ls="--", lw=0.8, color="#aaa")
    for i, (_, r) in enumerate(ps.iterrows()):
        d = r.get(fd, 0.0)
        c = sex_color.get(r["sex_label"], "#888")
        ax.hlines(i, 0, d, color=c, lw=1.2)
        ax.scatter(d, i, s=35, color=c, edgecolors="k", linewidths=0.4)
    ax.set_yticks(range(len(ps)))
    ax.set_yticklabels(ps["sample_id"], fontsize=6)
    ax.set_xlabel(f"X − autosomal {focal} (Δ proportion)")
    ax.set_title("Sex-biased signal per individual")
    ax.spines[["top", "right"]].set_visible(False)

    # c) stacked component bars, autosomes vs X, per individual (F | M ordered)
    ax = axes[1][0]
    psx = per_sample.sort_values(["sex_label", f"{focal}_auto"]).reset_index(drop=True)
    n = len(psx); xs = np.arange(n); bw = 0.38
    for off, suf, tag in [(-bw / 2, "_auto", "A"), (bw / 2, "_x", "X")]:
        bottom = np.zeros(n)
        for li, l in enumerate(labels):
            col = f"{l}{suf}"
            if col not in psx.columns:
                continue
            vals = psx[col].fillna(0).to_numpy()
            ax.bar(xs + off, vals, bw, bottom=bottom, color=palette[li % len(palette)],
                   edgecolor="white", linewidth=0.2,
                   label=l if suf == "_auto" else "_nolegend_")
            bottom += vals
    ax.set_xticks(xs); ax.set_xticklabels(psx["sample_id"], rotation=60, ha="right", fontsize=5)
    ax.set_ylabel("Component proportion"); ax.set_ylim(0, 1.05)
    ax.set_title("Ancestry: autosomes (A) vs X — per individual")
    ax.legend(fontsize=7, ncol=len(labels), frameon=False, loc="lower center")
    ax.spines[["top", "right"]].set_visible(False)

    # d) cohort mean auto vs X per ancestry (this K)
    ax = axes[1][1]
    cx = np.arange(len(csub)); bw = 0.38
    ax.bar(cx - bw / 2, csub["mean_auto"], bw, color="#aaaaaa", label="Autosomes")
    ax.bar(cx + bw / 2, csub["mean_x"], bw, color=COL_F, label="X chromosome")
    ax.set_xticks(cx); ax.set_xticklabels(csub["ancestry"], fontsize=8)
    ax.set_ylabel("Mean proportion (study)")
    ax.set_title("Cohort mean: autosomes vs X")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return focal


# --- main -------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--study-bed", type=Path, required=True, help="bvs cohort-bed prefix (all chr)")
    ap.add_argument("--sex-mapping", type=Path, required=True, help="participant_id<TAB>sex")
    ap.add_argument("--reference-dir", type=Path,
                    default=Path(os.environ.get("HGP1K_ADMIXTURE_REF",
                                                "/opt/biovault/reference/hgp1k_admixture")))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, default=None)
    ap.add_argument("--k-values", default=os.environ.get("BV_ADMIXTURE_K", "3,4,5"))
    ap.add_argument("--threads", default=os.environ.get("BV_THREADS", "8"))
    args = ap.parse_args()

    plink, plink2, admixture = tool("plink"), tool("plink2"), tool("admixture")
    ks = [int(x) for x in str(args.k_values).split(",") if x.strip()]
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    work = args.work_dir or (out / "_work")
    work.mkdir(parents=True, exist_ok=True)

    ref_auto = args.reference_dir / "reference_auto"
    ref_x = args.reference_dir / "reference_x"
    labels_df = pd.read_csv(args.reference_dir / "reference_labels.tsv", sep="\t",
                            dtype=str).rename(columns={"sample_id": "sample_id"})
    labels_df.columns = ["sample_id", "superpopulation"][:len(labels_df.columns)]

    # study sex (facet) -> plink sex file
    study_ids = read_fam_ids(args.study_bed)
    sex_map = load_sex_mapping(args.sex_mapping)
    sex_file = work / "study_sex.tsv"
    write_sex_file(study_ids, sex_map, sex_file)
    print(f"[hgp1k] study sex: {sum(v=='1' for v in sex_map.values())} male, "
          f"{sum(v=='2' for v in sex_map.values())} female, of {len(study_ids)} samples", flush=True)

    # 1) prep study, 2) split, 3) merge study+reference per compartment
    study_chrpos = prep_study(args.study_bed, sex_file, work, plink, plink2, args.threads)
    merged = {}
    study_n = {}
    for region, ref_pref in (("auto", ref_auto), ("x", ref_x)):
        s = split_compartment(study_chrpos, sex_file, region, work, plink2, args.threads)
        if s is None:
            print(f"[hgp1k] WARNING: no study {region} variants — skipping {region}", flush=True)
            continue
        merged[region] = merge_with_reference(s, ref_pref, region, work, plink, plink2, args.threads)
        study_n[region] = n_variants(s)

    if not merged:
        raise SystemExit("ERROR: no compartments produced (no autosome or X data)")

    # 4) combined = concat of the PRE-QC merged sets (identical sample sets, so
    #    --pmerge-list concatenates cleanly); QC+prune is then applied once each.
    to_qc = dict(merged)
    if "auto" in merged and "x" in merged:
        to_qc["combined"] = concat_compartments([merged["auto"], merged["x"]],
                                                work / "merged_combined", work, plink2, args.threads)

    # 5) QC + LD-prune each set
    runs = {}
    qc_lines = []
    for tag, m in to_qc.items():
        pr = qc_and_prune(m, tag, work, plink2, args.threads)
        runs[tag] = pr
        sn = study_n.get(tag, "—")
        qc_lines.append(f"{tag}: study={sn} merged={n_variants(m)} "
                        f"pruned={n_variants(pr)} samples={len(read_fam_ids(pr))}")

    pre_admixture_qc = "\n".join([
        "=== hgp1k sex-biased ADMIXTURE pre-run QC ===",
        f"K values: {ks}",
        f"reference: {args.reference_dir}",
        f"geno={GENO} mind={MIND} maf={MAF} hwe={HWE} "
        f"ld=({LD_WINDOW},{LD_STEP},{LD_R2})",
        *qc_lines,
        "NOTE: these pruned variant/sample counts are the exact PLINK BED inputs "
        "passed to ADMIXTURE.",
    ])
    print("[hgp1k] pre-ADMIXTURE QC summary:\n" + pre_admixture_qc, flush=True)

    # Combined sex map: study (from the facet) + reference (from the baked
    # frozen list, coded 1/2) — so reference rows also carry sex for the
    # negative-control figure.
    ref_sex = {}
    _rs_path = args.reference_dir / "reference_samples.tsv"
    if _rs_path.exists():
        _rs = pd.read_csv(_rs_path, sep="\t", dtype=str)
        if "sex" in _rs.columns:
            ref_sex = dict(zip(_rs["sample_id"], _rs["sex"].astype(str)))
    combined_sex = {**ref_sex, **sex_map}   # study facet wins on any overlap
    sex_df = pd.DataFrame({"sample_id": list(combined_sex), "sex": list(combined_sex.values())})

    # 6-7) ADMIXTURE + labelling for every run/K
    all_labeled = []
    comp_label_rows = []
    frames: dict[tuple[str, int], tuple[pd.DataFrame, dict]] = {}
    for tag, pref in runs.items():
        for k in ks:
            q_path, fam_ids = run_admixture(pref, tag, k, work, out, admixture, args.threads)
            df, rename = label_and_frame(q_path, fam_ids, k, labels_df, sex_df)
            frames[(tag, k)] = (df, rename)
            lab_out = out / f"admixture_{tag}_K{k}_labeled_Q.tsv"
            df.to_csv(lab_out, sep="\t", index=False, float_format="%.6g")
            for comp, lab in rename.items():
                comp_label_rows.append({"run": tag, "K": k, "component": comp, "ancestry_label": lab})
            all_labeled.append(lab_out)

    pd.DataFrame(comp_label_rows).to_csv(out / "component_labels.tsv", sep="\t", index=False)

    # 8) X-vs-autosome comparison (the sex-bias signal) — cohort table, the
    #    per-individual table, and a 4-panel figure (PNG + PDF) per K. Every
    #    plot's underlying data is dumped to TSV so figures are reproducible.
    comp_frames = []
    per_sample_by_k = {}
    for k in ks:
        if ("auto", k) in frames and ("x", k) in frames:
            a_df, a_map = frames[("auto", k)]
            x_df, x_map = frames[("x", k)]
            comp_frames.append(compare_x_vs_auto(a_df, x_df, k, a_map, x_map))
            ps, labels = build_per_sample(a_df, x_df, k, a_map, x_map)
            per_sample_by_k[k] = (ps, labels)
            ps.to_csv(out / f"sex_bias_per_sample_K{k}.tsv", sep="\t", index=False,
                      float_format="%.6g")
    sex_bias = pd.concat(comp_frames, ignore_index=True) if comp_frames else pd.DataFrame()
    sex_bias.to_csv(out / "sex_bias_x_vs_auto.tsv", sep="\t", index=False)
    if not sex_bias.empty:
        # cohort summary across all K (PNG + PDF)
        plot_x_vs_auto(sex_bias, out / "sex_bias_x_vs_auto.png")
        # per-K 4-panel figure (PNG + PDF), data behind it already in the TSVs above
        study_focal_by_k = {}
        for k, (ps, labels) in per_sample_by_k.items():
            study_focal_by_k[k] = plot_figure(
                ps, sex_bias, k, labels, out / f"figure_sex_biased_admixture_K{k}.png")
        print("[hgp1k] X-vs-autosome comparison:\n" + sex_bias.to_string(index=False), flush=True)
    else:
        study_focal_by_k = {}

    # 8b) reference negative control — the same figure/data for the 900 baked
    #     reference founders. They are unadmixed, so deltas sit at ~0: a built-in
    #     check that the pipeline manufactures no false sex-bias signal.
    ref_comp_frames = []
    ref_per_sample_by_k = {}
    for k in ks:
        if ("auto", k) in frames and ("x", k) in frames:
            a_df, a_map = frames[("auto", k)]
            x_df, x_map = frames[("x", k)]
            rc = compare_x_vs_auto(a_df, x_df, k, a_map, x_map, group="reference")
            if rc.empty:
                continue
            ref_comp_frames.append(rc)
            rps, rlabels = build_per_sample(a_df, x_df, k, a_map, x_map, group="reference")
            ref_per_sample_by_k[k] = (rps, rlabels)
            rps.to_csv(out / f"sex_bias_per_sample_K{k}_reference.tsv", sep="\t",
                       index=False, float_format="%.6g")
    ref_sex_bias = pd.concat(ref_comp_frames, ignore_index=True) if ref_comp_frames else pd.DataFrame()
    if not ref_sex_bias.empty:
        ref_sex_bias.to_csv(out / "sex_bias_x_vs_auto_reference.tsv", sep="\t", index=False)
        for k, (rps, rlabels) in ref_per_sample_by_k.items():
            plot_figure(rps, ref_sex_bias, k, rlabels,
                        out / f"figure_sex_biased_admixture_K{k}_reference.png",
                        focal_ancestry=study_focal_by_k.get(k))
        print("[hgp1k] reference negative-control (deltas should be ~0):\n"
              + ref_sex_bias.to_string(index=False), flush=True)

    (out / "qc_report.txt").write_text(
        "=== hgp1k sex-biased admixture QC ===\n"
        f"K values: {ks}\n"
        f"reference: {args.reference_dir}\n"
        f"geno={GENO} mind={MIND} maf={MAF} hwe={HWE} "
        f"ld=({LD_WINDOW},{LD_STEP},{LD_R2})\n"
        + "\n".join(qc_lines) + "\n"
        "NOTE: male X is haploid (dosage 0/1); X estimates include a per-sex "
        "breakdown (mean_x_female / mean_x_male_haploid) for that reason.\n")

    if os.environ.get("BV_KEEP_WORK", "0") != "1":
        shutil.rmtree(work, ignore_errors=True)
    print("[hgp1k] DONE", flush=True)


if __name__ == "__main__":
    main()

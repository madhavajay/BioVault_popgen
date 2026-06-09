#!/usr/bin/env python3
"""Generate a synthetic sex-biased admixture test cohort for the hgp1k flow.

Plants a *known* ancestry difference between autosomes and X so the pipeline's
output can be checked against ground truth:

  * autosomal AFR fraction  = --afr-auto   (e.g. 0.50)
  * X-chromosome AFR fraction = --afr-x     (e.g. 0.85)   -> female-biased AFR

Genotypes are drawn from the REAL 1KGP AFR/EUR allele frequencies (computed
from the baked reference BED), so ADMIXTURE — anchored to the same reference —
recovers the planted proportions. The remaining ancestry is EUR. Males are
hemizygous on X (one allele, written homozygous) so PLINK codes them haploid.

Output: per-participant DDNA-format genotype .txt under <out>/<id>/, plus
<out>/samplesheet.csv (participant_id, genotype_file, sex). Point the flow at
the samplesheet; expect AFR mean_x - mean_auto approx (afr_x - afr_auto).

By default this writes a compact file containing only the simulated HGP1K
overlap SNPs. Pass --template-ddna to instead write a full-size DDNA-style file:
every row from the template is copied, and the simulated sex-biased SNPs are
overlaid at matching chromosome/position rows.

Run inside an image with plink2 + numpy/pandas (the biovault-admixture image):
  python generate_sex_biased_testset.py --reference-dir .docker/reference/hgp1k_admixture \
    --locus-map tools/locus_map.tsv --out-dir <out> --n-samples 100
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


WORK = {}


def norm_chrom(value) -> str:
    c = str(value).strip()
    if c in ("23", "X", "x"):
        return "X"
    if c in ("24", "Y", "y"):
        return "Y"
    if c in ("25", "MT", "M", "mt", "m"):
        return "MT"
    return c.removeprefix("chr")


def load_template_ddna(path: Path):
    """Load a full DDNA file as reusable lines plus chr:pos lookup."""
    header: list[str] = []
    row_lines: list[str] = []
    row_fields: list[list[str]] = []
    index: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            raw = line.rstrip("\n\r")
            if not raw:
                continue
            if raw.startswith("#"):
                header.append(raw)
                continue
            fields = raw.split("\t")
            if len(fields) < 4:
                continue
            idx = len(row_lines)
            chrom = norm_chrom(fields[1])
            key = f"{chrom}:{fields[2]}"
            index.setdefault(key, idx)
            row_lines.append(raw)
            row_fields.append(fields)
    if not row_lines:
        raise ValueError(f"{path}: no DDNA genotype rows found")
    return header, row_lines, row_fields, index


def plink_freq(bed_prefix: Path, keep_ids: list[str], out: Path, tmp: Path) -> pd.DataFrame:
    """ALT-allele frequency per variant for a subset of reference samples."""
    keep = tmp / f"{out.name}.keep"
    with keep.open("w") as h:
        for s in keep_ids:
            h.write(f"0\t{s}\n")          # plink2 VCF-import FIDs are '0'
    subprocess.run(
        ["plink2", "--bfile", str(bed_prefix), "--keep", str(keep),
         "--freq", "--out", str(out), "--allow-no-sex"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    df = pd.read_csv(f"{out}.afreq", sep="\t")
    # columns: #CHROM ID REF ALT ALT_FREQS OBS_CT
    df = df.rename(columns={"#CHROM": "chrom", "ID": "id", "REF": "ref",
                            "ALT": "alt", "ALT_FREQS": "af"})
    return df[["id", "chrom", "ref", "alt", "af"]]


def build_panel(ref_dir: Path, locus_map: Path, n_auto: int, tmp: Path,
                rng: np.random.Generator, template_index: dict[str, int] | None = None):
    labels = pd.read_csv(ref_dir / "reference_labels.tsv", sep="\t")
    labels.columns = ["sample_id", "superpop"][:labels.shape[1]]
    afr = labels.loc[labels.superpop == "AFR", "sample_id"].tolist()
    eur = labels.loc[labels.superpop == "EUR", "sample_id"].tolist()

    panels = {}
    for region, bed in (("auto", ref_dir / "reference_auto"), ("x", ref_dir / "reference_x")):
        a = plink_freq(bed, afr, tmp / f"{region}_afr", tmp).rename(columns={"af": "afr_af"})
        e = plink_freq(bed, eur, tmp / f"{region}_eur", tmp)[["id", "af"]].rename(columns={"af": "eur_af"})
        m = a.merge(e, on="id")
        # informative SNPs only (some AFR/EUR difference, not fixed/zero)
        m = m[(m.afr_af > 0.0) | (m.eur_af > 0.0)]
        m = m[(m.afr_af < 1.0) | (m.eur_af < 1.0)]
        m["pos"] = m["id"].str.split(":").str[1].astype(int)
        m["template_key"] = m["chrom"].map(norm_chrom) + ":" + m["pos"].astype(str)
        if template_index is not None:
            m = m[m["template_key"].isin(template_index)].copy()
            m["template_idx"] = m["template_key"].map(template_index).astype(int)
        if region == "auto" and len(m) > n_auto:
            m = m.sample(n=n_auto, random_state=int(rng.integers(1 << 31))).copy()
        m = m.sort_values(["chrom", "pos"]).reset_index(drop=True)
        panels[region] = m
    # rsid lookup (chrom:pos -> rsid) from locus_map for the DDNA SNP-name col
    lm = pd.read_csv(locus_map, sep="\t", dtype={"chrom": str})
    lm["key"] = lm["chrom"].astype(str) + ":" + lm["pos"].astype(str)
    rs = dict(zip(lm["key"], lm["rsid"]))
    return panels, rs


def draw_region(af_mix: np.ndarray, n: int, haploid: bool, rng):
    """Draw genotype dosages (0/1/2 diploid, or 0/1 haploid) per variant."""
    if haploid:
        return rng.binomial(1, af_mix)        # one allele
    return rng.binomial(2, af_mix)            # two alleles


def geno_string(dosage: int, ref: str, alt: str, haploid: bool) -> str:
    """DDNA 2-char genotype call. ALT is the counted (A1) allele."""
    if haploid:                               # male X: write homozygous of the single allele
        return (alt + alt) if dosage == 1 else (ref + ref)
    return {0: ref + ref, 1: alt + ref, 2: alt + alt}[int(dosage)]


def write_ddna(path: Path, panel_auto, panel_x, g_auto, g_x, sex_is_male: bool, rs):
    rows = []
    for region, panel, g, hap in (("auto", panel_auto, g_auto, False),
                                  ("x", panel_x, g_x, sex_is_male)):
        ref = panel["ref"].to_numpy(); alt = panel["alt"].to_numpy()
        chrom = panel["chrom"].astype(str).to_numpy(); pos = panel["pos"].to_numpy()
        ids = panel["id"].to_numpy()
        for i in range(len(panel)):
            c = "X" if str(chrom[i]) in ("23", "X") else str(chrom[i])
            rsid = rs.get(f"{c}:{pos[i]}") or rs.get(ids[i]) or ids[i]
            gt = geno_string(g[i], ref[i], alt[i], hap)
            rows.append(f"{rsid}\t{c}\t{pos[i]}\t{gt}\t0.5000\t0.500\t0.0000")
    path.write_text("\n".join(rows) + "\n")


def write_full_ddna(path: Path, header: list[str], row_lines: list[str],
                    row_fields: list[list[str]], overlay: dict[int, str]):
    with path.open("w") as handle:
        for line in header:
            handle.write(line + "\n")
        for idx, line in enumerate(row_lines):
            gt = overlay.get(idx)
            if gt is None:
                handle.write(line + "\n")
                continue
            fields = row_fields[idx].copy()
            fields[3] = gt
            handle.write("\t".join(fields) + "\n")


def generate_sample(task: tuple[int, bool]) -> dict[str, str]:
    i, is_male = task
    pa = WORK["pa"]
    px = WORK["px"]
    af_auto = WORK["af_auto"]
    af_x = WORK["af_x"]
    out_dir = WORK["out_dir"]
    seed = WORK["seed"]
    template = WORK["template"]
    rs = WORK["rs"]

    rng = np.random.default_rng(seed + 1_000_003 * (i + 1))
    pid = f"SB{i:04d}"
    sex = "Male" if is_male else "Female"
    g_auto = draw_region(af_auto, len(pa), False, rng)
    g_x = draw_region(af_x, len(px), is_male, rng)
    d = out_dir / pid
    d.mkdir(parents=True, exist_ok=True)
    gfile = d / f"{pid}_X_X_GSAv3-DTC_GRCh38-synthetic.txt"
    if template is None:
        write_ddna(gfile, pa, px, g_auto, g_x, is_male, rs)
    else:
        header, row_lines, row_fields, _index = template
        overlay = {}
        for j, idx in enumerate(pa["template_idx"].to_numpy()):
            overlay[int(idx)] = geno_string(g_auto[j], pa["ref"].iat[j], pa["alt"].iat[j], False)
        for j, idx in enumerate(px["template_idx"].to_numpy()):
            overlay[int(idx)] = geno_string(g_x[j], px["ref"].iat[j], px["alt"].iat[j], is_male)
        write_full_ddna(gfile, header, row_lines, row_fields, overlay)
    return {"participant_id": pid, "genotype_file": str(gfile.resolve()), "sex": sex}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-dir", type=Path, required=True)
    ap.add_argument("--locus-map", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--n-samples", type=int, default=100)
    ap.add_argument("--n-auto-snps", type=int, default=60000)
    ap.add_argument("--template-ddna", type=Path, default=None,
                    help="Full DDNA file to use as the output template; sex-biased "
                         "SNPs are overlaid on matching chr:pos rows.")
    ap.add_argument("--afr-auto", type=float, default=0.50, help="autosomal AFR fraction")
    ap.add_argument("--afr-x", type=float, default=0.85, help="X AFR fraction (the planted bias)")
    ap.add_argument("--female-frac", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel sample writers. On Linux this uses fork so the "
                         "loaded full template is shared copy-on-write.")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    tmp = args.out_dir / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    template = None
    template_index = None
    if args.template_ddna:
        print(f"[gen] loading full DDNA template {args.template_ddna} ...", flush=True)
        template = load_template_ddna(args.template_ddna)
        template_index = template[3]
        print(f"[gen] template: {len(template[1])} genotype rows", flush=True)
    print(f"[gen] computing AFR/EUR frequencies from {args.reference_dir} ...", flush=True)
    panels, rs = build_panel(args.reference_dir, args.locus_map, args.n_auto_snps, tmp, rng,
                             template_index=template_index)
    pa, px = panels["auto"], panels["x"]
    if len(pa) == 0 or len(px) == 0:
        raise SystemExit(
            "template/reference intersection produced no usable SNPs "
            f"(auto={len(pa)} x={len(px)})"
        )
    print(f"[gen] panel: {len(pa)} autosomal + {len(px)} X SNPs", flush=True)

    af_auto = (args.afr_auto * pa["afr_af"] + (1 - args.afr_auto) * pa["eur_af"]).to_numpy()
    af_x = (args.afr_x * px["afr_af"] + (1 - args.afr_x) * px["eur_af"]).to_numpy()

    sexes = [rng.random() >= args.female_frac for _ in range(args.n_samples)]
    WORK.update({
        "pa": pa,
        "px": px,
        "af_auto": af_auto,
        "af_x": af_x,
        "out_dir": args.out_dir,
        "seed": args.seed,
        "template": template,
        "rs": rs,
    })

    tasks = list(enumerate(sexes))
    workers = max(1, args.workers)
    print(f"[gen] writing {args.n_samples} samples with workers={workers}", flush=True)
    if workers == 1:
        sheet = []
        for task in tasks:
            sheet.append(generate_sample(task))
            if len(sheet) % 25 == 0:
                print(f"[gen] {len(sheet)}/{args.n_samples} samples", flush=True)
    else:
        ctx = mp.get_context("fork")
        sheet = []
        with ctx.Pool(processes=workers) as pool:
            for row in pool.imap_unordered(generate_sample, tasks):
                sheet.append(row)
                if len(sheet) % 25 == 0:
                    print(f"[gen] {len(sheet)}/{args.n_samples} samples", flush=True)
        sheet.sort(key=lambda row: row["participant_id"])

    sheet_path = args.out_dir / "samplesheet.csv"
    with sheet_path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=["participant_id", "genotype_file", "sex"])
        w.writeheader(); w.writerows(sheet)

    nf = sum(1 for s in sheet if s["sex"] == "Female")
    (args.out_dir / "GROUND_TRUTH.txt").write_text(
        f"synthetic sex-biased admixture test set\n"
        f"n_samples={args.n_samples} (female={nf}, male={args.n_samples - nf})\n"
        f"autosomal AFR fraction = {args.afr_auto}\n"
        f"X AFR fraction         = {args.afr_x}\n"
        f"=> expected AFR delta_x_minus_auto ≈ {round(args.afr_x - args.afr_auto, 3)}\n"
        f"   (EUR delta ≈ {round(args.afr_auto - args.afr_x, 3)}; SAS ≈ 0)\n"
        f"auto SNPs={len(pa)} X SNPs={len(px)} seed={args.seed}\n"
        f"template_rows={len(template[1]) if template else 'compact'}\n")
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"[gen] DONE -> {sheet_path}\n[gen] ground truth: AFR X−auto ≈ {args.afr_x - args.afr_auto:+.2f}", flush=True)


if __name__ == "__main__":
    main()

"""
Build-time derivation of the AIMs gnomAD reference AF table.

Adapted from
04_population_level/aims_differential_snps/scripts/02_compute_gnomad_af_local.py
so it can run inside the biovault-popgen `tools` build stage (mirrors how
extract_loadings_variants.py derives loadings_variants.tsv from the HT). The
resulting small TSV is COPYed into the runtime image so the flow never touches
the ~80 GB HGDP+TGP VCFs at run time.

Usage:
  python derive_gnomad_aims_af.py <panel_tsv> <vcf_dir> <out_tsv>

  panel_tsv : panel_hgdp_tgp.tsv  (columns: sample, pop, super_pop, project)
  vcf_dir   : dir with gnomad.genomes.v3.1.2.hgdp_tgp.chr{1..22}.vcf.bgz
  out_tsv   : output: locus_key gnomAD_global gnomAD_AFR gnomAD_NFE gnomAD_SAS

Super-pop mapping (HGDP+TGP -> gnomAD-style label used downstream):
  AFR                         -> AFR
  EUR & pop != "FIN"          -> NFE
  CSA                         -> SAS
  all samples                 -> global
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CHROMS = list(range(1, 23))


def gt_to_alt(arr: np.ndarray) -> np.ndarray:
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    out[(arr == "0/0") | (arr == "0|0")] = 0.0
    out[(arr == "0/1") | (arr == "0|1") | (arr == "1/0") | (arr == "1|0")] = 1.0
    out[(arr == "1/1") | (arr == "1|1")] = 2.0
    return out


def parse_chrom_vcf(vcf: Path):
    samples = (
        subprocess.check_output(["bcftools", "query", "-l", str(vcf)], text=True)
        .strip()
        .split("\n")
    )
    raw = subprocess.check_output(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
        text=True,
    ).strip()
    if not raw:
        return None, None, samples
    rows = [line.split("\t") for line in raw.split("\n")]
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"] + samples)
    chrom_bare = df["CHROM"].str.replace("^chr", "", regex=True)
    df["locus_key"] = chrom_bare + "-" + df["POS"] + "-" + df["REF"] + "-" + df["ALT"]
    gts = df[samples].to_numpy()
    return df[["locus_key"]].copy(), gts, samples


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit("usage: derive_gnomad_aims_af.py <panel_tsv> <vcf_dir> <out_tsv>")
    panel_path, vcf_dir, out_path = (Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))

    panel = pd.read_csv(panel_path, sep="\t")
    groups = {
        "AFR": set(panel.loc[panel.super_pop == "AFR", "sample"]),
        "NFE": set(panel.loc[(panel.super_pop == "EUR") & (panel["pop"] != "FIN"), "sample"]),
        "SAS": set(panel.loc[panel.super_pop == "CSA", "sample"]),
        "global": set(panel["sample"]),
    }

    parts = []
    for chr_ in CHROMS:
        vcf = vcf_dir / f"gnomad.genomes.v3.1.2.hgdp_tgp.chr{chr_}.vcf.bgz"
        if not vcf.exists():
            print(f"[derive] chr{chr_}: VCF missing, skipping", flush=True)
            continue
        meta, gts, samples = parse_chrom_vcf(vcf)
        if meta is None:
            continue
        sample_idx = {s: i for i, s in enumerate(samples)}
        alt = gt_to_alt(gts)
        out_cols = {"locus_key": meta["locus_key"].values}
        for g, samp_set in groups.items():
            idxs = [sample_idx[s] for s in samp_set if s in sample_idx]
            if not idxs:
                af = np.full(len(meta), np.nan, dtype=np.float32)
            else:
                sub = alt[:, idxs]
                with np.errstate(invalid="ignore"):
                    ac = np.nansum(sub, axis=1)
                    an = (~np.isnan(sub)).sum(axis=1) * 2
                    af = np.where(an > 0, ac / an, np.nan).astype(np.float32)
            out_cols[f"gnomAD_{g}"] = af
        parts.append(pd.DataFrame(out_cols))
        print(f"[derive] chr{chr_}: {len(meta):,} variants", flush=True)

    final = (
        pd.concat(parts, ignore_index=True)
        .drop_duplicates(subset="locus_key", keep="first")
    )[["locus_key", "gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out_path, sep="\t", index=False, float_format="%.6f")
    print(f"[derive] wrote {len(final):,} loci -> {out_path}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Deterministically pick the disjoint biased-variant blocks used to inject
known-truth signal into the synthetic cohort, so each downstream analysis
has a verifiable positive control.

Blocks (all disjoint, drawn from the GSAv3 panel):

  island_structure[<island>]  ~400 autosomal SNPs, general panel.
      island K -> homozygous ALT; every other island -> homozygous REF.
      Drives: FST matrix, pca_qc clustering, AIMs differential SNPs.

  projection[<island>]         ~150 SNPs from GSA ∩ gnomAD loadings.
      Same island scheme. Drives: gnomad_projection_fast PC shift
      (only loadings variants move the projection).

  sex_block                    large chrX block (~10k SNPs), no autosomal.
      Mimics male X-hemizygosity: Male -> homozygous (low X het),
      Female -> heterozygous (normal diploid X het). Uniform across
      islands so invisible to island FST. Drives: sex_biased_admixture
      panel d (autosomal-vs-X heterozygosity by sex).

  singleton_island[<island>]   5 SNPs, extreme (in-island AF=1, else 0).
  singleton_sex                5 SNPs, extreme male/female.
      Crisp per-SNP binary asserts for the AIMs output.

Input:
  --illumina-ref   one bvs --format illumina file (gives rsid, [A/B]
                   allele pair, chr, pos for the whole panel).
  --loadings-tsv   loadings_variants.tsv (chrom, pos, ref, alt, id).

Output (JSON): {block_name: [{rsid, chromosome, position, reference,
                              alternate}, ...], ...} plus a "_meta" key.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

ISLANDS = ["Bahamas", "Barbados", "Bermuda", "BVI", "StLucia", "TT"]
AUTOSOMES = {str(c) for c in range(1, 23)}


def read_illumina_panel(path: Path):
    """Yield (rsid, chrom, pos, a1, a2) from a bvs Illumina GSGT file."""
    out = []
    with path.open() as f:
        in_data = False
        header_cols = None
        for line in f:
            line = line.rstrip("\n")
            if not in_data:
                if line.startswith("[Data]"):
                    in_data = True
                continue
            if header_cols is None:
                header_cols = line.split("\t")
                idx = {c: i for i, c in enumerate(header_cols)}
                i_snpname = idx.get("SNP Name")
                i_snp = idx.get("SNP")
                i_chr = idx.get("Chr")
                i_pos = idx.get("Position")
                continue
            parts = line.split("\t")
            if len(parts) <= max(i_snpname, i_snp, i_chr, i_pos):
                continue
            rsid = parts[i_snpname]
            snp = parts[i_snp]  # e.g. "[C/T]"
            chrom = parts[i_chr]
            pos = parts[i_pos]
            if not (snp.startswith("[") and "/" in snp and snp.endswith("]")):
                continue
            a1, a2 = snp[1:-1].split("/", 1)
            if a1 not in "ACGT" or a2 not in "ACGT" or a1 == a2:
                continue
            out.append((rsid, chrom, pos, a1, a2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--illumina-ref", required=True, type=Path)
    ap.add_argument("--loadings-tsv", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=100)
    ap.add_argument("--island-block", type=int, default=400)
    ap.add_argument("--proj-block", type=int, default=150)
    # Sex control mimics male X-hemizygosity: a large X-chromosome block
    # (Male->homozygous => low X het, Female->heterozygous => high X het).
    # Keep ~zero autosomal sex loci so autosomal het stays equal by sex
    # (clean panel-d contrast). Genome-mean het only moves if the X block
    # is a large fraction of the ~tens-of-thousands of X SNPs.
    ap.add_argument("--sex-auto", type=int, default=0)
    ap.add_argument("--sex-x", type=int, default=10000)
    ap.add_argument("--singletons", type=int, default=5)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    panel = read_illumina_panel(args.illumina_ref)
    if not panel:
        raise SystemExit(f"no variants parsed from {args.illumina_ref}")

    loadings = set()
    with args.loadings_tsv.open() as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                loadings.add((p[0], p[1]))

    autosomal = [v for v in panel if v[1] in AUTOSOMES]
    chrx = [v for v in panel if v[1] == "X"]
    in_load = [v for v in autosomal if (v[1], v[2]) in loadings]

    rng.shuffle(autosomal)
    rng.shuffle(chrx)
    rng.shuffle(in_load)

    used_rsids: set[str] = set()

    def take(pool, n, predicate=None):
        picked = []
        for v in pool:
            if len(picked) >= n:
                break
            if v[0] in used_rsids:
                continue
            if predicate and not predicate(v):
                continue
            picked.append(v)
            used_rsids.add(v[0])
        if len(picked) < n:
            raise SystemExit(
                f"pool exhausted: wanted {n}, got {len(picked)}")
        return picked

    def fmt(v):
        rsid, chrom, pos, a1, a2 = v
        return {"rsid": rsid, "chromosome": chrom, "position": int(pos),
                "reference": a1, "alternate": a2}

    blocks: dict = {}

    # Projection sub-blocks first (scarcest pool: GSA ∩ loadings).
    for isl in ISLANDS:
        blocks[f"projection[{isl}]"] = [fmt(v) for v in
                                        take(in_load, args.proj_block)]

    # Island structure blocks (general autosomal, excluding already used).
    for isl in ISLANDS:
        blocks[f"island_structure[{isl}]"] = [fmt(v) for v in
                                              take(autosomal, args.island_block)]

    # Sex block: autosomal + chrX.
    sex_auto = take(autosomal, args.sex_auto)
    sex_x = take(chrx, args.sex_x)
    blocks["sex_block"] = [fmt(v) for v in (sex_auto + sex_x)]

    # Singletons.
    for isl in ISLANDS:
        blocks[f"singleton_island[{isl}]"] = [fmt(v) for v in
                                              take(autosomal, args.singletons)]
    blocks["singleton_sex"] = [fmt(v) for v in
                               take(autosomal, args.singletons)]

    meta = {
        "seed": args.seed,
        "islands": ISLANDS,
        "counts": {k: len(v) for k, v in blocks.items()},
        "total_biased_rsids": len(used_rsids),
        "panel_size": len(panel),
        "gsa_loadings_intersection": len(in_load),
    }
    out = {"_meta": meta, **blocks}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))
    print(f"Wrote {args.out}")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

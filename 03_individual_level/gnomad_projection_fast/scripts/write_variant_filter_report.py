#!/usr/bin/env python3
"""Write per-variant filter diagnostics for the gnomAD fast projection step."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_bim(path: Path) -> list[dict[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                fields = line.strip().split()
            if len(fields) < 6:
                raise ValueError(f"{path}:{idx + 1}: expected 6 BIM columns")
            rows.append(
                {
                    "chrom": fields[0],
                    "variant_id": fields[1],
                    "cm": fields[2],
                    "pos": fields[3],
                    "a1": fields[4],
                    "a2": fields[5],
                }
            )
    return rows


def count_fam(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def hwe_exact(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
    """Exact HWE p-value, adapted from Wigginton et al. midpoint algorithm."""
    obs_homc = min(obs_hom1, obs_hom2)
    obs_homr = max(obs_hom1, obs_hom2)
    rare_copies = 2 * obs_homc + obs_hets
    genotypes = obs_hets + obs_homc + obs_homr
    if genotypes <= 0:
        return float("nan")

    probs = [0.0] * (rare_copies + 1)
    mid = int(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes))
    if (rare_copies & 1) ^ (mid & 1):
        mid += 1

    probs[mid] = 1.0
    total = probs[mid]

    curr_hets = mid
    curr_homr = (rare_copies - mid) // 2
    curr_homc = genotypes - curr_hets - curr_homr
    while curr_hets > 1:
        probs[curr_hets - 2] = (
            probs[curr_hets]
            * curr_hets
            * (curr_hets - 1)
            / (4.0 * (curr_homr + 1) * (curr_homc + 1))
        )
        total += probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = mid
    curr_homr = (rare_copies - mid) // 2
    curr_homc = genotypes - curr_hets - curr_homr
    while curr_hets <= rare_copies - 2:
        probs[curr_hets + 2] = (
            probs[curr_hets]
            * 4.0
            * curr_homr
            * curr_homc
            / ((curr_hets + 2) * (curr_hets + 1))
        )
        total += probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    if total <= 0:
        return float("nan")
    probs = [p / total for p in probs]
    observed = probs[obs_hets] if obs_hets < len(probs) else 0.0
    p_hwe = sum(p for p in probs if p <= observed)
    return min(1.0, max(0.0, p_hwe))


def iter_bed_variant_stats(bed_path: Path, n_samples: int, n_variants: int):
    bytes_per_variant = (n_samples + 3) // 4
    with bed_path.open("rb") as handle:
        magic = handle.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(f"{bed_path}: expected SNP-major PLINK BED magic bytes")
        for _variant_idx in range(n_variants):
            block = handle.read(bytes_per_variant)
            if len(block) != bytes_per_variant:
                raise ValueError(f"{bed_path}: truncated BED variant block")
            hom_a1 = het = hom_a2 = missing = 0
            seen = 0
            for byte in block:
                for shift in (0, 2, 4, 6):
                    if seen >= n_samples:
                        break
                    code = (byte >> shift) & 0b11
                    if code == 0b00:
                        hom_a1 += 1
                    elif code == 0b01:
                        missing += 1
                    elif code == 0b10:
                        het += 1
                    elif code == 0b11:
                        hom_a2 += 1
                    seen += 1
            called = hom_a1 + het + hom_a2
            allele1_count = 2 * hom_a1 + het
            allele2_count = 2 * hom_a2 + het
            allele_total = allele1_count + allele2_count
            maf = (
                min(allele1_count, allele2_count) / allele_total
                if allele_total
                else float("nan")
            )
            yield {
                "hom_a1": hom_a1,
                "het": het,
                "hom_a2": hom_a2,
                "missing_count": missing,
                "called_count": called,
                "missing_rate": missing / n_samples if n_samples else float("nan"),
                "allele1_count": allele1_count,
                "allele2_count": allele2_count,
                "maf": maf,
                "hwe_p": hwe_exact(het, hom_a1, hom_a2) if called else float("nan"),
            }


def fmt(value: float | int | str) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.8g}"
    return str(value)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pre-prefix", required=True, help="Prefix before QC filters, e.g. study_nodup")
    ap.add_argument("--post-prefix", required=True, help="Prefix after QC filters, e.g. study_qc")
    ap.add_argument("--out", required=True)
    ap.add_argument("--geno", type=float, required=True)
    ap.add_argument("--maf", type=float, required=True)
    ap.add_argument("--hwe", type=float, required=True)
    args = ap.parse_args()

    pre = Path(args.pre_prefix)
    post = Path(args.post_prefix)
    out = Path(args.out)
    variants = read_bim(pre.with_suffix(".bim"))
    kept = {row["variant_id"] for row in read_bim(post.with_suffix(".bim"))}
    n_samples = count_fam(pre.with_suffix(".fam"))

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "pos",
                "variant_id",
                "a1",
                "a2",
                "filter_reason",
                "missing_count",
                "sample_count",
                "missing_rate",
                "called_count",
                "hom_a1",
                "het",
                "hom_a2",
                "a1_count",
                "a2_count",
                "maf",
                "hwe_p",
            ]
        )
        for row, stats in zip(
            variants,
            iter_bed_variant_stats(pre.with_suffix(".bed"), n_samples, len(variants)),
        ):
            if row["variant_id"] in kept:
                continue
            reasons = []
            if stats["missing_rate"] > args.geno:
                reasons.append("GENO")
            if not math.isnan(stats["maf"]) and stats["maf"] < args.maf:
                reasons.append("MAF")
            if not math.isnan(stats["hwe_p"]) and stats["hwe_p"] < args.hwe:
                reasons.append("HWE")
            if not reasons:
                reasons.append("OTHER")
            writer.writerow(
                [
                    row["chrom"],
                    row["pos"],
                    row["variant_id"],
                    row["a1"],
                    row["a2"],
                    ",".join(reasons),
                    stats["missing_count"],
                    n_samples,
                    fmt(stats["missing_rate"]),
                    stats["called_count"],
                    stats["hom_a1"],
                    stats["het"],
                    stats["hom_a2"],
                    stats["allele1_count"],
                    stats["allele2_count"],
                    fmt(stats["maf"]),
                    fmt(stats["hwe_p"]),
                ]
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

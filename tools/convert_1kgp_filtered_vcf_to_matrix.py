#!/usr/bin/env python3
"""Convert filtered 1KGP VCFs into a compact dosage matrix.

The output is intended to be the reusable 1000 Genomes reference matrix for
BioVault PCA. It keeps only biallelic A/C/G/T SNP records, because a single
dosage value is only unambiguous when each locus has one REF and one ALT allele.

Outputs in --out-dir:
  hgp1k_dosage.npz       compressed numpy archive:
                           dosage uint8 [n_samples, n_variants]
                           0/1/2 = ALT allele dosage
                           missing_mask_packed stores missing calls, if any
                           samples, chrom, pos, rsid, ref, alt
  hgp1k_dosage.tsv       human-readable full matrix:
                           chrom/pos/rsid/ref/alt plus one dosage column per sample
  variants.tsv           chrom/pos/rsid/ref/alt per matrix column
  samples.tsv            sample_id per matrix row
  matrix_preview.tsv     small human-readable dosage preview; missing as blank
  matrix_report.txt      summary counts
"""

from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


BASES = {"A", "C", "G", "T"}
MISSING = np.uint8(255)
AUTOSOMES = [str(i) for i in range(1, 23)]


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str

    @property
    def key(self) -> tuple[str, int]:
        return self.chrom, self.pos


def setup_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")


def normalize_chrom(value: str) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return "MT" if text == "M" else text.upper()


def chromosome_list(raw: str) -> list[str]:
    if raw == "all":
        return AUTOSOMES
    out: list[str] = []
    for item in raw.split(","):
        chrom = normalize_chrom(item)
        if chrom not in AUTOSOMES + ["X", "Y", "MT"]:
            raise ValueError(f"unsupported chromosome: {item!r}")
        out.append(chrom)
    return out


def vcf_name_for_chr(chrom: str) -> str:
    if chrom == "X":
        return "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    return f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def find_vcfs(vcf_dir: Path, chromosomes: list[str]) -> list[Path]:
    vcfs: list[Path] = []
    for chrom in chromosomes:
        path = vcf_dir / vcf_name_for_chr(chrom)
        if not path.exists():
            logging.warning("missing VCF for chr%s: %s", chrom, path)
            continue
        if not Path(f"{path}.tbi").exists():
            logging.warning("missing index, skipping: %s.tbi", path)
            continue
        vcfs.append(path)
    if not vcfs:
        raise FileNotFoundError(f"no usable VCFs found in {vcf_dir}")
    return vcfs


def command_text(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)


def vcf_samples(vcf: Path) -> list[str]:
    text = command_text(["bcftools", "query", "-l", str(vcf)])
    return [line.strip() for line in text.splitlines() if line.strip()]


def load_locus_rsids(path: Path) -> dict[tuple[str, int], str]:
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str, "pos": int, "rsid": str})
    df["chrom"] = df["chrom"].map(normalize_chrom)
    out: dict[tuple[str, int], str] = {}
    conflicts: set[tuple[str, int]] = set()
    for chrom, pos, rsid in df[["chrom", "pos", "rsid"]].itertuples(index=False, name=None):
        key = (chrom, int(pos))
        if key in conflicts:
            continue
        old = out.get(key)
        if old is not None and old != rsid:
            out.pop(key, None)
            conflicts.add(key)
            continue
        out[key] = str(rsid)
    return out


def parse_gt_dosage(gt: str) -> int:
    gt = gt.split(":", 1)[0]
    if "." in gt:
        return int(MISSING)
    parts = gt.replace("|", "/").split("/")
    if len(parts) != 2:
        return int(MISSING)
    dosage = 0
    for allele in parts:
        if allele == "0":
            continue
        if allele == "1":
            dosage += 1
            continue
        return int(MISSING)
    return dosage


def stream_vcfs(
    vcfs: list[Path],
    locus_rsids: dict[tuple[str, int], str],
    max_variants: int | None,
) -> tuple[list[str], list[Variant], np.ndarray, dict[str, int]]:
    log = logging.getLogger(__name__)
    samples = vcf_samples(vcfs[0])
    sample_tuple = tuple(samples)
    rows: list[np.ndarray] = []
    variants: list[Variant] = []
    counts = {
        "vcf_rows": 0,
        "kept_biallelic_snps": 0,
        "dropped_non_biallelic_or_non_snp": 0,
    }

    for vcf in vcfs:
        if tuple(vcf_samples(vcf)) != sample_tuple:
            raise ValueError(f"{vcf}: sample list/order differs from first VCF")
        log.info("reading %s", vcf)
        cmd = [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n",
            str(vcf),
        ]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )
        assert proc.stdout is not None
        stopped_early = False
        chr_rows = 0
        chr_kept = 0
        for line in proc.stdout:
            counts["vcf_rows"] += 1
            chr_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5 + len(samples):
                counts["dropped_non_biallelic_or_non_snp"] += 1
                continue
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper()
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in BASES
                or alt not in BASES
                or "," in alt
            ):
                counts["dropped_non_biallelic_or_non_snp"] += 1
                continue
            dosage = np.fromiter(
                (parse_gt_dosage(gt) for gt in fields[5:]),
                dtype=np.uint8,
                count=len(samples),
            )
            variants.append(Variant(
                chrom=chrom,
                pos=pos,
                rsid=locus_rsids.get((chrom, pos), fields[2]),
                ref=ref,
                alt=alt,
            ))
            rows.append(dosage)
            counts["kept_biallelic_snps"] += 1
            chr_kept += 1
            if max_variants is not None and len(variants) >= max_variants:
                stopped_early = True
                proc.terminate()
                break
        if stopped_early:
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        rc = proc.wait()
        if rc != 0 and not stopped_early:
            raise RuntimeError(f"bcftools query failed for {vcf}:\n{stderr}")
        log.info("  rows=%s kept_biallelic_snps=%s", f"{chr_rows:,}", f"{chr_kept:,}")
        if max_variants is not None and len(variants) >= max_variants:
            break

    if not variants:
        raise ValueError("no biallelic SNP rows found")
    dosage = np.stack(rows, axis=1)  # samples x variants
    return samples, variants, dosage, counts


def drop_duplicate_positions(
    samples: list[str],
    variants: list[Variant],
    dosage: np.ndarray,
) -> tuple[list[str], list[Variant], np.ndarray, int]:
    key_counts: dict[tuple[str, int], int] = {}
    for variant in variants:
        key_counts[variant.key] = key_counts.get(variant.key, 0) + 1
    keep = np.array([key_counts[v.key] == 1 for v in variants], dtype=bool)
    return samples, [v for v, ok in zip(variants, keep) if ok], dosage[:, keep], int((~keep).sum())


def _format_dosage_row(values: np.ndarray) -> list[str | int]:
    return ["" if value == MISSING else int(value) for value in values]


def write_full_matrix_tsv(
    path: Path,
    samples: list[str],
    variants: list[Variant],
    dosage: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "rsid", "ref", "alt", *samples])
        for j, variant in enumerate(variants):
            writer.writerow([
                variant.chrom,
                variant.pos,
                variant.rsid,
                variant.ref,
                variant.alt,
                *_format_dosage_row(dosage[:, j]),
            ])


def write_tsvs(out_dir: Path, samples: list[str], variants: list[Variant], dosage: np.ndarray, preview_variants: int) -> None:
    with (out_dir / "samples.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id"])
        for sample in samples:
            writer.writerow([sample])

    with (out_dir / "variants.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "rsid", "ref", "alt"])
        for variant in variants:
            writer.writerow([variant.chrom, variant.pos, variant.rsid, variant.ref, variant.alt])

    n_preview = min(preview_variants, len(variants))
    with (out_dir / "matrix_preview.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "pos", "rsid", "ref", "alt", *samples])
        for j in range(n_preview):
            variant = variants[j]
            row = _format_dosage_row(dosage[:, j])
            writer.writerow([variant.chrom, variant.pos, variant.rsid, variant.ref, variant.alt, *row])


def main() -> None:
    setup_logging()
    log = logging.getLogger(__name__)
    repo_root = Path(__file__).resolve().parents[1]

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-dir", type=Path, default=repo_root / "data" / "1kgp_high_coverage" / "filtered")
    ap.add_argument("--locus-map", type=Path, default=repo_root / "tools" / "locus_map.tsv")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--chromosomes", default="all", help="all, 1, 1,2,3, etc.")
    ap.add_argument("--max-variants", type=int, default=None, help="testing cap")
    ap.add_argument("--preview-variants", type=int, default=25)
    ap.add_argument("--no-matrix-tsv", action="store_true", help="skip full hgp1k_dosage.tsv and write only the compact npz plus small TSVs")
    ap.add_argument("--keep-duplicate-positions", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    chromosomes = chromosome_list(args.chromosomes)
    vcfs = find_vcfs(args.vcf_dir, chromosomes)
    locus_rsids = load_locus_rsids(args.locus_map)
    samples, variants, dosage, counts = stream_vcfs(vcfs, locus_rsids, args.max_variants)
    duplicate_positions_dropped = 0
    if not args.keep_duplicate_positions:
        _ignored, variants, dosage, duplicate_positions_dropped = drop_duplicate_positions(samples, variants, dosage)

    log.info("writing compact matrix: %s samples x %s variants", f"{len(samples):,}", f"{len(variants):,}")
    missing_mask = dosage == MISSING
    dosage_out = dosage.copy()
    if missing_mask.any():
        dosage_out[missing_mask] = 0
    missing_mask_packed = np.packbits(missing_mask.reshape(-1))
    np.savez_compressed(
        args.out_dir / "hgp1k_dosage.npz",
        dosage=dosage_out,
        samples=np.array(samples, dtype=object),
        chrom=np.array([v.chrom for v in variants], dtype=object),
        pos=np.array([v.pos for v in variants], dtype=np.int64),
        rsid=np.array([v.rsid for v in variants], dtype=object),
        ref=np.array([v.ref for v in variants], dtype="S1"),
        alt=np.array([v.alt for v in variants], dtype="S1"),
        missing_mask_packed=missing_mask_packed,
        missing_mask_shape=np.array(missing_mask.shape, dtype=np.int64),
    )
    if not args.no_matrix_tsv:
        log.info("writing full human-readable dosage TSV")
        write_full_matrix_tsv(args.out_dir / "hgp1k_dosage.tsv", samples, variants, dosage)
    write_tsvs(args.out_dir, samples, variants, dosage, args.preview_variants)

    called = dosage != MISSING
    missing_rate = 1.0 - float(called.sum()) / float(dosage.size)
    dosage_values = dosage[called]
    hom_ref = int((dosage_values == 0).sum())
    het = int((dosage_values == 1).sum())
    hom_alt = int((dosage_values == 2).sum())

    report = args.out_dir / "matrix_report.txt"
    report.write_text(
        "=== 1KGP filtered VCF dosage matrix ===\n"
        f"VCF dir: {args.vcf_dir}\n"
        f"Locus map: {args.locus_map}\n"
        f"Chromosomes requested: {','.join(chromosomes)}\n"
        f"VCFs used: {len(vcfs)}\n"
        f"VCF rows read: {counts['vcf_rows']:,}\n"
        f"Biallelic SNP rows kept before duplicate-position drop: {counts['kept_biallelic_snps']:,}\n"
        f"Non-biallelic/non-SNP rows dropped: {counts['dropped_non_biallelic_or_non_snp']:,}\n"
        f"Duplicate-position rows dropped: {duplicate_positions_dropped:,}\n"
        f"Samples: {len(samples):,}\n"
        f"Variants: {len(variants):,}\n"
        f"Matrix shape: {len(samples):,} x {len(variants):,}\n"
        f"Full matrix TSV: {'no' if args.no_matrix_tsv else 'yes'}\n"
        f"Dosage encoding: 0/1/2 ALT allele dosage; missing calls stored in missing_mask_packed\n"
        f"Missing genotype rate: {missing_rate:.8f}\n"
        f"Dosage counts: hom_ref_0={hom_ref:,}, het_1={het:,}, hom_alt_2={hom_alt:,}\n"
        f"Elapsed seconds: {time.time() - t0:.1f}\n"
    )
    log.info("wrote %s", args.out_dir)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Joint PCA of BioVault samples and filtered 1000 Genomes high-coverage VCFs.

This pipeline deliberately does not project onto precomputed gnomAD loadings.
It builds a shared SNP matrix from:

  * 1000 Genomes high-coverage VCFs filtered to tools/locus_map.tsv positions.
  * BioVault/DDNA/Illumina genotype text files normalized by tools.genotype_normalizer.

It then computes allele frequencies for both groups, imputes missing genotypes
to the combined mean dosage, standardizes by the combined allele frequency, and
fits a PCA on the combined 1KGP + study matrix.

Outputs in out_dir:
  pca_scores.tsv              all 1KGP and study samples
  study_pca_projection.tsv    study-only legacy shape: s<TAB>scores
  allele_freqs.tsv            per-variant 1KGP and study allele frequencies
  qc_report.txt               run summary
  pca_projection.png          PC1/PC2 scatter
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

for _parent in Path(__file__).resolve().parents:
    if (_parent / "tools" / "genotype_normalizer.py").exists():
        sys.path.insert(0, str(_parent))
        break
sys.path.append("/opt/biovault")
from tools import genotype_normalizer as genoio  # noqa: E402


BASES = {"A", "C", "G", "T"}
AUTOSOMES = [str(i) for i in range(1, 23)]
MISSING = np.uint8(255)
REGION_CMAPS = {
    "AFR": "Oranges",
    "AMR": "Blues",
    "CSA": "YlGn",
    "EAS": "Greens",
    "EUR": "Reds",
    "MID": "Purples",
    "OCE": "PuRd",
}
REGION_ORDER = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]
REGION_ALIASES = {
    "SAS": "CSA",
    "AMR_HGDP": "AMR",
}


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


def chromosome_list(value: str) -> list[str]:
    if value == "all":
        return AUTOSOMES
    out: list[str] = []
    for raw in value.split(","):
        chrom = normalize_chrom(raw)
        if chrom not in AUTOSOMES:
            raise ValueError(f"unsupported chromosome for PCA: {raw!r}")
        out.append(chrom)
    return out


def vcf_name_for_chr(chrom: str) -> str:
    if chrom == "X":
        return "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    return f"1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def find_vcfs(vcf_dir: Path, chromosomes: list[str]) -> list[Path]:
    vcfs: list[Path] = []
    for chrom in chromosomes:
        candidate = vcf_dir / vcf_name_for_chr(chrom)
        if not candidate.exists():
            logging.warning("missing filtered 1KGP VCF for chr%s: %s", chrom, candidate)
            continue
        if (candidate.parent / f"{candidate.name}.aria2").exists():
            logging.warning("skipping active/incomplete aria2 download: %s", candidate)
            continue
        if not (candidate.parent / f"{candidate.name}.tbi").exists():
            logging.warning("skipping unindexed VCF: %s", candidate)
            continue
        vcfs.append(candidate)
    if not vcfs:
        raise FileNotFoundError(f"no usable filtered VCFs found in {vcf_dir}")
    return vcfs


def run_text(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)


def vcf_samples(vcf: Path) -> list[str]:
    text = run_text(["bcftools", "query", "-l", str(vcf)])
    return [line.strip() for line in text.splitlines() if line.strip()]


def parse_gt_dosage(gt: str) -> int:
    gt = gt.split(":", 1)[0]
    if "." in gt:
        return int(MISSING)
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2:
        return int(MISSING)
    dosage = 0
    for allele in alleles:
        if allele == "0":
            continue
        if allele == "1":
            dosage += 1
            continue
        return int(MISSING)
    return dosage


def load_locus_rsids(locus_map: Path) -> dict[tuple[str, int], str]:
    df = pd.read_csv(locus_map, sep="\t", dtype={"chrom": str, "pos": int, "rsid": str})
    df["chrom"] = df["chrom"].map(normalize_chrom)
    out: dict[tuple[str, int], str] = {}
    conflicts: set[tuple[str, int]] = set()
    for chrom, pos, rsid in df[["chrom", "pos", "rsid"]].itertuples(index=False, name=None):
        key = (chrom, int(pos))
        if key in conflicts:
            continue
        existing = out.get(key)
        if existing is not None and existing != rsid:
            out.pop(key, None)
            conflicts.add(key)
            continue
        out[key] = str(rsid)
    return out


def stream_vcf_matrix(
    vcfs: list[Path],
    locus_rsids: dict[tuple[str, int], str],
    max_variants: int | None = None,
) -> tuple[list[str], list[Variant], np.ndarray]:
    log = logging.getLogger(__name__)
    samples = vcf_samples(vcfs[0])
    rows: list[np.ndarray] = []
    variants: list[Variant] = []
    seen_sample_tuple = tuple(samples)

    for vcf in vcfs:
        current_samples = vcf_samples(vcf)
        if tuple(current_samples) != seen_sample_tuple:
            raise ValueError(f"{vcf}: sample list/order differs from first VCF")

        log.info("reading 1KGP VCF: %s", vcf)
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
        n_read = 0
        n_kept = 0
        stopped_early = False
        for line in proc.stdout:
            n_read += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 + len(samples):
                continue
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
                continue
            if "," in alt:
                continue
            rsid = locus_rsids.get((chrom, pos), fields[2])
            dosage = np.fromiter(
                (parse_gt_dosage(gt) for gt in fields[5:]),
                dtype=np.uint8,
                count=len(samples),
            )
            rows.append(dosage)
            variants.append(Variant(chrom=chrom, pos=pos, rsid=rsid, ref=ref, alt=alt))
            n_kept += 1
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
        if stderr.strip():
            log.warning("bcftools warning for %s: %s", vcf.name, stderr.strip().splitlines()[-1])
        log.info("  read %s VCF rows, kept %s biallelic SNP rows", f"{n_read:,}", f"{n_kept:,}")
        if max_variants is not None and len(variants) >= max_variants:
            break

    if not variants:
        raise ValueError("no biallelic SNP variants loaded from 1KGP VCFs")
    matrix = np.stack(rows, axis=1)  # samples x variants, uint8, 255 = missing
    return samples, variants, matrix


def _decode_array_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def load_reference_matrix_npz(path: Path) -> tuple[list[str], list[Variant], np.ndarray]:
    arc = np.load(path, allow_pickle=True)
    required = {"dosage", "samples", "chrom", "pos", "rsid", "ref", "alt"}
    missing = required - set(arc.files)
    if missing:
        raise ValueError(f"{path}: missing npz keys: {sorted(missing)}")
    dosage = arc["dosage"].astype(np.uint8, copy=True)
    if "missing_mask_packed" in arc.files and "missing_mask_shape" in arc.files:
        shape = tuple(int(x) for x in arc["missing_mask_shape"].tolist())
        if shape != dosage.shape:
            raise ValueError(f"{path}: missing mask shape {shape} does not match dosage shape {dosage.shape}")
        missing = np.unpackbits(arc["missing_mask_packed"], count=dosage.size).reshape(dosage.shape).astype(bool)
        if missing.any():
            dosage[missing] = MISSING
    samples = [_decode_array_value(value) for value in arc["samples"].tolist()]
    chrom = [_decode_array_value(value) for value in arc["chrom"].tolist()]
    pos = arc["pos"].astype(np.int64)
    rsid = [_decode_array_value(value) for value in arc["rsid"].tolist()]
    ref = [_decode_array_value(value) for value in arc["ref"].tolist()]
    alt = [_decode_array_value(value) for value in arc["alt"].tolist()]
    variants = [
        Variant(
            chrom=normalize_chrom(c),
            pos=int(p),
            rsid=r,
            ref=rf,
            alt=al,
        )
        for c, p, r, rf, al in zip(chrom, pos, rsid, ref, alt)
    ]
    if dosage.shape != (len(samples), len(variants)):
        raise ValueError(
            f"{path}: dosage shape {dosage.shape} does not match "
            f"{len(samples)} samples x {len(variants)} variants"
        )
    return samples, variants, dosage


def drop_duplicate_loci(variants: list[Variant], matrix: np.ndarray) -> tuple[list[Variant], np.ndarray, int]:
    counts: dict[tuple[str, int], int] = {}
    for variant in variants:
        counts[variant.key] = counts.get(variant.key, 0) + 1
    keep = np.array([counts[v.key] == 1 for v in variants], dtype=bool)
    return [v for v, ok in zip(variants, keep) if ok], matrix[:, keep], int((~keep).sum())


def discover_study_files(data_dir: Path) -> list[tuple[str, Path]]:
    sample_dirs = sorted(d for d in data_dir.iterdir() if d.is_dir())
    out: list[tuple[str, Path]] = []
    for sample_dir in sample_dirs:
        txts = sorted(sample_dir.glob("*.txt"))
        if txts:
            out.append((sample_dir.name, txts[0]))
        else:
            logging.warning("skip %s: no .txt genotype file", sample_dir.name)
    if not out:
        raise FileNotFoundError(f"no participant .txt files found under {data_dir}")
    return out


def load_study_matrix(
    sample_files: list[tuple[str, Path]],
    variants: list[Variant],
    min_gs: float,
) -> tuple[list[str], np.ndarray, list[dict[str, str]]]:
    key_to_idx = {variant.key: idx for idx, variant in enumerate(variants)}
    ref = np.array([v.ref for v in variants], dtype=object)
    alt = np.array([v.alt for v in variants], dtype=object)
    matrix = np.full((len(sample_files), len(variants)), MISSING, dtype=np.uint8)
    sample_ids: list[str] = []
    errors: list[dict[str, str]] = []

    for row_idx, (sample_id, path) in enumerate(sample_files):
        sample_ids.append(sample_id)
        try:
            df = genoio.read_pipeline_genotypes(path)
        except Exception as exc:
            errors.append({
                "participant_id": sample_id,
                "file": str(path),
                "code": "PARSE_FAILED",
                "message": str(exc),
            })
            continue
        if df.empty:
            continue
        df["chrom"] = df["chrom"].map(normalize_chrom)
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
        df["gs"] = pd.to_numeric(df["gs"], errors="coerce").fillna(1.0)
        df = df[df["pos"].notna() & df["chrom"].isin(AUTOSOMES)]
        df = df[df["gs"] >= min_gs]
        df = df[df["gt"].astype(str).str.fullmatch(r"[ACGT]{2}", na=False)]

        seen: set[int] = set()
        for chrom, pos_raw, gt in df[["chrom", "pos", "gt"]].itertuples(index=False, name=None):
            key = (chrom, int(pos_raw))
            idx = key_to_idx.get(key)
            if idx is None or idx in seen:
                continue
            seen.add(idx)
            a1, a2 = str(gt)[0], str(gt)[1]
            if (a1 not in {ref[idx], alt[idx]}) or (a2 not in {ref[idx], alt[idx]}):
                continue
            matrix[row_idx, idx] = np.uint8((a1 == alt[idx]) + (a2 == alt[idx]))
    return sample_ids, matrix, errors


def allele_stats(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    called = matrix != MISSING
    called_samples = called.sum(axis=0).astype(np.int64)
    alt_sum = np.where(called, matrix, 0).sum(axis=0).astype(np.float64)
    allele_n = called_samples * 2
    af = np.divide(alt_sum, allele_n, out=np.full_like(alt_sum, np.nan), where=allele_n > 0)
    missing_rate = 1.0 - (called_samples / max(matrix.shape[0], 1))
    return af, allele_n.astype(np.int64), missing_rate


def write_errors(path: Path, errors: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["participant_id", "file", "severity", "code", "message"])
        for row in errors:
            writer.writerow([
                row.get("participant_id", ""),
                row.get("file", ""),
                "ERROR",
                row.get("code", ""),
                row.get("message", ""),
            ])


def write_allele_freqs(
    path: Path,
    variants: list[Variant],
    ref_af: np.ndarray,
    ref_an: np.ndarray,
    study_af: np.ndarray,
    study_an: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chrom", "pos", "rsid", "ref", "alt",
            "hgp1k_af", "hgp1k_allele_number",
            "study_af", "study_allele_number",
        ])
        for i, variant in enumerate(variants):
            writer.writerow([
                variant.chrom,
                variant.pos,
                variant.rsid,
                variant.ref,
                variant.alt,
                "" if np.isnan(ref_af[i]) else f"{ref_af[i]:.8g}",
                int(ref_an[i]),
                "" if np.isnan(study_af[i]) else f"{study_af[i]:.8g}",
                int(study_an[i]),
            ])


def fill_standardized_block(
    out: np.ndarray,
    start_row: int,
    matrix: np.ndarray,
    keep_idx: np.ndarray,
    mean: np.ndarray,
    sd: np.ndarray,
    chunk_size: int,
) -> float:
    total_ss = 0.0
    n_cols = keep_idx.size
    for start in range(0, n_cols, chunk_size):
        end = min(start + chunk_size, n_cols)
        idx = keep_idx[start:end]
        block = matrix[:, idx].astype(np.float32, copy=True)
        miss = block == float(MISSING)
        block -= mean[start:end][None, :].astype(np.float32)
        block /= sd[start:end][None, :].astype(np.float32)
        block[miss] = 0.0
        out[start_row:start_row + matrix.shape[0], start:end] = block
        total_ss += float(np.square(block, dtype=np.float64).sum())
    return total_ss


def load_metadata(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".csv":
        sep = ","
    else:
        sep = r"\s+"
    df = pd.read_csv(path, sep=sep, dtype=str)
    cols = {c.lower(): c for c in df.columns}
    sample_col = (
        cols.get("sample_id")
        or cols.get("sample")
        or cols.get("sampleid")
        or cols.get("iid")
    )
    if not sample_col:
        raise ValueError(f"{path}: expected sample_id/sample/SampleID/IID column")
    df = df.rename(columns={sample_col: "sample_id"})
    rename = {}
    for canonical, options in {
        "population": ("population", "pop"),
        "superpopulation": ("superpopulation", "super_pop", "genetic_region"),
    }.items():
        for option in options:
            col = cols.get(option)
            if col and col != canonical:
                rename[col] = canonical
                break
    if rename:
        df = df.rename(columns=rename)
    if "superpopulation" in df.columns:
        df["region"] = df["superpopulation"].map(lambda x: REGION_ALIASES.get(str(x).strip(), str(x).strip()))
    elif "region" not in df.columns:
        df["region"] = ""
    if "population" not in df.columns:
        df["population"] = ""
    return df


def subpopulation_colors(ref: pd.DataFrame):
    from matplotlib import cm

    colors = {}
    for region in REGION_ORDER:
        pops = sorted(ref.loc[ref["region"] == region, "population"].dropna().unique())
        pops = [pop for pop in pops if str(pop).strip()]
        if not pops:
            continue
        cmap = cm.get_cmap(REGION_CMAPS.get(region, "Greys"))
        for i, pop in enumerate(pops):
            shade = 0.35 + 0.60 * (i / max(1, len(pops) - 1))
            colors[(region, pop)] = cmap(shade)
    return colors


def plot_scores(scores_df: pd.DataFrame, out_path: Path, metadata: pd.DataFrame | None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    df = scores_df.copy()
    if metadata is not None:
        df = df.merge(metadata, on="sample_id", how="left")

    fig, ax = plt.subplots(figsize=(12, 9))
    ref = df[df["source"] == "1KGP"]
    study = df[df["source"] == "BioVault"]

    handles = []
    if metadata is not None and {"region", "population"}.issubset(ref.columns):
        colors = subpopulation_colors(ref)
        for region in REGION_ORDER:
            region_ref = ref[ref["region"] == region]
            if region_ref.empty:
                continue
            for pop in sorted(region_ref["population"].dropna().unique()):
                pop_ref = region_ref[region_ref["population"] == pop]
                if pop_ref.empty:
                    continue
                color = colors.get((region, pop), "#8a8f98")
                ax.scatter(
                    pop_ref["PC1"],
                    pop_ref["PC2"],
                    s=13,
                    alpha=0.65,
                    color=color,
                    edgecolors="none",
                )
                handles.append(Line2D(
                    [0], [0],
                    marker="o",
                    linestyle="",
                    markerfacecolor=color,
                    markeredgecolor="none",
                    markersize=5,
                    label=f"[{region}] {pop} (n={len(pop_ref)})",
                ))
    else:
        ax.scatter(ref["PC1"], ref["PC2"], c="#8a8f98", s=9, alpha=0.35, label="1KGP")
        handles.append(Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor="#8a8f98",
            markeredgecolor="none", markersize=5, label=f"1KGP (n={len(ref)})",
        ))

    ax.scatter(
        study["PC1"], study["PC2"],
        c="#111111", marker="*", s=180,
        edgecolors="#f2c94c", linewidths=0.8,
    )
    handles.append(Line2D(
        [0], [0],
        marker="*",
        linestyle="",
        markerfacecolor="#111111",
        markeredgecolor="#f2c94c",
        markeredgewidth=0.8,
        markersize=13,
        label=f"BioVault (n={len(study)})",
    ))
    for _, row in study.iterrows():
        ax.annotate(
            str(row["sample_id"]),
            (row["PC1"], row["PC2"]),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.axhline(0, color="#b8bcc4", linewidth=0.6)
    ax.axvline(0, color="#b8bcc4", linewidth=0.6)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("BioVault samples and 1000 Genomes high-coverage PCA")
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
        ncol=2 if len(handles) > 30 else 1,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    setup_logging()
    log = logging.getLogger(__name__)

    repo_root = next(
        (p for p in Path(__file__).resolve().parents if (p / "tools" / "locus_map.tsv").exists()),
        Path.cwd(),
    )
    default_vcf_dir = repo_root / "data" / "1kgp_high_coverage" / "filtered"
    default_metadata = repo_root / "data" / "1kgp_high_coverage" / "20130606_g1k_3202_samples_ped_population.txt"

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--vcf-dir", type=Path, default=Path(os.environ.get("HGP1K_VCF_DIR", default_vcf_dir)))
    ap.add_argument("--matrix-npz", type=Path, default=Path(os.environ["HGP1K_MATRIX_NPZ"]) if os.environ.get("HGP1K_MATRIX_NPZ") else None)
    ap.add_argument("--locus-map", type=Path, default=Path(os.environ.get("LOCUS_MAP", repo_root / "tools" / "locus_map.tsv")))
    ap.add_argument("--chromosomes", default=os.environ.get("CHROMOSOMES", "all"))
    ap.add_argument(
        "--metadata-tsv",
        type=Path,
        default=Path(os.environ["HGP1K_METADATA_TSV"]) if os.environ.get("HGP1K_METADATA_TSV") else (
            default_metadata if default_metadata.exists() else None
        ),
    )
    ap.add_argument("--n-components", type=int, default=10)
    ap.add_argument("--min-gs", type=float, default=0.15)
    ap.add_argument("--min-af", type=float, default=0.01)
    ap.add_argument("--max-ref-missing", type=float, default=0.05)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--max-variants", type=int, default=None, help="debug/testing cap after loading VCF SNP rows")
    ap.add_argument("--random-state", type=int, default=7)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    chromosomes = chromosome_list(args.chromosomes)
    vcfs: list[Path] = []
    if args.matrix_npz:
        log.info("loading compact 1KGP matrix: %s", args.matrix_npz)
        ref_sample_ids, variants, ref_matrix = load_reference_matrix_npz(args.matrix_npz)
        n_dup_loci = 0
        if args.max_variants is not None:
            variants = variants[:args.max_variants]
            ref_matrix = ref_matrix[:, :args.max_variants]
        log.info("using compact matrix instead of VCF dir")
    else:
        vcfs = find_vcfs(args.vcf_dir, chromosomes)
        log.info("using chromosomes: %s", ",".join(chromosomes))
        log.info("using %d filtered 1KGP VCF(s)", len(vcfs))
        locus_rsids = load_locus_rsids(args.locus_map)
        ref_sample_ids, variants, ref_matrix = stream_vcf_matrix(vcfs, locus_rsids, args.max_variants)
        variants, ref_matrix, n_dup_loci = drop_duplicate_loci(variants, ref_matrix)
    log.info("1KGP matrix after duplicate-locus drop: %s samples x %s variants", ref_matrix.shape[0], ref_matrix.shape[1])

    sample_files = discover_study_files(args.data_dir)
    study_sample_ids, study_matrix, errors = load_study_matrix(sample_files, variants, args.min_gs)
    write_errors(args.out_dir / "errors.tsv", errors)
    log.info("study matrix: %s samples x %s variants", study_matrix.shape[0], study_matrix.shape[1])

    ref_af, ref_an, ref_missing = allele_stats(ref_matrix)
    study_af, study_an, study_missing = allele_stats(study_matrix)
    write_allele_freqs(args.out_dir / "allele_freqs.tsv", variants, ref_af, ref_an, study_af, study_an)

    combined = np.concatenate([ref_matrix, study_matrix], axis=0)
    combined_af, combined_an, _combined_missing = allele_stats(combined)
    af_ok = np.isfinite(combined_af) & (combined_af >= args.min_af) & (combined_af <= (1.0 - args.min_af))
    ref_missing_ok = ref_missing <= args.max_ref_missing
    keep = af_ok & ref_missing_ok
    keep_idx = np.flatnonzero(keep)
    if keep_idx.size < 2:
        raise SystemExit("ERROR: fewer than two variants passed PCA filters")

    n_components = min(args.n_components, combined.shape[0] - 1, keep_idx.size)
    if n_components < 2:
        raise SystemExit("ERROR: need at least two samples and two variants for PC1/PC2")

    mean = (2.0 * combined_af[keep_idx]).astype(np.float32)
    sd = np.sqrt(2.0 * combined_af[keep_idx] * (1.0 - combined_af[keep_idx])).astype(np.float32)
    sd = np.where(sd > 0, sd, 1.0).astype(np.float32)

    matrix_path = args.out_dir / "combined_standardized.float32.mmap"
    X = np.memmap(matrix_path, dtype=np.float32, mode="w+", shape=(combined.shape[0], keep_idx.size))
    ss_ref = fill_standardized_block(X, 0, ref_matrix, keep_idx, mean, sd, args.chunk_size)
    ss_study = fill_standardized_block(X, ref_matrix.shape[0], study_matrix, keep_idx, mean, sd, args.chunk_size)
    X.flush()

    from sklearn.utils.extmath import randomized_svd
    log.info("running randomized SVD PCA: %s samples x %s variants", X.shape[0], X.shape[1])
    U, S, _Vt = randomized_svd(
        X,
        n_components=n_components,
        n_iter=5,
        random_state=args.random_state,
    )
    scores = U * S[None, :]
    explained_var = (S ** 2) / max(X.shape[0] - 1, 1)
    total_var = (ss_ref + ss_study) / max(X.shape[0] - 1, 1)
    explained_ratio = explained_var / total_var if total_var > 0 else np.full_like(explained_var, np.nan)

    all_sample_ids = ref_sample_ids + study_sample_ids
    sources = ["1KGP"] * len(ref_sample_ids) + ["BioVault"] * len(study_sample_ids)
    score_cols = [f"PC{i}" for i in range(1, n_components + 1)]
    scores_df = pd.DataFrame(scores, columns=score_cols)
    scores_df.insert(0, "source", sources)
    scores_df.insert(0, "sample_id", all_sample_ids)
    scores_df.to_csv(args.out_dir / "pca_scores.tsv", sep="\t", index=False)

    with (args.out_dir / "study_pca_projection.tsv").open("w", encoding="utf-8") as handle:
        handle.write("s\tscores\n")
        for sid, row in zip(study_sample_ids, scores[len(ref_sample_ids):]):
            handle.write(f"{sid}\t[{','.join(repr(float(v)) for v in row)}]\n")

    metadata = load_metadata(args.metadata_tsv)
    plot_scores(scores_df, args.out_dir / "pca_projection.png", metadata)

    with (args.out_dir / "qc_report.txt").open("w", encoding="utf-8") as handle:
        handle.write("=== hgp1k_projection_fast QC ===\n")
        handle.write(f"VCF dir: {args.vcf_dir}\n")
        handle.write(f"Matrix NPZ: {args.matrix_npz or ''}\n")
        handle.write(f"Locus map: {args.locus_map}\n")
        handle.write(f"Chromosomes requested: {','.join(chromosomes)}\n")
        handle.write(f"VCFs used: {len(vcfs)}\n")
        handle.write(f"1KGP samples: {len(ref_sample_ids):,}\n")
        handle.write(f"Study samples: {len(study_sample_ids):,}\n")
        handle.write(f"Loaded biallelic SNP variants: {len(variants) + n_dup_loci:,}\n")
        handle.write(f"Duplicate-position variants dropped: {n_dup_loci:,}\n")
        handle.write(f"Variants after duplicate-position drop: {len(variants):,}\n")
        handle.write(f"PCA variants after AF/ref-missing filters: {keep_idx.size:,}\n")
        handle.write(f"min_af={args.min_af}, max_ref_missing={args.max_ref_missing}, min_gs={args.min_gs}\n")
        handle.write(f"Mean 1KGP missing rate across kept variants: {float(ref_missing[keep_idx].mean()):.6f}\n")
        handle.write(f"Mean study missing rate across kept variants: {float(study_missing[keep_idx].mean()):.6f}\n")
        handle.write("Explained variance ratio:\n")
        for i, value in enumerate(explained_ratio, start=1):
            handle.write(f"  PC{i}: {float(value):.8g}\n")
        handle.write(f"Elapsed seconds: {time.time() - t_start:.1f}\n")

    try:
        matrix_path.unlink()
    except OSError:
        pass

    log.info("done: outputs in %s", args.out_dir)


if __name__ == "__main__":
    main()

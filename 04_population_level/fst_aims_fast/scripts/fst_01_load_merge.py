"""
FST step 1 (flow fork) - load & merge per-country allele frequency files.

Forked from 04_population_level/fst_islands/scripts/01_load_merge.py. The only
behavioural change is that the population set is resolved from the `country`
facet (via popset.resolve_populations) instead of a hardcoded ISLAND_FILES
dict, and it explodes if any expected allele_freq_<pop>.tsv is missing/empty.

Paths are env-driven so main.nf can place inputs/outputs in the Nextflow
workdir:
  BV_RAW_DIR   dir containing allele_freq_<pop>.tsv   (required by the flow)
  BV_WORK_DIR  dir for data/ and logs/                (default: parents[1])
  BV_POPULATIONS  comma-separated normalized country labels (see popset.py)

Outputs (unchanged schema):
  data/merged/merged_allele_freq_annotated.tsv
  data/merged/merged_allele_freq.tsv
  data/merged/merged_allele_number.tsv
  data/merged/snp_overlap_summary.txt
"""

import logging
import os
from pathlib import Path

import pandas as pd

from popset import resolve_populations

BASE_DIR = Path(os.environ.get("BV_WORK_DIR", Path(__file__).resolve().parents[1]))
RAW_DIR = Path(os.environ.get("BV_RAW_DIR", BASE_DIR.parent / "raw_allele_freq_country"))
OUT_DIR = BASE_DIR / "data" / "merged"
LOG_DIR = BASE_DIR / "logs"

LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "01_load_merge.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


def load_population(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"locus_key": str, "rsid": str})
    log.info(
        f"  {label}: {len(df):,} SNPs  |  "
        f"unique locus_key: {df['locus_key'].nunique():,}  |  "
        f"unique rsid: {df['rsid'].nunique():,}"
    )

    mask_zero = df["allele_number"] == 0
    log.info(f"    {mask_zero.sum():,} rows with allele_number=0 -> NaN")
    df.loc[mask_zero, "allele_freq"] = float("nan")
    df.loc[mask_zero, "allele_number"] = float("nan")

    df = df.set_index("locus_key")
    out = df[["allele_freq", "allele_number"]].rename(
        columns={"allele_freq": f"{label}_freq", "allele_number": f"{label}_n"}
    )
    out["rsid"] = df["rsid"]
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    populations = resolve_populations(RAW_DIR)
    log.info(f"Populations (from country facet): {populations}")

    log.info("Loading per-country allele-frequency files ...")
    frames = []
    for label in populations:
        frames.append(load_population(RAW_DIR / f"allele_freq_{label}.tsv", label))

    log.info("Merging on locus_key (inner join) ...")
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.join(f.drop(columns="rsid"), how="inner")

    freq_cols = [c for c in merged.columns if c.endswith("_freq")]
    n_cols = [c for c in merged.columns if c.endswith("_n")]

    before = len(merged)
    merged = merged.dropna(subset=freq_cols + n_cols)
    after = len(merged)
    log.info(
        f"SNPs after dropping any-population missing: "
        f"{after:,}  (removed {before - after:,})"
    )

    freq_matrix = merged[["rsid"] + freq_cols].copy()
    freq_matrix.columns = ["rsid"] + [c.replace("_freq", "") for c in freq_cols]

    n_matrix = merged[n_cols].copy()
    n_matrix.columns = [c.replace("_n", "") for c in n_cols]

    freq_only = freq_matrix.drop(columns="rsid")

    freq_matrix.to_csv(OUT_DIR / "merged_allele_freq_annotated.tsv", sep="\t")
    freq_only.to_csv(OUT_DIR / "merged_allele_freq.tsv", sep="\t")
    n_matrix.to_csv(OUT_DIR / "merged_allele_number.tsv", sep="\t")

    log.info(
        f"Saved freq matrix (annotated) -> merged_allele_freq_annotated.tsv  {freq_matrix.shape}"
    )
    log.info(f"Saved freq matrix             -> merged_allele_freq.tsv  {freq_only.shape}")
    log.info(f"Saved allele-number matrix    -> merged_allele_number.tsv  {n_matrix.shape}")

    lines = [
        "Merge key      : locus_key (chr-pos-ref-alt)",
        f"Populations    : {populations}",
        f"Overlapping SNPs (all populations non-missing): {after:,}",
        "",
        "Per-population max allele_number (~ 2 x max individuals genotyped):",
    ]
    for col in n_matrix.columns:
        lines.append(f"  {col:12s}: {int(n_matrix[col].max())}")

    (OUT_DIR / "snp_overlap_summary.txt").write_text("\n".join(lines))
    for ln in lines:
        log.info(ln)

    log.info("FST step 1 complete.")


if __name__ == "__main__":
    main()

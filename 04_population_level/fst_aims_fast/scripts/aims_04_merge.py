"""
AIMs step 04 (flow fork) - build the master allele-frequency table.

Forked from 04_population_level/aims_differential_snps/scripts/04_merge_carib_gnomad.py.

Changes vs the original:
  * Population set comes from the country facet (popset), not a hardcoded
    ISLANDS list.
  * The gnomAD reference AF table is the small pre-baked TSV shipped in the
    image (BV_GNOMAD_AF, default /opt/biovault/reference/aims/
    gnomad_af_per_locus.tsv) so runtime never touches the 80 GB HGDP+TGP VCFs.
    This replaces the original AIMs step 02 entirely.

Env:
  BV_WORK_DIR   AIMs working tree (data/ output)         default: parents[1]
  BV_FST_DIR    FST working tree (reads merged_allele_freq_annotated.tsv)
  BV_GNOMAD_AF  pre-baked gnomAD AF-per-locus TSV
  BV_POPULATIONS  comma-separated normalized country labels

Output:
  data/master_af_table.tsv      locus_key rsid <pops...> gnomAD_global/_AFR/_NFE/_SAS
  data/master_af_table_summary.txt
"""

import os
from pathlib import Path

import pandas as pd

from popset import resolve_populations, require_columns

BASE = Path(os.environ.get("BV_WORK_DIR", Path(__file__).resolve().parents[1]))
FST_DIR = Path(os.environ.get("BV_FST_DIR", BASE.parent / "fst_islands"))
CARIB = FST_DIR / "data" / "merged" / "merged_allele_freq_annotated.tsv"
GNOMAD = Path(
    os.environ.get("BV_GNOMAD_AF", "/opt/biovault/reference/aims/gnomad_af_per_locus.tsv")
)
RAW_DIR = Path(os.environ.get("BV_RAW_DIR", BASE.parent / "raw_allele_freq_country"))
OUT = BASE / "data" / "master_af_table.tsv"
SUMM = BASE / "data" / "master_af_table_summary.txt"

GNOMAD_COLS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]


def main() -> None:
    populations = resolve_populations(RAW_DIR)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    if not GNOMAD.exists():
        raise SystemExit(
            f"gnomAD reference AF table not found at {GNOMAD}. The "
            f"biovault-popgen image must ship it (see build_docker.sh / "
            f"Dockerfile) or set BV_GNOMAD_AF."
        )

    carib = pd.read_csv(CARIB, sep="\t")
    require_columns(carib.columns, populations, f"merged AF table {CARIB}")
    g = pd.read_csv(GNOMAD, sep="\t")
    missing_g = [c for c in GNOMAD_COLS if c not in g.columns]
    if missing_g:
        raise SystemExit(
            f"gnomAD reference {GNOMAD} missing columns: {missing_g}. "
            f"Present: {list(g.columns)}"
        )

    n_carib_in = len(carib)
    df = carib.merge(g[["locus_key"] + GNOMAD_COLS], on="locus_key", how="left")

    mask_carib = df[populations].notna().all(axis=1)
    mask_g = df[GNOMAD_COLS].notna().all(axis=1)
    keep = mask_carib & mask_g
    out = df.loc[keep, ["locus_key", "rsid"] + populations + GNOMAD_COLS].copy()

    out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")

    lines = [
        f"Populations                  : {populations}",
        f"Caribbean SNPs in            : {n_carib_in:,}",
        f"After requiring all pops     : {mask_carib.sum():,}",
        f"After requiring gnomAD ref   : {keep.sum():,}",
        f"Final master_af_table rows   : {len(out):,}",
        "",
        "Allele-frequency distributions (mean +/- sd):",
    ]
    for c in populations + GNOMAD_COLS:
        lines.append(f"  {c:<18s} {out[c].mean():.4f} +/- {out[c].std():.4f}")
    SUMM.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[04] wrote {OUT}")
    print(f"[04] wrote {SUMM}")

    if len(out) == 0:
        raise SystemExit(
            "[04] master_af_table is empty - no SNPs survived the "
            "Caribbean-complete + gnomAD-covered intersection. Cannot run "
            "AIMs downstream."
        )


if __name__ == "__main__":
    main()

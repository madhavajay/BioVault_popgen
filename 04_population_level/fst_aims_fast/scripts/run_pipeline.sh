#!/usr/bin/env bash
# run_pipeline.sh - FST + AIMs pipeline (biovault-popgen container)
#
# Consumes the per-country allele_freq_<country>.tsv files produced by the
# split step and runs:
#   FST : 01 load/merge -> 02 pairwise WC84 -> 03 visualise
#   AIMs: 04 merge w/ bundled gnomAD ref -> 05 differential SNPs -> 06 AIMs
#
# Args:
#   $1  RAW_DIR : dir containing allele_freq_<country>.tsv
#   $2  WORK    : working tree root (fst/ and aims/ subtrees created here)
#   $3  RESULTS : dir to hoist published artefacts into
#   $4  BV_POPULATIONS : comma-separated normalized country labels
#
# Scripts live next to this file (bundled at /opt/biovault/scripts/
# bv_paper_fst_island_aims in the image).

set -euo pipefail

RAW_DIR="$1"
WORK="$2"
RESULTS="$3"
export BV_POPULATIONS="$4"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}:${PYTHONPATH:-}"
export BV_RAW_DIR="${RAW_DIR}"

FST_WORK="${WORK}/fst"
AIMS_WORK="${WORK}/aims"
mkdir -p "${FST_WORK}" "${AIMS_WORK}" "${RESULTS}"

source /opt/conda/etc/profile.d/conda.sh
conda activate biovault_popgen

echo "== FST 01 load/merge =="
BV_WORK_DIR="${FST_WORK}" python3 "${SCRIPT_DIR}/fst_01_load_merge.py"
echo "== FST 02 compute_fst =="
BV_WORK_DIR="${FST_WORK}" python3 "${SCRIPT_DIR}/fst_02_compute_fst.py"
echo "== FST 03 visualize =="
BV_WORK_DIR="${FST_WORK}" python3 "${SCRIPT_DIR}/fst_03_visualize.py"

echo "== AIMs 04 merge =="
BV_WORK_DIR="${AIMS_WORK}" BV_FST_DIR="${FST_WORK}" python3 "${SCRIPT_DIR}/aims_04_merge.py"
echo "== AIMs 05 differential SNPs =="
BV_WORK_DIR="${AIMS_WORK}" python3 "${SCRIPT_DIR}/aims_05_diff_snps.py"
echo "== AIMs 06 AIMs panels =="
BV_WORK_DIR="${AIMS_WORK}" python3 "${SCRIPT_DIR}/aims_06_dendrogram.py"

echo "== hoisting artefacts -> ${RESULTS} =="
cp "${FST_WORK}/data/merged/merged_allele_freq_annotated.tsv" "${RESULTS}/merged_allele_freq_annotated.tsv"
cp "${FST_WORK}/data/fst/fst_matrix.tsv"                       "${RESULTS}/fst_matrix.tsv"
cp "${AIMS_WORK}/data/master_af_table.tsv"                     "${RESULTS}/master_af_table.tsv"
cp "${AIMS_WORK}/data/differential_snps/all_outliers_long.tsv" "${RESULTS}/all_outliers_long.tsv"
cp "${AIMS_WORK}/data/aims/aims_combined.tsv"                  "${RESULTS}/aims_combined.tsv"

# country_map.tsv: participant -> normalized country (written by main.nf as
# country_map.tsv in RAW_DIR's parent staging dir; copy if present).
[ -f "${RAW_DIR}/country_map.tsv" ] && cp "${RAW_DIR}/country_map.tsv" "${RESULTS}/country_map.tsv" || true

{
    echo "Populations: ${BV_POPULATIONS}"
    echo ""
    echo "=== FST matrix ==="
    cat "${FST_WORK}/data/fst/fst_matrix.tsv"
    echo ""
    echo "=== master_af_table summary ==="
    cat "${AIMS_WORK}/data/master_af_table_summary.txt"
} > "${RESULTS}/population_level_summary.txt"

# plots (optional)
for p in "${FST_WORK}"/plots/*.png "${AIMS_WORK}"/plots/*.png "${AIMS_WORK}"/plots/*.pdf; do
    [ -f "${p}" ] && cp "${p}" "${RESULTS}/" || true
done

echo "Pipeline complete."

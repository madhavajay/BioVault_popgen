#!/usr/bin/env bash
# Build and package HGP1K ADMIXTURE reference presets.
#
# Presets:
#   balanced5  300 unrelated founders from each 1KGP superpopulation
#              (AFR, AMR, EAS, EUR, SAS) -> data/hgp1k_1500_sex_bias/
#   all        all 3202 1KGP samples from those superpopulations, including
#              related samples -> data/hgp1k_all_sex_bias/
#
# Requires bcftools and plink2 on PATH. If local plink2 is unavailable, run this
# inside the admixture tools image, e.g.:
#   docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work \
#     --entrypoint bash biovault-admixture:1.4.0-amd64-tools \
#     flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_reference_package.sh balanced5

set -euo pipefail

PRESET="${1:-balanced5}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"

case "${PRESET}" in
  balanced5)
    PACKAGE_NAME="${PACKAGE_NAME:-hgp1k_1500_sex_bias}"
    SELECT_ARGS=(--superpops all --per-pop 300 --seed "${SELECTION_SEED:-42}")
    SELECTION_POLICY="${SELECTION_POLICY:-300 unrelated founders per superpopulation (AFR, AMR, EAS, EUR, SAS)}"
    ;;
  all)
    PACKAGE_NAME="${PACKAGE_NAME:-hgp1k_all_sex_bias}"
    SELECT_ARGS=(--superpops all --all-samples)
    SELECTION_POLICY="${SELECTION_POLICY:-all 1KGP samples from AFR, AMR, EAS, EUR, SAS, including related samples}"
    ;;
  *)
    echo "ERROR: unknown preset '${PRESET}' (expected balanced5 or all)" >&2
    exit 1
    ;;
esac

BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/.docker/reference/${PACKAGE_NAME}}"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data/${PACKAGE_NAME}}"
SAMPLES_TSV="${SAMPLES_TSV:-${BUILD_DIR}/reference_samples.tsv}"

mkdir -p "${BUILD_DIR}"

python3 "${HERE}/select_reference_samples.py" \
  "${SELECT_ARGS[@]}" \
  --out "${SAMPLES_TSV}"

SAMPLES_TSV="${SAMPLES_TSV}" \
OUT_DIR="${BUILD_DIR}" \
"${HERE}/build_reference_bed.sh"

PACKAGE_NAME="${PACKAGE_NAME}" \
SRC_DIR="${BUILD_DIR}" \
OUT_DIR="${DATA_DIR}" \
SAMPLES_TSV="${SAMPLES_TSV}" \
SELECTION_POLICY="${SELECTION_POLICY}" \
SELECTION_SEED="${SELECTION_SEED:-42}" \
"${HERE}/pack_reference.sh"

echo "[package] ${PRESET} -> ${DATA_DIR}"

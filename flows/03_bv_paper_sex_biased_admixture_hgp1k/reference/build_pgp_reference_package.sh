#!/usr/bin/env bash
# Build HGP1K ADMIXTURE reference packages from the PGP-filtered 1KGP VCFs.
#
# This is intentionally separate from build_admixture_docker.sh: the outputs are
# test-data references and are not baked into the runtime container unless you
# explicitly point a run at them.
#
# Presets mirror build_reference_package.sh:
#   balanced5 -> data/hgp1k_1500_sex_bias_pgp/
#   all       -> data/hgp1k_all_sex_bias_pgp/
#
# Requires bcftools and plink2 on PATH. If local plink2 is unavailable, run in
# the existing tools image:
#   docker run --rm --platform linux/amd64 -v "$PWD:/work" -w /work \
#     --entrypoint bash biovault-admixture:1.4.0-amd64-tools \
#     flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_pgp_reference_package.sh balanced5

set -euo pipefail

PRESET="${1:-balanced5}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"

case "${PRESET}" in
  balanced5)
    DEFAULT_PACKAGE="hgp1k_1500_sex_bias_pgp"
    ;;
  all)
    DEFAULT_PACKAGE="hgp1k_all_sex_bias_pgp"
    ;;
  *)
    echo "ERROR: unknown preset '${PRESET}' (expected balanced5 or all)" >&2
    exit 1
    ;;
esac

PACKAGE_NAME="${PACKAGE_NAME:-${DEFAULT_PACKAGE}}"
VCF_DIR="${VCF_DIR:-${REPO_ROOT}/data/1kgp_high_coverage/filtered_pgp}"
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/.docker/reference/${PACKAGE_NAME}}"
DATA_DIR="${DATA_DIR:-${REPO_ROOT}/data/${PACKAGE_NAME}}"

[ -d "${VCF_DIR}" ] || { echo "ERROR: missing PGP-filtered VCF dir: ${VCF_DIR}" >&2; exit 1; }

VCF_DIR="${VCF_DIR}" \
PACKAGE_NAME="${PACKAGE_NAME}" \
BUILD_DIR="${BUILD_DIR}" \
DATA_DIR="${DATA_DIR}" \
"${HERE}/build_reference_package.sh" "${PRESET}"

echo "[pgp-package] ${PRESET} -> ${DATA_DIR}"

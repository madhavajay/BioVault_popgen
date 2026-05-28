#!/usr/bin/env bash
# Build the project Docker image and run any command inside it.
#
# Usage:
#   ./scripts/run_in_docker.sh
#   ./scripts/run_in_docker.sh bash 03_individual_level/pca_qc_fast/scripts/run_pipeline.sh
#   ./scripts/run_in_docker.sh bash 03_individual_level/gnomad_projection/scripts/run_pipeline_gnomad.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/image_versions.sh"
IMAGE="${IMAGE:-${BIOVAULT_IMAGE}}"
PLATFORM="${PLATFORM:-linux/amd64}"

docker build \
  --platform "${PLATFORM}" \
  -f "${ROOT_DIR}/Dockerfile" \
  -t "${IMAGE}" \
  "${ROOT_DIR}"

if [ "$#" -eq 0 ]; then
  set -- bash
fi

TTY_ARGS=()
if [ -t 0 ] && [ -t 1 ]; then
  TTY_ARGS=(-it)
fi

docker run --rm ${TTY_ARGS+"${TTY_ARGS[@]}"} \
  --platform "${PLATFORM}" \
  -v "${ROOT_DIR}:/work" \
  -w /work \
  -e CONDA_ENV="biovault_popgen" \
  -e THREADS="${THREADS:-}" \
  -e PARALLEL_CHRS="${PARALLEL_CHRS:-2}" \
  -e DATA_DIR="${DATA_DIR:-/work/03_individual_level}" \
  -e LOADINGS_HT="${LOADINGS_HT:-/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht}" \
  -e LOADINGS_VARIANTS_TSV="${LOADINGS_VARIANTS_TSV:-/opt/biovault/reference/pca_loadings/loadings_variants.tsv}" \
  "${IMAGE}" \
  "$@"

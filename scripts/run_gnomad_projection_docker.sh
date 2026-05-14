#!/usr/bin/env bash
# Run the gnomAD projection pipeline in the shared project Docker environment.
#
# Usage:
#   ./scripts/run_gnomad_projection_docker.sh
#   THREADS=8 PARALLEL_CHRS=1 ./scripts/run_gnomad_projection_docker.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "${ROOT_DIR}/scripts/run_in_docker.sh" \
  bash 03_individual_level/gnomad_projection/scripts/run_pipeline_gnomad.sh

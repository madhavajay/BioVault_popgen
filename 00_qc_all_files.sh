#!/usr/bin/env bash
# Run genotype input preflight QC.
#
# Examples:
#   ./00_qc_all_files.sh 00_qc_test_data/files
#   ./00_qc_all_files.sh --samplesheet 00_qc_test_data/samplesheet.csv
#   ./00_qc_all_files.sh --docker --samplesheet 00_qc_test_data/samplesheet.csv

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/results/qc_all_files}"
IMAGE="${IMAGE:-ghcr.io/madhavajay/biovault-popgen:0.2.6-fast}"
USE_DOCKER=0

ARGS=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --docker)
            USE_DOCKER=1
            shift
            ;;
        --output-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

if [ "${USE_DOCKER}" = "1" ]; then
    DOCKER_OUT="${OUT_DIR}"
    case "${OUT_DIR}" in
        "${ROOT_DIR}"/*) DOCKER_OUT="/work/${OUT_DIR#"${ROOT_DIR}/"}" ;;
    esac
    docker run --rm \
        -v "${ROOT_DIR}:/work" \
        -w /work \
        "${IMAGE}" \
        python /opt/biovault/scripts/qc_all_files/qc_all_files.py \
            --output-dir "${DOCKER_OUT}" "${ARGS[@]}"
else
    if command -v uv >/dev/null 2>&1; then
        uv run --with pandas --with numpy python "${ROOT_DIR}/00_qc_all_files/scripts/qc_all_files.py" \
            --output-dir "${OUT_DIR}" "${ARGS[@]}"
    else
        python3 "${ROOT_DIR}/00_qc_all_files/scripts/qc_all_files.py" \
            --output-dir "${OUT_DIR}" "${ARGS[@]}"
    fi
fi

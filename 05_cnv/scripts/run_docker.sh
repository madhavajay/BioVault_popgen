#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${IMAGE:-biovault-cnv:local}"
PLATFORM="${PLATFORM:-linux/amd64}"

: "${IDAT_DIR:?Set IDAT_DIR to the directory containing IDAT/BPM/EGT inputs}"
OUT_DIR="${OUT_DIR:-${IDAT_DIR}/results/cnv}"

docker build \
  --platform "${PLATFORM}" \
  -f "${ROOT_DIR}/05_cnv/Dockerfile" \
  -t "${IMAGE}" \
  "${ROOT_DIR}/05_cnv"

TTY_ARGS=()
if [ -t 0 ] && [ -t 1 ]; then
  TTY_ARGS=(-it)
fi

DOCKER_ARGS=(
  --rm
  --platform "${PLATFORM}"
  -v "${IDAT_DIR}:/input:ro"
  -v "${OUT_DIR}:/out"
  -e "STAGE=${STAGE:-all}"
  -e "FINAL_REPORT=${FINAL_REPORT:-}"
  -e "IAAP_COMMAND=${IAAP_COMMAND:-}"
  -e "CALL_CHRX=${CALL_CHRX:-0}"
  -e "MIN_SNP=${MIN_SNP:-3}"
  -e "MIN_LENGTH=${MIN_LENGTH:-}"
  -e "MIN_CONF=${MIN_CONF:-}"
  -e "SNP_POS_FILE=${SNP_POS_FILE:-}"
)

if [ -n "${IAAP_CLI:-}" ]; then
  DOCKER_ARGS+=(-v "${IAAP_CLI}:/opt/illumina/bin/iaap-cli:ro")
fi

if [ -n "${SEX_FILE:-}" ]; then
  DOCKER_ARGS+=(-v "${SEX_FILE}:/input/sex_file.tsv:ro" -e "SEX_FILE=/input/sex_file.tsv")
fi

docker run "${TTY_ARGS[@]}" "${DOCKER_ARGS[@]}" "${IMAGE}"

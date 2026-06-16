#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${IMAGE:-biovault-gtc2vcf:local}"
PLATFORM="${PLATFORM:-linux/amd64}"

: "${IDAT_DIR:?Set IDAT_DIR to the directory containing IDAT/BPM/EGT inputs}"
OUT_DIR="${OUT_DIR:-${IDAT_DIR}/results/cnv/gtc2vcf}"
BPM="${BPM:-${IDAT_DIR}/HumanOmniExpress-12v1_H.bpm}"
EGT="${EGT:-${IDAT_DIR}/HumanOmniExpress-12v1_H.egt}"

docker build \
  --platform "${PLATFORM}" \
  -f "${ROOT_DIR}/05_cnv/Dockerfile.gtc2vcf" \
  -t "${IMAGE}" \
  "${ROOT_DIR}/05_cnv"

mkdir -p "${OUT_DIR}"

TTY_ARGS=()
if [ -t 0 ] && [ -t 1 ]; then
  TTY_ARGS=(-it)
fi

docker run --rm "${TTY_ARGS[@]}" \
  --platform "${PLATFORM}" \
  -v "${IDAT_DIR}:/input:ro" \
  -v "${OUT_DIR}:/out" \
  "${IMAGE}" \
  bcftools +idat2gtc \
    --bpm "/input/$(basename "${BPM}")" \
    --egt "/input/$(basename "${EGT}")" \
    --idats /input \
    --output /out/gtc

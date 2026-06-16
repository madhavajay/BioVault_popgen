#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${IMAGE:-biovault-gtc2vcf:local}"
PLATFORM="${PLATFORM:-linux/amd64}"

: "${IDAT_DIR:?Set IDAT_DIR to the directory containing IDAT/BPM/EGT inputs}"
GTC_DIR="${GTC_DIR:-${IDAT_DIR}/results/cnv/gtc2vcf/gtc}"
OUT_DIR="${OUT_DIR:-${IDAT_DIR}/results/cnv/gtc2vcf}"
BPM="${BPM:-${IDAT_DIR}/HumanOmniExpress-12v1_H.bpm}"
EGT="${EGT:-${IDAT_DIR}/HumanOmniExpress-12v1_H.egt}"

MODE="${MODE:-all}"
ONE_GTC="${ONE_GTC:-}"
FINAL_REPORT="${FINAL_REPORT:-${OUT_DIR}/final_report.tsv}"

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

case "${MODE}" in
  one)
    if [ -z "${ONE_GTC}" ]; then
      ONE_GTC="$(find "${GTC_DIR}" -type f -name '*.gtc' | sort | head -n 1)"
    fi
    if [ -z "${ONE_GTC}" ] || [ ! -f "${ONE_GTC}" ]; then
      echo "No GTC file found for one-sample conversion" >&2
      exit 1
    fi
    docker run --rm "${TTY_ARGS[@]}" \
      --platform "${PLATFORM}" \
      -v "${IDAT_DIR}:/input:ro" \
      -v "${GTC_DIR}:/gtc:ro" \
      -v "${OUT_DIR}:/out" \
      "${IMAGE}" \
      bcftools +gtc2vcf \
        -b "/input/$(basename "${BPM}")" \
        -e "/input/$(basename "${EGT}")" \
        -O t \
        -o /out/one_sample_genomestudio.tsv \
        "/gtc/$(basename "${ONE_GTC}")"
    ;;
  all)
    docker run --rm "${TTY_ARGS[@]}" \
      --platform "${PLATFORM}" \
      -v "${IDAT_DIR}:/input:ro" \
      -v "${GTC_DIR}:/gtc:ro" \
      -v "${OUT_DIR}:/out" \
      "${IMAGE}" \
      bcftools +gtc2vcf \
        -b "/input/$(basename "${BPM}")" \
        -e "/input/$(basename "${EGT}")" \
        --gtcs /gtc \
        -O t \
        -o "/out/$(basename "${FINAL_REPORT}")"
    ;;
  *)
    echo "Unknown MODE=${MODE}; expected MODE=one or MODE=all" >&2
    exit 1
    ;;
esac

#!/usr/bin/env bash
set -euo pipefail

log() {
  printf '[05_cnv] %s\n' "$*" >&2
}

die() {
  printf '[05_cnv] ERROR: %s\n' "$*" >&2
  exit 1
}

first_match() {
  local pattern="$1"
  find "${IDAT_DIR}" -maxdepth 1 -type f -name "${pattern}" | sort | head -n 1
}

IDAT_DIR="${IDAT_DIR:-/input}"
OUT_DIR="${OUT_DIR:-/out}"
STAGE="${STAGE:-all}"

[ -d "${IDAT_DIR}" ] || die "IDAT_DIR does not exist: ${IDAT_DIR}"
mkdir -p "${OUT_DIR}" "${OUT_DIR}/iaap" "${OUT_DIR}/penncnv" "${OUT_DIR}/logs"

BPM="${BPM:-$(first_match '*.bpm')}"
EGT="${EGT:-$(first_match '*.egt')}"
SAMPLE_MAP="${SAMPLE_MAP:-$(find "${IDAT_DIR}" -maxdepth 1 -type f \( -iname '*sample*id*.txt' -o -iname '*chip*id*.txt' \) | sort | head -n 1)}"
SAMPLE_MANIFEST="${SAMPLE_MANIFEST:-${OUT_DIR}/iaap/idat_manifest.tsv}"
FINAL_REPORT="${FINAL_REPORT:-${OUT_DIR}/iaap/final_report.txt}"

[ -n "${BPM}" ] && [ -f "${BPM}" ] || die "Could not find BPM file. Set BPM explicitly."
[ -n "${EGT}" ] && [ -f "${EGT}" ] || die "Could not find EGT file. Set EGT explicitly."

log "IDAT_DIR=${IDAT_DIR}"
log "OUT_DIR=${OUT_DIR}"
log "BPM=${BPM}"
log "EGT=${EGT}"
log "SAMPLE_MAP=${SAMPLE_MAP:-none}"
log "FINAL_REPORT=${FINAL_REPORT}"

if [ -n "${SAMPLE_MAP:-}" ] && [ -f "${SAMPLE_MAP}" ]; then
  perl /opt/biovault-cnv/scripts/prepare_idat_manifest.pl \
    --idat-dir "${IDAT_DIR}" \
    --sample-map "${SAMPLE_MAP}" \
    --out "${SAMPLE_MANIFEST}" \
    --missing-out "${OUT_DIR}/iaap/idat_manifest_missing.tsv"
else
  log "No sample map detected; writing manifest directly from IDAT pairs."
  perl /opt/biovault-cnv/scripts/prepare_idat_manifest.pl \
    --idat-dir "${IDAT_DIR}" \
    --out "${SAMPLE_MANIFEST}" \
    --missing-out "${OUT_DIR}/iaap/idat_manifest_missing.tsv"
fi

case "${STAGE}" in
  all|iaap)
    if [ -s "${FINAL_REPORT}" ]; then
      log "Final Report already exists; skipping IAAP export: ${FINAL_REPORT}"
    else
      /opt/biovault-cnv/scripts/run_iaap_export.sh \
        "${IDAT_DIR}" "${BPM}" "${EGT}" "${SAMPLE_MANIFEST}" "${OUT_DIR}/iaap" "${FINAL_REPORT}"
    fi
    ;;
  penncnv)
    log "Skipping IAAP stage because STAGE=penncnv."
    ;;
  *)
    die "Invalid STAGE=${STAGE}. Use all, iaap, or penncnv."
    ;;
esac

case "${STAGE}" in
  all|penncnv)
    [ -s "${FINAL_REPORT}" ] || die "Final Report not found: ${FINAL_REPORT}. IAAP gencall produces GTC files; convert them to PennCNV signal text with LRR/BAF first, or set FINAL_REPORT to a GenomeStudio/converted report."
    /opt/biovault-cnv/scripts/run_penncnv.sh "${FINAL_REPORT}" "${OUT_DIR}/penncnv"
    ;;
  iaap)
    log "Skipping PennCNV stage because STAGE=iaap."
    ;;
esac

log "complete"

#!/usr/bin/env bash
set -euo pipefail

FINAL_REPORT="$1"
OUT_DIR="$2"

log() {
  printf '[05_cnv:penncnv] %s\n' "$*" >&2
}

die() {
  printf '[05_cnv:penncnv] ERROR: %s\n' "$*" >&2
  exit 1
}

[ -s "${FINAL_REPORT}" ] || die "Final Report is missing or empty: ${FINAL_REPORT}"

PENNCNV_HOME="${PENNCNV_HOME:-/opt/PennCNV-1.0.5}"
HMM="${HMM:-${PENNCNV_HOME}/lib/hhall.hmm}"
MIN_SNP="${MIN_SNP:-3}"
MIN_LENGTH="${MIN_LENGTH:-}"
MIN_CONF="${MIN_CONF:-}"
CALL_CHRX="${CALL_CHRX:-0}"
SEX_FILE="${SEX_FILE:-}"
SNP_POS_FILE="${SNP_POS_FILE:-}"

SIGNAL_DIR="${SIGNAL_DIR:-${OUT_DIR}/signal}"
LISTFILE="${LISTFILE:-${OUT_DIR}/signal_files.list}"
PFB="${PFB:-${OUT_DIR}/cohort.pfb}"
RAW_CNV="${RAW_CNV:-${OUT_DIR}/cnv_calls.rawcnv}"
LOGFILE="${LOGFILE:-${OUT_DIR}/detect_cnv.log}"

mkdir -p "${OUT_DIR}" "${SIGNAL_DIR}"

if [ -z "${SNP_POS_FILE}" ]; then
  AUTO_SNP_POS_FILE="${OUT_DIR}/snppos.tsv"
  if perl /opt/biovault-cnv/scripts/extract_snppos_from_final_report.pl \
      --final-report "${FINAL_REPORT}" \
      --out "${AUTO_SNP_POS_FILE}" \
      > "${OUT_DIR}/extract_snppos.log" 2>&1; then
    if [ -s "${AUTO_SNP_POS_FILE}" ] && [ "$(wc -l < "${AUTO_SNP_POS_FILE}" | tr -d ' ')" -gt 1 ]; then
      SNP_POS_FILE="${AUTO_SNP_POS_FILE}"
      log "extracted SNP positions: ${SNP_POS_FILE}"
    fi
  fi
fi

log "splitting Illumina Final Report into per-sample signal files"
if perl /opt/biovault-cnv/scripts/split_gtc2vcf_genomestudio_report.pl \
    --final-report "${FINAL_REPORT}" \
    --out-dir "${SIGNAL_DIR}" \
    > "${OUT_DIR}/split_gtc2vcf_genomestudio_report.log" 2>&1; then
  log "split gtc2vcf wide GenomeStudio report"
else
  log "gtc2vcf wide splitter did not apply; falling back to PennCNV split_illumina_report.pl"
  SPLIT_REPORT="${FINAL_REPORT}"
  SPLIT_WRITER_PID=""
  SPLIT_FIFO=""
  if ! head -n 1000 "${FINAL_REPORT}" | grep -qiE '^\[Data\][[:space:]]*$'; then
    SPLIT_FIFO="${OUT_DIR}/final_report.with_data.fifo"
    rm -f "${SPLIT_FIFO}"
    mkfifo "${SPLIT_FIFO}"
    ( printf '[Data]\n'; cat "${FINAL_REPORT}" ) > "${SPLIT_FIFO}" &
    SPLIT_WRITER_PID="$!"
    SPLIT_REPORT="${SPLIT_FIFO}"
    log "Final Report has no [Data] marker; streaming a PennCNV-compatible header"
  fi

  if ! split_illumina_report.pl \
    -prefix "${SIGNAL_DIR}/" \
    -suffix ".txt" \
    "${SPLIT_REPORT}" \
    > "${OUT_DIR}/split_illumina_report.log" 2>&1; then
    if [ -n "${SPLIT_WRITER_PID}" ]; then
      wait "${SPLIT_WRITER_PID}" 2>/dev/null || true
    fi
    [ -n "${SPLIT_FIFO}" ] && rm -f "${SPLIT_FIFO}"
    tail -n 40 "${OUT_DIR}/split_illumina_report.log" >&2 || true
    die "Illumina Final Report splitting failed"
  fi

  if [ -n "${SPLIT_WRITER_PID}" ]; then
    wait "${SPLIT_WRITER_PID}" || die "Final Report stream failed"
  fi
  [ -n "${SPLIT_FIFO}" ] && rm -f "${SPLIT_FIFO}"
fi

find "${SIGNAL_DIR}" -maxdepth 1 -type f -name '*.txt' | sort > "${LISTFILE}"
sample_count="$(wc -l < "${LISTFILE}" | tr -d ' ')"
[ "${sample_count}" -gt 0 ] || die "No signal files were produced in ${SIGNAL_DIR}"
log "signal files: ${sample_count}"

compile_args=(-listfile "${LISTFILE}" -output "${PFB}")
if [ -n "${SNP_POS_FILE}" ]; then
  [ -s "${SNP_POS_FILE}" ] || die "SNP_POS_FILE is set but missing or empty: ${SNP_POS_FILE}"
  compile_args+=(-snpposfile "${SNP_POS_FILE}")
fi

log "compiling cohort PFB"
if ! compile_pfb.pl "${compile_args[@]}" > "${OUT_DIR}/compile_pfb.log" 2>&1; then
  tail -n 40 "${OUT_DIR}/compile_pfb.log" >&2 || true
  die "PFB compilation failed. Ensure the Final Report has Chr/Position columns or set SNP_POS_FILE."
fi

[ -s "${PFB}" ] || die "PFB was not created: ${PFB}"

detect_args=(
  -test
  -hmm "${HMM}"
  -pfb "${PFB}"
  -listfile "${LISTFILE}"
  -logfile "${LOGFILE}"
  -output "${RAW_CNV}"
  -minsnp "${MIN_SNP}"
  -tabout
  -confidence
)

if [ -n "${MIN_LENGTH}" ]; then
  detect_args+=(-minlength "${MIN_LENGTH}")
fi

if [ -n "${MIN_CONF}" ]; then
  detect_args+=(-minconf "${MIN_CONF}")
fi

log "calling autosomal CNVs"
detect_cnv.pl "${detect_args[@]}"

if [ ! -e "${RAW_CNV}" ]; then
  log "PennCNV completed without CNV calls; creating empty output ${RAW_CNV}"
  : > "${RAW_CNV}"
fi

filter_cnv.pl "${RAW_CNV}" > "${OUT_DIR}/cnv_calls.filteredcnv" 2> "${OUT_DIR}/filter_cnv.log" || true

if [ "${CALL_CHRX}" = "1" ]; then
  [ -n "${SEX_FILE}" ] && [ -s "${SEX_FILE}" ] || die "CALL_CHRX=1 requires SEX_FILE with two columns: signal filename and sex"
  log "calling chrX CNVs"
  detect_cnv.pl "${detect_args[@]}" \
    -chrx \
    -sexfile "${SEX_FILE}" \
    -output "${OUT_DIR}/cnv_calls.chrx.rawcnv" \
    -logfile "${OUT_DIR}/detect_cnv.chrx.log"
fi

{
  printf 'final_report\t%s\n' "${FINAL_REPORT}"
  printf 'signal_dir\t%s\n' "${SIGNAL_DIR}"
  printf 'signal_files\t%s\n' "${sample_count}"
  printf 'pfb\t%s\n' "${PFB}"
  printf 'hmm\t%s\n' "${HMM}"
  printf 'snp_pos_file\t%s\n' "${SNP_POS_FILE}"
  printf 'raw_cnv\t%s\n' "${RAW_CNV}"
  printf 'min_snp\t%s\n' "${MIN_SNP}"
  printf 'min_length\t%s\n' "${MIN_LENGTH}"
  printf 'min_conf\t%s\n' "${MIN_CONF}"
} > "${OUT_DIR}/cnv_run_summary.tsv"

log "wrote ${RAW_CNV}"

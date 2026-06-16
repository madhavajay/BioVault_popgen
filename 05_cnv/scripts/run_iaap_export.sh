#!/usr/bin/env bash
set -euo pipefail

IDAT_DIR="$1"
BPM="$2"
EGT="$3"
SAMPLE_MANIFEST="$4"
OUT_DIR="$5"
FINAL_REPORT="$6"

log() {
  printf '[05_cnv:iaap] %s\n' "$*" >&2
}

die() {
  printf '[05_cnv:iaap] ERROR: %s\n' "$*" >&2
  exit 1
}

mkdir -p "${OUT_DIR}"

IAAP_BIN="${IAAP_BIN:-/opt/illumina/bin/iaap-cli}"
IAAP_COMMAND="${IAAP_COMMAND:-}"

write_template() {
  cat > "${OUT_DIR}/IAAP_COMMAND.template.txt" <<EOF
Set IAAP_COMMAND to override the default Illumina IAAP CLI command if needed.

Default IAAP command for GTC generation:

  CLR_ICU_VERSION_OVERRIDE="\$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" {IAAP_BIN} \\
    gencall \\
    "{BPM}" \\
    "{EGT}" \\
    "{OUT_DIR}/gtc" \\
    --idat-folder "{IDAT_DIR}" \\
    --output-gtc \\
    --gender-estimate-call-rate-threshold 0.0

Available placeholders:
  {IAAP_BIN}          ${IAAP_BIN}
  {IDAT_DIR}          ${IDAT_DIR}
  {BPM}               ${BPM}
  {EGT}               ${EGT}
  {SAMPLE_MANIFEST}   ${SAMPLE_MANIFEST}
  {OUT_DIR}           ${OUT_DIR}
  {FINAL_REPORT}      ${FINAL_REPORT}

PennCNV still needs signal text with LRR/BAF. IAAP gencall writes GTC files;
convert those with bcftools +gtc2vcf or GenomeStudio before running PennCNV.
If you choose to have IAAP_COMMAND write a Final Report, it must contain:
  SNP Name
  Chr
  Position
  Sample ID
  B Allele Freq
  Log R Ratio
EOF
}

if [ -z "${IAAP_COMMAND}" ]; then
  write_template
  IAAP_COMMAND='export LANG="en_US.UTF-8"; if command -v uconv >/dev/null 2>&1; then export CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed '\''s/.* //g'\'')"; fi; mkdir -p "{OUT_DIR}/gtc"; "{IAAP_BIN}" gencall "{BPM}" "{EGT}" "{OUT_DIR}/gtc" --idat-folder "{IDAT_DIR}" --output-gtc --gender-estimate-call-rate-threshold 0.0'
fi

if [ ! -x "${IAAP_BIN}" ] && ! command -v "${IAAP_BIN}" >/dev/null 2>&1; then
  die "IAAP_COMMAND is not set. Wrote template: ${OUT_DIR}/IAAP_COMMAND.template.txt"
fi

cmd="${IAAP_COMMAND}"
cmd="${cmd//\{IAAP_BIN\}/${IAAP_BIN}}"
cmd="${cmd//\{IDAT_DIR\}/${IDAT_DIR}}"
cmd="${cmd//\{BPM\}/${BPM}}"
cmd="${cmd//\{EGT\}/${EGT}}"
cmd="${cmd//\{SAMPLE_MANIFEST\}/${SAMPLE_MANIFEST}}"
cmd="${cmd//\{OUT_DIR\}/${OUT_DIR}}"
cmd="${cmd//\{FINAL_REPORT\}/${FINAL_REPORT}}"

log "running IAAP command"
printf '%s\n' "${cmd}" > "${OUT_DIR}/iaap_command.expanded.sh"
bash -lc "${cmd}" 2>&1 | tee "${OUT_DIR}/iaap.log"

if [ -s "${FINAL_REPORT}" ]; then
  log "wrote ${FINAL_REPORT}"
elif find "${OUT_DIR}/gtc" -maxdepth 1 -type f -name '*.gtc' 2>/dev/null | grep -q .; then
  log "wrote GTC files under ${OUT_DIR}/gtc"
else
  die "IAAP command completed but produced neither Final Report nor GTC files"
fi

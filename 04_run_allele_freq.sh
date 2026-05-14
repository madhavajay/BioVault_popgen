#!/usr/bin/env bash
# Compute per-island allele frequencies using biosynth's bvs tool
# (emit-long + aggregate-long), driven by an island_mapping.tsv.
#
# Stage 1 (per participant): bvs emit-long  <txt> -> <pid>.bvlr  (idempotent)
# Stage 2 (per island):      bvs aggregate-long  <subset of .bvlr> -> allele_freq_<label>.tsv
#
# Usage:
#   bash 04_run_allele_freq.sh                       # use all participants in mapping
#   bash 04_run_allele_freq.sh --limit 50            # cap per-island participants
#
# Overrides (env):
#   MAPPING        island_mapping.tsv path
#   DATA_DIR       dir containing <pid>/*.txt
#   OUT_DIR        where allele_freq_<label>.tsv files land
#   BVLR_DIR       where per-participant .bvlr files cache
#   IMAGE          biosynth docker image
#   PARALLEL       emit-long fan-out (default 8)
#
# Output filenames match what 04_population_level/fst_islands/scripts/
# 01_load_merge.py expects: BVI / TT keep case; the rest lowercase.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPPING="${MAPPING:-${ROOT_DIR}/01_mock_data_generation/output/island_mapping.tsv}"
DATA_DIR="${DATA_DIR:-${ROOT_DIR}/01_mock_data_generation/output}"
OUT_DIR="${OUT_DIR:-${ROOT_DIR}/04_population_level/raw_allele_freq_country}"
BVLR_DIR="${BVLR_DIR:-${ROOT_DIR}/04_population_level/.bvlr_cache}"
IMAGE="${IMAGE:-ghcr.io/openmined/biosynth:latest}"
PARALLEL="${PARALLEL:-8}"

LIMIT=""
while [ $# -gt 0 ]; do
    case "$1" in
        --limit) LIMIT="$2"; shift 2 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
    esac
done

[ -f "${MAPPING}" ] || { echo "ERROR: mapping not found: ${MAPPING}" >&2; exit 1; }

mkdir -p "${OUT_DIR}" "${BVLR_DIR}"

ISLANDS=(BVI TT Bahamas Barbados Bermuda StLucia)
island_label() {
    case "$1" in
        BVI)      echo BVI ;;
        TT)       echo TT ;;
        Bahamas)  echo bahamas ;;
        Barbados) echo barbados ;;
        Bermuda)  echo bermuda ;;
        StLucia)  echo stlucia ;;
        *)        return 1 ;;
    esac
}

# ── Stage 1: emit-long per participant ────────────────────────────────────────
TODO=$(mktemp)
trap 'rm -f "${TODO}"' EXIT

awk -F'\t' -v limit="${LIMIT:-0}" '
    NR == 1 { next }                                          # header
    limit > 0 {
        c[$2]++
        if (c[$2] > limit) next
    }
    { print }
' "${MAPPING}" > "${TODO}.mapping"

cut -f1 "${TODO}.mapping" | while read -r pid; do
    BVLR="${BVLR_DIR}/${pid}.bvlr"
    [ -s "${BVLR}" ] && continue
    TXT="$(find "${DATA_DIR}/${pid}" -maxdepth 1 -name '*.txt' -print -quit 2>/dev/null || true)"
    [ -n "${TXT}" ] || { echo "skip ${pid}: no .txt" >&2; continue; }
    printf '%s\t%s\n' "${pid}" "${TXT}"
done > "${TODO}"

N_TODO=$(wc -l < "${TODO}" | tr -d ' ')
echo "Stage 1: emit-long pending: ${N_TODO} (cached: $(( $(wc -l < "${TODO}.mapping") - N_TODO )))"

if [ "${N_TODO}" -gt 0 ]; then
    < "${TODO}" xargs -L 1 -P "${PARALLEL}" \
      sh -c '
        pid=$1; txt=$2
        docker run --rm --platform linux/amd64 \
            --entrypoint "" \
            -v "'"${ROOT_DIR}"':'"${ROOT_DIR}"'" \
            -w "'"${ROOT_DIR}"'" \
            "'"${IMAGE}"'" \
            bvs emit-long \
                --input "$txt" \
                --output "'"${BVLR_DIR}"'/$pid.bvlr" \
                --participant "$pid" >/dev/null
        ' _
fi

# ── Stage 2: aggregate per island ─────────────────────────────────────────────
echo "Stage 2: aggregate per island ..."
for ISLAND in "${ISLANDS[@]}"; do
    LABEL="$(island_label "${ISLAND}")"
    LIST="${BVLR_DIR}/.list_${ISLAND}.txt"

    awk -F'\t' -v island="${ISLAND}" '$2 == island { print $1 }' "${TODO}.mapping" \
      | while read -r pid; do
          BVLR="${BVLR_DIR}/${pid}.bvlr"
          [ -s "${BVLR}" ] && printf '%s\n' "${BVLR}"
        done > "${LIST}"

    COUNT=$(wc -l < "${LIST}" | tr -d ' ')
    if [ "${COUNT}" -eq 0 ]; then
        echo "  ${ISLAND}: 0 .bvlr files, skipping"
        continue
    fi

    OUT="${OUT_DIR}/allele_freq_${LABEL}.tsv"
    echo "  ${ISLAND} (${COUNT} participants) -> ${OUT}"
    docker run --rm --platform linux/amd64 \
        --entrypoint "" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${ROOT_DIR}" \
        "${IMAGE}" \
        bvs aggregate-long \
            --input-list "${LIST}" \
            --allele-freq-tsv "${OUT}" >/dev/null
done

rm -f "${TODO}.mapping"

echo
echo "=== output ==="
for f in "${OUT_DIR}"/allele_freq_*.tsv; do
    [ -f "${f}" ] || continue
    printf '  %s  (%s, %s rows)\n' "${f}" "$(du -h "${f}" | cut -f1)" "$(( $(wc -l < "${f}") - 1 ))"
done

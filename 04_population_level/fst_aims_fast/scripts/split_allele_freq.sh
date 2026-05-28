#!/usr/bin/env bash
# split_allele_freq.sh - per-country allele-frequency split (biosynth container)
#
# Runs inside the pinned BIOSYNTH_IMAGE. Mirrors what the repo-level
# 04_run_allele_freq.sh does, but driven by the `country` participant facet
# instead of a static island_mapping.tsv:
#
#   Stage 1 (per participant): bvs emit-long  <genotype> -> <pid>.bvlr
#   Stage 2 (per country):     bvs aggregate-long  <country .bvlr files>
#                                 -> allele_freq_<country_norm>.tsv
#
# Args:
#   $1  mapping TSV: participant_id <TAB> country_norm <TAB> genotype_filename
#   $2  dir containing the staged genotype files (flat, by filename)
#   $3  output dir for allele_freq_<country_norm>.tsv
#
# Fail-loud: any country whose participants produce no usable .bvlr aborts the
# whole step (the downstream popgen scripts also re-check, belt and braces).

set -euo pipefail

MAPPING="$1"
GENO_DIR="$2"
OUT_DIR="$3"
BVLR_DIR="$(pwd)/bvlr"

mkdir -p "${OUT_DIR}" "${BVLR_DIR}"

[ -s "${MAPPING}" ] || { echo "ERROR: empty mapping: ${MAPPING}" >&2; exit 1; }

echo "Stage 1: emit-long per participant ..."
while IFS=$'\t' read -r pid country fname; do
    [ -n "${pid}" ] || continue
    src="${GENO_DIR}/${fname}"
    [ -s "${src}" ] || { echo "ERROR: missing genotype for ${pid}: ${src}" >&2; exit 1; }
    bvs emit-long --input "${src}" --output "${BVLR_DIR}/${pid}.bvlr" --participant "${pid}" >/dev/null
done < "${MAPPING}"

echo "Stage 2: aggregate-long per country ..."
COUNTRIES="$(cut -f2 "${MAPPING}" | sort -u)"
FAIL=0
for country in ${COUNTRIES}; do
    CDIR="$(pwd)/agg_${country}"
    mkdir -p "${CDIR}"
    n=0
    while IFS=$'\t' read -r pid c _; do
        [ "${c}" = "${country}" ] || continue
        b="${BVLR_DIR}/${pid}.bvlr"
        if [ -s "${b}" ]; then
            cp "${b}" "${CDIR}/"
            n=$((n + 1))
        fi
    done < "${MAPPING}"

    if [ "${n}" -eq 0 ]; then
        echo "ERROR: country '${country}' produced 0 usable .bvlr files" >&2
        FAIL=1
        continue
    fi

    OUT="${OUT_DIR}/allele_freq_${country}.tsv"
    echo "  ${country} (${n} participants) -> ${OUT}"
    bvs aggregate-long \
        --input "${CDIR}" \
        --matrix-tsv "$(pwd)/matrix_${country}.tsv" \
        --allele-freq-tsv "${OUT}" >/dev/null
    [ -s "${OUT}" ] || { echo "ERROR: empty AF output for ${country}: ${OUT}" >&2; FAIL=1; }
done

if [ "${FAIL}" -ne 0 ]; then
    echo "ERROR: per-country allele-frequency split failed (see above)" >&2
    exit 1
fi

echo "=== output ==="
for f in "${OUT_DIR}"/allele_freq_*.tsv; do
    [ -f "${f}" ] || continue
    printf '  %s  (%s rows)\n' "${f}" "$(( $(wc -l < "${f}") - 1 ))"
done

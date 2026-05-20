#!/usr/bin/env bash
# Step 0: generate a synthetic cohort with deliberately injected,
# known-truth biases (BIOVAULT.md). Self-contained chain:
#
#   1. one Illumina-format reference panel (cached) — gives the rsid /
#      allele-pair / chr / pos universe for block selection.
#   2. make_bias_blocks.py  -> .cache/biasgen/bias_blocks.json
#        island_structure / projection / sex (X-hemizygosity) / singletons
#   3. generate_biased_cohort.py -> 01_mock_data_generation/output/
#        per-participant DDNA+Illumina dirs + cohort_spec.tsv +
#        island_mapping.tsv + sex_mapping.tsv
#
# Usage:
#   bash 01_generate_biased.sh                 # 100 participants, seed 100
#   bash 01_generate_biased.sh --count 1000    # scale up
#   bash 01_generate_biased.sh --count 100 --seed 7
#   bash 01_generate_biased.sh --slow          # accepted; no fast variant
#
# bvs: host (0.1.22, == pinned container) by default; fast, no docker.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COUNT=100
SEED=100
while [ $# -gt 0 ]; do
    case "$1" in
        --count) COUNT="$2"; shift 2 ;;
        --seed)  SEED="$2";  shift 2 ;;
        --slow)  shift ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 2 ;;
    esac
done

command -v bvs >/dev/null 2>&1 || {
    echo "ERROR: host bvs not found (need 0.1.22). Install or run via docker." >&2
    exit 1
}

BIASGEN="${ROOT_DIR}/.cache/biasgen"
REF_DIR="${BIASGEN}/REF"
BLOCKS="${BIASGEN}/bias_blocks.json"
LOADINGS_TSV="${ROOT_DIR}/.docker/reference/pca_loadings/loadings_variants.tsv"
SCRIPTS="${ROOT_DIR}/01_mock_data_generation/scripts"
OUT_DIR="${ROOT_DIR}/01_mock_data_generation/output"

[ -s "${LOADINGS_TSV}" ] || { echo "ERROR: missing ${LOADINGS_TSV} (run build_docker.sh / populate the PCA cache)" >&2; exit 1; }

mkdir -p "${REF_DIR}"

# 1. Illumina reference panel (cached, one file is enough).
REF_FILE="$(find "${REF_DIR}" -name '*.txt' 2>/dev/null | head -1 || true)"
if [ -z "${REF_FILE}" ]; then
    echo "[0/3] generating Illumina reference panel (one-time, cached) ..."
    bvs synthetic --format illumina \
        --output "${REF_DIR}/{id}_GSAv3.txt" \
        --count 1 --seed 1 --alt-frequency 0.0
    REF_FILE="$(find "${REF_DIR}" -name '*.txt' | head -1)"
else
    echo "[0/3] using cached Illumina reference panel: $(basename "${REF_FILE}")"
fi

# 2. Bias blocks (deterministic for a given seed).
echo "[1/3] selecting bias blocks (seed=${SEED}) ..."
python3 "${SCRIPTS}/make_bias_blocks.py" \
    --illumina-ref "${REF_FILE}" \
    --loadings-tsv "${LOADINGS_TSV}" \
    --out "${BLOCKS}" \
    --seed "${SEED}" | tail -3

# 3. Cohort.
echo "[2/3] generating ${COUNT} biased participants (seed=${SEED}) ..."
rm -rf "${OUT_DIR}"
python3 "${SCRIPTS}/generate_biased_cohort.py" \
    --count "${COUNT}" \
    --bias-blocks "${BLOCKS}" \
    --out-dir "${OUT_DIR}" \
    --seed "${SEED}"

echo
echo "[3/3] done. cohort -> ${OUT_DIR}"
echo "  cohort_spec.tsv / island_mapping.tsv / sex_mapping.tsv written"

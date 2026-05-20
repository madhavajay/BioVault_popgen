#!/usr/bin/env bash
# Generate synthetic GSA genotypes for N participants and assign each to a
# Caribbean island. Wipes the existing mock output dir before regenerating.
#
# Usage:
#   bash 01_mock_data_generation.sh                       # default count
#   bash 01_mock_data_generation.sh --count 1000
#   bash 01_mock_data_generation.sh --count 500 --seed 42 --min 50
#   bash 01_mock_data_generation.sh --count 1000 --no-clean   # keep existing dirs
#
# Flags:
#   --count N       Number of synthetic participants (default 1000)
#   --seed N        Biosynth + island-split RNG seed (default 100)
#   --min N         Minimum participants per island (default 100)
#   --threads N     Biosynth thread count (default: host CPU count)
#   --no-clean      Don't wipe ${OUTPUT_DIR} before regenerating

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${ROOT_DIR}/01_mock_data_generation"
OUTPUT_DIR="${BASE_DIR}/output"
GENERATE_SCRIPT="${BASE_DIR}/scripts/generate_mock_genotypes.sh"

COUNT=1000
SEED=100
MIN_PER_ISLAND=100
THREADS=""
CLEAN=1

while [ $# -gt 0 ]; do
    case "$1" in
        --count)     COUNT="$2"; shift 2 ;;
        --seed)      SEED="$2"; shift 2 ;;
        --min)       MIN_PER_ISLAND="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        --no-clean)  CLEAN=0; shift ;;
        -h|--help)
            sed -n '2,18p' "$0"
            exit 0
            ;;
        *)
            echo "ERROR: unknown arg: $1" >&2
            exit 2
            ;;
    esac
done

if ! [[ "${COUNT}" =~ ^[0-9]+$ ]] || [ "${COUNT}" -lt 1 ]; then
    echo "ERROR: --count must be a positive integer (got '${COUNT}')" >&2
    exit 2
fi

if [ "${CLEAN}" = "1" ] && [ -d "${OUTPUT_DIR}" ]; then
    echo "Cleaning ${OUTPUT_DIR} ..."
    rm -rf "${OUTPUT_DIR}"
fi

export COUNT SEED MIN_PER_ISLAND
[ -n "${THREADS}" ] && export THREADS

bash "${GENERATE_SCRIPT}"

echo
echo "=== mock data ready ==="
echo "  participants:    $(find "${OUTPUT_DIR}" -maxdepth 1 -mindepth 1 -type d -name '[0-9]*' | wc -l | tr -d ' ')"
echo "  island mapping:  ${OUTPUT_DIR}/island_mapping.tsv"
if [ -f "${OUTPUT_DIR}/island_mapping.tsv" ]; then
    echo "  island counts:"
    awk -F'\t' 'NR>1 {c[$2]++} END {for (k in c) printf "    %-10s %d\n", k, c[k]}' \
        "${OUTPUT_DIR}/island_mapping.tsv" | sort
fi

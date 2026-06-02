#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="${1:?missing data_dir}"
WORKING="${2:?missing working_dir}"
OUT_DIR="${3:?missing output_dir}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

mkdir -p "${WORKING}" "${OUT_DIR}"

if [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    source /opt/conda/etc/profile.d/conda.sh
    conda activate biovault_popgen
fi

VCF_DIR="${HGP1K_VCF_DIR:-${REPO_ROOT}/data/1kgp_high_coverage/filtered}"
LOCUS_MAP="${LOCUS_MAP:-${REPO_ROOT}/tools/locus_map.tsv}"
CHROMOSOMES="${CHROMOSOMES:-all}"
THREADS="${THREADS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"
MIN_GS="${MIN_GS:-0.15}"
MIN_AF="${MIN_AF:-0.01}"
MAX_REF_MISSING="${MAX_REF_MISSING:-0.05}"
N_COMPONENTS="${N_COMPONENTS:-10}"

metadata_args=()
if [ -n "${HGP1K_METADATA_TSV:-}" ]; then
    metadata_args=(--metadata-tsv "${HGP1K_METADATA_TSV}")
fi

matrix_args=()
if [ -n "${HGP1K_MATRIX_NPZ:-}" ]; then
    matrix_args=(--matrix-npz "${HGP1K_MATRIX_NPZ}")
fi

max_variant_args=()
if [ -n "${MAX_VARIANTS:-}" ]; then
    max_variant_args=(--max-variants "${MAX_VARIANTS}")
fi

"${PYTHON:-python3}" "${SCRIPT_DIR}/hgp1k_pipeline.py" \
    "${DATA_DIR}" \
    "${OUT_DIR}" \
    --vcf-dir "${VCF_DIR}" \
    --locus-map "${LOCUS_MAP}" \
    --chromosomes "${CHROMOSOMES}" \
    --n-components "${N_COMPONENTS}" \
    --min-gs "${MIN_GS}" \
    --min-af "${MIN_AF}" \
    --max-ref-missing "${MAX_REF_MISSING}" \
    --chunk-size "${CHUNK_SIZE:-5000}" \
    "${metadata_args[@]}" \
    "${matrix_args[@]}" \
    "${max_variant_args[@]}"

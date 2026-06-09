#!/usr/bin/env bash
# Generate a synthetic sex-biased admixture test cohort, runnable straight
# through the flow (flow.yaml / desktop app) to verify the pipeline recovers a
# KNOWN planted signal: autosomal AFR fraction vs X AFR fraction.
#
# Usage:
#   ./generate_testset.sh [OUT_DIR] [N_SAMPLES]
# Env knobs:
#   AFR_AUTO=0.5  AFR_X=0.85   the planted ancestry (expect AFR Δ ≈ AFR_X-AFR_AUTO)
#   N_AUTO_SNPS=60000  FEMALE_FRAC=0.5  SEED=42
#   WORKERS=8  write participant files in parallel
#   TEMPLATE_DDNA=/path/to/full_ddna.txt  write full-size files and overlay the
#     planted sex-biased SNPs at matching chr:pos rows
#   IMAGE=ghcr.io/madhavajay/biovault-admixture:0.2.5-fast
#
# Needs the baked reference at .docker/reference/hgp1k_admixture (reassemble via
# flows/.../reference/reassemble_reference.sh if missing). Writes <out>/
# samplesheet.csv + per-participant genotype dirs + GROUND_TRUTH.txt.

set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_DIR="${1:-${ROOT_DIR}/testdata_sexbias}"
N="${2:-100}"
IMAGE="${IMAGE:-ghcr.io/madhavajay/biovault-admixture:0.2.5-fast}"
REF="${ROOT_DIR}/.docker/reference/hgp1k_admixture"

if [ ! -s "${REF}/reference_auto.bed" ]; then
  echo "Reference missing — reassembling from committed shards..."
  bash "${ROOT_DIR}/flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/reassemble_reference.sh"
fi
mkdir -p "${OUT_DIR}"

DOCKER_MOUNTS=(-v "${ROOT_DIR}:/work" -v "${OUT_DIR}:/out")
TEMPLATE_ARG=""
if [ -n "${TEMPLATE_DDNA:-}" ]; then
  if [ ! -f "${TEMPLATE_DDNA}" ]; then
    echo "ERROR: TEMPLATE_DDNA does not exist: ${TEMPLATE_DDNA}" >&2
    exit 1
  fi
  TEMPLATE_DIR="$(cd "$(dirname "${TEMPLATE_DDNA}")" && pwd)"
  TEMPLATE_FILE="$(basename "${TEMPLATE_DDNA}")"
  DOCKER_MOUNTS+=(-v "${TEMPLATE_DIR}:/template:ro")
  TEMPLATE_ARG="--template-ddna /template/${TEMPLATE_FILE}"
fi

docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" \
  "${DOCKER_MOUNTS[@]}" --entrypoint bash "${IMAGE}" -lc "
    export HOME=/tmp MPLCONFIGDIR=/tmp
    python3 /work/flows/03_bv_paper_sex_biased_admixture_hgp1k/testdata/generate_sex_biased_testset.py \
      --reference-dir /work/.docker/reference/hgp1k_admixture \
      --locus-map /work/tools/locus_map.tsv \
      --out-dir /out --n-samples '${N}' \
      --n-auto-snps '${N_AUTO_SNPS:-60000}' \
      --afr-auto '${AFR_AUTO:-0.5}' --afr-x '${AFR_X:-0.85}' \
      --female-frac '${FEMALE_FRAC:-0.5}' --seed '${SEED:-42}' \
      --workers '${WORKERS:-1}' \
      ${TEMPLATE_ARG}
  "

# The generator runs inside the container with --out-dir /out, so genotype_file
# paths in the samplesheet are container paths. Rewrite them to the host OUT_DIR
# so the samplesheet imports/resolves on the host.
if [ -f "${OUT_DIR}/samplesheet.csv" ]; then
  tmp="${OUT_DIR}/samplesheet.csv.tmp"
  sed "s|/out/|${OUT_DIR}/|g" "${OUT_DIR}/samplesheet.csv" > "${tmp}" && mv "${tmp}" "${OUT_DIR}/samplesheet.csv"
fi

echo
echo "Test set -> ${OUT_DIR}/samplesheet.csv"
cat "${OUT_DIR}/GROUND_TRUTH.txt"

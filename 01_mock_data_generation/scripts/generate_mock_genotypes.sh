#!/usr/bin/env bash
# generate_mock_genotypes.sh
# ---------------------------
# Generate synthetic GSA microarray genotype files using OpenMined biosynth,
# then assign each participant to a Caribbean island.
#
# Requirements: Docker (linux/amd64 platform), python3 on host.
#
# Usage:        bash generate_mock_genotypes.sh
#
# Overrides:
#   COUNT=100 bash generate_mock_genotypes.sh           # number of participants (default 1000)
#   SEED=42 bash generate_mock_genotypes.sh             # biosynth + island-split RNG seed
#   THREADS=8 bash generate_mock_genotypes.sh           # biosynth thread count
#   MIN_PER_ISLAND=80 bash generate_mock_genotypes.sh   # island-split floor
#
# Output:
#   ./output/{id}/{id}_X_X_GSAv3-DTC_GRCh38-...txt       (one per participant)
#   ./output/vcf/...vcf.gz                               (after VCF conversion)
#   ./output/island_mapping.tsv                          (participant_id → island)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "${SCRIPT_DIR}")"
cd "${BASE_DIR}"

COUNT="${COUNT:-1000}"
SEED="${SEED:-100}"
THREADS="${THREADS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)}"
MIN_PER_ISLAND="${MIN_PER_ISLAND:-100}"
PYTHON="${PYTHON:-python3}"

echo "Generating ${COUNT} synthetic genotypes (seed=${SEED}, threads=${THREADS}) ..."
docker run --platform linux/amd64 --rm \
  -v "${BASE_DIR}:/work" \
  -w /work \
  ghcr.io/openmined/biosynth:0.1.22 \
  synthetic \
    --output "output/{id}/{id}_X_X_GSAv3-DTC_GRCh38-{month}-{day}-{year}.txt" \
    --count "${COUNT}" \
    --threads "${THREADS}" \
    --alt-frequency 0.5 \
    --seed "${SEED}"

echo "Converting genotypes to VCF ..."
docker run --platform linux/amd64 --rm \
  -v "${BASE_DIR}:/work" \
  -w /work \
  ghcr.io/openmined/biosynth:0.1.22 \
  genotype-to-vcf \
    --input output \
    --outdir output/vcf \
    --gzip

echo "Assigning participants to Caribbean islands ..."
"${PYTHON}" "${SCRIPT_DIR}/make_island_mapping.py" \
    --seed "${SEED}" \
    --min "${MIN_PER_ISLAND}" \
    "${BASE_DIR}/output"

echo
echo "Done."
echo "  Genotypes:       ${BASE_DIR}/output/{id}/"
echo "  VCFs:            ${BASE_DIR}/output/vcf/"
echo "  Island mapping:  ${BASE_DIR}/output/island_mapping.tsv"

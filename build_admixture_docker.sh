#!/usr/bin/env bash
# Build a local amd64 ADMIXTURE 1.4 image for the sex-biased admixture prototype.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="${VERSION:-1.4.0}"
IMAGE_NAME="${IMAGE_NAME:-ghcr.io/madhavajay/biovault-admixture}"
PLATFORM="${PLATFORM:-linux/amd64}"
LOCAL_TAG="${LOCAL_TAG:-biovault-admixture:${VERSION}-amd64}"
PUSH="${PUSH:-0}"

IMAGE_VERSIONED="${IMAGE_NAME}:${VERSION}-amd64"
IMAGE_LATEST="${IMAGE_NAME}:latest-amd64"

REF_DIR="${ROOT_DIR}/.docker/reference/hgp1k_admixture"
FORCE_REFERENCE="${FORCE_REFERENCE:-0}"

# 1) Build the `tools` stage (plink/plink2/bcftools/admixture + scripts), used
#    both to run the reference pre-gen and as the base of the final image.
docker build \
  --platform "${PLATFORM}" \
  --build-arg ADMIXTURE_VERSION="${VERSION}" \
  --target tools \
  -f "${ROOT_DIR}/docker/admixture/Dockerfile" \
  -t "${LOCAL_TAG}-tools" \
  "${ROOT_DIR}"

# 2) Obtain the baked HGP1K ADMIXTURE reference BED. Preference order:
#      a. reuse an existing BED in .docker/
#      b. reassemble from the committed shards (data/hgp1k_900_sex_bias) — works
#         in CI without the raw 1KGP VCFs
#      c. regenerate from the filtered 1KGP VCFs (needs data/ + bcftools+plink2)
REF_GEN="${ROOT_DIR}/flows/03_bv_paper_sex_biased_admixture_hgp1k/reference"
SHARD_DIR="${ROOT_DIR}/data/hgp1k_900_sex_bias"
if [ "${FORCE_REFERENCE}" = "1" ] || [ ! -s "${REF_DIR}/reference_auto.bed" ] || [ ! -s "${REF_DIR}/reference_x.bed" ]; then
  if [ "${FORCE_REFERENCE}" != "1" ] && ls "${SHARD_DIR}/hgp1k_900_sex_bias.tar.gz".* >/dev/null 2>&1; then
    echo "==> Reassembling HGP1K ADMIXTURE reference BED from committed shards…"
    bash "${REF_GEN}/reassemble_reference.sh"
  else
    echo "==> Generating HGP1K ADMIXTURE reference BED (reads data/1kgp_high_coverage/filtered)…"
    docker run --rm \
      --platform "${PLATFORM}" \
      -v "${ROOT_DIR}:/work" -w /work \
      -e CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}" \
      -e JOBS="${REF_JOBS:-4}" -e THREADS="${REF_THREADS:-4}" \
      --entrypoint bash \
      "${LOCAL_TAG}-tools" \
      flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_reference_bed.sh
  fi
else
  echo "==> Reusing existing reference BED at ${REF_DIR} (FORCE_REFERENCE=1 to rebuild)"
fi

[ -s "${REF_DIR}/reference_auto.bed" ] || { echo "ERROR: reference_auto.bed missing after pre-gen" >&2; exit 1; }
[ -s "${REF_DIR}/reference_x.bed" ]    || { echo "ERROR: reference_x.bed missing after pre-gen" >&2; exit 1; }

# 3) Build the final runtime image with the reference baked in.
docker build \
  --platform "${PLATFORM}" \
  --build-arg ADMIXTURE_VERSION="${VERSION}" \
  --target runtime \
  -f "${ROOT_DIR}/docker/admixture/Dockerfile" \
  -t "${LOCAL_TAG}" \
  -t "${IMAGE_VERSIONED}" \
  -t "${IMAGE_LATEST}" \
  "${ROOT_DIR}"

echo
echo "Built:"
echo "  ${LOCAL_TAG}"
echo "  ${IMAGE_VERSIONED}"
echo "  ${IMAGE_LATEST}"

if [ "${PUSH}" = "1" ]; then
  docker push "${IMAGE_VERSIONED}"
  docker push "${IMAGE_LATEST}"
fi

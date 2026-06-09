#!/usr/bin/env bash
# Build a BioVault OpenADMIXTURE image with an HGP1K PLINK reference baked in.
#
# Default mode is for the real all-sample DDNA-locus reference package:
#   data/hgp1k_all_sex_bias/
#
# PGP-filtered data is test data and is never selected implicitly. To bake the
# PGP package, pass --pgp explicitly:
#   ./build_openadmixture_docker.sh --pgp

set -euo pipefail

usage() {
  cat <<'EOF'
usage: ./build_openadmixture_docker.sh [--pgp] [--push]

Builds a derived OpenADMIXTURE image with the HGP1K reference copied to:
  /opt/biovault/reference/hgp1k_admixture

Modes:
  default   Uses data/hgp1k_all_sex_bias/        (real all-sample DDNA-locus package)
  --pgp     Uses data/hgp1k_1500_sex_bias_pgp/   (PGP-filtered test package)

Environment:
  BASE_IMAGE       Base OpenADMIXTURE image. Defaults to
                   ghcr.io/madhavajay/openadmixture.jl:madhava-update.
  IMAGE_NAME       Output repository
  PLATFORM         Docker platform. Defaults to the local base image platform
                   for local builds, otherwise linux/amd64.
  LOCAL_TAG        Local tag override
  LOCAL_LATEST     Local latest tag override
  PUSH=1           Push remote tags
  ALLOW_PGP_PUSH=1 Permit pushing a --pgp image
EOF
}

PGP=0
PUSH="${PUSH:-0}"
while [ "$#" -gt 0 ]; do
  case "$1" in
    --pgp)
      PGP=1
      shift
      ;;
    --push)
      PUSH=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_IMAGE="${BASE_IMAGE:-ghcr.io/madhavajay/openadmixture.jl:madhava-update}"
IMAGE_NAME="${IMAGE_NAME:-ghcr.io/madhavajay/biovault-openadmixture}"
PLATFORM="${PLATFORM:-}"

if [ -z "${PLATFORM}" ]; then
  if [ "${PUSH}" != "1" ] && docker image inspect "${BASE_IMAGE}" >/dev/null 2>&1; then
    os="$(docker image inspect "${BASE_IMAGE}" --format '{{.Os}}')"
    arch="$(docker image inspect "${BASE_IMAGE}" --format '{{.Architecture}}')"
    PLATFORM="${os}/${arch}"
  else
    PLATFORM="linux/amd64"
  fi
fi

REF_GEN="${ROOT_DIR}/flows/03_bv_paper_sex_biased_admixture_hgp1k/reference"
REF_DIR="${ROOT_DIR}/.docker/reference/openadmixture_hgp1k"

if [ "${PGP}" = "1" ]; then
  PACKAGE_NAME="${PACKAGE_NAME:-hgp1k_1500_sex_bias_pgp}"
  FLAVOR="${FLAVOR:-pgp1500}"
  if [ "${PUSH}" = "1" ] && [ "${ALLOW_PGP_PUSH:-0}" != "1" ]; then
    echo "ERROR: refusing to push PGP test-data image without ALLOW_PGP_PUSH=1" >&2
    exit 1
  fi
else
  PACKAGE_NAME="${PACKAGE_NAME:-hgp1k_all_sex_bias}"
  FLAVOR="${FLAVOR:-ddnaall}"
fi

SHARD_DIR="${SHARD_DIR:-${ROOT_DIR}/data/${PACKAGE_NAME}}"
LOCAL_TAG="${LOCAL_TAG:-biovault-openadmixture:${FLAVOR}-amd64}"
LOCAL_LATEST="${LOCAL_LATEST:-biovault-openadmixture:latest}"
IMAGE_FLAVOR="${IMAGE_NAME}:${FLAVOR}-amd64"
IMAGE_FLAVOR_PLAIN="${IMAGE_NAME}:${FLAVOR}"
IMAGE_LATEST="${IMAGE_NAME}:latest-amd64"
IMAGE_LATEST_PLAIN="${IMAGE_NAME}:latest"

if ! ls "${SHARD_DIR}/${PACKAGE_NAME}.tar.gz".* >/dev/null 2>&1; then
  echo "ERROR: missing reference shards in ${SHARD_DIR}" >&2
  echo "Expected files matching: ${SHARD_DIR}/${PACKAGE_NAME}.tar.gz.*" >&2
  if [ "${PGP}" = "1" ]; then
    echo "Generate them with: flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_pgp_reference_package.sh balanced5" >&2
  else
    echo "Generate/commit the real DDNA-locus package before running CI." >&2
  fi
  exit 1
fi

echo "==> Reassembling ${PACKAGE_NAME} into ${REF_DIR}"
rm -rf "${REF_DIR}"
SHARD_DIR="${SHARD_DIR}" \
DEST_DIR="${REF_DIR}" \
PACKAGE_NAME="${PACKAGE_NAME}" \
bash "${REF_GEN}/reassemble_reference.sh"

[ -s "${REF_DIR}/reference_auto.bed" ] || { echo "ERROR: reference_auto.bed missing" >&2; exit 1; }
[ -s "${REF_DIR}/reference_x.bed" ]    || { echo "ERROR: reference_x.bed missing" >&2; exit 1; }

echo "==> Building OpenADMIXTURE image (${FLAVOR}) from ${BASE_IMAGE}"
tags=(-t "${LOCAL_TAG}" -t "${LOCAL_LATEST}" -t "${IMAGE_FLAVOR}" -t "${IMAGE_FLAVOR_PLAIN}")
if [ "${PGP}" != "1" ] || [ "${PUSH}" != "1" ]; then
  tags+=(-t "${IMAGE_LATEST}" -t "${IMAGE_LATEST_PLAIN}")
fi

docker build \
  --platform "${PLATFORM}" \
  --build-arg "BASE_IMAGE=${BASE_IMAGE}" \
  -f "${ROOT_DIR}/docker/openadmixture/Dockerfile" \
  "${tags[@]}" \
  "${ROOT_DIR}"

echo
echo "Built:"
echo "  ${LOCAL_TAG}"
echo "  ${LOCAL_LATEST}"
echo "  ${IMAGE_FLAVOR}"
echo "  ${IMAGE_FLAVOR_PLAIN}"
if [ "${PGP}" != "1" ] || [ "${PUSH}" != "1" ]; then
  echo "  ${IMAGE_LATEST}"
  echo "  ${IMAGE_LATEST_PLAIN}"
fi

if [ "${PUSH}" = "1" ]; then
  docker push "${IMAGE_FLAVOR}"
  docker push "${IMAGE_FLAVOR_PLAIN}"
  if [ "${PGP}" != "1" ]; then
    docker push "${IMAGE_LATEST}"
    docker push "${IMAGE_LATEST_PLAIN}"
  fi
fi

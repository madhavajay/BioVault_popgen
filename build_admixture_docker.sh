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

docker build \
  --platform "${PLATFORM}" \
  --build-arg ADMIXTURE_VERSION="${VERSION}" \
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

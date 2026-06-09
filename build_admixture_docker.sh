#!/usr/bin/env bash
# Build a BioVault ADMIXTURE 1.4 image with an HGP1K PLINK reference baked in.
#
# Default mode is for the real all-sample DDNA-locus reference package:
#   data/hgp1k_all_sex_bias/
#
# PGP-filtered data is test data and is never selected implicitly. To bake the
# PGP package, pass --pgp explicitly:
#   ./build_admixture_docker.sh --pgp

set -euo pipefail

usage() {
  cat <<'EOF'
usage: ./build_admixture_docker.sh [--pgp] [--push]

Builds a derived ADMIXTURE image with the HGP1K reference copied to:
  /opt/biovault/reference/hgp1k_admixture

Modes:
  default   Uses data/hgp1k_all_sex_bias/        (real all-sample DDNA-locus package)
  --pgp     Uses data/hgp1k_1500_sex_bias_pgp/   (PGP-filtered test package)

Examples:
  ./build_admixture_docker.sh
  ./build_admixture_docker.sh --pgp
  PACKAGE_NAME=hgp1k_1500_sex_bias FLAVOR=ddna1500 ./build_admixture_docker.sh
  PACKAGE_NAME=hgp1k_all_sex_bias_pgp FLAVOR=pgpall ./build_admixture_docker.sh --pgp

Environment:
  VERSION          Image version tag. Defaults to 0.2.5-fast.
  ADMIXTURE_VERSION ADMIXTURE binary version. Defaults to 1.4.0.
  ADMIXTURE_URL    ADMIXTURE binary archive URL forwarded to Docker.
  IMAGE_NAME       Output repository.
  PLATFORM         Docker platform. Defaults to linux/amd64.
  PACKAGE_NAME     Reference package directory/archive basename override.
  FLAVOR           Dataset flavor tag override.
  LOCAL_TAG        Local version tag override.
  LOCAL_LATEST     Local latest tag override.
  PUSH=1           Push remote tags.
  ALLOW_PGP_PUSH=1 Permit pushing a --pgp image.
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
VERSION="${VERSION:-0.2.5-fast}"
ADMIXTURE_VERSION="${ADMIXTURE_VERSION:-1.4.0}"
IMAGE_NAME="${IMAGE_NAME:-ghcr.io/madhavajay/biovault-admixture}"
PLATFORM="${PLATFORM:-linux/amd64}"

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

REF_GEN="${ROOT_DIR}/flows/03_bv_paper_sex_biased_admixture_hgp1k/reference"
REF_DIR="${ROOT_DIR}/.docker/reference/hgp1k_admixture"
SHARD_DIR="${SHARD_DIR:-${ROOT_DIR}/data/${PACKAGE_NAME}}"

LOCAL_TAG="${LOCAL_TAG:-biovault-admixture:${VERSION}-amd64}"
LOCAL_LATEST="${LOCAL_LATEST:-biovault-admixture:latest}"
LOCAL_FLAVOR="${LOCAL_FLAVOR:-biovault-admixture:${FLAVOR}-amd64}"
IMAGE_VERSIONED="${IMAGE_NAME}:${VERSION}-amd64"
IMAGE_VERSIONED_PLAIN="${IMAGE_NAME}:${VERSION}"
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

echo "==> Building ADMIXTURE image ${VERSION} (${FLAVOR}) with ADMIXTURE ${ADMIXTURE_VERSION}"
tags=(
  -t "${LOCAL_TAG}"
  -t "${LOCAL_LATEST}"
  -t "${LOCAL_FLAVOR}"
  -t "${IMAGE_VERSIONED}"
  -t "${IMAGE_VERSIONED_PLAIN}"
  -t "${IMAGE_FLAVOR}"
  -t "${IMAGE_FLAVOR_PLAIN}"
)
if [ "${PGP}" != "1" ] || [ "${PUSH}" != "1" ]; then
  tags+=(-t "${IMAGE_LATEST}" -t "${IMAGE_LATEST_PLAIN}")
fi
build_args=(--build-arg "ADMIXTURE_VERSION=${ADMIXTURE_VERSION}")
if [ -n "${ADMIXTURE_URL:-}" ]; then
  build_args+=(--build-arg "ADMIXTURE_URL=${ADMIXTURE_URL}")
fi

docker build \
  --platform "${PLATFORM}" \
  "${build_args[@]}" \
  --target runtime \
  -f "${ROOT_DIR}/docker/admixture/Dockerfile" \
  "${tags[@]}" \
  "${ROOT_DIR}"

echo
echo "Built:"
echo "  ${LOCAL_TAG}"
echo "  ${LOCAL_LATEST}"
echo "  ${LOCAL_FLAVOR}"
echo "  ${IMAGE_VERSIONED}"
echo "  ${IMAGE_VERSIONED_PLAIN}"
echo "  ${IMAGE_FLAVOR}"
echo "  ${IMAGE_FLAVOR_PLAIN}"
if [ "${PGP}" != "1" ] || [ "${PUSH}" != "1" ]; then
  echo "  ${IMAGE_LATEST}"
  echo "  ${IMAGE_LATEST_PLAIN}"
fi

if [ "${PUSH}" = "1" ]; then
  docker push "${IMAGE_VERSIONED}"
  docker push "${IMAGE_VERSIONED_PLAIN}"
  docker push "${IMAGE_FLAVOR}"
  docker push "${IMAGE_FLAVOR_PLAIN}"
  if [ "${PGP}" != "1" ]; then
    docker push "${IMAGE_LATEST}"
    docker push "${IMAGE_LATEST_PLAIN}"
  fi
fi

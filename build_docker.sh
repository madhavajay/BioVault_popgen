#!/usr/bin/env bash
# Build the BioVault_popgen Docker image.
#
# Overrides:
#   VERSION=0.1.3 ./build_docker.sh                   # tag ghcr.io/madhavajay/biovault-popgen:0.1.3 + :latest
#   IMAGE_NAME=ghcr.io/foo/biovault-popgen ./build_docker.sh
#   PLATFORM=linux/arm64 ./build_docker.sh
#   BUILD_FAST=0 ./build_docker.sh                    # skip :<version>-fast + :fast
#   FORCE_REFERENCE_CACHE=1 ./build_docker.sh         # re-mirror loadings cache from GCS

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="${VERSION:-0.1.3}"
IMAGE_NAME="${IMAGE_NAME:-ghcr.io/madhavajay/biovault-popgen}"
IMAGE_VERSIONED="${IMAGE_NAME}:${VERSION}"
IMAGE_LATEST="${IMAGE_NAME}:latest"
FAST_VERSION="${FAST_VERSION:-${VERSION}-fast}"
IMAGE_FAST_VERSIONED="${IMAGE_NAME}:${FAST_VERSION}"
IMAGE_FAST_LATEST="${IMAGE_NAME}:fast"
PLATFORM="${PLATFORM:-linux/amd64}"
BUILD_FAST="${BUILD_FAST:-1}"
TOOLS_IMAGE="${TOOLS_IMAGE:-${IMAGE_NAME}:tools}"
LOADINGS_HT_SOURCE="${LOADINGS_HT_SOURCE:-gs://gcp-public-data--gnomad/release/3.1/pca/gnomad.v3.1.pca_loadings.ht}"
CACHE_DIR="${ROOT_DIR}/.docker/reference/pca_loadings"
CACHED_HT="${CACHE_DIR}/gnomad.v3.1.pca_loadings.ht"
CACHED_HT_TAR="${CACHE_DIR}/gnomad.v3.1.pca_loadings.ht.tar.gz"
CACHED_TSV="${CACHE_DIR}/loadings_variants.tsv"
AIMS_CACHE_DIR="${ROOT_DIR}/.docker/reference/aims"
AIMS_AF_TSV="${AIMS_CACHE_DIR}/gnomad_af_per_locus.tsv"
HGDP_TGP_VCF_DIR="${HGDP_TGP_VCF_DIR:-${ROOT_DIR}/.docker/reference/hgdp_tgp_vcf}"
HGDP_TGP_PANEL="${ROOT_DIR}/03_individual_level/gnomad_projection/reference/panel_hgdp_tgp.tsv"

docker build \
  --platform "${PLATFORM}" \
  --target tools \
  -f "${ROOT_DIR}/Dockerfile" \
  -t "${TOOLS_IMAGE}" \
  "${ROOT_DIR}"

# A fresh clone may carry the compact committed archive instead of the
# expanded 11k-file Hail Table directory. Expand it before the cache check so
# the rest of the build path stays identical for local, CI, and old checkouts.
if [ ! -f "${CACHED_HT}/_SUCCESS" ] && [ -s "${CACHED_HT_TAR}" ]; then
  echo "Expanding cached HT archive ${CACHED_HT_TAR}"
  mkdir -p "${CACHE_DIR}"
  tar -xzf "${CACHED_HT_TAR}" -C "${CACHE_DIR}"
fi

# 1. Mirror gnomAD HT to local cache via gsutil. Resumable: rerun to fill in
#    any missing files. Force a full re-mirror with FORCE_REFERENCE_CACHE=1.
if [ "${FORCE_REFERENCE_CACHE:-0}" = "1" ] || [ ! -f "${CACHED_HT}/_SUCCESS" ]; then
  echo "Mirroring ${LOADINGS_HT_SOURCE} -> ${CACHED_HT}"
  mkdir -p "${CACHE_DIR}"
  if command -v gsutil >/dev/null 2>&1; then
    GSUTIL=(gsutil)
  elif command -v uvx >/dev/null 2>&1; then
    GSUTIL=(uvx --from gsutil gsutil)
  else
    echo "ERROR: cache fetch needs gsutil or uvx on the host." >&2
    echo "       Install one of:" >&2
    echo "         pip install gsutil" >&2
    echo "         curl -LsSf https://astral.sh/uv/install.sh | sh   # provides uvx" >&2
    echo "       Or pre-populate the cache yourself:" >&2
    echo "         gsutil -m cp -r ${LOADINGS_HT_SOURCE} ${CACHE_DIR}/" >&2
    exit 1
  fi
  "${GSUTIL[@]}" -m cp -r "${LOADINGS_HT_SOURCE}" "${CACHE_DIR}/"
else
  echo "Using cached HT at ${CACHED_HT}"
fi

# 2. Derive the variant TSV from the local HT via Hail in-container. Local-only
#    Hail read, so no GCS in the loop.
if [ "${FORCE_REFERENCE_CACHE:-0}" = "1" ] || [ ! -s "${CACHED_TSV}" ]; then
  echo "Deriving ${CACHED_TSV} from local HT"
  docker run --rm \
    --platform "${PLATFORM}" \
    -v "${ROOT_DIR}:/work" \
    -w /work \
    "${TOOLS_IMAGE}" \
    python /tmp/extract_loadings_variants.py \
      "/work/.docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht" \
      "/work/.docker/reference/pca_loadings/loadings_variants.tsv"
else
  echo "Using cached TSV at ${CACHED_TSV}"
fi

# 2b. Derive the float64 loadings.npz from the local HT via Hail in-container
#     (same pattern as the TSV). Baked so gnomad_projection_fast skips the
#     ~45 s Hail extraction on every flow run. Skip if cached.
CACHED_NPZ="${CACHE_DIR}/loadings.npz"
if [ "${FORCE_REFERENCE_CACHE:-0}" = "1" ] || [ ! -s "${CACHED_NPZ}" ]; then
  echo "Deriving ${CACHED_NPZ} from local HT"
  docker run --rm \
    --platform "${PLATFORM}" \
    -v "${ROOT_DIR}:/work" \
    -w /work \
    "${TOOLS_IMAGE}" \
    python /work/03_individual_level/gnomad_projection_fast/scripts/extract_loadings_matrix.py \
      --ht  "/work/.docker/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht" \
      --out "/work/.docker/reference/pca_loadings/loadings.npz"
else
  echo "Using cached loadings.npz at ${CACHED_NPZ}"
fi

# 3. Derive the small AIMs gnomAD AF table from the HGDP+TGP VCFs. Same shape
#    as the loadings cache: skip if already present, otherwise mirror the
#    ~80 GB VCFs locally and derive the few-hundred-KB TSV in-container. Only
#    this TSV is COPYed into the runtime image (the VCFs never are).
if [ "${FORCE_REFERENCE_CACHE:-0}" = "1" ] || [ ! -s "${AIMS_AF_TSV}" ]; then
  echo "Deriving ${AIMS_AF_TSV} from HGDP+TGP VCFs"
  mkdir -p "${AIMS_CACHE_DIR}" "${HGDP_TGP_VCF_DIR}"
  if [ ! -s "${HGDP_TGP_VCF_DIR}/gnomad.genomes.v3.1.2.hgdp_tgp.chr22.vcf.bgz" ]; then
    echo "Mirroring HGDP+TGP VCFs -> ${HGDP_TGP_VCF_DIR} (~80 GB, resumable)"
    bash "${ROOT_DIR}/02_reference_panels/scripts/download_gnomad_v3_hgdp_tgp.sh" \
      "${HGDP_TGP_VCF_DIR}"
  fi
  docker run --rm \
    --platform "${PLATFORM}" \
    -v "${ROOT_DIR}:/work" \
    -v "${HGDP_TGP_VCF_DIR}:/hgdp_tgp_vcf" \
    -w /work \
    "${TOOLS_IMAGE}" \
    python /work/04_population_level/fst_aims_fast/build/derive_gnomad_aims_af.py \
      "/work/03_individual_level/gnomad_projection/reference/panel_hgdp_tgp.tsv" \
      "/hgdp_tgp_vcf" \
      "/work/.docker/reference/aims/gnomad_af_per_locus.tsv"
else
  echo "Using cached AIMs AF table at ${AIMS_AF_TSV}"
fi

docker build \
  --platform "${PLATFORM}" \
  --target runtime \
  -f "${ROOT_DIR}/Dockerfile" \
  -t "${IMAGE_VERSIONED}" \
  -t "${IMAGE_LATEST}" \
  "${ROOT_DIR}"

if [ "${BUILD_FAST}" = "1" ]; then
  docker build \
    --platform "${PLATFORM}" \
    --target fast-runtime \
    -f "${ROOT_DIR}/Dockerfile" \
    -t "${IMAGE_FAST_VERSIONED}" \
    -t "${IMAGE_FAST_LATEST}" \
    "${ROOT_DIR}"
fi

echo
echo "Built:"
echo "  ${IMAGE_VERSIONED}"
echo "  ${IMAGE_LATEST}"
if [ "${BUILD_FAST}" = "1" ]; then
  echo "  ${IMAGE_FAST_VERSIONED}"
  echo "  ${IMAGE_FAST_LATEST}"
fi

#!/usr/bin/env bash
# Reassemble the baked HGP1K ADMIXTURE reference BED from the committed shards
# in data/hgp1k_900_sex_bias/ into .docker/reference/hgp1k_admixture/, so the
# image can be built without the raw 1KGP VCFs (which are not in CI).
#
# Requires: tar, cat. Optional b3sum/shasum for verification.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"

SHARD_DIR="${SHARD_DIR:-${REPO_ROOT}/data/hgp1k_900_sex_bias}"
DEST_DIR="${DEST_DIR:-${REPO_ROOT}/.docker/reference/hgp1k_admixture}"
PACKAGE_NAME="${PACKAGE_NAME:-$(basename "${SHARD_DIR}")}"
ARCHIVE="${ARCHIVE:-${PACKAGE_NAME}.tar.gz}"

shards=( "${SHARD_DIR}/${ARCHIVE}".* )
[ -e "${shards[0]}" ] || { echo "ERROR: no shards in ${SHARD_DIR}" >&2; exit 1; }

mkdir -p "${DEST_DIR}"
tmp="$(mktemp)"
cat "${SHARD_DIR}/${ARCHIVE}".* > "${tmp}"
tar -xzf "${tmp}" -C "${DEST_DIR}"
rm -f "${tmp}"

for f in reference_auto.bed reference_x.bed reference_labels.tsv; do
    [ -s "${DEST_DIR}/${f}" ] || { echo "ERROR: ${f} missing after extract" >&2; exit 1; }
done
echo "[reassemble] reference BED -> ${DEST_DIR}"
ls -1 "${DEST_DIR}" | grep -E '^reference_'

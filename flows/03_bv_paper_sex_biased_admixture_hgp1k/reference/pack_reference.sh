#!/usr/bin/env bash
# Package the baked HGP1K ADMIXTURE reference BED into committable, <100MB
# GitHub-pushable shards under data/hgp1k_900_sex_bias/ — same scheme as
# data/hgp1k_dosage_split (tar.gz | split -b 95m + a .yaml sidecar).
#
# The shards carry the runtime files (reference_auto/reference_x BED + labels +
# the frozen sample-selection list). Reassemble with reassemble_reference.sh.
#
# Requires: tar, split, b3sum (or shasum fallback).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"

SRC_DIR="${SRC_DIR:-${REPO_ROOT}/.docker/reference/hgp1k_admixture}"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/data/hgp1k_900_sex_bias}"
ARCHIVE="hgp1k_900_sex_bias.tar.gz"
SHARD_SIZE="${SHARD_SIZE:-95m}"
SAMPLES_TSV="${HERE}/reference_samples.tsv"
PED="${REPO_ROOT}/data/1kgp_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

CONTENTS=(reference_auto.bed reference_auto.bim reference_auto.fam
          reference_x.bed reference_x.bim reference_x.fam
          reference_labels.tsv reference_samples.tsv build_manifest.txt)

for f in "${CONTENTS[@]}"; do
    [ -s "${SRC_DIR}/${f}" ] || { echo "ERROR: missing ${SRC_DIR}/${f}" >&2; exit 1; }
done

b3() { if command -v b3sum >/dev/null 2>&1; then b3sum "$1" | awk '{print $1}';
       else shasum -a 256 "$1" | awk '{print $1}'; fi; }
HASH_ALGO=$(command -v b3sum >/dev/null 2>&1 && echo b3sum || echo sha256)

rm -rf "${OUT_DIR}"; mkdir -p "${OUT_DIR}"

# 1) tar.gz the runtime files (flat paths) then split into shards.
tar -C "${SRC_DIR}" -czf - "${CONTENTS[@]}" \
    | split -b "${SHARD_SIZE}" - "${OUT_DIR}/${ARCHIVE}."

# reconstruct the full archive once to record its hash/size, then drop it.
cat "${OUT_DIR}/${ARCHIVE}".* > "${OUT_DIR}/${ARCHIVE}"
ARCHIVE_SIZE=$(wc -c < "${OUT_DIR}/${ARCHIVE}" | tr -d ' ')
ARCHIVE_HASH=$(b3 "${OUT_DIR}/${ARCHIVE}")
rm -f "${OUT_DIR}/${ARCHIVE}"

# 2) human-readable selection list outside the archive (browsable on GitHub).
cp "${SAMPLES_TSV}" "${OUT_DIR}/reference_samples.tsv"

# 3) per-superpopulation / per-population selection summary.
echo "[pack] selection summary:"
sel_summary=$(awk -F'\t' 'NR>1{sp[$2]++; pop[$2"/"$3]++; sx[$2"/"$4]++}
    END{for(k in sp) printf "  %s total=%d\n", k, sp[k]}' "${SAMPLES_TSV}" | sort)
echo "${sel_summary}"

# 4) yaml sidecar (mirrors data/hgp1k_dosage_split/*.yaml).
YAML="${OUT_DIR}/hgp1k_900_sex_bias.yaml"
{
  echo "version: 1"
  echo "description: >-"
  echo "  Baked HGP1K ADMIXTURE reference for the sex-biased admixture flow:"
  echo "  900 unrelated 1000 Genomes founders (300 AFR + 300 EUR + 300 SAS),"
  echo "  PLINK BED keyed by chr:pos, pre-split into autosomes (reference_auto)"
  echo "  and X (reference_x, sex-coded for male haploid handling). Reassemble"
  echo "  into .docker/reference/hgp1k_admixture/ to bake the image without the"
  echo "  raw 1KGP VCFs (not in CI)."
  echo "hash_algo: ${HASH_ALGO}"
  echo "shard_size: \"${SHARD_SIZE}\""
  echo "archive_format: tar.gz"
  echo "selection:"
  echo "  source_ped: $(basename "${PED}")"
  echo "  policy: unrelated founders (FatherID==0 && MotherID==0)"
  echo "  seed: 42"
  echo "  per_superpop: { AFR: 300, EUR: 300, SAS: 300 }"
  echo "  sample_list: reference_samples.tsv  # sample_id, superpopulation, population, sex"
  echo "  per_population:"
  awk -F'\t' 'NR>1{c[$2"\t"$3]++} END{for(k in c){split(k,a,"\t"); printf "    %s_%s: %d\n", a[1], a[2], c[k]}}' \
      "${SAMPLES_TSV}" | sort
  echo "reference:"
  grep -E '^reference_(auto|x):' "${SRC_DIR}/build_manifest.txt" | sed 's/^/  /'
  echo "archive:"
  echo "  name: ${ARCHIVE}"
  echo "  size_bytes: ${ARCHIVE_SIZE}"
  echo "  ${HASH_ALGO}: ${ARCHIVE_HASH}"
  echo "  contents:"
  for f in "${CONTENTS[@]}"; do echo "    - ${f}"; done
  echo "  shards:"
  for s in "${OUT_DIR}/${ARCHIVE}".*; do
    [ -f "$s" ] || continue
    echo "    - name: $(basename "$s")"
    echo "      size_bytes: $(wc -c < "$s" | tr -d ' ')"
    echo "      ${HASH_ALGO}: $(b3 "$s")"
  done
  echo "instructions:"
  echo "  reassemble: cat ${ARCHIVE}.* > ${ARCHIVE}"
  echo "  extract: tar xzf ${ARCHIVE} -C <dest_dir>"
  echo "  helper: flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/reassemble_reference.sh"
} > "${YAML}"

echo "[pack] wrote shards + ${YAML}"
ls -la "${OUT_DIR}"
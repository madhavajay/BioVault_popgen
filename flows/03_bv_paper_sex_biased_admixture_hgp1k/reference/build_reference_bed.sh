#!/usr/bin/env bash
# Build the baked HGP1K ADMIXTURE reference: a PLINK BED of the frozen 1000
# Genomes subsample (reference_samples.tsv), for ALL chromosomes incl. X,
# keyed by chr:pos, pre-split into autosomes and X.
#
# Why this is NOT the projection dosage matrix (hgp1k_dosage.npz):
#   * ADMIXTURE needs hard genotype CALLS in PLINK .bed/.bim/.fam, not floats.
#   * the dosage matrix is autosomes-only; this analysis needs chrX.
# So we go back to the filtered 1KGP VCFs (data/1kgp_high_coverage/filtered),
# already restricted to the array (locus_map) sites.
#
# Variant keying: study (bvs cohort-bed) emits a MIX of rsID and chr:pos IDs,
# so neither side can be keyed by rsID. Both sides are re-keyed to chr:pos
# (`plink2 --set-all-var-ids @:#`), which is unambiguous in GRCh38 and the same
# scheme the runtime analysis applies to the study BED. The filtered-VCF IDs
# ("." on chrX, chr:pos:ref:alt on autosomes) are therefore irrelevant.
#
# chrX needs sex (males hemizygous) + --split-par so PAR is diploid and non-PAR
# males are haploid; sex comes from reference_samples.tsv col 4.
#
# Requires: bcftools, plink2 on PATH. Run inside an image that has both, e.g.
#   ./scripts/run_in_docker.sh bash \
#     flows/03_bv_paper_sex_biased_admixture_hgp1k/reference/build_reference_bed.sh
#
# Tunables (env): CHROMS, JOBS, THREADS, OUT_DIR, KEEP_WORK.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"

VCF_DIR="${VCF_DIR:-${REPO_ROOT}/data/1kgp_high_coverage/filtered}"
SAMPLES_TSV="${SAMPLES_TSV:-${HERE}/reference_samples.tsv}"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/.docker/reference/hgp1k_admixture}"
WORK_DIR="${WORK_DIR:-${OUT_DIR}/_work}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}"
JOBS="${JOBS:-4}"
THREADS="${THREADS:-4}"
BUILD="${BUILD:-hg38}"

for t in bcftools plink2; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found on PATH" >&2; exit 1; }
done
[ -s "${SAMPLES_TSV}" ] || { echo "ERROR: missing ${SAMPLES_TSV}" >&2; exit 1; }
[ -d "${VCF_DIR}" ]     || { echo "ERROR: missing VCF dir ${VCF_DIR}" >&2; exit 1; }

mkdir -p "${OUT_DIR}" "${WORK_DIR}"

# --- sample id list + sex map (IID SEX; 1=male 2=female) ---
SAMPLE_IDS="${WORK_DIR}/sample_ids.txt"
tail -n +2 "${SAMPLES_TSV}" | cut -f1 | sort -u > "${SAMPLE_IDS}"
N_SAMPLES=$(wc -l < "${SAMPLE_IDS}" | tr -d ' ')
SEX_FILE="${WORK_DIR}/sex.tsv"
{ printf '#IID\tSEX\n'; tail -n +2 "${SAMPLES_TSV}" | awk -F'\t' '{print $1"\t"$4}'; } > "${SEX_FILE}"
echo "[ref] ${N_SAMPLES} reference samples"

find_vcf() { ls "${VCF_DIR}"/*Illumina.chr${1}.filtered.*vcf.gz 2>/dev/null | head -1 || true; }

# --- build one chromosome -> ${WORK_DIR}/ref_chr<chr>.{bed,bim,fam} ---
build_one() {
    local chr="$1" vcf out renamemap
    vcf="$(find_vcf "${chr}")"
    out="${WORK_DIR}/ref_chr${chr}"
    if [ -z "${vcf}" ]; then echo "[ref][chr${chr}] no VCF, skip" >&2; return 0; fi
    echo "[ref][chr${chr}] $(basename "${vcf}")"

    renamemap="${WORK_DIR}/rename_chr${chr}.txt"
    printf 'chr%s\t%s\n' "${chr}" "${chr}" > "${renamemap}"

    # subset 900 samples, biallelic SNPs only, strip chr prefix
    bcftools view -S "${SAMPLE_IDS}" --force-samples -m2 -M2 -v snps \
            --threads "${THREADS}" "${vcf}" -Ou \
        | bcftools annotate --rename-chrs "${renamemap}" --threads "${THREADS}" \
            -Oz -o "${out}.vcf.gz"

    local x_args=()
    if [ "${chr}" = "X" ]; then
        x_args=(--update-sex "${SEX_FILE}" --split-par "${BUILD}")
    fi
    # chr:pos var IDs; drop same-position duplicates; biallelic SNPs.
    plink2 --vcf "${out}.vcf.gz" "${x_args[@]}" \
           --set-all-var-ids '@:#' --rm-dup exclude-all \
           --max-alleles 2 --snps-only \
           --make-bed --out "${out}" \
           --threads "${THREADS}" --memory 4000 >"${out}.plink.log" 2>&1
    rm -f "${out}.vcf.gz"
    echo "[ref][chr${chr}] $(wc -l < "${out}.bim" 2>/dev/null || echo 0) SNPs"
}
export -f build_one find_vcf
export WORK_DIR VCF_DIR SAMPLE_IDS SEX_FILE THREADS BUILD

printf '%s\n' ${CHROMS} | xargs -P "${JOBS}" -I{} bash -c 'build_one "$@"' _ {}

# --- merge per-chromosome filesets (same samples, disjoint variants) ---
MERGE_LIST="${WORK_DIR}/merge_list.txt"; : > "${MERGE_LIST}"
for chr in ${CHROMS}; do
    p="${WORK_DIR}/ref_chr${chr}"
    [ -s "${p}.bed" ] && [ -s "${p}.bim" ] && echo "${p}" >> "${MERGE_LIST}"
done
[ -s "${MERGE_LIST}" ] || { echo "ERROR: no chromosomes built" >&2; exit 1; }

FIRST=$(head -1 "${MERGE_LIST}")
tail -n +2 "${MERGE_LIST}" > "${WORK_DIR}/merge_tail.txt"
if [ -s "${WORK_DIR}/merge_tail.txt" ]; then
    plink2 --bfile "${FIRST}" --pmerge-list "${WORK_DIR}/merge_tail.txt" bfile \
           --make-bed --out "${OUT_DIR}/reference" \
           --threads "${THREADS}" --memory 6000 >"${WORK_DIR}/merge.log" 2>&1 || {
        echo "ERROR: merge failed — see ${WORK_DIR}/merge.log" >&2; exit 1; }
else
    for e in bed bim fam; do cp "${FIRST}.${e}" "${OUT_DIR}/reference.${e}"; done
fi

# --- pre-split autosomes / X (cheap row filter on the binary BED) ---
# re-stamp sex on X so male hemizygous handling is guaranteed downstream.
plink2 --bfile "${OUT_DIR}/reference" --chr 1-22 --make-bed \
       --out "${OUT_DIR}/reference_auto" --threads "${THREADS}" >/dev/null 2>&1
plink2 --bfile "${OUT_DIR}/reference" --chr X --update-sex "${SEX_FILE}" --make-bed \
       --out "${OUT_DIR}/reference_x" --threads "${THREADS}" >/dev/null 2>&1

# --- labels: sample_id -> superpopulation (steering / anchor labels) ---
cut -f1,2 "${SAMPLES_TSV}" > "${OUT_DIR}/reference_labels.tsv"

{
    echo "# HGP1K ADMIXTURE reference (chr:pos keyed)"
    echo "samples_tsv: ${SAMPLES_TSV}"
    echo "n_samples: ${N_SAMPLES}"
    echo "chromosomes: ${CHROMS}"
    echo "reference (all): $(wc -l < "${OUT_DIR}/reference.bim") SNPs, $(wc -l < "${OUT_DIR}/reference.fam") samples"
    echo "reference_auto:  $(wc -l < "${OUT_DIR}/reference_auto.bim") SNPs"
    echo "reference_x:     $(wc -l < "${OUT_DIR}/reference_x.bim") SNPs"
    echo "superpop_counts:"; tail -n +2 "${OUT_DIR}/reference_labels.tsv" | cut -f2 | sort | uniq -c
} > "${OUT_DIR}/build_manifest.txt"

# keep the frozen selection list alongside the BED (which 1KGP samples chosen).
if [ "$(cd "$(dirname "${SAMPLES_TSV}")" && pwd)/$(basename "${SAMPLES_TSV}")" != \
     "$(cd "${OUT_DIR}" && pwd)/reference_samples.tsv" ]; then
    cp "${SAMPLES_TSV}" "${OUT_DIR}/reference_samples.tsv"
fi

# the all-chr reference.* (bed/bim/fam + pgen/psam/pvar merge intermediates) is
# only used to derive the auto/X split; runtime reads reference_auto /
# reference_x only. Drop it to roughly halve the baked size.
if [ "${KEEP_ALL_CHR:-0}" != "1" ]; then
    rm -f "${OUT_DIR}"/reference.bed "${OUT_DIR}"/reference.bim \
          "${OUT_DIR}"/reference.fam "${OUT_DIR}"/reference.log \
          "${OUT_DIR}"/reference.pgen "${OUT_DIR}"/reference.psam \
          "${OUT_DIR}"/reference.pvar
fi

echo "[ref] DONE"; cat "${OUT_DIR}/build_manifest.txt"
[ "${KEEP_WORK:-0}" = "1" ] || rm -rf "${WORK_DIR}"

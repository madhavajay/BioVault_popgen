#!/usr/bin/env bash
# Run an individual-level pipeline inside the biovault-popgen docker image
# against a variable subset of mock participant genotypes.
#
# Modes:
#   (default)        gnomad_projection      — slow Hail projection (baked HT)
#   --fast           gnomad_projection_fast — numpy projection (bit-identical)
#   --qc             pca_qc_fast            — within-cohort QC + cohort PCA
#   --sex            sex_biased_admixture_fast — X-hemizygosity / sex facet
#   --slow           use the non-_fast implementation for the selected mode
#
# Usage:
#   bash 03_individual_level.sh              # all participants, projection
#   bash 03_individual_level.sh --fast 50    # 50 participants, fast projection
#   bash 03_individual_level.sh --qc 3       # 3 participants, qc
#   bash 03_individual_level.sh --qc --slow  # original pca_qc -> results/pca_qc
#   bash 03_individual_level.sh --sex 100    # 100 participants, sex-bias
#   bash 03_individual_level.sh --sex --slow # original sex_biased_admixture
#
# Overrides:
#   SAMPLES_SRC=/path/to/dir
#   IMAGE=ghcr.io/madhavajay/biovault-popgen:0.1.2
#   RESULTS_ROOT=/path/to/results

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT_DIR}/scripts/image_versions.sh"
SAMPLES_SRC="${SAMPLES_SRC:-${ROOT_DIR}/01_mock_data_generation/output}"
IMAGE="${IMAGE:-${BIOVAULT_IMAGE}}"
SUBSET_DIR="${ROOT_DIR}/03_individual_level/.samples"
RESULTS_ROOT="${RESULTS_ROOT:-${ROOT_DIR}/results}"

MODE="projection"
SLOW=0
N=""
for arg in "$@"; do
    case "${arg}" in
        --qc)        MODE="qc" ;;
        --fast)      MODE="projection_fast" ;;
        --sex)       MODE="sex" ;;
        --slow)      SLOW=1 ;;
        --proj*)     MODE="projection" ;;
        -h|--help)
            sed -n '2,21p' "$0"
            exit 0
            ;;
        *)
            if [[ "${arg}" =~ ^[0-9]+$ ]]; then
                N="${arg}"
            else
                echo "ERROR: unknown arg: ${arg}" >&2
                exit 2
            fi
            ;;
    esac
done

if [ "${SLOW}" = "1" ]; then
    case "${MODE}" in
        qc)              MODE="qc_slow" ;;
        sex)             MODE="sex_slow" ;;
        projection_fast) MODE="projection" ;;
        projection)      MODE="projection" ;;
    esac
fi

[ -d "${SAMPLES_SRC}" ] || { echo "ERROR: samples src missing: ${SAMPLES_SRC}" >&2; exit 1; }

# Collect numeric participant ids in the source, sorted. (Avoiding `mapfile`
# so this also runs on macOS's bash 3.2.)
ALL_PIDS=()
while IFS= read -r pid; do
    ALL_PIDS+=("${pid}")
done < <(
    find "${SAMPLES_SRC}" -maxdepth 1 -mindepth 1 -type d -name '[0-9]*' -exec basename {} \; \
        | sort
)
TOTAL=${#ALL_PIDS[@]}
[ "${TOTAL}" -gt 0 ] || { echo "ERROR: no participant dirs found in ${SAMPLES_SRC}" >&2; exit 1; }

if [ -z "${N}" ] || ! [[ "${N}" =~ ^[0-9]+$ ]] || [ "${N}" -gt "${TOTAL}" ]; then
    N="${TOTAL}"
fi

echo "Linking ${N} of ${TOTAL} participants from ${SAMPLES_SRC} -> ${SUBSET_DIR}"
rm -rf "${SUBSET_DIR}"
mkdir -p "${SUBSET_DIR}"
for pid in "${ALL_PIDS[@]:0:${N}}"; do
    ln -s "${SAMPLES_SRC}/${pid}" "${SUBSET_DIR}/${pid}"
done

# --- mode dispatch ----------------------------------------------------------

case "${MODE}" in
qc_slow)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/pca_qc"
    FAST_DIR="${ROOT_DIR}/03_individual_level/pca_qc_fast"
    TASK_DIR="${RESULTS_ROOT}/pca_qc"
    WORK_BASE="${TASK_DIR}/work/pca_qc"
    rm -rf "${TASK_DIR}"
    mkdir -p "${WORK_BASE}/scripts" "${WORK_BASE}/data" "${WORK_BASE}/plots" "${WORK_BASE}/logs"
    cp "${PIPELINE_DIR}/scripts/"*.py "${WORK_BASE}/scripts/"
    cp "${PIPELINE_DIR}/scripts/"*.sh "${WORK_BASE}/scripts/"
    echo "Running pca_qc (original Python QC + PCA) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${WORK_BASE}" \
        -e BIOVAULT_DATA_DIR="${SUBSET_DIR}" \
        "${IMAGE}" \
        bash -c '
            set -euo pipefail
            USER_ID="$(id -u)"; GROUP_ID="$(id -g)"
            if ! getent passwd "${USER_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
            fi
            if ! getent group "${GROUP_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${GROUP_ID}:" >> /etc/group
            fi
            export HOME=/tmp
            source /opt/conda/etc/profile.d/conda.sh
            conda activate biovault_popgen
            bash '"${WORK_BASE}"'/scripts/run_pipeline.sh
        '
    mkdir -p "${TASK_DIR}"
    cp "${WORK_BASE}/data/pca/pca.eigenvec" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/data/pca/pca.eigenval" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/data/merged/snp_info.tsv" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/pca_pc1_pc2.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/pca_pc3_pc4.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/01_merge.log" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/02_encode.log" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/03b_python_pca.log" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/04_plot_pca.log" "${TASK_DIR}/" 2>/dev/null || true
    echo
    echo "=== pca_qc outputs (N=${N}) -> ${TASK_DIR} ==="
    for f in \
        "${TASK_DIR}/snp_info.tsv" \
        "${TASK_DIR}/pca.eigenvec" \
        "${TASK_DIR}/pca.eigenval" \
        "${TASK_DIR}/pca_pc1_pc2.png" \
        "${TASK_DIR}/pca_pc3_pc4.png" \
        "${TASK_DIR}/01_merge.log" \
        "${TASK_DIR}/02_encode.log" \
        "${TASK_DIR}/03b_python_pca.log" \
        "${TASK_DIR}/04_plot_pca.log" \
        "${WORK_BASE}/data/merged/genotype_matrix_raw.tsv" \
        "${WORK_BASE}/data/merged/genotype_matrix_numeric.tsv" \
        "${WORK_BASE}/data/plink/genotypes.ped" \
        "${WORK_BASE}/data/plink/genotypes.map"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

qc)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/pca_qc_fast"
    TASK_DIR="${RESULTS_ROOT}/pca_qc_fast"
    WORK_BASE="${TASK_DIR}/work/pca_qc_fast"
    rm -rf "${TASK_DIR}"
    mkdir -p "${WORK_BASE}/scripts" "${WORK_BASE}/data" "${WORK_BASE}/plots" "${WORK_BASE}/logs"
    cp "${PIPELINE_DIR}/scripts/"*.py "${WORK_BASE}/scripts/"
    echo "Running pca_qc_fast (within-cohort QC + PCA) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${WORK_BASE}" \
        -e BIOVAULT_DATA_DIR="${SUBSET_DIR}" \
        "${IMAGE}" \
        bash -c '
            set -euo pipefail
            USER_ID="$(id -u)"; GROUP_ID="$(id -g)"
            if ! getent passwd "${USER_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
            fi
            if ! getent group "${GROUP_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${GROUP_ID}:" >> /etc/group
            fi
            export HOME=/tmp
            source /opt/conda/etc/profile.d/conda.sh
            conda activate biovault_popgen
            python3 '"${WORK_BASE}"'/scripts/fast_pipeline.py
        '
    mkdir -p "${TASK_DIR}"
    cp "${WORK_BASE}/data/pca/pca.eigenvec" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/data/pca/pca.eigenval" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/data/merged/snp_info.tsv" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/pca_pc1_pc2.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/pca_pc3_pc4.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/fast_pipeline.log" "${TASK_DIR}/" 2>/dev/null || true
    echo
    echo "=== pca_qc_fast outputs (N=${N}) -> ${TASK_DIR} ==="
    for f in \
        "${TASK_DIR}/snp_info.tsv" \
        "${TASK_DIR}/pca.eigenvec" \
        "${TASK_DIR}/pca.eigenval" \
        "${TASK_DIR}/pca_pc1_pc2.png" \
        "${TASK_DIR}/pca_pc3_pc4.png" \
        "${TASK_DIR}/fast_pipeline.log" \
        "${WORK_BASE}/data/merged/genotype_matrix_raw.tsv" \
        "${WORK_BASE}/data/merged/genotype_matrix_numeric.tsv" \
        "${WORK_BASE}/data/plink/genotypes.ped" \
        "${WORK_BASE}/data/plink/genotypes.map"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

sex)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/sex_biased_admixture_fast"
    ORIG_DIR="${ROOT_DIR}/03_individual_level/sex_biased_admixture"
    TASK_DIR="${RESULTS_ROOT}/sex_biased_admixture_fast"
    WORK_ROOT="${TASK_DIR}/work"
    WORK_FAST="${WORK_ROOT}/sex_biased_admixture_fast"
    WORK_ORIG="${WORK_ROOT}/sex_biased_admixture"
    rm -rf "${TASK_DIR}"
    mkdir -p "${WORK_FAST}/scripts" "${WORK_ORIG}/scripts"
    cp "${PIPELINE_DIR}/scripts/"*.py "${WORK_FAST}/scripts/"
    cp "${ORIG_DIR}/scripts/"*.py "${WORK_ORIG}/scripts/"
    echo "Running sex_biased_admixture_fast (X-hemizygosity, sex facet) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${WORK_FAST}" \
        -e BIOVAULT_DATA_DIR="${SUBSET_DIR}" \
        -e BIOVAULT_SEX_MAPPING="${SAMPLES_SRC}/sex_mapping.tsv" \
        "${IMAGE}" \
        bash -c '
            set -euo pipefail
            USER_ID="$(id -u)"; GROUP_ID="$(id -g)"
            if ! getent passwd "${USER_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
            fi
            if ! getent group "${GROUP_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${GROUP_ID}:" >> /etc/group
            fi
            export HOME=/tmp
            source /opt/conda/etc/profile.d/conda.sh
            conda activate biovault_popgen
            python3 '"${WORK_FAST}"'/scripts/fast_sex_biased_admixture.py
        '
    mkdir -p "${TASK_DIR}"
    cp "${WORK_FAST}/results/sex_bias_results.tsv" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_FAST}/results"/nmf_variant_filter_*.tsv "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_FAST}/plots/figure4_sex_biased_admixture.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_FAST}/plots/figure4_sex_biased_admixture.pdf" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_FAST}/logs/sex_biased_admixture.log" "${TASK_DIR}/" 2>/dev/null || \
        cp "${WORK_ORIG}/logs/sex_biased_admixture.log" "${TASK_DIR}/" 2>/dev/null || true
    echo
    echo "=== sex_biased_admixture_fast outputs (N=${N}) -> ${TASK_DIR} ==="
    for f in \
        "${TASK_DIR}/sex_bias_results.tsv" \
        "${TASK_DIR}/nmf_variant_filter_autosomes.tsv" \
        "${TASK_DIR}/nmf_variant_filter_x.tsv" \
        "${TASK_DIR}/figure4_sex_biased_admixture.png" \
        "${TASK_DIR}/figure4_sex_biased_admixture.pdf" \
        "${TASK_DIR}/sex_biased_admixture.log"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

sex_slow)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/sex_biased_admixture"
    TASK_DIR="${RESULTS_ROOT}/sex_biased_admixture"
    WORK_BASE="${TASK_DIR}/work/sex_biased_admixture"
    rm -rf "${TASK_DIR}"
    mkdir -p "${WORK_BASE}/scripts"
    cp "${PIPELINE_DIR}/scripts/"*.py "${WORK_BASE}/scripts/"
    echo "Running sex_biased_admixture (original analysis) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${WORK_BASE}" \
        -e BIOVAULT_DATA_DIR="${SUBSET_DIR}" \
        -e BIOVAULT_SEX_MAPPING="${SAMPLES_SRC}/sex_mapping.tsv" \
        "${IMAGE}" \
        bash -c '
            set -euo pipefail
            USER_ID="$(id -u)"; GROUP_ID="$(id -g)"
            if ! getent passwd "${USER_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:/tmp:/bin/bash" >> /etc/passwd
            fi
            if ! getent group "${GROUP_ID}" >/dev/null 2>&1; then
                echo "biovault:x:${GROUP_ID}:" >> /etc/group
            fi
            export HOME=/tmp
            source /opt/conda/etc/profile.d/conda.sh
            conda activate biovault_popgen
            python3 '"${WORK_BASE}"'/scripts/sex_biased_admixture.py
        '
    mkdir -p "${TASK_DIR}"
    cp "${WORK_BASE}/results/sex_bias_results.tsv" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/results"/nmf_variant_filter_*.tsv "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/figure4_sex_biased_admixture.png" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/plots/figure4_sex_biased_admixture.pdf" "${TASK_DIR}/" 2>/dev/null || true
    cp "${WORK_BASE}/logs/sex_biased_admixture.log" "${TASK_DIR}/" 2>/dev/null || true
    echo
    echo "=== sex_biased_admixture outputs (N=${N}) -> ${TASK_DIR} ==="
    for f in \
        "${TASK_DIR}/sex_bias_results.tsv" \
        "${TASK_DIR}/nmf_variant_filter_autosomes.tsv" \
        "${TASK_DIR}/nmf_variant_filter_x.tsv" \
        "${TASK_DIR}/figure4_sex_biased_admixture.png" \
        "${TASK_DIR}/figure4_sex_biased_admixture.pdf" \
        "${TASK_DIR}/sex_biased_admixture.log"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

projection_fast)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/gnomad_projection_fast"
    OUT_DIR="${RESULTS_ROOT}/gnomad_projection_fast"
    WORK_DIR="${OUT_DIR}/work"
    rm -rf "${OUT_DIR}"
    mkdir -p "${OUT_DIR}" "${WORK_DIR}"

    echo "Running gnomad_projection_fast (same QC as slow, numpy PCA project) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${PIPELINE_DIR}" \
        -e GENO="${GENO:-}" \
        -e MIND="${MIND:-}" \
        -e MAF="${MAF:-}" \
        -e HWE="${HWE:-}" \
        -e MIN_GS="${MIN_GS:-}" \
        "${IMAGE}" \
        bash "${PIPELINE_DIR}/scripts/run_fast_projection.sh" \
            "${SUBSET_DIR}" "${WORK_DIR}" "${OUT_DIR}"

    echo
    echo "=== gnomad_projection_fast outputs (N=${N}) ==="
    for f in \
        "${OUT_DIR}/study_pca_projection.tsv" \
        "${OUT_DIR}/qc_report.txt" \
        "${OUT_DIR}/pca_projection.png"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

projection)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/gnomad_projection"
    OUT_DIR="${RESULTS_ROOT}/gnomad_projection"
    WORK_DIR="${OUT_DIR}/work"
    rm -rf "${OUT_DIR}"
    mkdir -p "${OUT_DIR}" "${WORK_DIR}"

    echo "Running gnomad_projection (project onto baked HGDP+1kGP PCA space) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${PIPELINE_DIR}" \
        -e GENO="${GENO:-}" \
        -e MIND="${MIND:-}" \
        -e MAF="${MAF:-}" \
        -e HWE="${HWE:-}" \
        -e MIN_GS="${MIN_GS:-}" \
        "${IMAGE}" \
        bash /opt/biovault/scripts/gnomad_projection/run_flow_projection.sh \
            "${SUBSET_DIR}" "${WORK_DIR}" "${OUT_DIR}"

    echo
    echo "=== gnomad_projection outputs (N=${N}) ==="
    for f in \
        "${OUT_DIR}/study_pca_projection.tsv" \
        "${OUT_DIR}/qc_report.txt" \
        "${OUT_DIR}/pca_projection.png" \
        "${OUT_DIR}/hail.log"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;
esac

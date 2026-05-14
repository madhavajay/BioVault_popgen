#!/usr/bin/env bash
# Run an individual-level pipeline inside the biovault-popgen docker image
# against a variable subset of mock participant genotypes.
#
# Modes:
#   (default)        gnomad_projection  — project samples onto the gnomAD
#                                          HGDP+1kGP PCA space (uses baked HT)
#   --qc             pca_qc_fast        — within-cohort QC + cohort PCA
#
# Usage:
#   bash 03_individual_level.sh              # all participants, projection
#   bash 03_individual_level.sh 5            # 5 participants, projection
#   bash 03_individual_level.sh --qc         # all participants, qc
#   bash 03_individual_level.sh --qc 3       # 3 participants, qc
#
# Overrides:
#   SAMPLES_SRC=/path/to/dir
#   IMAGE=biovault-popgen:0.1.0

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES_SRC="${SAMPLES_SRC:-${ROOT_DIR}/01_mock_data_generation/output}"
IMAGE="${IMAGE:-biovault-popgen:0.1.0}"
SUBSET_DIR="${ROOT_DIR}/03_individual_level/.samples"

MODE="projection"
N=""
for arg in "$@"; do
    case "${arg}" in
        --qc)        MODE="qc" ;;
        --fast)      MODE="projection_fast" ;;
        --proj*)     MODE="projection" ;;
        -h|--help)
            sed -n '2,18p' "$0"
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
qc)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/pca_qc_fast"
    echo "Running pca_qc_fast (within-cohort QC + PCA) in ${IMAGE} ..."
    docker run --rm \
        --platform linux/amd64 \
        -u "$(id -u):$(id -g)" \
        -v "${ROOT_DIR}:${ROOT_DIR}" \
        -w "${PIPELINE_DIR}" \
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
            python3 '"${PIPELINE_DIR}"'/scripts/fast_pipeline.py
        '
    echo
    echo "=== pca_qc_fast outputs (N=${N}) ==="
    for f in \
        "${PIPELINE_DIR}/data/merged/genotype_matrix_raw.tsv" \
        "${PIPELINE_DIR}/data/merged/genotype_matrix_numeric.tsv" \
        "${PIPELINE_DIR}/data/merged/snp_info.tsv" \
        "${PIPELINE_DIR}/data/plink/genotypes.ped" \
        "${PIPELINE_DIR}/data/plink/genotypes.map" \
        "${PIPELINE_DIR}/data/pca/pca.eigenvec" \
        "${PIPELINE_DIR}/data/pca/pca.eigenval" \
        "${PIPELINE_DIR}/plots/pca_pc1_pc2.png" \
        "${PIPELINE_DIR}/plots/pca_pc3_pc4.png" \
        "${PIPELINE_DIR}/logs/fast_pipeline.log"; do
        if [ -e "${f}" ]; then
            printf '  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
        else
            printf '  %s  (missing)\n' "${f}"
        fi
    done
    ;;

projection_fast)
    PIPELINE_DIR="${ROOT_DIR}/03_individual_level/gnomad_projection_fast"
    OUT_DIR="${PIPELINE_DIR}/results/local"
    WORK_DIR="${PIPELINE_DIR}/working/local"
    rm -rf "${OUT_DIR}" "${WORK_DIR}"
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
        -e HWE="${HWE:-0}" \
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
    OUT_DIR="${PIPELINE_DIR}/results/local"
    WORK_DIR="${PIPELINE_DIR}/working/local"
    rm -rf "${OUT_DIR}" "${WORK_DIR}"
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
        -e HWE="${HWE:-0}" \
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

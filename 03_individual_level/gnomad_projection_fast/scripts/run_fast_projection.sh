#!/usr/bin/env bash
# Fast equivalent of ../gnomad_projection/scripts/run_flow_projection.sh.
#
# Same DDNA -> PLINK -> QC pipeline, but the Hail-based PC projection is
# replaced by a numpy implementation that reproduces hl.experimental.pc_project
# bit-for-bit at float64 precision. Avoids ~45 s JVM warm-up.
#
# Usage:
#   run_fast_projection.sh <data_dir> <working_dir> <output_dir>
#
# Same I/O contract as run_flow_projection.sh:
#   <output_dir>/study_pca_projection.tsv
#   <output_dir>/qc_report.txt
#   <output_dir>/pca_projection.png

set -euo pipefail

DATA_DIR="${1:?missing data_dir}"
WORKING="${2:?missing working_dir}"
OUT_DIR="${3:?missing output_dir}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Mirror the passwd/group + HOME setup from the slow runner so any subprocess
# that touches the JVM (e.g. extract_loadings_matrix.py during first-time
# loadings export) also works.
USER_ID="$(id -u)"
GROUP_ID="$(id -g)"
if ! getent passwd "${USER_ID}" >/dev/null 2>&1 && [ -w /etc/passwd ]; then
    echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:${OUT_DIR}/.home:/bin/bash" >> /etc/passwd
fi
if ! getent group "${GROUP_ID}" >/dev/null 2>&1 && [ -w /etc/group ]; then
    echo "biovault:x:${GROUP_ID}:" >> /etc/group
fi
HOME_DIR="${OUT_DIR}/.home"
mkdir -p "${HOME_DIR}" /tmp/.ivy2
export HOME="${HOME_DIR}"
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Duser.home=${HOME_DIR} -Divy.default.ivy.user.dir=/tmp/.ivy2"
export PYSPARK_SUBMIT_ARGS="--conf spark.jars.ivy=/tmp/.ivy2 pyspark-shell"

source /opt/conda/etc/profile.d/conda.sh
conda activate biovault_popgen

mkdir -p "${WORKING}" "${OUT_DIR}"

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
GENO="${GENO:-0.05}"
# Sample missingness is handled by fast_convert_ddna_to_plink.py as a
# projection-specific loading-overlap count. Do not apply --mind over all
# 76,399 gnomAD loading rows: SNP arrays only overlap ~9k of them, so the
# old --mind 0.1 denominator incorrectly removes every healthy array sample.
MIND="${MIND:-1}"
MAF="${MAF:-0.01}"
HWE="${HWE:-1e-4}"
MIN_GS="${MIN_GS:-0.15}"
EXPECTED_LOADINGS_OVERLAP="${EXPECTED_LOADINGS_OVERLAP:-8971}"
MIN_LOADINGS_RATIO="${MIN_LOADINGS_RATIO:-0.95}"

export LOADINGS_HT="${LOADINGS_HT:-/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht}"
LOADINGS_NPZ="${LOADINGS_NPZ:-/opt/biovault/reference/pca_loadings/loadings.npz}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# One-time numpy export of the loadings table. Idempotent — only runs when
# the .npz is missing.
if [ ! -s "${LOADINGS_NPZ}" ]; then
    log "Extracting loadings matrix (one-time) ..."
    NPZ_DIR="$(dirname "${LOADINGS_NPZ}")"
    if [ ! -w "${NPZ_DIR}" ]; then
        LOADINGS_NPZ="${OUT_DIR}/loadings.npz"
    fi
    python3 "${SCRIPT_DIR}/extract_loadings_matrix.py" \
        --ht "${LOADINGS_HT}" \
        --out "${LOADINGS_NPZ}"
fi
export LOADINGS_NPZ

if [ -n "${BV_PREBUILT_RAW:-}" ]; then
    # Fast path: study_raw bed/bim/fam already built by `bvs project-bed` (byte-for-byte
    # identical to fast_convert --loadings-npz). Skip the Python parse; QC/projection
    # below are unchanged so outputs match exactly.
    log "STEP 1: using prebuilt study_raw from bvs project-bed: ${BV_PREBUILT_RAW}.{bed,bim,fam}"
    cp "${BV_PREBUILT_RAW}.bed" "${WORKING}/study_raw.bed"
    cp "${BV_PREBUILT_RAW}.bim" "${WORKING}/study_raw.bim"
    cp "${BV_PREBUILT_RAW}.fam" "${WORKING}/study_raw.fam"
    _raw_dir="$(dirname "${BV_PREBUILT_RAW}")"
    [ -f "${_raw_dir}/errors.tsv" ] && cp "${_raw_dir}/errors.tsv" "${OUT_DIR}/errors.tsv"
else
    log "STEP 1: DDNA -> PLINK bed/bim/fam (vectorized, no tped intermediary)"
    python3 "${SCRIPT_DIR}/fast_convert_ddna_to_plink.py" \
        "${DATA_DIR}" "${WORKING}/study_raw" --min-gs "${MIN_GS}" --workers "${THREADS}" \
        --loadings-npz "${LOADINGS_NPZ}" \
        --expected-loadings-overlap "${EXPECTED_LOADINGS_OVERLAP}" \
        --min-loadings-ratio "${MIN_LOADINGS_RATIO}"
    [ -f "${WORKING}/errors.tsv" ] && cp "${WORKING}/errors.tsv" "${OUT_DIR}/errors.tsv"
    [ -f "${WORKING}/warnings.tsv" ] && cp "${WORKING}/warnings.tsv" "${OUT_DIR}/warnings.tsv"
fi

log "STEP 2: QC"
plink2 --bfile "${WORKING}/study_raw" --rm-dup exclude-all \
       --make-bed --out "${WORKING}/study_nodup" --threads "${THREADS}"
plink2 --bfile "${WORKING}/study_nodup" \
       --geno "${GENO}" --mind "${MIND}" --maf "${MAF}" --hwe "${HWE}" \
       --make-bed --out "${WORKING}/study_qc" --threads "${THREADS}"

python3 "${SCRIPT_DIR}/write_variant_filter_report.py" \
    --pre-prefix "${WORKING}/study_nodup" \
    --post-prefix "${WORKING}/study_qc" \
    --out "${OUT_DIR}/filtered.tsv" \
    --geno "${GENO}" \
    --maf "${MAF}" \
    --hwe "${HWE}"

{
    echo "=== QC Report ==="
    echo "Input SNPs: $(wc -l < "${WORKING}/study_raw.bim")"
    echo "Input samples: $(wc -l < "${WORKING}/study_raw.fam")"
    echo "Final SNPs: $(wc -l < "${WORKING}/study_qc.bim")"
    echo "Final samples: $(wc -l < "${WORKING}/study_qc.fam")"
    echo "geno=${GENO}, mind=${MIND}, maf=${MAF}, hwe=${HWE}"
    echo "expected_loadings_overlap=${EXPECTED_LOADINGS_OVERLAP}, min_loadings_ratio=${MIN_LOADINGS_RATIO}"
} > "${OUT_DIR}/qc_report.txt"

log "STEP 3: Project onto gnomAD PCA space (numpy)"
python3 "${SCRIPT_DIR}/fast_pca_project.py" "${WORKING}/study_qc" "${OUT_DIR}"

log "Done. Outputs in ${OUT_DIR}/"

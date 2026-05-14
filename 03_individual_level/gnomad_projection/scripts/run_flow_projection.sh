#!/usr/bin/env bash
# Streamlined gnomAD-projection-only pipeline for the BioVault popgen flow.
# Skips per-chromosome reference VCF downloads and ADMIXTURE; reads the baked
# gnomAD v3.1 PCA loadings HT to project study samples directly.
#
# Usage:
#   run_flow_projection.sh <data_dir> <working_dir> <output_dir>
#
# Inputs:
#   data_dir: directory whose immediate subdirs are per-participant folders,
#             each containing one DDNA *_GSAv3-DTC_GRCh38*.txt file.
#
# Outputs (in output_dir):
#   study_pca_projection.tsv
#   qc_report.txt
#   hail.log

set -euo pipefail

DATA_DIR="${1:?missing data_dir}"
WORKING="${2:?missing working_dir}"
OUT_DIR="${3:?missing output_dir}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Nextflow runs docker with `-u <host_uid>:<host_gid>`. That host UID has no
# entry in the container's /etc/passwd, which makes JAAS's UnixLoginModule
# return a null username and Spark's Hadoop UGI then throws a Kerberos NPE on
# init. Add a passwd/group entry for the runtime UID before anything JVM runs.
USER_ID="$(id -u)"
GROUP_ID="$(id -g)"
if ! getent passwd "${USER_ID}" >/dev/null 2>&1; then
    echo "biovault:x:${USER_ID}:${GROUP_ID}:biovault:${OUT_DIR}/.home:/bin/bash" >> /etc/passwd
fi
if ! getent group "${GROUP_ID}" >/dev/null 2>&1; then
    echo "biovault:x:${GROUP_ID}:" >> /etc/group
fi
export HADOOP_USER_NAME="${HADOOP_USER_NAME:-biovault}"
export USER="${USER:-biovault}"

source /opt/conda/etc/profile.d/conda.sh
conda activate biovault_popgen

# Nextflow runs docker with `-u <host_uid>:<host_gid>`. The host UID typically
# has no entry in the container's /etc/passwd, so the JVM resolves
# `user.home` to "?", and Spark/Ivy then dies with
#   basedir must be absolute: ?/.ivy2/local
# Anchor a writable home, point Ivy at /tmp, and override user.home for every
# JVM Spark spawns via JAVA_TOOL_OPTIONS.
HOME_DIR="${OUT_DIR}/.home"
mkdir -p "${HOME_DIR}" /tmp/.ivy2
export HOME="${HOME_DIR}"
export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Duser.home=${HOME_DIR} -Divy.default.ivy.user.dir=/tmp/.ivy2"
export PYSPARK_SUBMIT_ARGS="--conf spark.jars.ivy=/tmp/.ivy2 pyspark-shell"

mkdir -p "${WORKING}" "${OUT_DIR}"

THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"
GENO="${GENO:-0.05}"
MIND="${MIND:-0.1}"
MAF="${MAF:-0.01}"
HWE="${HWE:-1e-4}"
MIN_GS="${MIN_GS:-0.15}"

export LOADINGS_HT="${LOADINGS_HT:-/opt/biovault/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht}"
export LOADINGS_VARIANTS_TSV="${LOADINGS_VARIANTS_TSV:-/opt/biovault/reference/pca_loadings/loadings_variants.tsv}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "STEP 1: DDNA -> PLINK"
python3 "${SCRIPT_DIR}/convert_ddna_to_plink.py" \
    "${DATA_DIR}" "${WORKING}/study_raw" "${MIN_GS}"

WORKING="${WORKING}" python3 - <<'FILTER'
import os
from pathlib import Path
working = Path(os.environ["WORKING"])
inp = working / "study_raw.tped"
out = working / "study_raw.biallelic.tped"
n_total = n_keep = 0
with inp.open() as f, out.open("w") as fo:
    for line in f:
        n_total += 1
        parts = line.rstrip("\n").split()
        observed = {a for a in parts[4:] if a not in {"0", "N", "-", "."}}
        if len(observed) <= 2:
            fo.write(line)
            n_keep += 1
print(f"Filtered TPED: {n_keep}/{n_total} biallelic SNPs kept")
FILTER

plink2 --tped "${WORKING}/study_raw.biallelic.tped" \
       --tfam "${WORKING}/study_raw.tfam" \
       --make-bed --out "${WORKING}/study_raw" \
       --threads "${THREADS}" --allow-extra-chr
rm -f "${WORKING}/study_raw.tped"

log "STEP 2: QC"
plink2 --bfile "${WORKING}/study_raw" --rm-dup exclude-all \
       --make-bed --out "${WORKING}/study_nodup" --threads "${THREADS}"
plink2 --bfile "${WORKING}/study_nodup" \
       --geno "${GENO}" --mind "${MIND}" --maf "${MAF}" --hwe "${HWE}" \
       --make-bed --out "${WORKING}/study_qc" --threads "${THREADS}"

{
    echo "=== QC Report ==="
    echo "Input SNPs: $(wc -l < "${WORKING}/study_raw.bim")"
    echo "Input samples: $(wc -l < "${WORKING}/study_raw.fam")"
    echo "Final SNPs: $(wc -l < "${WORKING}/study_qc.bim")"
    echo "Final samples: $(wc -l < "${WORKING}/study_qc.fam")"
    echo "geno=${GENO}, mind=${MIND}, maf=${MAF}, hwe=${HWE}"
} > "${OUT_DIR}/qc_report.txt"

log "STEP 3: Project onto gnomAD PCA space"
python3 "${SCRIPT_DIR}/pca_project.py" "${WORKING}/study_qc" "${OUT_DIR}"

log "Done. Outputs in ${OUT_DIR}/"

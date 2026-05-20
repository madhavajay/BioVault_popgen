#!/usr/bin/env bash
# Population-level pipeline (BIOVAULT.md Step 3): per-country allele
# frequencies -> pairwise Weir & Cockerham 1984 FST -> AIMs differential
# SNPs, driven by island_mapping.tsv (the `country`/`sex`-style facet).
#
# Stage A : 04_run_allele_freq.sh  (host bvs 0.1.22, DDNA + Illumina,
#           .bvlr-cached) -> raw_allele_freq_country/allele_freq_<pop>.tsv
# Stage B : FST 01-03 + AIMs 04-06 inside biovault-popgen (BV_POPULATIONS
#           unset -> auto-discovers the per-pop AF files).
#
# Usage:
#   bash 04_population_level.sh                 # all participants in mapping
#   bash 04_population_level.sh --limit 50      # cap per-island participants
#
# Overrides: IMAGE, MAPPING, DATA_DIR, RESULTS_ROOT (see 04_run_allele_freq.sh).

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${IMAGE:-biovault-popgen:0.1.1}"
RESULTS_ROOT="${RESULTS_ROOT:-${ROOT_DIR}/results}"
TASK_DIR="${RESULTS_ROOT}/population_level"
RAW_DIR="${TASK_DIR}/raw_allele_freq_country"
WORK="${TASK_DIR}/work"
RES="${TASK_DIR}"
BVLR_DIR="${TASK_DIR}/.bvlr_cache"

# `clean` wipes the per-participant .bvlr cache + per-country AF + FST/AIMs
# work/results, forcing a full recompute. Use after a regenerated cohort or
# a suspected-corrupt cache.
if [ "${1:-}" = "clean" ]; then
    echo "== clean: removing ${TASK_DIR} =="
    rm -rf "${TASK_DIR}"
    echo "cleaned. re-run without 'clean' to recompute."
    exit 0
fi

echo "== Stage A: per-country allele frequencies =="
OUT_DIR="${RAW_DIR}" BVLR_DIR="${BVLR_DIR}" bash "${ROOT_DIR}/04_run_allele_freq.sh" "$@"
python3 - "$ROOT_DIR" "$RAW_DIR" <<'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
raw = Path(sys.argv[2])
mapping = Path(__import__("os").environ.get(
    "MAPPING", root / "01_mock_data_generation/output/island_mapping.tsv"
))

def norm(raw_label: str) -> str:
    table = {
        "BVI": "bvi",
        "TT": "tt",
        "Bahamas": "bahamas",
        "Barbados": "barbados",
        "Bermuda": "bermuda",
        "StLucia": "stlucia",
    }
    if raw_label in table:
        return table[raw_label]
    import re
    return re.sub(r"[^a-z0-9]+", "_", raw_label.strip().lower()).strip("_")

lines = ["participant_id\tcountry"]
for ln in mapping.read_text().splitlines()[1:]:
    if not ln.strip():
        continue
    pid, country = ln.split("\t")[:2]
    lines.append(f"{pid}\t{norm(country)}")
(raw / "country_map.tsv").write_text("\n".join(lines) + "\n")
PY

echo
echo "== Stage B: FST + AIMs (${IMAGE}) =="
rm -rf "${WORK}"; mkdir -p "${WORK}" "${RES}"
docker run --rm \
    --platform linux/amd64 \
    -u "$(id -u):$(id -g)" \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    "${IMAGE}" bash -c '
        set -euo pipefail
        if ! getent passwd "$(id -u)" >/dev/null 2>&1; then
            echo "biovault:x:$(id -u):$(id -g):biovault:/tmp:/bin/bash" >> /etc/passwd
        fi
        export HOME=/tmp
        SD="'"${ROOT_DIR}"'/04_population_level/fst_aims_fast/scripts"
        export PYTHONPATH="$SD"
        export BV_RAW_DIR="'"${RAW_DIR}"'"
        unset BV_POPULATIONS                       # auto-discover allele_freq_*.tsv
        FST="'"${WORK}"'/fst"; AIM="'"${WORK}"'/aims"
        mkdir -p "$FST" "$AIM"
        source /opt/conda/etc/profile.d/conda.sh
        conda activate biovault_popgen
        echo "-- FST 01 load/merge --";  BV_WORK_DIR="$FST" python3 "$SD/fst_01_load_merge.py"
        echo "-- FST 02 compute_fst --"; BV_WORK_DIR="$FST" python3 "$SD/fst_02_compute_fst.py"
        echo "-- FST 03 visualize --";   BV_WORK_DIR="$FST" python3 "$SD/fst_03_visualize.py"
        echo "-- AIMs 04 merge --";      BV_WORK_DIR="$AIM" BV_FST_DIR="$FST" python3 "$SD/aims_04_merge.py"
        echo "-- AIMs 05 diff snps --";  BV_WORK_DIR="$AIM" python3 "$SD/aims_05_diff_snps.py"
        echo "-- AIMs 06 panels --";     BV_WORK_DIR="$AIM" python3 "$SD/aims_06_dendrogram.py"
    '

# Hoist the headline artefacts.
cp "${WORK}/fst/data/fst/fst_matrix.tsv"                       "${RES}/" 2>/dev/null || true
cp "${WORK}/fst/data/merged/merged_allele_freq_annotated.tsv"  "${RES}/" 2>/dev/null || true
cp "${WORK}/aims/data/master_af_table.tsv"                     "${RES}/" 2>/dev/null || true
cp "${WORK}/aims/data/differential_snps/all_outliers_long.tsv" "${RES}/" 2>/dev/null || true
cp "${WORK}/aims/data/aims/aims_combined.tsv"                  "${RES}/" 2>/dev/null || true
cp "${RAW_DIR}/country_map.tsv"                                "${RES}/" 2>/dev/null || true
{
    POPS="$(find "${RAW_DIR}" -maxdepth 1 -name 'allele_freq_*.tsv' -exec basename {} .tsv \; \
        | sed 's/^allele_freq_//' | sort | paste -sd, -)"
    echo "Populations: ${POPS}"
    echo ""
    echo "=== FST matrix ==="
    cat "${WORK}/fst/data/fst/fst_matrix.tsv"
    echo ""
    echo "=== master_af_table summary ==="
    cat "${WORK}/aims/data/master_af_table_summary.txt"
} > "${RES}/population_level_summary.txt" 2>/dev/null || true
cp "${WORK}"/fst/plots/*.png "${WORK}"/aims/plots/*.png "${WORK}"/aims/plots/*.pdf "${RES}/" 2>/dev/null || true

echo
echo "=== population-level outputs -> ${RES} ==="
for f in fst_matrix.tsv merged_allele_freq_annotated.tsv master_af_table.tsv \
         all_outliers_long.tsv aims_combined.tsv; do
    p="${RES}/${f}"
    if [ -e "${p}" ]; then printf '  %s  (%s)\n' "${f}" "$(du -h "${p}" | cut -f1)"
    else printf '  %s  (missing)\n' "${f}"; fi
done

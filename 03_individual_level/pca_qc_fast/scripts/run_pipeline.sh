#!/usr/bin/env bash
# =============================================================================
# Fast PCA/QC runner.
#
# Usage:
#   bash run_pipeline.sh
#   PYTHON=/path/to/python bash run_pipeline.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON="${PYTHON:-python3}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$BASE_DIR/.cache/matplotlib}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$BASE_DIR/.cache}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

echo "============================================="
echo "  Fast Ancestry PCA Pipeline"
echo "============================================="
echo ""

"$PYTHON" "$SCRIPT_DIR/fast_pipeline.py"

echo ""
echo "============================================="
echo "  Pipeline complete."
echo "  PCA plot: $BASE_DIR/plots/pca_pc1_pc2.png"
echo "============================================="

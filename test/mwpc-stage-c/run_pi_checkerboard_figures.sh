#!/usr/bin/env bash
# PI deliverables from the existing full checkerboard Stage-C projection.
# Run from anywhere:
#   bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh
# Alternative position of MS2:
#   bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh --ms2-z 550 --out test_runs/mwpc_stage_c/pi_checkerboard_figures_ms2_550
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO"

python3 "$SCRIPT_DIR/test_checkerboard_projection.py"
python3 "$SCRIPT_DIR/test_pi_checkerboard_figures.py"

python3 -c "import matplotlib, numpy" || {
  echo "ERROR: matplotlib and numpy required. Install: python3 -m pip install matplotlib numpy" >&2
  exit 1
}

python3 "$SCRIPT_DIR/make_pi_checkerboard_figures.py" "$@"

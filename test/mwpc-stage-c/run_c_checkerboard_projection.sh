#!/usr/bin/env bash
# Project saved GEANT muon reference states onto the 298-chamber checkerboard.
# Usage: bash test/mwpc-stage-c/run_c_checkerboard_projection.sh [smoke|full]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO"

mode="smoke"
if [[ $# -gt 0 ]]; then mode="$1"; shift; fi

python3 "$SCRIPT_DIR/test_checkerboard_projection.py"

case "$mode" in
  smoke)
    echo "Stage-C checkerboard: 10,000 reference rows/source (smoke)"
    python3 "$SCRIPT_DIR/mwpc_checkerboard_projection.py" \
      --limit-per-source 10000 --no-plots \
      --out "$REPO/test_runs/mwpc_stage_c/checkerboard_smoke" "$@"
    ;;
  full)
    echo "Stage-C checkerboard: all historical muon reference rows"
    python3 "$SCRIPT_DIR/mwpc_checkerboard_projection.py" \
      --write-crossings \
      --out "$REPO/test_runs/mwpc_stage_c/checkerboard_full" "$@"
    ;;
  *)
    echo "Usage: bash test/mwpc-stage-c/run_c_checkerboard_projection.sh [smoke|full]" >&2
    exit 2
    ;;
esac

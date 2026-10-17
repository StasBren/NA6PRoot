#!/usr/bin/env bash
# PI: large, resumable Heed/Garfield wire-collection response per chamber.
#
# Usage:
#   bash test/mwpc-stage-c/run_pi_per_chamber_response.sh smoke
#   bash test/mwpc-stage-c/run_pi_per_chamber_response.sh full
#   bash test/mwpc-stage-c/run_pi_per_chamber_response.sh full
#
# Repeated "full" calls continue from checkpoint. No Geant4 is run.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi

python3 -c "import numpy" || {
  echo "ERROR: numpy is required in the active Python environment" >&2
  exit 1
}

echo "Testing deterministic Heed/Garfield sampler and checkpoint/resume..."
python3 "$HERE/test_pi_per_chamber_response.py"

if [[ "$MODE" == "smoke" ]]; then
  echo "Smoke mode: up to 100 crossings, separate output directory."
  python3 "$HERE/run_pi_per_chamber_response.py" \
    --hours 0.1 --checkpoint-every 100 --max-per-source 100 \
    --out "$REPO/test_runs/mwpc_stage_c/pi_per_chamber_response_smoke" "$@"
elif [[ "$MODE" == "full" ]]; then
  echo "Full mode: up to 2.5 hours; rerun to resume if unfinished."
  python3 "$HERE/run_pi_per_chamber_response.py" \
    --hours 2.5 --checkpoint-every 10000 \
    --out "$REPO/test_runs/mwpc_stage_c/pi_per_chamber_response" "$@"
else
  echo "Usage: bash $0 [smoke|full]" >&2
  exit 2
fi

#!/usr/bin/env bash
# Produce PI module layout with all numbered anodes from existing CSV.
# Usage: bash test/mwpc-stage-c/run_pi_station_numbered_anodes.sh [det_id] [source]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
DETECTOR="17"
SOURCE="jpsi40-bore"
if [[ $# -gt 0 ]]; then DETECTOR="$1"; shift; fi
if [[ $# -gt 0 ]]; then SOURCE="$1"; shift; fi
COUNT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_fixed_wire_counts.csv"
OUTPUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_explanation"
if [[ ! -s "$COUNT" ]]; then
  echo "Missing finished fixed-wire counts file: $COUNT" >&2
  exit 2
fi
python3 "$HERE/plot_pi_station_numbered_anodes.py" \
  --counts "$COUNT" --out "$OUTPUT" \
  --det-id "$DETECTOR" --source "$SOURCE"
echo "Done: $OUTPUT"

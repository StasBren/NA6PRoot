#!/usr/bin/env bash
# Recreate 20x30 Garfield Stage-C drift-time heatmaps from saved electron trials.
# Reuses hybrid boundary map where available, otherwise complete baseline map.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
HYBRID="$REPO/test_runs/mwpc_stage_c/c_transport_box_map_hybrid_boundary_2p4/c_transport_box_electrons.csv"
BASELINE="$REPO/test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/c_transport_box_electrons.csv"
if [[ -s "$HYBRID" ]]; then
  INPUT="$HYBRID"
elif [[ -s "$BASELINE" ]]; then
  INPUT="$BASELINE"
else
  echo "ERROR: missing both saved Garfield maps:" >&2
  echo "$HYBRID" >&2
  echo "$BASELINE" >&2
  exit 1
fi
python3 "$HERE/plot_c_drift_time_grid.py" \
  --input "$INPUT" \
  --out "$REPO/test_runs/mwpc_stage_c/c_drift_grid_figures" "$@"

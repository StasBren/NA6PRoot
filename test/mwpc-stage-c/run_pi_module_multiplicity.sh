#!/usr/bin/env bash
# PI fixed-wire per-DETECTOR 0/1/2/3/4/5+ probability report and refreshed UI.
# Reuses old GEANT-reference projections and precomputed Heed/Garfield maps.
# Does NOT rerun Geant4 or microscopic Garfield electron trajectories.
# Run after the long Geant4 campaign if competing for CPU resources.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
python3 -c "import numpy" >/dev/null
python3 "$HERE/test_pi_fixed_wire_maps.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps"
python3 "$HERE/pi_fixed_wire_occupancy.py" --out "$OUT"
python3 "$HERE/build_pi_station_dashboard.py" \
  --out "$OUT/pi_station_interactive.html"
echo "Per-module multiplicity:"
echo "  $OUT/pi_fixed_wire_module_multiplicity.csv"
echo "Interactive station application updated:"
echo "  $OUT/pi_station_interactive.html"

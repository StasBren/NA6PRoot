#!/usr/bin/env bash
# Rebuild PI fixed-anode occupancy atlas from existing checkerboard crossings
# + saved Heed and Garfield kernels. No Geant4 or Garfield drift simulation.
#
#   bash test/mwpc-stage-c/run_pi_fixed_wire_atlas.sh smoke
#   bash test/mwpc-stage-c/run_pi_fixed_wire_atlas.sh full
#
# full: 18 global station maps, 18 contact sheets, 894 chamber details,
#       9 physics comparisons, 3 raw CSVs and settings/limitations.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi

python3 -c "import numpy, matplotlib" || {
  echo "ERROR: numpy and matplotlib required in active Python" >&2
  exit 1
}
python3 "$HERE/test_pi_fixed_wire_maps.py"

OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps"
case "$MODE" in
  smoke)
    OUT="$OUT""_smoke"
    echo "SMOKE: process 200 chamber crossings per channel; skip detail sheets"
    python3 "$HERE/pi_fixed_wire_occupancy.py" \
      --max-per-source 200 --out "$OUT" "$@"
    python3 "$HERE/plot_pi_fixed_wire_atlas.py" \
      --input "$OUT" --no-details
    ;;
  full)
    echo "FULL: all 2.49 million projected chamber crossings, 298 chambers"
    python3 "$HERE/pi_fixed_wire_occupancy.py" \
      --out "$OUT" "$@"
    python3 "$HERE/plot_pi_fixed_wire_atlas.py" --input "$OUT"
    ;;
  *)
    echo "Usage: bash $0 [smoke|full]" >&2
    exit 2
    ;;
esac

echo
echo "COMPLETE: $OUT"
echo "  station overview : $OUT/station_overview/"
echo "  per chamber      : $OUT/per_chamber/"
echo "  contact sheets   : $OUT/contact_sheets/"
echo "  physical tests   : $OUT/physics_summary/"
echo "  per-wire CSV     : $OUT/pi_fixed_wire_counts.csv"

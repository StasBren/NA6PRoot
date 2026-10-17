#!/usr/bin/env bash
# Improved PI physics validation: actual XY crossing map + geometric-anode
# footprint + 3 Garfield interpretations, station and individual chamber.
# Reuses completed ~2.49M wire-response CSV. Does NOT re-run Heed/Garfield.
#
#  bash test/mwpc-stage-c/run_pi_wire_validation_v2.sh preview
#  bash test/mwpc-stage-c/run_pi_wire_validation_v2.sh full
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

MODE="preview"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
INPUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps"
CONTROLS="$INPUT/physical_validation"
PROJECTION="$REPO/test_runs/mwpc_stage_c/checkerboard_full"

if [[ "$MODE" != preview && "$MODE" != full ]]; then
  echo "Usage: bash $0 [preview|full]" >&2
  exit 2
fi
if [[ ! -s "$INPUT/pi_fixed_wire_counts.csv" ]]; then
  echo "ERROR: completed physical-wire sampling CSV missing:" >&2
  echo "  $INPUT/pi_fixed_wire_counts.csv" >&2
  echo "Run the prior full fixed-wire atlas once first." >&2
  exit 3
fi
for channel in jpsi40-bore phi40-bore omega40-bore; do
  if [[ ! -s "$PROJECTION/crossings_$channel""_ms2_530.csv" ]]; then
    echo "ERROR: missing source crossing data for $channel" >&2
    exit 4
  fi
done

mkdir -p "$CONTROLS"
echo "=== STEP 1: regression test with synthetic Heed/Garfield and PNG ==="
python3 "$HERE/test_pi_fixed_wire_maps.py"

echo "=== STEP 2: physical XY + geometric wire controls (no Heed/Garfield) ==="
python3 "$HERE/prepare_pi_wire_geometry_controls.py" \
  --crossings "$PROJECTION" --out "$CONTROLS"

echo "=== STEP 3: physically aligned comparison maps ==="
if [[ "$MODE" == preview ]]; then
  python3 "$HERE/plot_pi_wire_physical_validation_v2.py" \
    --input "$INPUT" --controls "$CONTROLS" --no-details
else
  python3 "$HERE/plot_pi_wire_physical_validation_v2.py" \
    --input "$INPUT" --controls "$CONTROLS"
fi

echo
echo "COMPLETED: $MODE"
echo "  $CONTROLS/atlas_v2/station_validation/"
echo "  $CONTROLS/atlas_v2/contact_sheets/"
echo "  $CONTROLS/atlas_v2/chamber_profiles/"
echo "  $CONTROLS/pi_geom_wire_counts.csv"
echo "  $CONTROLS/pi_crossing_density_xy.npz"

#!/usr/bin/env bash
# Run *new* physically centered 3+3mm Garfield anode timing-angle PILOT on
# real native Geant4 gas segments (the FIRST CHUNK only for each channel).
#
# First stage smoke: all model/source checks + 100 native gas segments/source.
# Second stage full: all 3,000 generated events/source, ~75k native gas hits.
#
# No new Geant4 and no new Garfield: reads already VALIDATED centered map.
# No strip output; no old 2+4 field interpolation.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
mode="${1:-smoke}"
python3 -m py_compile \
 "$HERE/pi_centered_native_timing_pilot.py" \
 "$HERE/test_pi_centered_native_timing_pilot.py"
python3 "$HERE/test_pi_centered_native_timing_pilot.py"
case "$mode" in
 smoke)
   echo "RUN centered pilot smoke: 100 TRUE Geant4 gas segments/source"
   python3 "$HERE/pi_centered_native_timing_pilot.py" \
    --max-per-source 100 \
    --out "$REPO/test_runs/mwpc_stage_c/pi_centered_timing_native_smoke" \
    --html "$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_native_timing_smoke.html"
 ;;
 all)
   echo "RUN centered pilot: all gas hits in FIRST 3,000 native generated events/source"
   python3 "$HERE/pi_centered_native_timing_pilot.py"
 ;;
 *)
   echo "Usage: bash $0 [smoke|all]" >&2
   exit 2
 ;;
esac
echo "PASS: centered 3+3mm Geant4+Heed+Garfield first-chunk timing-angle pilot $mode"

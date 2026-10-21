#!/usr/bin/env bash
# PI first native 10ns electron drift-time histograms for all 6 Garfield maps.
# Reads frozen Geant4 first-chunk TGeo-local exports and Heed calibration.
# Does NOT start ROOT/Geant4/Garfield microscopic transport.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

python3 -m py_compile \
  "$HERE/pi_six_geometry_native_timing.py" \
  "$HERE/pi_six_geometry_native_timing_html.py" \
  "$HERE/test_pi_six_geometry_native_timing.py"
python3 "$HERE/test_pi_six_geometry_native_timing.py"

if [[ $# -ge 1 && "$1" == "smoke" ]]; then
  if [[ $# -ne 1 ]]; then
    echo "Usage: $0 smoke | $0 full" >&2; exit 2
  fi
  echo "FIRST SMALL NATIVE TEST: first 200 real muon gas segments/channel."
  python3 "$HERE/pi_six_geometry_native_timing.py" \
      --limit-per-source 200 \
      --out "$REPO/test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/smoke"
  python3 "$HERE/pi_six_geometry_native_timing_html.py" \
      --input "$REPO/test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/smoke/pi_six_geometry_native_timing.json" \
      --out "$REPO/test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/smoke/pi_six_geometry_native_timing.html"
elif [[ $# -eq 0 || ($# -eq 1 && "$1" == "full") ]]; then
  echo "FULL NATIVE FIRST CHUNK: 3000 parent events PER Jpsi/Phi/Omega channel."
  echo "All real gas-hit segments; NO new Geant4/Garfield transport."
  python3 "$HERE/pi_six_geometry_native_timing.py"
  python3 "$HERE/pi_six_geometry_native_timing_html.py"
else
  echo "Usage: $0 smoke | $0 full" >&2
  exit 2
fi
echo "PASS: native 10ns drift histograms for 6 Garfield fields, station/chamber/anode."

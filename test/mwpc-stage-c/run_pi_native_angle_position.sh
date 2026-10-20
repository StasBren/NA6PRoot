#!/usr/bin/env bash
# Instant PI angle-vs-module XY map from already completed native Geant4
# all-201-chunk quality CSVs. No Geant4, Garfield or ROOT invocation.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile "$HERE/pi_native_angle_position.py" \
    "$HERE/test_pi_native_angle_position.py"
python3 "$HERE/test_pi_native_angle_position.py"
python3 "$HERE/pi_native_angle_position.py"
python3 "$HERE/refresh_pi_dashboard_navigation.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_angle_position.html"
echo "PASS: full Geant4 muon angle vs module position:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows path: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

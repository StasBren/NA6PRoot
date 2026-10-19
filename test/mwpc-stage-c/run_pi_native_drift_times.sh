#!/usr/bin/env bash
# Plot unbiased FIRST-CHUNK native Geant4->Heed->saved Garfield primary
# electron arrival time distributions, station and chamber by chamber.
# This is a pilot timing-window dependence, NOT the 200k/source full fold.
# No ROOT export, Geant4 transport or microscopic Garfield reruns.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile "$HERE/pi_native_drift_times.py" \
 "$HERE/test_pi_native_drift_times.py"
python3 "$HERE/test_pi_native_drift_times.py"
python3 "$HERE/pi_native_drift_times.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_drift_times.html"
echo "First-chunk timing analysis ready:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows path: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

#!/usr/bin/env bash
# PI dashboard UX upgrade using EXISTING successful Geant4/Heed/Garfield
# centered 3+3mm pilot. NO rerun of any transport or Heed resampling.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/build_pi_centered_timing_dashboard_v2.py" \
 "$HERE/test_pi_centered_timing_dashboard_v2.py"
echo "Checking all original 3+3 mm native physics totals and PI plot labels..."
python3 "$HERE/test_pi_centered_timing_dashboard_v2.py"
echo "Producing improved dashboard from PRE-EXISTING native 3+3 HTML JSON..."
python3 "$HERE/build_pi_centered_timing_dashboard_v2.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_timing_dashboard_v2.html"
echo "PASS: clickable XY, filled angular bins, 1D birth-coordinate time slices,"
echo "      explicit axes and physically explained wire multiplicity."
echo "OPEN:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows file: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

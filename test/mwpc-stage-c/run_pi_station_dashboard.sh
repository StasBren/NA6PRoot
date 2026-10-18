#!/usr/bin/env bash
# Generate the full Stage-C offline interactive station browser from the
# existing 3-channel/298-module Garfield and XY projection CSVs.
# Does not rerun Geant4, Heed, Garfield or the physics folding.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

python3 -c "import numpy" || {
 echo "ERROR: need Python 3 and numpy" >&2
 exit 1
}
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_station_interactive.html"
python3 "$HERE/build_pi_station_dashboard.py" --out "$OUT"
echo
echo "Open this HTML in a browser:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
 echo "Windows path: $(wslpath -w "$OUT")"
 if command -v explorer.exe >/dev/null 2>&1; then
  explorer.exe "$(wslpath -w "$OUT")" || true
 fi
fi

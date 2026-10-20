#!/usr/bin/env bash
# Update navigation across already-existing offline PI dashboards without
# re-running Geant4, ROOT, Heed, Garfield or changing underlying JSON data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/refresh_pi_dashboard_navigation.py" \
 "$HERE/test_refresh_pi_dashboard_navigation.py"
python3 "$HERE/test_refresh_pi_dashboard_navigation.py"
python3 "$HERE/refresh_pi_dashboard_navigation.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_latest_dashboards.html"
echo "PASS: current PI dashboards now crosslink latest valid HTML only"
echo "START PAGE: $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows path: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

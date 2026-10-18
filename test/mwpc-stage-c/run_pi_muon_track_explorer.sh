#!/usr/bin/env bash
# Offline real-event Muon Track Explorer, with saved Heed/Garfield single
# electron transport only for deterministically selected projected gas hits.
#
# Does NOT rerun Geant4 or Garfield.
# While the large native Geant4 campaign is running:
#   bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh geometry-only
#
# Once enough spare CPU/RAM is available, use the richer collected-charge view:
#   bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
python3 -m py_compile \
  "$HERE/build_pi_muon_track_explorer.py" \
  "$HERE/test_pi_muon_track_explorer.py"
python3 "$HERE/test_pi_muon_track_explorer.py"
MODE="${1:-charges}"
if [[ "$MODE" == "geometry-only" && $# -le 1 ]]; then
  python3 "$HERE/build_pi_muon_track_explorer.py" --geometry-only
elif [[ "$MODE" == "charges" && $# -le 1 ]]; then
  python3 "$HERE/build_pi_muon_track_explorer.py"
else
  echo "Usage: bash $0 [charges|geometry-only]" >&2
  exit 2
fi
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_muon_track_explorer.html"
echo "Muon Track Explorer ready:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows path: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

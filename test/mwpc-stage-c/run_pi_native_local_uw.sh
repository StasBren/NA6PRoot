#!/usr/bin/env bash
# first: true native first-chunk CSVs, 3000 Geant4 events per source.
# all: each original ROOT chunk READ-ONLY, ROOT exporter + independent
#      quality checks with resumable checkpoints. NO new Geant4 or Garfield.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
source "$HERE/activate_pi_stagec_env.sh"
python3 -m py_compile \
 "$HERE/pi_native_local_uw_density.py" \
 "$HERE/pi_native_angle_position.py" \
 "$HERE/test_pi_native_local_uw_density.py" \
 "$HERE/test_pi_native_angle_position.py"
python3 "$HERE/test_pi_native_local_uw_density.py"
python3 "$HERE/test_pi_native_angle_position.py"
case "${1:-first}" in
 first|all)
   if [[ $# -gt 1 ]]; then echo "Usage: $0 [first|all]" >&2; exit 2; fi
   python3 "$HERE/pi_native_local_uw_density.py" "${1:-first}"
   ;;
 *)
   echo "Usage: $0 [first|all]" >&2; exit 2
   ;;
esac
python3 "$HERE/pi_native_angle_position.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_angle_position.html"
echo "PASS: global XY angular PDFs + selected chamber local u-w geometry,"
echo "  true native gas-hit midpoint 2D density and angle-cut 1D projections."
echo "Dashboard: $OUT"
if command -v wslpath >/dev/null 2>&1; then
 echo "Windows path: $(wslpath -w "$OUT")"
 if command -v explorer.exe >/dev/null 2>&1; then
  explorer.exe "$(wslpath -w "$OUT")" || true
 fi
fi

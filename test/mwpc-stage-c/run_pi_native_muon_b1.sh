#!/usr/bin/env bash
# True native Geant4 muons -> actual Heed clusters -> full microscopic
# Garfield ion births for four symmetric ±3mm chamber geometries.
# No original ROOT, drift LUT, archived gain runs or historical 2+4mm
# physics files are modified.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(build|prepare2|prepare8|check2|check8|smoke2|smoke8|status)$ ]];then
  echo "Usage: bash $0 build|prepare2|prepare8|check2|check8|smoke2|smoke8|status" >&2
  exit 2
fi
PY=python3
DRIVER="$HERE/pi_run_native_muon_b1.py"
"$PY" -m py_compile \
    "$HERE/pi_prepare_muon_avalanche_b1.py" \
    "$DRIVER" \
    "$HERE/test_pi_native_muon_b1.py"
"$PY" "$HERE/test_pi_native_muon_b1.py"

if [[ "$MODE" == "status" ]];then
  "$PY" "$DRIVER" status
  exit
fi
if [[ "$MODE" == "prepare2" || "$MODE" == "prepare8" ]];then
  "$PY" "$DRIVER" prepare
  exit
fi

if [[ -f "$HOME/na6p/env.sh" ]];then
  # shellcheck source=/dev/null
  source "$HOME/na6p/env.sh"
fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  # shellcheck source=/dev/null
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi

if [[ "$MODE" == "build" ]];then
  cmake -S "$REPO/test/mwpc-avalanche" \
        -B "$REPO/test/mwpc-avalanche/build" \
        -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
  cmake --build "$REPO/test/mwpc-avalanche/build" \
        --target mwpc_phase_b_native_muon_3p3 -j2
  echo "PASS target compiled; select original Geant4 muons next with 'prepare2'."
  exit
fi
if [[ "$MODE" == "check2" || "$MODE" == "check8" ]];then
  "$PY" "$DRIVER" "$MODE"
  exit
fi
"$PY" "$DRIVER" "check${MODE#smoke}"
nice -n 15 "$PY" "$DRIVER" "$MODE"

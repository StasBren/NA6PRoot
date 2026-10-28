#!/usr/bin/env bash
# Reuse the original Garfield++ Phase-A avalanche engine with symmetric
# 3+3 mm cathode gaps for the four pitch/FW geometries.
# NO old 2+4 outputs overwritten, NO full muon physics in this pilot.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(build|check|smoke5|pilot100|status)$ ]];then
  echo "Usage: bash $0 build|check|smoke5|pilot100|status" >&2
  exit 2
fi
PY=python3
SCRIPT="$HERE/pi_avalanche_3p3_gain_campaign.py"
"$PY" -m py_compile "$SCRIPT"
if [[ "$MODE" == status ]];then
  "$PY" "$SCRIPT" status
  exit
fi
# Ensure Garfield is available in WSL for standalone Phase A.
if [[ -f "$HOME/na6p/env.sh" ]];then
  # shellcheck source=/dev/null
  source "$HOME/na6p/env.sh"
fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  # shellcheck source=/dev/null
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
if [[ "$MODE" == build ]];then
  cmake -S "$REPO/test/mwpc-avalanche" \
        -B "$REPO/test/mwpc-avalanche/build" \
        -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
  cmake --build "$REPO/test/mwpc-avalanche/build" \
        --target mwpc_phase_a_gain -j2
  "$PY" "$SCRIPT" check
  exit
fi
"$PY" "$SCRIPT" check
if [[ "$MODE" == check ]];then exit;fi
# The source/rebuilt binary and raw original runs remain untouched.
nice -n 15 "$PY" "$SCRIPT" "$MODE"

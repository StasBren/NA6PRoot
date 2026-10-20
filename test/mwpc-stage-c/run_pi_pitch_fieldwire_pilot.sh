#!/usr/bin/env bash
# Stage-C NEW independent Garfield-only B=0 +3/-3 mm pitch/field-wire pilot.
# Does not modify frozen Geant4 or existing verified 4mm LUT.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

MODE=smoke
FIELD_DIAM_UM=""
SAMPLES=""
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 ]]; then FIELD_DIAM_UM="$1"; shift; fi
if [[ $# -gt 0 ]]; then SAMPLES="$1"; shift; fi
if [[ $# -gt 0 ]]; then echo "Too many arguments" >&2; exit 2; fi
if [[ "$MODE" != "smoke" && "$MODE" != "pilot" ]]; then
  echo "Usage: bash $0 smoke|pilot FIELD_WIRE_DIAM_UM [SAMPLES]" >&2
  exit 2
fi
if [[ -z "$FIELD_DIAM_UM" ]]; then
  echo "An explicit field-wire diameter in micrometers is required." >&2
  echo "Example exploratory parameter only: 100 um (not verified hardware)." >&2
  exit 2
fi
if [[ -z "$SAMPLES" ]]; then
  if [[ "$MODE" == "smoke" ]]; then SAMPLES=2
  else SAMPLES=20
  fi
fi

python3 -m py_compile "$HERE/pi_pitch_fieldwire_pilot.py" \
 "$HERE/test_pi_pitch_fieldwire_pilot.py"
python3 "$HERE/test_pi_pitch_fieldwire_pilot.py"

source "$HERE/activate_pi_stagec_env.sh"
cmake -S "$HERE" -B "$HERE/build" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$HERE/build" -j4 \
  --target mwpc_stage_c_transport_box_map

echo "NEW Garfield physics calculation: six independent configs."
echo "p=3,3.5,4mm; with/without grounded midpoint field wires; B=0."
echo "Field-wire diameter is a provisional, explicit input: $FIELD_DIAM_UM um."
echo "Mode=$MODE; trials/box=$SAMPLES."
python3 "$HERE/pi_pitch_fieldwire_pilot.py" \
 --mode "$MODE" --samples "$SAMPLES" \
 --field-wire-diam-um "$FIELD_DIAM_UM"
echo "PASS: new independent pitch/field-wire $MODE finished."

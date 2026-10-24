#!/usr/bin/env bash
# Stage-C NEW independent Garfield-only B=0 +3/-3 mm pitch/field-wire pilot.
# Does not modify frozen Geant4 or existing verified 4mm LUT.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

MODE=smoke
FIELD_DIAM_UM="75"
SAMPLES=""
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 ]]; then FIELD_DIAM_UM="$1"; shift; fi
if [[ $# -gt 0 ]]; then SAMPLES="$1"; shift; fi
if [[ $# -gt 0 ]]; then echo "Too many arguments" >&2; exit 2; fi
if [[ "$MODE" != "smoke" && "$MODE" != "pilot" ]]; then
  echo "Usage: bash $0 smoke|pilot FIELD_WIRE_DIAM_UM [SAMPLES]" >&2
  exit 2
fi
# Study default: 75 um field-wire diameter, requested by the PI workflow.
# Users may still supply an explicit replacement as the second argument.
if [[ -z "$SAMPLES" ]]; then
  if [[ "$MODE" == "smoke" ]]; then SAMPLES=2
  else SAMPLES=20
  fi
fi

python3 -m py_compile "$HERE/pi_pitch_fieldwire_pilot.py" \
 "$HERE/test_pi_pitch_fieldwire_pilot.py"
python3 "$HERE/test_pi_pitch_fieldwire_pilot.py"

source "$HERE/activate_pi_stagec_env.sh"
# Geant4 activation is NOT Garfield/ROOT activation. The binary builds even
# when libGarfield.so cannot be found by the runtime ELF loader (exit 127).
GARFIELD_PREFIX="$HOME/na6p/install/garfieldpp"
GARFIELD_SETUP="$GARFIELD_PREFIX/share/Garfield/setupGarfield.sh"
if [[ -f "$GARFIELD_SETUP" ]]; then
  # Garfield's generated setup may reference optional variables with no
  # default; retain this script's normal nounset checks afterward.
  set +u
  source "$GARFIELD_SETUP"
  set -u
fi
if [[ -d "$GARFIELD_PREFIX/lib" ]]; then
  export LD_LIBRARY_PATH="$GARFIELD_PREFIX/lib:${LD_LIBRARY_PATH:-}"
fi
if [[ -d "$GARFIELD_PREFIX/lib64" ]]; then
  export LD_LIBRARY_PATH="$GARFIELD_PREFIX/lib64:${LD_LIBRARY_PATH:-}"
fi
# ROOT's root-config may not be in PATH even when ROOT is installed.
# On this Stage-C workstation, ROOT is under ~/na6p/install/root/lib.
# Resolve the real libraries directly rather than trusting PATH alone.
ROOT_LIBDIR=""
ROOT_CANDIDATES=()
if [[ -n "${NA6P_ROOT_LIBDIR:-}" ]]; then
  ROOT_CANDIDATES+=("$NA6P_ROOT_LIBDIR")
fi
if [[ -n "${ROOTSYS:-}" ]]; then
  ROOT_CANDIDATES+=("$ROOTSYS/lib" "$ROOTSYS/lib64")
fi
ROOT_CANDIDATES+=(
  "$HOME/na6p/install/root/lib"
  "$HOME/na6p/install/root/lib64"
)
if command -v root-config >/dev/null 2>&1; then
  ROOT_CANDIDATES+=("$(root-config --libdir)")
elif [[ -x "$HOME/na6p/install/root/bin/root-config" ]]; then
  ROOT_CANDIDATES+=("$("$HOME/na6p/install/root/bin/root-config" --libdir)")
fi
for CANDIDATE in "${ROOT_CANDIDATES[@]}"; do
  if [[ -f "$CANDIDATE/libRIO.so" &&
        -f "$CANDIDATE/libTree.so" &&
        -f "$CANDIDATE/libCore.so" ]]; then
    ROOT_LIBDIR="$(cd "$CANDIDATE" && pwd -P)"
    break
  fi
done
if [[ -z "$ROOT_LIBDIR" ]]; then
  echo "ERROR: ROOT shared libraries (libRIO/libTree/libCore) not found." >&2
  echo "Searched ROOT library dirs:" >&2
  printf '  %s\\n' "${ROOT_CANDIDATES[@]}" >&2
  echo "Set NA6P_ROOT_LIBDIR=/absolute/path/to/root/lib if needed." >&2
  exit 127
fi
export LD_LIBRARY_PATH="$ROOT_LIBDIR:${LD_LIBRARY_PATH:-}"
echo "ROOT runtime libraries: $ROOT_LIBDIR"
cmake -S "$HERE" -B "$HERE/build" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$HERE/build" -j4 \
  --target mwpc_stage_c_transport_box_map

BINARY="$HERE/build/mwpc_stage_c_transport_box_map"
if [[ ! -x "$BINARY" ]]; then
  echo "ERROR: Garfield executable is missing: $BINARY" >&2
  exit 127
fi
echo "Garfield dynamic linker preflight:"
ldd "$BINARY"
UNRESOLVED="$(ldd "$BINARY" | grep 'not found' || true)"
if [[ -n "$UNRESOLVED" ]]; then
  echo "ERROR: Garfield shared libraries unresolved:" >&2
  echo "$UNRESOLVED" >&2
  echo "Garfield setup: $GARFIELD_SETUP" >&2
  echo "ROOT runtime libdir: $ROOT_LIBDIR" >&2
  exit 127
fi
echo "NEW Garfield physics calculation: six independent configs."
echo "p=3,3.5,4mm; with/without grounded midpoint field wires; B=0."
echo "Field-wire diameter: $FIELD_DIAM_UM um (requested study setting)."
echo "Mode=$MODE; trials/box=$SAMPLES."
# Keep historical defaults unchanged; new high-stat LUTs require an explicit
# namespace and seed to avoid modifying/overwriting the frozen 20-trial maps.
EXTRA_ARGS=()
if [[ -n "${NA60_PI_PITCH_LUT_OUT:-}" ]]; then
  EXTRA_ARGS+=(--out "$NA60_PI_PITCH_LUT_OUT")
fi
if [[ -n "${NA60_PI_PITCH_LUT_SEED:-}" ]]; then
  EXTRA_ARGS+=(--seed "$NA60_PI_PITCH_LUT_SEED")
fi
python3 "$HERE/pi_pitch_fieldwire_pilot.py" \
 --mode "$MODE" --samples "$SAMPLES" \
 --field-wire-diam-um "$FIELD_DIAM_UM" \
 "${EXTRA_ARGS[@]}"
echo "PASS: new independent pitch/field-wire $MODE finished."

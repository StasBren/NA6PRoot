#!/usr/bin/env bash
# Read-only frozen 2100V K256 muon -> provisional U0/V0 ion response.
# Does NOT touch original microscopic campaign, FEM meshes or past output.
# Usage: bash .../run_p3_2100_stereo_ion_response.sh preflight|pilot|muon
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
MODE="${1:-preflight}"
MUON="${MWPC_STEREO_EVENT:-0301}"
DATA="${MWPC_STEREO_DATA:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
MAP_RUN="${MWPC_STEREO_MAP_RUN:-$REPO/test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1}"
MOBILITY="${MWPC_STEREO_MOBILITY:-$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt}"
BUILD="$REPO/test/mwpc-stage-c/build-fem-strip-probe"
BIN="$BUILD/mwpc_stage_c_stereo_ion_response_2100"
OUTPUT="${MWPC_STEREO_OUTPUT:-$REPO/test_runs/mwpc_stage_c/p3_2100_stereo_ion_response_v1}"
SOURCE="$DATA/muon_$MUON"
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
fail() { echo "STEREO U/V STOP: $*" >&2; exit 2; }
[[ "$MUON" =~ ^[0-9]{4}$ ]] || fail "event must be four-digit original index"
case "$MODE" in preflight|pilot|muon) ;; *) fail "usage preflight|pilot|muon";; esac

preflight() {
  [[ -d "$SOURCE" && -s "$SOURCE/ion_births.root" &&
     -s "$SOURCE/complete.txt" && -s "$SOURCE/SHA256SUMS" ]] ||
      fail "source muon_$MUON is not SHA-sealed"
  grep -Fxq 'PASS CATHODE ROOT PHYSICS AUDIT' "$SOURCE/root_audit.log" ||
      fail "microscopic audit not passed"
  grep -Fxq 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$SOURCE/physics.log" ||
      fail "microscopic physics not verified uncapped"
  (cd "$SOURCE" && sha256sum -c SHA256SUMS >/dev/null) ||
      fail "frozen source muon checksums failed"
  [[ -s "$MAP_RUN/SHA256SUMS" ]] || fail "FEM run not completed"
  (cd "$MAP_RUN" && sha256sum -c SHA256SUMS >/dev/null) ||
      fail "FEM map checksum failure"
  [[ -s "$MOBILITY" ]] || fail "ion mobility file absent"
  python3 - "$MAP_RUN/input_config.json" <<'PY'
import json, sys
c=json.load(open(sys.argv[1]))
assert c["geometry"]["half_wires"]==4
assert c["geometry"]["gap_plus_mm"]==c["geometry"]["gap_minus_mm"]==3
assert c["geometry"]["wire_pitch_mm"]==3
assert c["geometry"]["wire_diameter_um"]==30
assert c["electrostatics"]["anode_voltage_V"]==2100
assert c["weighting"]["stereo_sign_convention"]=="plus_negative_u"
assert abs(c["geometry"]["strip_pitch_mm"]-1.6915632233569817)<1.e-10
assert abs(c["geometry"]["strip_width_mm"]-1.6915632233569817)<1.e-10
assert abs(c["geometry"]["tan_alpha"]-.1)<1.e-12
print("PASS exact candidate P3 2100V U/V geometry config")
PY
  echo "PASS READ-ONLY SOURCE+MAP preflight"
  echo "source=$SOURCE/ion_births.root"
  echo "map=$MAP_RUN/map"
  echo "PROVISIONAL: FEM near-wire geometry/convergence not yet certified"
}

preflight
if [[ "$MODE" == preflight ]]; then
  exit 0
fi
[[ -x "$BIN" ]] || fail "compile $BIN first (stage-C CMake target)"
mkdir -p "$OUTPUT"
PREFIX="$OUTPUT/muon_${MUON}_${MODE}_central_uv"
for f in _summary.csv _channels.csv _waveforms.csv _pilot.txt _complete.txt; do
  [[ ! -e "$PREFIX$f" ]] || fail "output exists ($PREFIX$f); refusing overwrite"
done
echo "BEGIN $MODE $(date -Is) source=muon_$MUON"
args=(--ion-file "$SOURCE/ion_births.root"
      --map-dir "$MAP_RUN/map"
      --mobility "$MOBILITY"
      --out-prefix "$PREFIX"
      --dt-ns 10 --window-us 100 --step-mm .05 --half-strips 0)
if [[ "$MODE" == pilot ]]; then args+=(--ion-limit 1);fi
nice -n 10 "$BIN" "${args[@]}"
echo "RESULT=$PREFIX"
echo "END $MODE $(date -Is)"

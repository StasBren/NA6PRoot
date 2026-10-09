#!/usr/bin/env bash
set -euo pipefail

# Full Stage-C2b event with numerical Elmer strip weighting.
#
# Usage:
#   bash test/mwpc-stage-c/run_c2_full_numerical.sh \
#     [momentum_GeV] [theta_u_deg] [theta_w_deg] [seed] [config_json] [events]

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
BUILD_DIR="${STAGE_DIR}/build"

P_GEV="${1:-5}"
THETA_U="${2:-10}"
THETA_W="${3:-25}"
SEED="${4:-120001}"
CONFIG="${5:-${NUM_DIR}/weighting_config_c2b.json}"
EVENTS="${6:-1}"

TAG="p${P_GEV}_tu${THETA_U}_tw${THETA_W}_seed${SEED}"
TAG="${TAG//./p}"
TAG="${TAG//-/m}"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_full_numerical/${TAG}"
MAP_DIR="${OUT_DIR}/map"
EVENT_DIR="${OUT_DIR}/event"

missing=0
for cmd in python3 gmsh ElmerGrid ElmerSolver cmake; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required command: ${cmd}"
    missing=1
  fi
done
if [[ "${missing}" -ne 0 ]]; then
  exit 2
fi

MOBILITY="${MWPC_ION_MOBILITY:-}"
if [[ -z "${MOBILITY}" ]]; then
  MOBILITY="$(find "${HOME}/na6p/install/garfieldpp" \
    -type f -name IonMobility_Ar+_Ar.txt -print -quit 2>/dev/null || true)"
fi
if [[ -z "${MOBILITY}" || ! -f "${MOBILITY}" ]]; then
  echo "Could not find IonMobility_Ar+_Ar.txt."
  echo "Set MWPC_ION_MOBILITY=/full/path/IonMobility_Ar+_Ar.txt"
  exit 2
fi

read_json() {
  python3 -c 'import json,sys
obj=json.load(open(sys.argv[1]))
for key in sys.argv[2].split("."):
    obj=obj[key]
print(obj)' "${CONFIG}" "$1"
}

GAP_MINUS="$(read_json geometry.gap_minus_mm)"
GAP_PLUS="$(read_json geometry.gap_plus_mm)"
WIRE_PITCH="$(read_json geometry.wire_pitch_mm)"
WIRE_DIAM="$(read_json geometry.wire_diameter_um)"
MAP_HALF_WIRES="$(read_json geometry.half_wires)"
STRIP_PITCH="$(read_json geometry.strip_pitch_mm)"
STRIP_WIDTH="$(read_json geometry.strip_width_mm)"
TAN_ALPHA="$(read_json geometry.tan_alpha)"
MAP_U_HALF="$(read_json geometry.u_half_span_mm)"
MAP_W_HALF="$(read_json geometry.w_half_span_mm)"
HALF_STRIPS="$(read_json weighting.translated_half_strips)"
HV="$(read_json electrostatics.anode_voltage_V)"

rm -rf "${OUT_DIR}"
mkdir -p "${MAP_DIR}" "${EVENT_DIR}"

echo "=== C2b numerical full-event setup ==="
echo "config       : ${CONFIG}"
echo "p            : ${P_GEV} GeV/c"
echo "theta_u      : ${THETA_U} deg"
echo "theta_w      : ${THETA_W} deg"
echo "seed         : ${SEED}"
echo "events       : ${EVENTS}"
echo "gap          : ${GAP_MINUS} + ${GAP_PLUS} mm"
echo "wire pitch   : ${WIRE_PITCH} mm"
echo "strip pitch  : ${STRIP_PITCH} mm"
echo "tan(alpha)   : ${TAN_ALPHA}"
echo "strip bank   : -${HALF_STRIPS} .. +${HALF_STRIPS}"
echo "FEM patch    : +/-${MAP_U_HALF} mm (u), +/-${MAP_W_HALF} mm (w)"
echo

python3 "${NUM_DIR}/generate_weighting_model.py" \
  --config "${CONFIG}" \
  --output "${MAP_DIR}"

cd "${MAP_DIR}"

echo
echo "=== Gmsh: C2b 3-D quadratic mesh ==="
gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh

echo
echo "=== ElmerGrid ==="
ElmerGrid 14 2 chamber.msh -autoclean

echo
echo "=== Elmer: reference physical map ==="
ElmerSolver drift.sif

echo
echo "=== Elmer: central minus/plus weighting maps ==="
ElmerSolver weight_minus_0.sif
ElmerSolver weight_plus_0.sif

echo
echo "=== Build C2 targets ==="
cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 \
  --target mwpc_stage_c2_weighting_probe mwpc_stage_c2_full_muon_numerical

echo
echo "=== Map sanity probe on the enlarged C2b patch ==="
"${BUILD_DIR}/mwpc_stage_c2_weighting_probe" \
  --map-dir "${MAP_DIR}" \
  --gap-minus-mm "${GAP_MINUS}" \
  --gap-plus-mm "${GAP_PLUS}" \
  --profile-half-w-mm 4.0 \
  --profile-step-mm 0.1

cd "${EVENT_DIR}"

echo
echo "=== Full muon with direct numerical strip weighting ==="
"${BUILD_DIR}/mwpc_stage_c2_full_muon_numerical" \
  --numerical-map-dir "${MAP_DIR}" \
  --map-u-half-mm "${MAP_U_HALF}" \
  --map-w-half-mm "${MAP_W_HALF}" \
  --map-wire-pitch-mm "${WIRE_PITCH}" \
  --map-wire-diam-um "${WIRE_DIAM}" \
  --map-half-wires "${MAP_HALF_WIRES}" \
  --gap-minus-mm "${GAP_MINUS}" \
  --gap-plus-mm "${GAP_PLUS}" \
  --wire-pitch-mm "${WIRE_PITCH}" \
  --wire-diam-um "${WIRE_DIAM}" \
  --hv "${HV}" \
  --strip-pitch-mm "${STRIP_PITCH}" \
  --strip-width-mm "${STRIP_WIDTH}" \
  --tan-alpha "${TAN_ALPHA}" \
  --half-strips "${HALF_STRIPS}" \
  --momentum-gev "${P_GEV}" \
  --u0-mm 1.0 \
  --w0-mm 0.0 \
  --theta-u-deg "${THETA_U}" \
  --theta-w-deg "${THETA_W}" \
  --b-tesla 0 \
  --observation-ns 100 \
  --ion-rk-dt-ns 5 \
  --half-wires 6 \
  --events "${EVENTS}" \
  --max-seeds 0 \
  --base-seed "${SEED}" \
  --ion-mobility "${MOBILITY}" \
  --output-prefix c2_numerical

echo
echo "C2b event outputs:"
echo "  ${EVENT_DIR}/c2_numerical_event_summary.csv"
echo "  ${EVENT_DIR}/c2_numerical_strip_summary.csv"
echo "  ${EVENT_DIR}/c2_numerical_seed_summary.csv"
echo "  ${EVENT_DIR}/c2_numerical_wire_summary.csv"

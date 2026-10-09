#!/usr/bin/env bash
set -euo pipefail

# C2c: scan the track intercept across one 4-mm wire cell.
#
# This deliberately starts at normal incidence to isolate the detector's
# across-wire transfer function from finite-gap angular centroid shifts.
#
# Usage:
#   bash test/mwpc-stage-c/run_c2_u_cell_scan.sh [events_per_point] [base_seed] [config_json]
#
# Default scan points (mm at v=0):
#   -1.5, -0.75, 0, +0.75, +1.5
#
# The same random-seed block is reused at each u0 point to make the scan as
# nearly paired as possible.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
BUILD_DIR="${STAGE_DIR}/build"

EVENTS="${1:-5}"
BASE_SEED="${2:-130001}"
CONFIG="${3:-${NUM_DIR}/weighting_config_c2b.json}"

U_POINTS=(-1.5 -0.75 0.0 0.75 1.5)

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_u_cell_scan"
MAP_DIR="${OUT_DIR}/map"
EVENT_ROOT="${OUT_DIR}/events"
COMBINED="${OUT_DIR}/c2_u_cell_scan_event_summary.csv"

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
mkdir -p "${MAP_DIR}" "${EVENT_ROOT}"

echo "=== C2c across-wire cell scan ==="
echo "config            : ${CONFIG}"
echo "events / point    : ${EVENTS}"
echo "base seed         : ${BASE_SEED}"
echo "u0 points [mm]    : ${U_POINTS[*]}"
echo "theta_u/theta_w   : 0 / 0 deg"
echo "wire pitch        : ${WIRE_PITCH} mm"
echo

# One numerical map for the whole scan.
python3 "${NUM_DIR}/generate_weighting_model.py" \
  --config "${CONFIG}" \
  --output "${MAP_DIR}"

cd "${MAP_DIR}"
gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh
ElmerGrid 14 2 chamber.msh -autoclean
ElmerSolver drift.sif
ElmerSolver weight_minus_0.sif
ElmerSolver weight_plus_0.sif

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 \
  --target mwpc_stage_c2_weighting_probe mwpc_stage_c2_full_muon_numerical

"${BUILD_DIR}/mwpc_stage_c2_weighting_probe" \
  --map-dir "${MAP_DIR}" \
  --gap-minus-mm "${GAP_MINUS}" \
  --gap-plus-mm "${GAP_PLUS}" \
  --profile-half-w-mm 4.0 \
  --profile-step-mm 0.1

rm -f "${COMBINED}"
first=1

for U0 in "${U_POINTS[@]}"; do
  tag="${U0//-/m}"
  tag="${tag//./p}"
  EVENT_DIR="${EVENT_ROOT}/u_${tag}"
  mkdir -p "${EVENT_DIR}"
  cd "${EVENT_DIR}"

  echo
  echo "=== u0 = ${U0} mm ==="

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
    --momentum-gev 5 \
    --u0-mm "${U0}" \
    --w0-mm 0 \
    --theta-u-deg 0 \
    --theta-w-deg 0 \
    --b-tesla 0 \
    --observation-ns 100 \
    --ion-rk-dt-ns 5 \
    --half-wires 6 \
    --events "${EVENTS}" \
    --max-seeds 0 \
    --base-seed "${BASE_SEED}" \
    --ion-mobility "${MOBILITY}" \
    --output-prefix c2_u_scan

  if [[ "${first}" -eq 1 ]]; then
    cat c2_u_scan_event_summary.csv > "${COMBINED}"
    first=0
  else
    tail -n +2 c2_u_scan_event_summary.csv >> "${COMBINED}"
  fi
done

echo
echo "C2c combined event summary:"
echo "  ${COMBINED}"
echo
echo "Per-position outputs:"
echo "  ${EVENT_ROOT}/u_*/c2_u_scan_*.csv"

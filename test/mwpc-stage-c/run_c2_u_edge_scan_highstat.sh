#!/usr/bin/env bash
set -euo pipefail

# High-statistics C2 scan focused on the four wire-cell boundary points where
# the 30-event dense scan showed rare wire-sharing tails and apparent +u/-u
# asymmetry.
#
# Points:
#   u0 = -1.9, -1.7, +1.7, +1.9 mm
#
# Default:
#   250 events / point = 1000 full microscopic events
#   independent seed block beginning at 150001
#
# Usage:
#   bash test/mwpc-stage-c/run_c2_u_edge_scan_highstat.sh \
#     [events_per_point] [base_seed] [config_json]
#
# The run is resumable by complete u0 point. An interrupted point is rerun
# from the beginning; completed points are skipped. The same seed block is
# reused at all four positions for lower-variance distribution comparisons.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
BUILD_DIR="${STAGE_DIR}/build"

EVENTS="${1:-250}"
BASE_SEED="${2:-150001}"
CONFIG="${3:-${NUM_DIR}/weighting_config_c2b.json}"

U_POINTS=(-1.9 -1.7 1.7 1.9)

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_u_edge_scan_highstat"
EVENT_ROOT="${OUT_DIR}/events"
RUN_LOG="${OUT_DIR}/scan_progress.log"
COMBINED_EVENT="${OUT_DIR}/c2_u_edge_event_summary.csv"
COMBINED_STRIP="${OUT_DIR}/c2_u_edge_strip_summary.csv"
COMBINED_WIRE="${OUT_DIR}/c2_u_edge_wire_summary.csv"

DENSE_MAP_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_u_cell_scan_dense/map"
PRIVATE_MAP_DIR="${OUT_DIR}/map"

mkdir -p "${OUT_DIR}" "${EVENT_ROOT}"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo
echo "================================================================"
echo "C2 high-stat boundary scan"
echo "started           : $(date --iso-8601=seconds)"
echo "config            : ${CONFIG}"
echo "points [mm]       : ${U_POINTS[*]}"
echo "events / point    : ${EVENTS}"
echo "base seed         : ${BASE_SEED}"
echo "total events      : $(( EVENTS * ${#U_POINTS[@]} ))"
echo "output            : ${OUT_DIR}"
echo "================================================================"

if ! [[ "${EVENTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "events_per_point must be a positive integer."
  exit 2
fi

missing=0
for cmd in python3 gmsh ElmerGrid ElmerSolver cmake sha256sum; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required command: ${cmd}"
    missing=1
  fi
done
[[ "${missing}" -eq 0 ]] || exit 2

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

CONFIG_HASH="$(sha256sum "${CONFIG}" | awk '{print $1}')"

map_matches() {
  local map_dir="$1"
  [[ -f "${map_dir}/chamber/drift.result" ]] || return 1
  [[ -f "${map_dir}/chamber/weight_minus_0.result" ]] || return 1
  [[ -f "${map_dir}/chamber/weight_plus_0.result" ]] || return 1
  [[ -f "${map_dir}/config.sha256" ]] || return 1
  [[ "$(cat "${map_dir}/config.sha256")" == "${CONFIG_HASH}" ]]
}

if map_matches "${DENSE_MAP_DIR}"; then
  MAP_DIR="${DENSE_MAP_DIR}"
  echo "Reusing dense-scan FEM map: ${MAP_DIR}"
elif map_matches "${PRIVATE_MAP_DIR}"; then
  MAP_DIR="${PRIVATE_MAP_DIR}"
  echo "Reusing edge-scan FEM map: ${MAP_DIR}"
else
  MAP_DIR="${PRIVATE_MAP_DIR}"
  echo "Building FEM map for edge scan..."
  rm -rf "${MAP_DIR}"
  mkdir -p "${MAP_DIR}"

  python3 "${NUM_DIR}/generate_weighting_model.py" \
    --config "${CONFIG}" \
    --output "${MAP_DIR}"

  cd "${MAP_DIR}"
  gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh
  ElmerGrid 14 2 chamber.msh -autoclean
  ElmerSolver drift.sif
  ElmerSolver weight_minus_0.sif
  ElmerSolver weight_plus_0.sif
  echo "${CONFIG_HASH}" > "${MAP_DIR}/config.sha256"
fi

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 \
  --target mwpc_stage_c2_weighting_probe mwpc_stage_c2_full_muon_numerical

"${BUILD_DIR}/mwpc_stage_c2_weighting_probe" \
  --map-dir "${MAP_DIR}" \
  --gap-minus-mm "${GAP_MINUS}" \
  --gap-plus-mm "${GAP_PLUS}" \
  --profile-half-w-mm 4.0 \
  --profile-step-mm 0.2

point_complete() {
  local csv="$1"
  local expected="$2"
  local u0="$3"

  [[ -f "${csv}" ]] || return 1

  python3 - "${csv}" "${expected}" "${u0}" "${BASE_SEED}" <<'PY'
import csv, sys
path, expected, u0, seed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
if len(rows) != expected or not rows:
    raise SystemExit(1)
if any(abs(float(r["u0_mm"]) - u0) > 1e-9 for r in rows):
    raise SystemExit(1)
seeds = [int(r["random_seed"]) for r in rows]
if seeds != list(range(seed, seed + expected)):
    raise SystemExit(1)
raise SystemExit(0)
PY
}

rebuild_combined_one() {
  local suffix="$1"
  local destination="$2"
  local first=1
  rm -f "${destination}"

  for U0 in "${U_POINTS[@]}"; do
    local tag="${U0//-/m}"
    tag="${tag//./p}"
    local event_dir="${EVENT_ROOT}/u_${tag}"
    local event_csv="${event_dir}/c2_u_dense_event_summary.csv"
    local src="${event_dir}/c2_u_dense_${suffix}_summary.csv"

    if point_complete "${event_csv}" "${EVENTS}" "${U0}" && [[ -f "${src}" ]]; then
      if [[ "${first}" -eq 1 ]]; then
        cat "${src}" > "${destination}"
        first=0
      else
        tail -n +2 "${src}" >> "${destination}"
      fi
    fi
  done
}

rebuild_combined() {
  rebuild_combined_one event "${COMBINED_EVENT}"
  rebuild_combined_one strip "${COMBINED_STRIP}"
  rebuild_combined_one wire "${COMBINED_WIRE}"
}

trap 'echo; echo "Edge scan interrupted at $(date --iso-8601=seconds). Re-run the same command to resume."; rebuild_combined || true' INT TERM EXIT

completed=0
for U0 in "${U_POINTS[@]}"; do
  tag="${U0//-/m}"
  tag="${tag//./p}"
  EVENT_DIR="${EVENT_ROOT}/u_${tag}"
  CSV="${EVENT_DIR}/c2_u_dense_event_summary.csv"

  if point_complete "${CSV}" "${EVENTS}" "${U0}"; then
    completed=$((completed + 1))
    echo
    echo "[${completed}/${#U_POINTS[@]}] u0=${U0} mm already complete -> skip"
    continue
  fi

  echo
  echo "----------------------------------------------------------------"
  echo "u0 = ${U0} mm"
  echo "started: $(date --iso-8601=seconds)"
  echo "----------------------------------------------------------------"

  rm -rf "${EVENT_DIR}"
  mkdir -p "${EVENT_DIR}"
  cd "${EVENT_DIR}"

  stdbuf -oL -eL "${BUILD_DIR}/mwpc_stage_c2_full_muon_numerical" \
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
    --output-prefix c2_u_dense

  if ! point_complete "${CSV}" "${EVENTS}" "${U0}"; then
    echo "Point u0=${U0} did not produce the expected complete CSV."
    exit 6
  fi

  completed=$((completed + 1))
  echo "completed: $(date --iso-8601=seconds)"
  echo "progress : ${completed}/${#U_POINTS[@]} points"
  rebuild_combined
done

rebuild_combined

python3 "${STAGE_DIR}/analyze_c2_edge_symmetry.py" \
  --scan-dir "${OUT_DIR}"

trap - INT TERM EXIT

echo
echo "================================================================"
echo "C2 high-stat boundary scan COMPLETE"
echo "finished          : $(date --iso-8601=seconds)"
echo "event summary     : ${COMBINED_EVENT}"
echo "strip summary     : ${COMBINED_STRIP}"
echo "wire summary      : ${COMBINED_WIRE}"
echo "progress log      : ${RUN_LOG}"
echo "================================================================"

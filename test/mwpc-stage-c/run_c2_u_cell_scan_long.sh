#!/usr/bin/env bash
set -euo pipefail

# Long/resumable C2 dense scan across one 4-mm wire cell.
#
# Default workload:
#   u0 = -1.9 .. +1.9 mm in 0.2-mm steps  -> 20 positions
#   30 full microscopic events / position -> 600 events total
#
# Usage:
#   bash test/mwpc-stage-c/run_c2_u_cell_scan_long.sh #     [events_per_point] [base_seed] [step_mm] [config_json]
#
# Optional environment variables:
#   C2_U_MIN_MM   (default -1.9)
#   C2_U_MAX_MM   (default +1.9)
#
# The script is resumable. A completed point is skipped on the next run.
# An interrupted/incomplete point is rerun from scratch. The same random-seed
# block is reused at every u0, which is useful for paired comparisons.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
BUILD_DIR="${STAGE_DIR}/build"

EVENTS="${1:-30}"
BASE_SEED="${2:-140001}"
STEP_MM="${3:-0.2}"
CONFIG="${4:-${NUM_DIR}/weighting_config_c2b.json}"

U_MIN_MM="${C2_U_MIN_MM:--1.9}"
U_MAX_MM="${C2_U_MAX_MM:-1.9}"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_u_cell_scan_dense"
MAP_DIR="${OUT_DIR}/map"
EVENT_ROOT="${OUT_DIR}/events"
COMBINED="${OUT_DIR}/c2_u_dense_event_summary.csv"
POINT_SUMMARY="${OUT_DIR}/c2_u_dense_point_summary.csv"
RUN_LOG="${OUT_DIR}/scan_progress.log"
CONFIG_HASH_FILE="${MAP_DIR}/config.sha256"

mkdir -p "${OUT_DIR}" "${EVENT_ROOT}"

exec > >(tee -a "${RUN_LOG}") 2>&1

echo
echo "================================================================"
echo "C2 dense across-wire scan"
echo "started           : $(date --iso-8601=seconds)"
echo "config            : ${CONFIG}"
echo "events / point    : ${EVENTS}"
echo "base seed         : ${BASE_SEED}"
echo "u0 range [mm]     : ${U_MIN_MM} .. ${U_MAX_MM}"
echo "step [mm]         : ${STEP_MM}"
echo "output            : ${OUT_DIR}"
echo "================================================================"

missing=0
for cmd in python3 gmsh ElmerGrid ElmerSolver cmake sha256sum; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required command: ${cmd}"
    missing=1
  fi
done
if [[ "${missing}" -ne 0 ]]; then
  exit 2
fi

if ! [[ "${EVENTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "events_per_point must be a positive integer."
  exit 2
fi

python3 - "${U_MIN_MM}" "${U_MAX_MM}" "${STEP_MM}" <<'PY'
from decimal import Decimal
import sys
a, b, s = map(Decimal, sys.argv[1:])
if s <= 0 or b <= a:
    raise SystemExit("Require step > 0 and u_max > u_min.")
n = (b - a) / s
if n != n.to_integral_value():
    raise SystemExit("(u_max-u_min)/step must be an integer.")
PY

MOBILITY="${MWPC_ION_MOBILITY:-}"
if [[ -z "${MOBILITY}" ]]; then
  MOBILITY="$(find "${HOME}/na6p/install/garfieldpp"     -type f -name IonMobility_Ar+_Ar.txt -print -quit 2>/dev/null || true)"
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

mapfile -t U_POINTS < <(
  python3 - "${U_MIN_MM}" "${U_MAX_MM}" "${STEP_MM}" <<'PY'
from decimal import Decimal
import sys
a, b, s = map(Decimal, sys.argv[1:])
n = int((b - a) / s)
for i in range(n + 1):
    x = a + i * s
    # Avoid -0.0 and preserve a compact decimal label.
    if x == 0:
        x = Decimal("0")
    print(format(x, "f"))
PY
)

echo "u0 points         : ${#U_POINTS[@]}"
echo "points            : ${U_POINTS[*]}"
echo "total events      : $(( EVENTS * ${#U_POINTS[@]} ))"
echo

# ---------------------------------------------------------------------------
# Numerical map: build once and reuse on resume if the config is unchanged.
# ---------------------------------------------------------------------------
CONFIG_HASH="$(sha256sum "${CONFIG}" | awk '{print $1}')"
MAP_READY=0
if [[ -f "${MAP_DIR}/chamber/drift.result" &&       -f "${MAP_DIR}/chamber/weight_minus_0.result" &&       -f "${MAP_DIR}/chamber/weight_plus_0.result" &&       -f "${CONFIG_HASH_FILE}" ]]; then
  OLD_HASH="$(cat "${CONFIG_HASH_FILE}")"
  if [[ "${OLD_HASH}" == "${CONFIG_HASH}" ]]; then
    MAP_READY=1
  fi
fi

if [[ "${MAP_READY}" -eq 1 ]]; then
  echo "Reusing existing FEM map."
else
  echo "Building FEM map once for the dense scan..."
  rm -rf "${MAP_DIR}"
  mkdir -p "${MAP_DIR}"

  python3 "${NUM_DIR}/generate_weighting_model.py"     --config "${CONFIG}"     --output "${MAP_DIR}"

  cd "${MAP_DIR}"
  gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh
  ElmerGrid 14 2 chamber.msh -autoclean
  ElmerSolver drift.sif
  ElmerSolver weight_minus_0.sif
  ElmerSolver weight_plus_0.sif

  echo "${CONFIG_HASH}" > "${CONFIG_HASH_FILE}"
fi

# Build once.
cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}"   -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4   --target mwpc_stage_c2_weighting_probe mwpc_stage_c2_full_muon_numerical

# Quick map sanity only once per invocation.
"${BUILD_DIR}/mwpc_stage_c2_weighting_probe"   --map-dir "${MAP_DIR}"   --gap-minus-mm "${GAP_MINUS}"   --gap-plus-mm "${GAP_PLUS}"   --profile-half-w-mm 4.0   --profile-step-mm 0.2

point_complete() {
  local csv="$1"
  local expected="$2"
  local u0="$3"

  [[ -f "${csv}" ]] || return 1

  python3 - "${csv}" "${expected}" "${u0}" "${BASE_SEED}" <<'PY'
import csv, math, sys
path, expected, u0, seed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
if len(rows) != expected:
    raise SystemExit(1)
if not rows:
    raise SystemExit(1)
if any(abs(float(r["u0_mm"]) - u0) > 1e-9 for r in rows):
    raise SystemExit(1)
seeds = [int(r["random_seed"]) for r in rows]
if seeds != list(range(seed, seed + expected)):
    raise SystemExit(1)
raise SystemExit(0)
PY
}

rebuild_combined() {
  local first=1
  rm -f "${COMBINED}"

  for U0 in "${U_POINTS[@]}"; do
    local tag="${U0//-/m}"
    tag="${tag//./p}"
    local csv="${EVENT_ROOT}/u_${tag}/c2_u_dense_event_summary.csv"
    if point_complete "${csv}" "${EVENTS}" "${U0}"; then
      if [[ "${first}" -eq 1 ]]; then
        cat "${csv}" > "${COMBINED}"
        first=0
      else
        tail -n +2 "${csv}" >> "${COMBINED}"
      fi
    fi
  done

  if [[ -f "${COMBINED}" ]]; then
    python3 "${STAGE_DIR}/analyze_c2_u_dense_scan.py"       --input "${COMBINED}"       --output "${POINT_SUMMARY}"
  fi
}

trap 'echo; echo "Dense scan interrupted at $(date --iso-8601=seconds). Re-run the same command to resume."; rebuild_combined || true' INT TERM EXIT

completed=0
for U0 in "${U_POINTS[@]}"; do
  tag="${U0//-/m}"
  tag="${tag//./p}"
  EVENT_DIR="${EVENT_ROOT}/u_${tag}"
  CSV="${EVENT_DIR}/c2_u_dense_event_summary.csv"

  if point_complete "${CSV}" "${EVENTS}" "${U0}"; then
    completed=$((completed + 1))
    echo
    echo "[$completed/${#U_POINTS[@]}] u0=${U0} mm already complete -> skip"
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

  stdbuf -oL -eL "${BUILD_DIR}/mwpc_stage_c2_full_muon_numerical"     --numerical-map-dir "${MAP_DIR}"     --map-u-half-mm "${MAP_U_HALF}"     --map-w-half-mm "${MAP_W_HALF}"     --map-wire-pitch-mm "${WIRE_PITCH}"     --map-wire-diam-um "${WIRE_DIAM}"     --map-half-wires "${MAP_HALF_WIRES}"     --gap-minus-mm "${GAP_MINUS}"     --gap-plus-mm "${GAP_PLUS}"     --wire-pitch-mm "${WIRE_PITCH}"     --wire-diam-um "${WIRE_DIAM}"     --hv "${HV}"     --strip-pitch-mm "${STRIP_PITCH}"     --strip-width-mm "${STRIP_WIDTH}"     --tan-alpha "${TAN_ALPHA}"     --half-strips "${HALF_STRIPS}"     --momentum-gev 5     --u0-mm "${U0}"     --w0-mm 0     --theta-u-deg 0     --theta-w-deg 0     --b-tesla 0     --observation-ns 100     --ion-rk-dt-ns 5     --half-wires 6     --events "${EVENTS}"     --max-seeds 0     --base-seed "${BASE_SEED}"     --ion-mobility "${MOBILITY}"     --output-prefix c2_u_dense

  if ! point_complete "${CSV}" "${EVENTS}" "${U0}"; then
    echo "Point u0=${U0} did not produce the expected complete CSV."
    exit 6
  fi

  completed=$((completed + 1))
  echo "completed: $(date --iso-8601=seconds)"
  echo "progress : ${completed}/${#U_POINTS[@]} points"

  # Keep combined/summary files useful even while the scan is still running.
  rebuild_combined
done

rebuild_combined

trap - INT TERM EXIT

echo
echo "================================================================"
echo "C2 dense scan COMPLETE"
echo "finished          : $(date --iso-8601=seconds)"
echo "event summary     : ${COMBINED}"
echo "point summary     : ${POINT_SUMMARY}"
echo "progress log      : ${RUN_LOG}"
echo "================================================================"

#!/usr/bin/env bash
set -euo pipefail

# Resumable Stage-C wire-only scan across the two boundaries of the central
# 4-mm anode-wire cell.
#
# Default scan:
#   u0 = +/- (1.5, 1.6, ..., 2.5) mm
#   100 events / point = 2200 microscopic events
#
# It intentionally does NOT load Elmer maps or transport ions for cathode
# signals.  Only Heed -> microscopic electron avalanche -> anode-wire
# collection is simulated, so this isolates the coordinate associated with the
# wire plane and is substantially lighter than the full C2 cathode response.
#
# Usage:
#   bash test/mwpc-stage-c/run_c_wire_boundary_scan.sh [events] [base_seed]
#
# Environment overrides:
#   C_WIRE_ABS_MIN_MM  default 1.5
#   C_WIRE_ABS_MAX_MM  default 2.5
#   C_WIRE_STEP_MM     default 0.1

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"

EVENTS="${1:-100}"
BASE_SEED="${2:-160001}"

ABS_MIN="${C_WIRE_ABS_MIN_MM:-1.5}"
ABS_MAX="${C_WIRE_ABS_MAX_MM:-2.5}"
STEP="${C_WIRE_STEP_MM:-0.1}"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_wire_boundary_scan"
EVENT_ROOT="${OUT_DIR}/events"
RUN_LOG="${OUT_DIR}/scan_progress.log"
COMBINED_EVENT="${OUT_DIR}/c_wire_event_summary.csv"
COMBINED_WIRE="${OUT_DIR}/c_wire_wire_summary.csv"
POINT_SUMMARY="${OUT_DIR}/c_wire_point_summary.csv"
TRANSITION_SUMMARY="${OUT_DIR}/c_wire_transition_summary.csv"

mkdir -p "${OUT_DIR}" "${EVENT_ROOT}"
exec > >(tee -a "${RUN_LOG}") 2>&1

if ! [[ "${EVENTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "events must be a positive integer."
  exit 2
fi

mapfile -t U_POINTS < <(
python3 - "${ABS_MIN}" "${ABS_MAX}" "${STEP}" <<'PY'
from decimal import Decimal
import sys
a, b, s = map(Decimal, sys.argv[1:])
if a <= 0 or b <= a or s <= 0:
    raise SystemExit("Require 0 < abs_min < abs_max and step > 0.")
n = (b - a) / s
if n != n.to_integral_value():
    raise SystemExit("(abs_max-abs_min)/step must be an integer.")
vals = [a + i*s for i in range(int(n)+1)]
for x in reversed(vals):
    print(format(-x, "f"))
for x in vals:
    print(format(x, "f"))
PY
)

echo
echo "================================================================"
echo "Stage-C wire-only boundary scan"
echo "started           : $(date --iso-8601=seconds)"
echo "points            : ${U_POINTS[*]}"
echo "events / point    : ${EVENTS}"
echo "base seed         : ${BASE_SEED}"
echo "total events      : $(( EVENTS * ${#U_POINTS[@]} ))"
echo "output            : ${OUT_DIR}"
echo "================================================================"

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 --target mwpc_stage_c_wire_response

point_complete() {
  local csv="$1"
  local expected="$2"
  local u0="$3"

  [[ -f "${csv}" ]] || return 1

  python3 - "${csv}" "${expected}" "${u0}" "${BASE_SEED}" <<'PY'
import csv, sys
path, expected, u0, base = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
if len(rows) != expected or not rows:
    raise SystemExit(1)
if any(abs(float(r["u0_mm"]) - u0) > 1e-9 for r in rows):
    raise SystemExit(1)
seeds = [int(r["random_seed"]) for r in rows]
if seeds != list(range(base, base + expected)):
    raise SystemExit(1)
raise SystemExit(0)
PY
}

rebuild_combined() {
  local first_event=1
  local first_wire=1
  rm -f "${COMBINED_EVENT}" "${COMBINED_WIRE}"

  for U0 in "${U_POINTS[@]}"; do
    local tag="${U0//-/m}"
    tag="${tag//./p}"
    local d="${EVENT_ROOT}/u_${tag}"
    local ev="${d}/c_wire_event_summary.csv"
    local wi="${d}/c_wire_wire_summary.csv"

    if point_complete "${ev}" "${EVENTS}" "${U0}"; then
      if [[ "${first_event}" -eq 1 ]]; then
        cat "${ev}" > "${COMBINED_EVENT}"
        first_event=0
      else
        tail -n +2 "${ev}" >> "${COMBINED_EVENT}"
      fi

      if [[ -f "${wi}" ]]; then
        if [[ "${first_wire}" -eq 1 ]]; then
          cat "${wi}" > "${COMBINED_WIRE}"
          first_wire=0
        else
          tail -n +2 "${wi}" >> "${COMBINED_WIRE}"
        fi
      fi
    fi
  done

  if [[ -f "${COMBINED_EVENT}" ]]; then
    python3 "${STAGE_DIR}/analyze_c_wire_coordinate.py" \
      --input "${COMBINED_EVENT}" \
      --output "${POINT_SUMMARY}" \
      --transition-output "${TRANSITION_SUMMARY}" \
      --wire-pitch-mm 4.0 || true
  fi
}

trap 'echo; echo "Wire scan interrupted at $(date --iso-8601=seconds). Re-run the same command to resume."; rebuild_combined || true' INT TERM EXIT

completed=0
for U0 in "${U_POINTS[@]}"; do
  tag="${U0//-/m}"
  tag="${tag//./p}"
  D="${EVENT_ROOT}/u_${tag}"
  CSV="${D}/c_wire_event_summary.csv"

  if point_complete "${CSV}" "${EVENTS}" "${U0}"; then
    completed=$((completed + 1))
    echo "[${completed}/${#U_POINTS[@]}] u0=${U0} mm complete -> skip"
    continue
  fi

  echo
  echo "----------------------------------------------------------------"
  echo "u0 = ${U0} mm"
  echo "started: $(date --iso-8601=seconds)"
  echo "----------------------------------------------------------------"

  rm -rf "${D}"
  mkdir -p "${D}"
  cd "${D}"

  stdbuf -oL -eL "${BUILD_DIR}/mwpc_stage_c_wire_response" \
    --gap-minus-mm 2.5 \
    --gap-plus-mm 2.5 \
    --wire-pitch-mm 4.0 \
    --wire-diam-um 30.0 \
    --hv 1800 \
    --b-tesla 0 \
    --half-wires 6 \
    --momentum-gev 5 \
    --u0-mm "${U0}" \
    --w0-mm 0 \
    --theta-u-deg 0 \
    --theta-w-deg 0 \
    --events "${EVENTS}" \
    --base-seed "${BASE_SEED}" \
    --max-seeds 0 \
    --output-prefix c_wire

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
trap - INT TERM EXIT

echo
echo "================================================================"
echo "Stage-C wire-only boundary scan COMPLETE"
echo "finished          : $(date --iso-8601=seconds)"
echo "event summary     : ${COMBINED_EVENT}"
echo "wire summary      : ${COMBINED_WIRE}"
echo "point summary     : ${POINT_SUMMARY}"
echo "transition summary: ${TRANSITION_SUMMARY}"
echo "================================================================"

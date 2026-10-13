#!/usr/bin/env bash
set -euo pipefail

# Resumable Stage-C single-electron transport response map.
#
# Default baseline:
#   one 4-mm wire cell: u in [-2,+2] mm
#   full 5-mm gas gap:  v in [-2.5,+2.5] mm
#   0.2 x 0.2 mm boxes -> 20 x 25 = 500 boxes
#   50 electrons / box -> 25,000 microscopic drift electrons
#
# The scan runs one u-column at a time, so completed columns survive an
# interruption. No Elmer map and no ion transport are used.
#
# Usage:
#   bash test/mwpc-stage-c/run_c_transport_box_map.sh \
#     [electrons_per_box] [base_seed]
#
# Optional environment variables:
#   C_TRANSPORT_DU_MM       default 0.2
#   C_TRANSPORT_DV_MM       default 0.2
#   C_TRANSPORT_B_TESLA     default 0
#   C_TRANSPORT_START_E_EV  default 0.1

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"

ELECTRONS="${1:-50}"
BASE_SEED="${2:-170001}"

DU_MM="${C_TRANSPORT_DU_MM:-0.2}"
DV_MM="${C_TRANSPORT_DV_MM:-0.2}"
B_TESLA="${C_TRANSPORT_B_TESLA:-0}"
START_E_EV="${C_TRANSPORT_START_E_EV:-0.1}"

U_MIN_MM="-2.0"
U_MAX_MM="+2.0"
V_MIN_MM="-2.5"
V_MAX_MM="+2.5"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_box_map"
COL_DIR="${OUT_DIR}/columns"
RUN_LOG="${OUT_DIR}/scan_progress.log"
COMBINED="${OUT_DIR}/c_transport_box_electrons.csv"
SUMMARY="${OUT_DIR}/c_transport_box_summary.csv"

mkdir -p "${OUT_DIR}" "${COL_DIR}"
exec > >(tee -a "${RUN_LOG}") 2>&1

if ! [[ "${ELECTRONS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "electrons_per_box must be a positive integer."
  exit 2
fi

read -r N_U N_V < <(
python3 - "${U_MIN_MM}" "${U_MAX_MM}" "${V_MIN_MM}" "${V_MAX_MM}" "${DU_MM}" "${DV_MM}" <<'PY'
from decimal import Decimal
import sys
umin, umax, vmin, vmax, du, dv = map(Decimal, sys.argv[1:])
nu = (umax - umin) / du
nv = (vmax - vmin) / dv
if du <= 0 or dv <= 0 or nu != nu.to_integral_value() or nv != nv.to_integral_value():
    raise SystemExit("Grid ranges must be integer multiples of du/dv.")
print(int(nu), int(nv))
PY
)

EXPECTED_PER_COLUMN=$(( N_V * ELECTRONS ))

echo
echo "================================================================"
echo "Stage-C transport box map"
echo "started               : $(date --iso-8601=seconds)"
echo "grid                  : ${N_U} x ${N_V}"
echo "du,dv [mm]            : ${DU_MM}, ${DV_MM}"
echo "electrons / box       : ${ELECTRONS}"
echo "electrons / u-column  : ${EXPECTED_PER_COLUMN}"
echo "total electrons       : $(( N_U * EXPECTED_PER_COLUMN ))"
echo "B [T]                 : ${B_TESLA}"
echo "start energy [eV]     : ${START_E_EV}"
echo "output                 : ${OUT_DIR}"
echo "================================================================"

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 \
  --target mwpc_stage_c_transport_box_map

column_complete() {
  local csv="$1"
  local expected="$2"
  local iu="$3"
  [[ -f "${csv}" ]] || return 1

  python3 - "${csv}" "${expected}" "${iu}" <<'PY'
import csv, sys
path, expected, iu = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with open(path, newline="") as f:
    rows = list(csv.DictReader(f))
if len(rows) != expected:
    raise SystemExit(1)
if any(int(r["iu"]) != iu for r in rows):
    raise SystemExit(1)
raise SystemExit(0)
PY
}

rebuild_combined() {
  local first=1
  rm -f "${COMBINED}"

  for ((iu=0; iu<N_U; ++iu)); do
    local csv="${COL_DIR}/iu_$(printf '%03d' "${iu}").csv"
    if column_complete "${csv}" "${EXPECTED_PER_COLUMN}" "${iu}"; then
      if [[ "${first}" -eq 1 ]]; then
        cat "${csv}" > "${COMBINED}"
        first=0
      else
        tail -n +2 "${csv}" >> "${COMBINED}"
      fi
    fi
  done

  if [[ -f "${COMBINED}" ]]; then
    python3 "${STAGE_DIR}/analyze_c_transport_box_map.py" \
      --input "${COMBINED}" \
      --output "${SUMMARY}"
  fi
}

trap 'echo; echo "Transport-map scan interrupted at $(date --iso-8601=seconds). Re-run the same command to resume."; rebuild_combined || true' INT TERM EXIT

completed=0
for ((iu=0; iu<N_U; ++iu)); do
  CSV="${COL_DIR}/iu_$(printf '%03d' "${iu}").csv"

  if column_complete "${CSV}" "${EXPECTED_PER_COLUMN}" "${iu}"; then
    completed=$((completed + 1))
    echo "[${completed}/${N_U}] u-column ${iu} complete -> skip"
    continue
  fi

  echo
  echo "----------------------------------------------------------------"
  echo "u-column ${iu}/$((N_U - 1))"
  echo "started: $(date --iso-8601=seconds)"
  echo "----------------------------------------------------------------"

  rm -f "${CSV}"

  stdbuf -oL -eL "${BUILD_DIR}/mwpc_stage_c_transport_box_map" \
    --gap-minus-mm 2.5 \
    --gap-plus-mm 2.5 \
    --wire-pitch-mm 4.0 \
    --wire-diam-um 30.0 \
    --hv 1800 \
    --b-tesla "${B_TESLA}" \
    --half-wires 4 \
    --u-min-mm "${U_MIN_MM}" \
    --u-max-mm "${U_MAX_MM}" \
    --v-min-mm "${V_MIN_MM}" \
    --v-max-mm "${V_MAX_MM}" \
    --du-mm "${DU_MM}" \
    --dv-mm "${DV_MM}" \
    --electrons-per-box "${ELECTRONS}" \
    --base-seed "${BASE_SEED}" \
    --start-energy-ev "${START_E_EV}" \
    --w0-mm 0 \
    --w-half-span-mm 20 \
    --iu-min "${iu}" \
    --iu-max "${iu}" \
    --output "${CSV}"

  if ! column_complete "${CSV}" "${EXPECTED_PER_COLUMN}" "${iu}"; then
    echo "u-column ${iu} did not produce the expected complete CSV."
    exit 6
  fi

  completed=$((completed + 1))
  echo "completed: $(date --iso-8601=seconds)"
  echo "progress : ${completed}/${N_U} u-columns"
  rebuild_combined
done

rebuild_combined
trap - INT TERM EXIT

echo
echo "================================================================"
echo "Stage-C transport box map COMPLETE"
echo "finished       : $(date --iso-8601=seconds)"
echo "raw electrons  : ${COMBINED}"
echo "box summary    : ${SUMMARY}"
echo "progress log   : ${RUN_LOG}"
echo "================================================================"

#!/usr/bin/env bash
set -euo pipefail

# Targeted high-statistics Stage-C transport run for the two u-bins nearest
# each wire-cell boundary:
#   iu=0  -> u=-1.9 mm
#   iu=1  -> u=-1.7 mm
#   iu=18 -> u=+1.7 mm
#   iu=19 -> u=+1.9 mm
#
# All 30 v bins are rerun. Default 2000 electrons/box:
#   4 * 30 * 2000 = 240,000 microscopic drifts.
# The run is resumable column-by-column and produces a hybrid full map.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"

ELECTRONS="${1:-2000}"
BASE_SEED="${C_BOUNDARY_BASE_SEED:-210001}"
U_COLUMNS=(0 1 18 19)
N_V=30

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_boundary_highstat_2p4"
COL_DIR="${OUT_DIR}/columns"
COMBINED="${OUT_DIR}/c_transport_boundary_highstat_electrons.csv"
ANALYSIS_PREFIX="${OUT_DIR}/c_transport_boundary_highstat"

BASELINE="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/c_transport_box_electrons.csv"
HYBRID_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_box_map_hybrid_boundary_2p4"
HYBRID="${HYBRID_DIR}/c_transport_box_electrons.csv"

if ! [[ "${ELECTRONS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "electrons_per_box must be a positive integer." >&2
  exit 2
fi

EXPECTED_PER_COLUMN=$((N_V * ELECTRONS))
TOTAL_TARGETED=$((4 * EXPECTED_PER_COLUMN))
mkdir -p "${COL_DIR}" "${HYBRID_DIR}"

echo
echo "================================================================"
echo "Targeted MWPC boundary transport"
echo "u columns             : ${U_COLUMNS[*]}"
echo "u centres [mm]        : -1.9 -1.7 +1.7 +1.9"
echo "v bins / column       : ${N_V}"
echo "electrons / box       : ${ELECTRONS}"
echo "total microscopic e-  : ${TOTAL_TARGETED}"
echo "baseline map          : ${BASELINE}"
echo "================================================================"

if [[ ! -f "${BASELINE}" ]]; then
  echo "Baseline map not found: ${BASELINE}" >&2
  exit 3
fi

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 --target mwpc_stage_c_transport_box_map

column_complete() {
  local csv="$1"
  local iu="$2"
  python3 - "$csv" "${EXPECTED_PER_COLUMN}" "$iu" <<'PY'
import csv
import os
import sys

path, expected, iu = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
if not os.path.isfile(path):
    raise SystemExit(1)
with open(path, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
if len(rows) != expected:
    raise SystemExit(1)
if any(int(float(r["iu"])) != iu for r in rows):
    raise SystemExit(1)
ivs = {int(float(r["iv"])) for r in rows}
if ivs != set(range(30)):
    raise SystemExit(1)
raise SystemExit(0)
PY
}

for iu in "${U_COLUMNS[@]}"; do
  csv="${COL_DIR}/iu_$(printf '%03d' "${iu}").csv"
  if column_complete "${csv}" "${iu}"; then
    echo "iu=${iu}: complete -> skip"
    continue
  fi
  echo
  echo "----------------------------------------------------------------"
  echo "Running high-stat u-column iu=${iu}"
  echo "----------------------------------------------------------------"
  rm -f "${csv}"

  "${BUILD_DIR}/mwpc_stage_c_transport_box_map" \
    --gap-minus-mm 2 \
    --gap-plus-mm 4 \
    --wire-pitch-mm 4 \
    --wire-diam-um 30 \
    --hv 1800 \
    --b-tesla 0 \
    --half-wires 4 \
    --u-min-mm -2 \
    --u-max-mm 2 \
    --v-min-mm -2 \
    --v-max-mm 4 \
    --du-mm 0.2 \
    --dv-mm 0.2 \
    --electrons-per-box "${ELECTRONS}" \
    --base-seed "$((BASE_SEED + iu * 100000))" \
    --start-energy-ev 0.1 \
    --w0-mm 0 \
    --w-half-span-mm 20 \
    --iu-min "${iu}" \
    --iu-max "${iu}" \
    --output "${csv}"

  if ! column_complete "${csv}" "${iu}"; then
    echo "Incomplete output for iu=${iu}: ${csv}" >&2
    exit 4
  fi
done

rm -f "${COMBINED}"
first=1
for iu in "${U_COLUMNS[@]}"; do
  csv="${COL_DIR}/iu_$(printf '%03d' "${iu}").csv"
  if [[ "${first}" -eq 1 ]]; then
    cat "${csv}" > "${COMBINED}"
    first=0
  else
    tail -n +2 "${csv}" >> "${COMBINED}"
  fi
done

python3 "${STAGE_DIR}/analyze_c_transport_timing_shapes.py" \
  --input "${COMBINED}" \
  --output-prefix "${ANALYSIS_PREFIX}"

python3 "${STAGE_DIR}/merge_c_transport_maps.py" \
  --baseline "${BASELINE}" \
  --override "${COMBINED}" \
  --output "${HYBRID}"

echo
echo "================================================================"
echo "Boundary high-stat study COMPLETE"
echo "targeted raw map      : ${COMBINED}"
echo "timing boxes          : ${ANALYSIS_PREFIX}_boxes.csv"
echo "timing box/wires      : ${ANALYSIS_PREFIX}_box_wires.csv"
echo "hybrid full map       : ${HYBRID}"
echo "================================================================"

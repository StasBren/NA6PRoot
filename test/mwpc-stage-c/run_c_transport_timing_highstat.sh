#!/usr/bin/env bash
set -euo pipefail

# Targeted high-statistics timing study for the PI box-response idea.
#
# The full transport map uses 50 electrons/box, which is enough for mean
# response and wire probabilities but too small to judge Gaussianity of the
# drift-time distribution.  This script reruns only four informative u-columns
# with much higher statistics:
#
#   iu=0   -> u=-1.9 mm  (left cell boundary)
#   iu=10  -> u=+0.1 mm  (near central wire)
#   iu=17  -> u=+1.5 mm  (approaching separatrix)
#   iu=19  -> u=+1.9 mm  (right cell boundary)
#
# All 30 v bins of the 2+4 mm chamber are retained.  With the default
# 500 electrons/box this is 4*30*500 = 60,000 microscopic drifts.
#
# Usage:
#   bash test/mwpc-stage-c/run_c_transport_timing_highstat.sh [electrons_per_box]
#
# Default: 500 electrons/box.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"
OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_timing_highstat_2p4"

ELECTRONS="${1:-500}"
BASE_SEED="${C_TIMING_BASE_SEED:-190001}"
U_COLUMNS=(0 10 17 19)

mkdir -p "${OUT_DIR}/columns"

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}"   -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4   --target mwpc_stage_c_transport_box_map

COMBINED="${OUT_DIR}/c_transport_timing_highstat_electrons.csv"
rm -f "${COMBINED}"
first=1

for iu in "${U_COLUMNS[@]}"; do
  out="${OUT_DIR}/columns/iu_$(printf '%03d' "${iu}").csv"
  echo
  echo "=== high-stat timing u-column iu=${iu} ==="

  "${BUILD_DIR}/mwpc_stage_c_transport_box_map"     --gap-minus-mm 2     --gap-plus-mm 4     --wire-pitch-mm 4     --wire-diam-um 30     --hv 1800     --b-tesla 0     --half-wires 4     --u-min-mm -2     --u-max-mm 2     --v-min-mm -2     --v-max-mm 4     --du-mm 0.2     --dv-mm 0.2     --electrons-per-box "${ELECTRONS}"     --base-seed "$((BASE_SEED + iu * 100000))"     --start-energy-ev 0.1     --w0-mm 0     --w-half-span-mm 20     --iu-min "${iu}"     --iu-max "${iu}"     --output "${out}"

  if [[ "${first}" -eq 1 ]]; then
    cat "${out}" > "${COMBINED}"
    first=0
  else
    tail -n +2 "${out}" >> "${COMBINED}"
  fi
done

python3 "${STAGE_DIR}/analyze_c_transport_timing_shapes.py"   --input "${COMBINED}"   --output-prefix "${OUT_DIR}/c_transport_timing_highstat"

echo
echo "Targeted PI timing study COMPLETE"
echo "  raw            : ${COMBINED}"
echo "  boxes          : ${OUT_DIR}/c_transport_timing_highstat_boxes.csv"
echo "  box/wire       : ${OUT_DIR}/c_transport_timing_highstat_box_wires.csv"
echo "  representative : ${OUT_DIR}/c_transport_timing_highstat_representative.csv"

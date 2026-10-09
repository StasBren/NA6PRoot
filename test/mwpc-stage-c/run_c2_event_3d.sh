#!/usr/bin/env bash
set -euo pipefail

# C2 visualization V1:
#   FEM mesh + cathodes/wires + muon track + Heed clusters/primaries +
#   selected microscopic avalanche electron paths + ion-birth cloud +
#   numerical E/weighting-potential slice.
#
# Usage:
#   bash test/mwpc-stage-c/run_c2_event_3d.sh [seed] [detailed_seeds] [config_json]
#
# The script reuses the C2 cell-scan FEM map if available. Otherwise it builds
# one numerical map from the supplied config.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
VIS_DIR="${STAGE_DIR}/visualization"
BUILD_DIR="${STAGE_DIR}/build"

SEED="${1:-120001}"
DETAILED_SEEDS="${2:-2}"
CONFIG="${3:-${NUM_DIR}/weighting_config_c2b.json}"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_event_3d/seed_${SEED}"
EVENT_DIR="${OUT_DIR}/event"
VTK_DIR="${OUT_DIR}/vtk"

REUSE_MAP="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_u_cell_scan/map"
if [[ -f "${REUSE_MAP}/chamber/drift.result" &&       -f "${REUSE_MAP}/chamber/weight_minus_0.result" &&       -f "${REUSE_MAP}/chamber/weight_plus_0.result" ]]; then
  MAP_DIR="${REUSE_MAP}"
  BUILD_MAP=0
else
  MAP_DIR="${OUT_DIR}/map"
  BUILD_MAP=1
fi

missing=0
for cmd in python3 cmake; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required command: ${cmd}"
    missing=1
  fi
done
if [[ "${BUILD_MAP}" -eq 1 ]]; then
  for cmd in gmsh ElmerGrid ElmerSolver; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
      echo "Missing required command: ${cmd}"
      missing=1
    fi
  done
fi
if [[ "${missing}" -ne 0 ]]; then
  exit 2
fi

VIZ_PY="${HOME}/na6p/venvs/mwpc-viz/bin/python"
if [[ ! -x "${VIZ_PY}" ]]; then
  VIZ_PY="$(command -v python3)"
fi
if ! "${VIZ_PY}" -c 'import meshio' >/dev/null 2>&1; then
  echo "Python environment does not contain meshio."
  echo "Expected visualization venv:"
  echo "  source ~/na6p/venvs/mwpc-viz/bin/activate"
  echo "  python -m pip install meshio"
  exit 2
fi

MOBILITY="${MWPC_ION_MOBILITY:-}"
if [[ -z "${MOBILITY}" ]]; then
  MOBILITY="$(find "${HOME}/na6p/install/garfieldpp"     -type f -name IonMobility_Ar+_Ar.txt -print -quit 2>/dev/null || true)"
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
HV="$(read_json electrostatics.anode_voltage_V)"

mkdir -p "${EVENT_DIR}" "${VTK_DIR}"

echo "=== C2 3D event visualization ==="
echo "seed              : ${SEED}"
echo "detailed seeds    : ${DETAILED_SEEDS}"
echo "config            : ${CONFIG}"
echo "FEM map           : ${MAP_DIR}"
echo "reuse FEM map     : $([[ "${BUILD_MAP}" -eq 0 ]] && echo yes || echo no)"
echo "output            : ${OUT_DIR}"
echo

if [[ "${BUILD_MAP}" -eq 1 ]]; then
  rm -rf "${MAP_DIR}"
  mkdir -p "${MAP_DIR}"

  python3 "${NUM_DIR}/generate_weighting_model.py"     --config "${CONFIG}"     --output "${MAP_DIR}"

  cd "${MAP_DIR}"
  gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh
  ElmerGrid 14 2 chamber.msh -autoclean
  ElmerSolver drift.sif
  ElmerSolver weight_minus_0.sif
  ElmerSolver weight_plus_0.sif
fi

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}"   -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 --target mwpc_stage_c2_event_3d_dump

cd "${EVENT_DIR}"

"${BUILD_DIR}/mwpc_stage_c2_event_3d_dump"   --map-dir "${MAP_DIR}"   --gap-minus-mm "${GAP_MINUS}"   --gap-plus-mm "${GAP_PLUS}"   --wire-pitch-mm "${WIRE_PITCH}"   --wire-diam-um "${WIRE_DIAM}"   --hv "${HV}"   --half-wires 6   --momentum-gev 5   --u0-mm 1.0   --w0-mm 0.0   --theta-u-deg 10   --theta-w-deg 25   --b-tesla 0   --seed "${SEED}"   --detailed-seeds "${DETAILED_SEEDS}"   --max-electron-paths-per-seed 160   --max-ion-births-per-seed 5000   --field-u-half-mm 8   --field-nu 161   --field-nv 101   --ion-mobility "${MOBILITY}"   --output-prefix c2_event3d

"${VIZ_PY}" "${VIS_DIR}/export_event_to_vtk.py"   --event-dir "${EVENT_DIR}"   --map-dir "${MAP_DIR}"   --config "${CONFIG}"   --prefix c2_event3d   --output "${VTK_DIR}"

echo
echo "=== 3D scene ready ==="
echo "VTK directory:"
echo "  ${VTK_DIR}"
echo
echo "If ParaView is installed in WSL:"
echo "  paraview --script=\"${VTK_DIR}/open_scene.py\""
echo
echo "Or open ${VTK_DIR}/chamber.vtu manually in ParaView and add the *.vtk layers."

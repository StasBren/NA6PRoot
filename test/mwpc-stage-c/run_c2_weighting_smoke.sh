#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
NUM_DIR="${STAGE_DIR}/numerical-weighting"
BUILD_DIR="${STAGE_DIR}/build"
OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c2_weighting_smoke"
CONFIG="${1:-${NUM_DIR}/weighting_config.json}"

missing=0
for cmd in python3 gmsh ElmerGrid ElmerSolver cmake; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Missing required command: ${cmd}"
    missing=1
  fi
done

if [[ "${missing}" -ne 0 ]]; then
  echo
  echo "C2 numerical weighting requires Gmsh and Elmer in addition to the"
  echo "existing Garfield++ environment. Install/activate those tools, then rerun."
  exit 2
fi

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

echo "=== C2a numerical weighting preflight ==="
echo "config      : ${CONFIG}"
echo "output      : ${OUT_DIR}"
echo "gmsh        : $(command -v gmsh)"
echo "ElmerGrid   : $(command -v ElmerGrid)"
echo "ElmerSolver : $(command -v ElmerSolver)"
echo

python3 "${NUM_DIR}/generate_weighting_model.py"   --config "${CONFIG}"   --output "${OUT_DIR}"

cd "${OUT_DIR}"

echo
echo "=== Gmsh: quadratic 3-D mesh ==="
gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh

echo
echo "=== ElmerGrid: Gmsh -> Elmer mesh ==="
ElmerGrid 14 2 chamber.msh -autoclean

if [[ ! -f chamber/mesh.header ]]; then
  echo "Expected Elmer mesh chamber/mesh.header was not produced."
  exit 3
fi

echo
echo "=== Elmer: physical drift field ==="
ElmerSolver drift.sif

echo
echo "=== Elmer: weighting fields ==="
shopt -s nullglob
weight_sifs=(weight_*.sif)
if [[ "${#weight_sifs[@]}" -eq 0 ]]; then
  echo "No weighting SIF files were generated."
  exit 4
fi
for sif in "${weight_sifs[@]}"; do
  echo "--- ${sif}"
  ElmerSolver "${sif}"
done

if [[ ! -f chamber/drift.result ]]; then
  echo "Expected chamber/drift.result was not produced."
  exit 5
fi
if [[ ! -f chamber/weight_minus_0.result ||       ! -f chamber/weight_plus_0.result ]]; then
  echo "Central-strip weighting results are missing."
  exit 5
fi

echo
echo "=== Garfield++: import and boundary sanity probe ==="
cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}"   -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 --target mwpc_stage_c2_weighting_probe

"${BUILD_DIR}/mwpc_stage_c2_weighting_probe"   --map-dir "${OUT_DIR}"   --gap-minus-mm 2.5   --gap-plus-mm 2.5   --profile-half-w-mm 4.0   --profile-step-mm 0.1

echo
echo "C2a outputs:"
echo "  ${OUT_DIR}/manifest.json"
echo "  ${OUT_DIR}/chamber.msh"
echo "  ${OUT_DIR}/chamber/drift.result"
echo "  ${OUT_DIR}/chamber/weight_minus_0.result"
echo "  ${OUT_DIR}/chamber/weight_plus_0.result"
echo "  ${OUT_DIR}/central_strip_weighting_profile.csv"

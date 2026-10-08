#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"
OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c1_smoke"

mkdir -p "${OUT_DIR}"

MOBILITY="${MWPC_ION_MOBILITY:-}"
if [[ -z "${MOBILITY}" ]]; then
  MOBILITY="$(find "${HOME}/na6p/install/garfieldpp"     -type f -name IonMobility_Ar+_Ar.txt -print -quit 2>/dev/null || true)"
fi

if [[ -z "${MOBILITY}" || ! -f "${MOBILITY}" ]]; then
  echo "Could not find IonMobility_Ar+_Ar.txt."
  echo "Set MWPC_ION_MOBILITY=/full/path/IonMobility_Ar+_Ar.txt"
  exit 2
fi

cmake -S "${STAGE_DIR}" -B "${BUILD_DIR}"   -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4

cd "${OUT_DIR}"

"${BUILD_DIR}/mwpc_stage_c1_full_muon"   --momentum-gev 5   --u0-mm 1.0   --w0-mm 0.0   --theta-u-deg 10   --theta-w-deg 25   --b-tesla 0   --observation-ns 100   --ion-rk-dt-ns 5   --half-wires 6   --half-strips 12   --events 1   --max-seeds 8   --base-seed 120001   --ion-mobility "${MOBILITY}"   --output-prefix c1_smoke

echo
echo "Smoke outputs:"
echo "  ${OUT_DIR}/c1_smoke_event_summary.csv"
echo "  ${OUT_DIR}/c1_smoke_strip_summary.csv"
echo "  ${OUT_DIR}/c1_smoke_seed_summary.csv"
echo "  ${OUT_DIR}/c1_smoke_wire_summary.csv"

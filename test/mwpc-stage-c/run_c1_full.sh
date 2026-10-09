#!/usr/bin/env bash
set -euo pipefail

# Full (all-Heed-seeds) Stage-C1 event.
#
# Usage:
#   bash test/mwpc-stage-c/run_c1_full.sh [momentum_GeV] [theta_u_deg] [theta_w_deg] [seed]
#
# Defaults are chosen near the low-mass downstream core identified in C0.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"
BUILD_DIR="${STAGE_DIR}/build"

P_GEV="${1:-5}"
THETA_U="${2:-10}"
THETA_W="${3:-25}"
SEED="${4:-130001}"

TAG="p${P_GEV}_tu${THETA_U}_tw${THETA_W}_seed${SEED}"
TAG="${TAG//./p}"
TAG="${TAG//-/m}"

OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c1_full/${TAG}"
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

echo "Running full Stage-C1 event:"
echo "  p       = ${P_GEV} GeV/c"
echo "  theta_u = ${THETA_U} deg"
echo "  theta_w = ${THETA_W} deg"
echo "  seed    = ${SEED}"
echo "  output  = ${OUT_DIR}"
echo

"${BUILD_DIR}/mwpc_stage_c1_full_muon"   --gap-minus-mm 2.5   --gap-plus-mm 2.5   --momentum-gev "${P_GEV}"   --u0-mm 1.0   --w0-mm 0.0   --theta-u-deg "${THETA_U}"   --theta-w-deg "${THETA_W}"   --b-tesla 0   --observation-ns 100   --ion-rk-dt-ns 5   --half-wires 6   --half-strips 12   --events 1   --max-seeds 0   --base-seed "${SEED}"   --ion-mobility "${MOBILITY}"   --output-prefix c1_full

echo
echo "Full-event outputs:"
echo "  ${OUT_DIR}/c1_full_event_summary.csv"
echo "  ${OUT_DIR}/c1_full_strip_summary.csv"
echo "  ${OUT_DIR}/c1_full_seed_summary.csv"
echo "  ${OUT_DIR}/c1_full_wire_summary.csv"

#!/usr/bin/env bash
set -euo pipefail

# Build a high-statistics Heed ionisation library for the PI MWPC study.
#
# Geometry:
#   gas gap = 2 mm + 4 mm = 6 mm
#   Ar/CO2 = 70/30
#   wire pitch = 4 mm
#   wire diameter = 30 um
#   HV = 1800 V
#
# This is Heed-only: no drift, avalanche, weighting field or electronics.
# It is therefore cheap compared with Stage-C microscopic transport.
#
# phaseA_heed_statistics starts the muon 0.1 mm inside the +v cathode to avoid
# the exact boundary.  The sampled Heed path is therefore 5.9 mm, even though
# the physical chamber gap is 6.0 mm.  The analysis uses 5.9 mm when extracting
# the cluster density and then the PI generator extrapolates that density over
# the full physical path.
#
# Usage:
#   bash test/mwpc-stage-c/run_c_heed_library_2p4.sh [tracks] [momentum_GeV]
#
# Default:
#   20000 tracks, 10 GeV/c mu-

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
AVAL_DIR="${ROOT_DIR}/test/mwpc-avalanche"
BUILD_DIR="${AVAL_DIR}/build"
OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_heed_library_2p4"

TRACKS="${1:-20000}"
MOMENTUM_GEV="${2:-10}"

mkdir -p "${OUT_DIR}"

cmake -S "${AVAL_DIR}" -B "${BUILD_DIR}" \
  -DGarfield_DIR="${HOME}/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "${BUILD_DIR}" -j4 --target mwpc_phase_a_heed_stats

PREFIX="${OUT_DIR}/heed_2p4"

"${BUILD_DIR}/mwpc_phase_a_heed_stats" \
  --pitch-mm 4 \
  --wire-diam-um 30 \
  --gap-minus-mm 2 \
  --gap-plus-mm 4 \
  --hv 1800 \
  --u0-mm 1 \
  --momentum-gev "${MOMENTUM_GEV}" \
  --theta-deg 0 \
  --phi-deg 0 \
  --tracks "${TRACKS}" \
  --output-prefix "${PREFIX}"

python3 "${ROOT_DIR}/test/mwpc-stage-c/analyze_c_heed_library.py" \
  --tracks "${PREFIX}_tracks.csv" \
  --clusters "${PREFIX}_clusters.csv" \
  --reference-length-mm 5.9 \
  --output "${OUT_DIR}/heed_2p4_summary.txt"

echo
echo "Wrote:"
echo "  ${PREFIX}_tracks.csv"
echo "  ${PREFIX}_clusters.csv"
echo "  ${OUT_DIR}/heed_2p4_summary.txt"

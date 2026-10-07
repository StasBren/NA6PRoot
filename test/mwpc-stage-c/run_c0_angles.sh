#!/usr/bin/env bash
set -euo pipefail

# Canonical Stage-C0 runner.
#
# Usage:
#   bash test/mwpc-stage-c/run_c0_angles.sh [SUITE]
#
# SUITE is the saved high-statistics MNP33-bore production directory.
# No Geant4 transport is run here.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FULL_DEFAULT="${HOME}/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup"
KNOWN_SUITE="${FULL_DEFAULT}/mnp33-bore-runs/20260912-202112-0wrxhdg_"

SUITE="${1:-${MWPC_STAGE_C_BORE_SUITE:-}}"
if [[ -z "${SUITE}" ]]; then
  POINTER="${FULL_DEFAULT}/latest-mnp33-bore-production.txt"
  if [[ -f "${POINTER}" ]]; then
    SUITE="$(cat "${POINTER}")"
  elif [[ -d "${KNOWN_SUITE}" ]]; then
    SUITE="${KNOWN_SUITE}"
  else
    echo "Could not locate the saved high-statistics MNP33-bore suite."
    echo "Pass it explicitly:"
    echo "  bash test/mwpc-stage-c/run_c0_angles.sh /path/to/suite"
    exit 2
  fi
fi

if [[ ! -f "${SUITE}/manifest.json" ]]; then
  echo "Missing manifest.json in suite: ${SUITE}"
  exit 2
fi

OUT="${MWPC_STAGE_C_C0_OUT:-${ROOT_DIR}/test_runs/mwpc_stage_c/c0_angles}"

echo "Stage C0 suite : ${SUITE}"
echo "Stage C0 output: ${OUT}"

python3 "${ROOT_DIR}/test/mwpc-stage-c/analyze_saved_bore_angles.py"   "${SUITE}"   --config "${ROOT_DIR}/test/mwpc-stage-c/stage_c_reference.json"   --output "${OUT}"

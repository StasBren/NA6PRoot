#!/usr/bin/env bash
# NA60+/DiCE: reproducible LUT500 vs midline, from verified frozen WSL inputs.
# Does not reset, overwrite, or re-run any historical Geant4/Garfield campaign.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
case "${1:-}" in
  test)
    "$PYTHON_BIN" -m unittest discover -s "$HERE" -p 'test_*.py' -v
    ;;
  export)
    "$PYTHON_BIN" "$REPO/test/mwpc-stage-c/pi_2m_export_completed_chunks.py" \
      --max-chunks-per-source 40 --root "${ROOT_BIN:-$HOME/na6p/install/root/bin/root}"
    ;;
  select)
    "$PYTHON_BIN" "$HERE/run_comparison.py" select \
      --repo "$REPO" --chunks 40 --n 100000
    ;;
  pilot|full)
    if [[ "$1" == pilot ]]; then n=10000; else n=100000; fi
    lut="$(bash "$HERE/ensure_frozen_inputs.sh" lut)"
    calibration="$(bash "$HERE/ensure_frozen_inputs.sh" calibration)"
    "$PYTHON_BIN" "$HERE/run_comparison.py" compare \
       --repo "$REPO" --selected-n 100000 --limit "$n" \
       --lut "$lut" --calibration "$calibration" --weights donor
    ;;
  *)
    cat <<'USAGE'
Usage: bash test/mwpc-stage-c/lut500-vs-geometry-100k/run.sh COMMAND
  test    -- seven synthetic tests (no physics data)
  export  -- verified/restartable existing 2M Geant4 ROOT-to-CSV export
  select  -- nine hash-ranked 100k crossing samples (3 channels x 3 stations)
  pilot   -- run first 10k in each sample with frozen LUT500 + J/psi donors
  full    -- run all 100k in each sample with frozen LUT500 + J/psi donors
Optional env vars: ROOT_BIN, LUT_ZIP, CALIBRATION_NPZ, PYTHON_BIN.
Requires the pre-existing 2M Geant4 parent chunks and Heed calibration files in WSL.
USAGE
    exit 2
    ;;
esac
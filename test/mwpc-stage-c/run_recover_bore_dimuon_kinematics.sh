#!/usr/bin/env bash
# Rebuild high-statistics J/psi, Phi and Omega parent/daughter kinematics
# from existing MCKine.root files, WITHOUT new Geant4.
#
# Usage:
#   bash test/mwpc-stage-c/run_recover_bore_dimuon_kinematics.sh smoke
#   bash test/mwpc-stage-c/run_recover_bore_dimuon_kinematics.sh full

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
BASE="$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/mnp33-bore-runs/20260912-202112-0wrxhdg_"
OUT="$REPO/test_runs/mwpc_stage_c/recovered_dimuon_kinematics"
MACRO="$HERE/recover_bore_dimuon_kinematics.C"

case "$MODE" in
  smoke) MAX_EVENTS=200; OUT="$OUT""_smoke" ;;
  full) MAX_EVENTS=0 ;;
  *) echo "Usage: bash $0 [smoke|full]" >&2; exit 2 ;;
esac

for channel in jpsi40-bore phi40-bore omega40-bore; do
  for chunk in 001 002; do
    input="$BASE/$channel-chunk$chunk/MCKine.root"
    if [[ ! -s "$input" ]]; then
      echo "MISSING/EMPTY: $input" >&2
      exit 1
    fi
  done
done

if ! command -v root >/dev/null 2>&1; then
  if [[ -f "$HOME/na6p/env.sh" ]]; then
    # Source the existing NA6PRoot environment; do not reinstall anything.
    # shellcheck disable=SC1091
    source "$HOME/na6p/env.sh"
  fi
fi
if ! command -v root >/dev/null 2>&1; then
  echo "ROOT was not found. Source your usual NA6PRoot environment first." >&2
  exit 1
fi

mkdir -p "$OUT"
echo "Reading original ROOT kinematics (mode: $MODE)"
echo "Source: $BASE"
echo "Output: $OUT"
root -l -b -q "$MACRO+(\"$BASE\",\"$OUT\",$MAX_EVENTS)"

echo
echo "Output files:"
ls -lh "$OUT"/kinematics_*.png "$OUT"/kinematics_truth_summary.csv
echo
echo "No new Geant4/Garfield jobs were started."

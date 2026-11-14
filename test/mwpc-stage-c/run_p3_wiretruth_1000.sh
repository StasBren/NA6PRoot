#!/usr/bin/env bash
# Isolated, reproducible 1000-real-Geant4-muon FULL Garfield++ wire-truth study.
# Run from a separate Git worktree; never touches frozen native binaries/results.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
SETUP="$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
if [[ -f "$SETUP" ]]; then source "$SETUP"; fi
if ! command -v cmake >/dev/null || ! command -v python3 >/dev/null; then
  echo 'ERROR: cmake and python3 are required.' >&2;exit 2
fi
SCRIPT="$HERE/p3_wiretruth_1000_v1.py"
BUILD="$REPO/test/mwpc-avalanche/build"
TARGET='mwpc_p3_noFW_penning057_wiretruth_v1'
cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$BUILD" --target "$TARGET" -j2
python3 "$SCRIPT" selftest
python3 "$SCRIPT" prepare
python3 "$SCRIPT" status

# Native first-chunk gate. Never carry on if C++ compiling or the exact
# per-anode charge/ion conservation checks fail.
FIRST="$REPO/test_runs/mwpc_stage_c/p3_noFW_penning057_wiretruth_1000_v1/inputs/chunk_0000.csv"
"$BUILD/$TARGET" --validate-input "$FIRST"
"$BUILD/$TARGET" --validate-heed "$FIRST"
python3 -u "$SCRIPT" run --max-chunks 1
python3 "$SCRIPT" audit

if [[ "${NA60_WIRETRUTH_PILOT_ONLY:-0}" == '1' ]]; then
  echo 'PILOT-ONLY: 20 native muons verified. Remove NA60_WIRETRUTH_PILOT_ONLY=1 for full 1000.'
  exit 0
fi
python3 -u "$SCRIPT" run --max-chunks 50
python3 "$SCRIPT" audit
python3 "$SCRIPT" analyse
python3 "$SCRIPT" export
python3 "$SCRIPT" status
printf '\nRESULTS: %s\n' "$REPO/test_runs/mwpc_stage_c/p3_noFW_penning057_wiretruth_1000_v1"

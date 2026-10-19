#!/usr/bin/env bash
# Read existing complete full-campaign native QA outputs; make a shareable,
# conservation-checked text report. Does not invoke ROOT or Geant4.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile "$HERE/summarize_pi_native_full_qa.py"
python3 "$HERE/summarize_pi_native_full_qa.py"
echo
echo "Send this small text file to ChatGPT for PI physics interpretation:"
echo "  $REPO/test_runs/mwpc_stage_c/pi_native_full_quality/pi_native_QA_readout.txt"

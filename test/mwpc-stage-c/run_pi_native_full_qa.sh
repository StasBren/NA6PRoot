#!/usr/bin/env bash
# PI actual-Geant4 full 298-chamber validation, no cathode strips/electronics,
# NO new Geant4 and NO new microscopic Garfield.
#
# Start with one validated native chunk per channel:
#   bash test/mwpc-stage-c/run_pi_native_full_qa.sh first
#
# Then stream the full frozen 201-chunk campaign, resuming finished chunks:
#   bash test/mwpc-stage-c/run_pi_native_full_qa.sh all
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
source "$HERE/activate_pi_stagec_env.sh"
python3 -m py_compile \
 "$HERE/pi_native_full_campaign_qa.py" \
 "$HERE/test_pi_native_full_campaign_qa.py"
python3 "$HERE/test_pi_native_full_campaign_qa.py"
case "${1:-first}" in
 first)
   if [[ $# -gt 1 ]]; then echo "Usage: $0 first|all" >&2; exit 2; fi
   python3 "$HERE/pi_native_full_campaign_qa.py" --first-chunks 1
   ;;
 all)
   if [[ $# -gt 1 ]]; then echo "Usage: $0 first|all" >&2; exit 2; fi
   python3 "$HERE/pi_native_full_campaign_qa.py"
   ;;
 *)
   echo "Usage: $0 first|all" >&2
   exit 2
   ;;
esac

#!/usr/bin/env bash
# Analyze existing already-merged 100k-parent/source MWPC/Heed/Garfield data.
# No Geant4, Garfield++, Heed, ROOT processes or physics refolding.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
  "$HERE/analyze_pi_100k_anode_timing.py" \
  "$HERE/test_pi_100k_anode_timing.py"
echo "TEST: read-only station/anode moment+time-bin conservation (synthetic)"
python3 "$HERE/test_pi_100k_anode_timing.py"
echo "ANALYZE: 18 already-merged real source/station JSON files;"
echo "no new physics simulation and no changes to existing V3 HTML."
python3 "$HERE/analyze_pi_100k_anode_timing.py"
echo "ANALYSIS COMPLETE. To open the independent English anode-time report:"
echo '  explorer.exe "$(wslpath -w "$HOME/na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/anode_timing_report_en.html")"'

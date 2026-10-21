#!/usr/bin/env bash
# READ-ONLY analysis of six completed full-pitch Garfield field-wire pilots.
# No ROOT/Geant4/Garfield process or new physics simulation is launched.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
  "$HERE/pi_pitch_fieldwire_paired_analysis.py" \
  "$HERE/pi_pitch_fieldwire_paired_html.py" \
  "$HERE/test_pi_pitch_fieldwire_paired_analysis.py"
python3 "$HERE/test_pi_pitch_fieldwire_paired_analysis.py"
python3 "$HERE/pi_pitch_fieldwire_paired_analysis.py" "$@"
echo "PASS: paired Garfield pitch/field-wire comparison and offline dashboard."

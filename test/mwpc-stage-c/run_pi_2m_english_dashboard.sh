#!/usr/bin/env bash
# Read-only HTML localisation AFTER native 100k/channel merge has passed.
# NO GEANT4 / HEED / GARFIELD / ROOT processes, NO new physical sampling.
# Creates English *_en.html siblings of all original 18 XY V3 maps,
# all 18 detailed per-anode spectra pages, and the new English portal.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE=build
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 || ! "$MODE" =~ ^(build|check)$ ]]; then
  echo "Usage: bash $0 build|check" >&2
  exit 2
fi
python3 -m py_compile \
  "$HERE/build_pi_100k_english_dashboard.py" \
  "$HERE/test_pi_100k_english_dashboard.py"
python3 "$HERE/test_pi_100k_english_dashboard.py"
if [[ "$MODE" == check ]]; then
  python3 "$HERE/build_pi_100k_english_dashboard.py" --dry-run
else
  python3 "$HERE/build_pi_100k_english_dashboard.py"
fi

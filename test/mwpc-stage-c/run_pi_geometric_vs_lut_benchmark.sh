#!/usr/bin/env bash
# Matched geometric Poisson-Heed binomial vs six frozen Garfield LUT.
# Does NOT run Geant4 or microscopic Garfield physics, so no C++ rebuild.
# Common Poisson clusters and empirical Heed marks are reused exactly.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -ne 0 || ! "$MODE" =~ ^(smoke|benchmark)$ ]]; then
  echo "Usage: bash $0 smoke|benchmark" >&2
  exit 2
fi
python3 -m py_compile \
 "$HERE/pi_geometric_vs_lut_benchmark.py" \
 "$HERE/test_pi_geometric_vs_lut_benchmark.py"
python3 "$HERE/test_pi_geometric_vs_lut_benchmark.py"
if [[ "$MODE" == "smoke" ]]; then
  N=200
else
  N=2000
fi
python3 "$HERE/pi_geometric_vs_lut_benchmark.py" \
  --channel Jpsi --chunk 1 --station 0 --sample-muon-hits "$N"

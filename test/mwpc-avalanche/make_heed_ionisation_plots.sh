#!/usr/bin/env bash
set -euo pipefail

EXE="./build/mwpc_phase_a_heed_stats"
PLOTTER="plot_heed_ionisation_statistics.py"

if [[ ! -x "$EXE" ]]; then
  echo "Missing $EXE. Build test/mwpc-avalanche first."
  exit 1
fi

PREFIX="heed_ionisation_5000"

"$EXE" \
  --pitch-mm 4 \
  --gap-minus-mm 2 \
  --gap-plus-mm 4 \
  --wire-diam-um 30 \
  --hv 1800 \
  --u0-mm 1.0 \
  --momentum-gev 10 \
  --theta-deg 0 \
  --phi-deg 0 \
  --tracks 5000 \
  --output-prefix "$PREFIX"

python3 "$PLOTTER" \
  --tracks-csv "${PREFIX}_tracks.csv" \
  --clusters-csv "${PREFIX}_clusters.csv" \
  --output-prefix "$PREFIX"

echo
echo "Finished Heed ionisation-statistics ensemble."

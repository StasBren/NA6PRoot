#!/usr/bin/env bash
set -euo pipefail

# Stage-B ensemble study of one-seed intrinsic position reconstruction.
#
# Default: 50 independent microscopic avalanches at each of 9 true positions.
# Use e.g.
#   EVENTS=10 bash run_stageB_resolution_ensemble.sh
# for a quick pilot, or
#   EVENTS=100 bash run_stageB_resolution_ensemble.sh
# for higher statistics.
#
# Zero/small-gain events are NOT resampled. They are retained as
# unreconstructable events so the ensemble is not gain-conditioned.

EXE="${EXE:-./build/mwpc_phase_b_resolution_ensemble}"
EVENTS="${EVENTS:-50}"
BASE_SEED="${BASE_SEED:-30000}"
ION_DT_NS="${ION_DT_NS:-1.0}"

"${EXE}" \
  --gap-minus-mm 2 \
  --gap-plus-mm 4 \
  --wire-pitch-mm 4 \
  --wire-diam-um 30 \
  --hv 1800 \
  --b-tesla 0 \
  --side minus \
  --strip-pitch-mm 1.7 \
  --strip-width-mm 1.7 \
  --half-strips 5 \
  --seed-u-mm 0 \
  --seed-distance-mm 1.0 \
  --seed-energy-ev 0.1 \
  --positions-mm=-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8 \
  --windows-ns=25,50,100,200 \
  --reco-strips=-1,0,1 \
  --events-per-position "${EVENTS}" \
  --base-seed "${BASE_SEED}" \
  --avalanche-limit 50000 \
  --min-ions-for-reco 1 \
  --ion-rk-dt-ns "${ION_DT_NS}" \
  --output stageB_resolution_ensemble.csv

python3 analyze_stageB_resolution_ensemble.py \
  --input stageB_resolution_ensemble.csv \
  --representative-window-ns 100 \
  --output-prefix stageB_resolution_ensemble

#!/usr/bin/env bash
set -euo pipefail

# Phase-B1b controlled sigma_w scan.
# Default: +4 mm cathode side, 200 ions/run, converged B1a numerics.

EXE="${EXE:-./build/mwpc_phase_b_ion_cloud}"

common=(
  --gap-minus-mm 2
  --gap-plus-mm 4
  --wire-pitch-mm 4
  --wire-diam-um 30
  --hv 1800
  --side plus
  --strip-pitch-mm 1.7
  --strip-width-mm 1.7
  --half-strips 5
  --ions 200
  --ion-weight 1
  --sigma-u-um 5
  --sigma-v-um 5
  --dt-ns 5
  --tmax-us 300
  --max-step-mm 0.05
  --seed 12345
)

"${EXE}" "${common[@]}" --sigma-w-um   50 --output-prefix phaseB_sigmaw_005
"${EXE}" "${common[@]}" --sigma-w-um  200 --output-prefix phaseB_sigmaw_020
"${EXE}" "${common[@]}" --sigma-w-um  500 --output-prefix phaseB_sigmaw_050
"${EXE}" "${common[@]}" --sigma-w-um 1000 --output-prefix phaseB_sigmaw_100

python3 plot_phaseB_sigmaw_scan.py \
  --runs \
  "0.05=phaseB_sigmaw_005" \
  "0.20=phaseB_sigmaw_020" \
  "0.50=phaseB_sigmaw_050" \
  "1.00=phaseB_sigmaw_100" \
  --output-prefix phaseB_sigmaw_scan

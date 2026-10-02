#!/usr/bin/env bash
set -euo pipefail

# Stage-B controlled synthetic avalanche-ion cloud scan.
# Same -2 mm cathode side as the single-pair study.
#
# Synthetic cloud:
#   u ~ N(0, sigma_u)
#   v ~ N(v0, sigma_v)
#   w ~ N(0, sigma_w)
#
# We scan sigma_w because the ideal strips are infinite in u and segmented in w,
# so sigma_w is the cloud width most directly visible in strip sharing.

EXE="${EXE:-./build/mwpc_phase_b_ion_cloud}"

common=(
  --gap-minus-mm 2
  --gap-plus-mm 4
  --wire-pitch-mm 4
  --wire-diam-um 30
  --hv 1800
  --side minus
  --strip-pitch-mm 1.7
  --strip-width-mm 1.7
  --half-strips 5
  --ions 500
  --ion-weight 1
  --sigma-u-um 5
  --sigma-v-um 5
  --dt-ns 5
  --tmax-us 30
  --max-step-mm 0.05
  --seed 12345
)

"${EXE}" "${common[@]}" --sigma-w-um   50 --output-prefix stageB_cloud_sw005
"${EXE}" "${common[@]}" --sigma-w-um  200 --output-prefix stageB_cloud_sw020
"${EXE}" "${common[@]}" --sigma-w-um  500 --output-prefix stageB_cloud_sw050
"${EXE}" "${common[@]}" --sigma-w-um 1000 --output-prefix stageB_cloud_sw100

python3 plot_phaseB_sigmaw_scan.py \
  --runs \
  "0.05=stageB_cloud_sw005" \
  "0.20=stageB_cloud_sw020" \
  "0.50=stageB_cloud_sw050" \
  "1.00=stageB_cloud_sw100" \
  --output-prefix stageB_cloud_width_scan

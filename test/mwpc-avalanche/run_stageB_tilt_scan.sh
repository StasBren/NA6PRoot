#!/usr/bin/env bash
set -euo pipefail

# Stage B3a: one ideal tilted strip family.
#
# Fixed physical chamber:
#   Ar/CO2 70:30, 30 um wires, 4 mm wire pitch,
#   2+4 mm cathode gaps, 1800 V, B=0.
#
# Fixed readout sandbox point:
#   p_strip = 1.7 mm, s_strip = 1.7 mm,
#   readout on the 2-mm cathode, T_obs = 100 ns.
#
# Scan:
#   tan(alpha) = 0, 0.05, 0.10
#   anode-wire index m = -2 .. +2
#   along-wire position w0 across one strip pitch.
#
# The SAME central-wire microscopic avalanche ensemble is translated by integer
# wire pitches to isolate the ideal periodic readout geometry.
#
# Key test:
#   x_alpha = w0 cos(alpha) - u_wire sin(alpha)
# Curves from different wires should collapse when plotted against x_alpha.

EXE="${EXE:-./build/mwpc_phase_b_tilt_scan}"
EVENTS="${EVENTS:-30}"

"${EXE}" \
  --gap-minus-mm 2 \
  --gap-plus-mm 4 \
  --wire-pitch-mm 4 \
  --wire-diam-um 30 \
  --hv 1800 \
  --b-tesla 0 \
  --side minus \
  --seed-distance-mm 1.0 \
  --seed-energy-ev 0.1 \
  --strip-pitch-mm 1.7 \
  --strip-width-mm 1.7 \
  --tan-alphas=0,0.05,0.10 \
  --wire-indices=-2,-1,0,1,2 \
  --xi-w-values=-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5 \
  --half-strips 7 \
  --events "${EVENTS}" \
  --base-seed 52000 \
  --avalanche-limit 50000 \
  --observation-ns 100 \
  --ion-rk-dt-ns 1 \
  --output stageB_tilt_scan.csv

python3 analyze_stageB_tilt_scan.py \
  --input stageB_tilt_scan.csv \
  --target-tan-alpha 0.10 \
  --output-prefix stageB_tilt_scan

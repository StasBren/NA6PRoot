#!/usr/bin/env bash
set -euo pipefail

# Stage B2: first readout-geometry study.
#
# Hold the physical chamber and observation time fixed:
#   Ar/CO2 70:30, 2+4 mm gaps, 30 um wires, 4 mm wire pitch,
#   1800 V, B=0, readout on the 2-mm side, alpha=0, T=100 ns.
#
# Vary only strip pitch p and strip width s. The same physical avalanche
# ensemble is reused for every geometry and every normalized sub-strip position.
#
# This is deliberately a SIGNAL-FORMATION scan, not a final resolution scan.

EXE="${EXE:-./build/mwpc_phase_b_strip_geometry_scan}"
EVENTS="${EVENTS:-30}"

"${EXE}" \
  --gap-minus-mm 2 \
  --gap-plus-mm 4 \
  --wire-pitch-mm 4 \
  --wire-diam-um 30 \
  --hv 1800 \
  --b-tesla 0 \
  --side minus \
  --seed-u-mm 0 \
  --seed-distance-mm 1.0 \
  --seed-energy-ev 0.1 \
  --strip-pitches-mm=1.5,1.7,1.9,2.1,2.3 \
  --strip-widths-mm=1.0,1.3,1.5,1.7,2.0 \
  --xi-values=-0.45,-0.30,-0.15,0,0.15,0.30,0.45 \
  --half-strips 6 \
  --events "${EVENTS}" \
  --base-seed 41000 \
  --avalanche-limit 50000 \
  --observation-ns 100 \
  --ion-rk-dt-ns 1 \
  --output stageB_strip_geometry_scan.csv

python3 analyze_stageB_strip_geometry_scan.py \
  --input stageB_strip_geometry_scan.csv \
  --baseline-pitch-mm 1.7 \
  --baseline-width-mm 1.7 \
  --output-prefix stageB_strip_geometry_scan

#!/usr/bin/env bash
set -euo pipefail

# Controlled Stage-B position scan across one 1.7 mm strip pitch.
#
# Key idea:
#   use the SAME Garfield random seed at every w0.
# In the ideal geometry the physical drift field is translationally invariant
# along w, so this reproduces (to numerical precision) the same microscopic
# avalanche translated relative to the segmented readout strips.
#
# This isolates the readout-position response from avalanche-to-avalanche
# fluctuations. It is NOT yet a resolution measurement.

EXE="${EXE:-./build/mwpc_phase_b_real_avalanche}"
SEED="${SEED:-24680}"

positions=(-0.80 -0.60 -0.40 -0.20 0.00 0.20 0.40 0.60 0.80)

for w in "${positions[@]}"; do
  tag=$(printf "%+.2f" "${w}" | sed 's/+/_p/; s/-/_m/; s/\./p/')
  prefix="stageB_posscan${tag}"

  echo
  echo "=== w0 = ${w} mm -> ${prefix} ==="

  "${EXE}" \
    --gap-minus-mm 2 \
    --gap-plus-mm 4 \
    --wire-pitch-mm 4 \
    --wire-diam-um 30 \
    --hv 1800 \
    --side minus \
    --strip-pitch-mm 1.7 \
    --strip-width-mm 1.7 \
    --half-strips 5 \
    --seed-u-mm 0 \
    --seed-w-mm "${w}" \
    --seed-distance-mm 1.0 \
    --seed-energy-ev 0.1 \
    --avalanche-limit 50000 \
    --min-avalanche-ions 500 \
    --max-avalanche-attempts 20 \
    --random-seed "${SEED}" \
    --electron-dt-ns 0.02 \
    --electron-tmax-ns 200 \
    --ion-dt-ns 5 \
    --ion-tmax-us 50 \
    --max-ion-step-mm 0.05 \
    --b-tesla 0 \
    --output-prefix "${prefix}"
done

python3 analyze_stageB_position_scan.py \
  --prefix-glob "stageB_posscan*_event_summary.csv" \
  --windows-ns 25,50,100,200 \
  --cluster-strips=-1,0,1 \
  --output-prefix stageB_position_scan

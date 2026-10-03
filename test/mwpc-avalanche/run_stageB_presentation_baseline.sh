#!/usr/bin/env bash
set -euo pipefail

# Reproducible Stage-B signal-formation baseline for presentation.
#
# The random seed is selected from the existing ensemble study at w0=0:
# among successful non-zero avalanches, choose the event closest to the
# median positive gain.  Then rerun exactly that one avalanche with the full
# microscopic electron signal + explicit DriftLineRKF positive-ion transport.
#
# This is an illustrative median-scale avalanche, not an ensemble average.

ENSEMBLE_CSV="${ENSEMBLE_CSV:-stageB_resolution_ensemble.csv}"
PREFIX="${PREFIX:-stageB_presentation_baseline}"

if [[ ! -f "${ENSEMBLE_CSV}" ]]; then
  echo "Missing ${ENSEMBLE_CSV}."
  echo "Finish the Stage-B ensemble run first, or set ENSEMBLE_CSV=..."
  exit 2
fi

SEED=$(python3 select_stageB_presentation_seed.py \
  --input "${ENSEMBLE_CSV}" \
  --true-w-mm 0 \
  --seed-only)

echo
python3 select_stageB_presentation_seed.py \
  --input "${ENSEMBLE_CSV}" \
  --true-w-mm 0
echo
echo "Re-running presentation baseline with random seed ${SEED}"

./build/mwpc_phase_b_real_avalanche \
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
  --seed-w-mm 0 \
  --seed-distance-mm 1.0 \
  --seed-energy-ev 0.1 \
  --avalanche-limit 50000 \
  --min-avalanche-ions 0 \
  --max-avalanche-attempts 1 \
  --random-seed "${SEED}" \
  --electron-dt-ns 0.02 \
  --electron-tmax-ns 200 \
  --ion-dt-ns 5 \
  --ion-tmax-us 50 \
  --max-ion-step-mm 0.05 \
  --b-tesla 0 \
  --output-prefix "${PREFIX}"

python3 analyze_stageB_real_avalanche.py \
  --prefix "${PREFIX}" \
  --output-prefix "${PREFIX}_analysis"

echo
echo "Presentation baseline outputs:"
echo "  ${PREFIX}_event_summary.csv"
echo "  ${PREFIX}_analysis_current_components.png"
echo "  ${PREFIX}_analysis_signal_composition_vs_time.png"
echo "  ${PREFIX}_analysis_early_spatial_sharing.png"

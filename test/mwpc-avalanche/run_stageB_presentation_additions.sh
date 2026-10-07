#!/usr/bin/env bash
set -euo pipefail

# Generate the two synthesis plots requested for the revised Stage-B story.
#
# Prerequisites:
#   1) run_stageB_position_scan.sh has produced
#      stageB_position_scan_summary.csv and stageB_posscan_p0p40_* waveforms.
#
# No new Garfield simulation is required here; these plots reuse existing
# Stage-B outputs.

POSITION_SUMMARY="${POSITION_SUMMARY:-stageB_position_scan_summary.csv}"
KERNEL_PREFIX="${KERNEL_PREFIX:-stageB_posscan_p0p40}"

if [[ ! -f "${POSITION_SUMMARY}" ]]; then
  echo "Missing ${POSITION_SUMMARY}"
  echo "Run ./run_stageB_position_scan.sh first."
  exit 2
fi

for suffix in _electron_waveforms.csv _ion_waveforms.csv _strip_summary.csv _event_summary.csv; do
  if [[ ! -f "${KERNEL_PREFIX}${suffix}" ]]; then
    echo "Missing ${KERNEL_PREFIX}${suffix}"
    echo "Run ./run_stageB_position_scan.sh first, or set KERNEL_PREFIX=..."
    exit 2
  fi
done

python3 plot_stageB_cog_time_stability.py   --input "${POSITION_SUMMARY}"   --reference-window-ns 200   --output-prefix stageB_cog_time_stability

python3 plot_stageB_response_kernel.py   --prefix "${KERNEL_PREFIX}"   --cluster-strips=-2,-1,0,1,2   --max-window-ns 5000   --output stageB_response_kernel.png

echo
echo "Wrote:"
echo "  stageB_cog_time_stability.png"
echo "  stageB_cog_time_stability_metrics.csv"
echo "  stageB_response_kernel.png"


# Optional microscopic diagnostic for the multiple peaks in the prompt
# electron-current pulse.  Reuses the already exported electron endpoints.
if [[ -f "${KERNEL_PREFIX}_electron_endpoints.csv" ]]; then
  python3 plot_stageB_avalanche_burst_structure.py \
    --prefix "${KERNEL_PREFIX}" \
    --bin-ns 0.02 \
    --output stageB_avalanche_burst_structure.png
else
  echo "Skipping avalanche-burst diagnostic: missing ${KERNEL_PREFIX}_electron_endpoints.csv"
fi

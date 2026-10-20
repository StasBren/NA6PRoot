#!/usr/bin/env bash
# Read-only physics preflight and CPU-side LUT microbenchmark.
# DOES NOT start Geant4, Heed or a new Garfield transport scan.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

python3 -m py_compile \
 "$HERE/pi_quadrant_symmetry_audit.py" \
 "$HERE/test_pi_quadrant_symmetry_audit.py" \
 "$HERE/pi_stagec_performance_report.py" \
 "$HERE/test_pi_stagec_performance_report.py" \
 "$HERE/pi_centered_lut_quality_benchmark.py" \
 "$HERE/test_pi_centered_lut_quality_benchmark.py" \
 "$HERE/run_pi_geant4_campaign.py"

python3 "$HERE/test_pi_quadrant_symmetry_audit.py"
python3 "$HERE/test_pi_stagec_performance_report.py"
python3 "$HERE/test_pi_centered_lut_quality_benchmark.py"

QA="$REPO/test_runs/mwpc_stage_c/pi_native_full_quality"
if [[ -f "$QA/pi_native_all_QA_report.json" ]]; then
  python3 "$HERE/pi_quadrant_symmetry_audit.py"
else
  python3 "$HERE/pi_quadrant_symmetry_audit.py" --geometry-only
  echo "NOTE: full Geant4 QA missing, quadrant count/angle comparison not run."
fi

SUITE="$REPO/test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1"
if [[ -f "$SUITE/campaign.json" ]]; then
  python3 "$HERE/pi_stagec_performance_report.py" --suite "$SUITE"
else
  echo "NOTE: original Geant4 campaign manifest not present."
fi

LUT="$REPO/test_runs/mwpc_stage_c/c_transport_box_map_gm3p0_gp3p0"
if [[ -f "$LUT/c_transport_box_electrons.csv" && \
      -f "$LUT/pi_centered_transport_manifest.json" ]]; then
  python3 "$HERE/pi_centered_lut_quality_benchmark.py"
else
  echo "NOTE: centered p=4 mm Garfield LUT/manifest missing; "
  echo "      no lookup-precision or timing claims until it is available."
fi
echo "PREFLIGHT COMPLETE. Geant4 quadrant-only generation remains DISABLED."

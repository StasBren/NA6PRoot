#!/usr/bin/env bash
# Refresh BOTH offline MWPC dashboards and Geant4 event viewer HTML from
# EXISTING local CSV libraries. No native ROOT export, no Geant4 transport,
# no microscopic Garfield electron simulation, no cathode strip model.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
OUTPUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps"
NATIVE="$REPO/test_runs/mwpc_stage_c/pi_native_anode_first_chunk/crossings"
for source in jpsi40-bore phi40-bore omega40-bore; do
  file="$NATIVE/crossings_${source}_ms2_530.csv"
  if [[ ! -s "$file" ]]; then
    echo "Missing existing native muon gas-hit CSV: $file" >&2
    echo "First run the existing validated native exporter:" >&2
    echo " bash test/mwpc-stage-c/run_pi_native_anode_first_chunk.sh" >&2
    exit 2
  fi
done
echo "Refreshing station statistics dashboard..."
bash "$HERE/run_pi_station_dashboard.sh"
echo "Refreshing native-Geant4 Muon Track Explorer from preserved gas CSVs..."
python3 "$HERE/build_pi_muon_track_explorer.py" \
  --crossings "$NATIVE" --data-kind native \
  --sample-modulo 1 --max-events-per-source 32 \
  --out "$OUTPUT/pi_muon_track_explorer_native.html"
echo "PASS: both HTMLs updated, no new ROOT export or Geant4"
echo "  $OUTPUT/pi_station_interactive.html"
echo "  $OUTPUT/pi_muon_track_explorer_native.html"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows browser path (native Geant4 viewer):"
  echo "  $(wslpath -w "$OUTPUT/pi_muon_track_explorer_native.html")"
fi

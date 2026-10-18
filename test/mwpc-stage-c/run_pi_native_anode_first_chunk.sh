#!/usr/bin/env bash
# Sample actual current-layout MWPC Geant4 gas hits (not projected history!)
# from first completed chunk of each source and open anode-only explorer.
#
# One-time native ROOT export links EACH saved hit to the corresponding
# remapped MCKine TParticle PDG in the SAME merged chunk.
# No new Geant4 transport, Garfield microscopic simulations or detector
# installation changes. The frozen 200k/channel results are read only.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
source "$HERE/activate_pi_stagec_env.sh"

SUITE="$REPO/test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1"
STUDY="$REPO/test_runs/mwpc_stage_c/pi_native_anode_first_chunk"
EXPORT="$STUDY/native_exports"
CROSSINGS="$STUDY/crossings"
mkdir -p "$EXPORT" "$CROSSINGS"
python3 -m py_compile \
 "$HERE/prepare_pi_native_hit_explorer.py" \
 "$HERE/build_pi_muon_track_explorer.py" \
 "$HERE/pi_anode_digitizer.py" \
 "$HERE/test_pi_anode_digitizer.py"
python3 "$HERE/test_pi_anode_digitizer.py"
python3 "$HERE/test_pi_native_anode_pipeline.py"
for channel in jpsi phi omega; do
  CHUNK="$SUITE/$channel/chunk0001"
  if [[ ! -f "$CHUNK/completed.json" ]]; then
    echo "ERROR: no validated $channel first full-campaign chunk: $CHUNK" >&2
    exit 2
  fi
  for file in HitsMuonSpec.root MCKine.root geometry.root; do
    if [[ ! -s "$CHUNK/$file" ]]; then
      echo "ERROR: missing native real-hit input $CHUNK/$file" >&2
      exit 2
    fi
  done
  CSV="$EXPORT/${channel}_geant_gas_crossings.csv"
  LOG="$EXPORT/${channel}_native_export.log"
  echo "Exporting $channel true Geant4 MWPC muon hits from $(basename "$CHUNK") ..."
  MACRO="$HERE/export_c_geant_muon_crossings.C+(\"$CHUNK/HitsMuonSpec.root\",\"$CHUNK/geometry.root\",\"$CSV\",\"$CHUNK/MCKine.root\")"
  root -l -b -q "$MACRO" > "$LOG" 2>&1 || {
    tail -n 55 "$LOG" >&2
    exit 3
  }
  if [[ ! -s "$CSV" ]] || ! grep -q 'GEANT MWPC crossing export COMPLETE' "$LOG"; then
    echo "ERROR: native ROOT exporter did not complete for $channel" >&2
    tail -n 55 "$LOG" >&2
    exit 4
  fi
  tail -n 6 "$LOG"
done

python3 "$HERE/prepare_pi_native_hit_explorer.py" \
 --export-dir "$EXPORT" --out "$CROSSINGS" --chunk chunk0001

# Use each actual Geant4 gas crossing and PDG, then resample ionisation
# with existing Heed and Garfield single-electron lookups.
OUTPUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_muon_track_explorer_native.html"
python3 "$HERE/build_pi_muon_track_explorer.py" \
 --crossings "$CROSSINGS" --data-kind native \
 --sample-modulo 1 --max-events-per-source 32 --out "$OUTPUT"

echo "New Geant4 true-gas-hit Muon Track Explorer:"
echo "  $OUTPUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows path: $(wslpath -w "$OUTPUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUTPUT")" || true
  fi
fi

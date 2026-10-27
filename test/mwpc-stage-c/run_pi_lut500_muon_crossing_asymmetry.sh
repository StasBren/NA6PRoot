#!/usr/bin/env bash
# Re-fold frozen Jpsi Geant4/Heed with p3 and optional p4 LUT500
# at true single-muon physical gas-crossing level.
# No Garfield, Geant4, or new Heed library simulation.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(check|pilot3|full3|pilot34|full34|status)$ ]];then
  echo "Usage: bash $0 check|pilot3|full3|pilot34|full34|status" >&2
  exit 2
fi
PY="$REPO/test_runs/mwpc_stage_c/.venv_pi_slides/bin/python"
BASE="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1"
if [[ "$MODE" == status ]];then
  for folder in "$BASE"/pilot3 "$BASE"/full3 "$BASE"/pilot34 "$BASE"/full34;do
    if [[ -f "$folder/manifest.json" ]];then
      echo "PRESENT: $folder"
    else
      echo "MISSING: $folder"
    fi
  done
  exit
fi
if [[ ! -x "$PY" ]];then
  echo "Missing expected scientific Python venv: $PY" >&2;exit 2
fi
"$PY" -m py_compile \
 "$HERE/pi_lut500_muon_crossing_asymmetry.py" \
 "$HERE/test_pi_lut500_muon_crossing_asymmetry.py"
"$PY" "$HERE/test_pi_lut500_muon_crossing_asymmetry.py"
if [[ "$MODE" == check ]];then
  echo "PASS: Python syntax and crossing-level synthetic conservation tests"
  exit
fi
case "$MODE" in
 pilot3) OUTPUT="$BASE/pilot3"; PITCH=(3);CHUNKS=1;MAX=200;;
 full3) OUTPUT="$BASE/full3"; PITCH=(3);CHUNKS=2;MAX=0;;
 pilot34) OUTPUT="$BASE/pilot34"; PITCH=(3 4);CHUNKS=1;MAX=200;;
 full34) OUTPUT="$BASE/full34"; PITCH=(3 4);CHUNKS=2;MAX=0;;
esac
if [[ "$MODE" == *34 ]];then
  echo "Requiring the separate p4 LUT500 to be COMPLETE and SHA-verified:"
  bash "$HERE/run_pi_pitch_fieldwire_500_p4.sh" verify
fi
if [[ -e "$OUTPUT/manifest.json" ||
      -e "$OUTPUT/summary.json" ||
      -e "$OUTPUT/crossings.csv.gz" ||
      -e "$OUTPUT/crossings.csv.gz.incomplete" ]];then
  echo "SAFE STOP: output already exists, inspect rather than overwrite: $OUTPUT" >&2
  exit 2
fi
if [[ "$MAX" == 0 && ! -f "$BASE/pilot3/manifest.json" ]];then
  echo "Run pilot3 first before any full response analysis." >&2
  exit 2
fi
if [[ "$MODE" == full34 && ! -f "$BASE/pilot34/manifest.json" ]];then
  echo "Run pilot34 first." >&2
  exit 2
fi
echo "One actual gas crossing = one data point (NOT pooled J/psi parent A)."
echo "Mode=$MODE, parent slots=$((CHUNKS*50000)), max-crossings/station/chunk=$MAX"
echo "Pitches: ${PITCH[*]} mm. Observables: hemisphere A, top/bottom sector A."
echo "Starting READ-ONLY Geant4/Heed/LUT500 fold; no new microscopic Garfield."
"$PY" "$HERE/pi_lut500_muon_crossing_asymmetry.py" \
 --chunks "$CHUNKS" --stations 0 1 2 3 4 5 \
 --pitches "${PITCH[@]}" --max-crossings "$MAX" --out "$OUTPUT"
echo "DONE: $OUTPUT"
if command -v wslpath >/dev/null 2>&1;then
  echo "Windows folder: $(wslpath -w "$OUTPUT")"
fi

#!/usr/bin/env bash
# p4mm LUT20 per-physical-muon-crossing response, independent of ongoing p4 LUT500.
# Pure READ ONLY native Geant4 crossing + stable Heed birth + frozen LUT20 lookup.
# Do not overwrite old 2M LUT20, 100k p3 LUT500 or 1.2M Garfield p4 LUT500.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(check|pilot|full|status)$ ]];then
  echo "Usage: bash $0 check|pilot|full|status" >&2
  exit 2
fi
BASE="$REPO/test_runs/mwpc_stage_c/pi_lut20_p4_muon_crossing_asymmetry_v1"
PY="$REPO/test_runs/mwpc_stage_c/.venv_pi_slides/bin/python"
SCRIPT="$HERE/pi_lut500_muon_crossing_asymmetry.py"
TEST="$HERE/test_pi_lut500_muon_crossing_asymmetry.py"

if [[ "$MODE" == "status" ]];then
  for folder in "$BASE/pilot200" "$BASE/full100k";do
    if [[ -f "$folder/manifest.json" ]];then
      echo "AVAILABLE: $folder"
    elif [[ -f "$folder/crossings.csv.gz.incomplete" ]];then
      echo "INCOMPLETE, preserve/inspect before restarting: $folder"
    else
      echo "NOT YET RUN: $folder"
    fi
  done
  exit 0
fi
if [[ ! -x "$PY" ]];then
  echo "Missing Python environment: $PY" >&2
  exit 2
fi
"$PY" -m py_compile "$SCRIPT" "$TEST"
"$PY" "$TEST"
if [[ "$MODE" == "check" ]];then
  echo "PASS: syntax and synthetic LUT20/LUT500 crossing-conservation tests."
  echo "No Geant4/Garfield/Heed physics run."
  exit 0
fi

if [[ "$MODE" == "pilot" ]];then
  OUTPUT="$BASE/pilot200"
  CHUNKS=1; MAX=200
else
  OUTPUT="$BASE/full100k"
  CHUNKS=2; MAX=0
  if [[ ! -f "$BASE/pilot200/manifest.json" ]];then
    echo "Run pilot and inspect the figures before full." >&2
    exit 2
  fi
fi
if [[ -e "$OUTPUT/summary.json" ||
      -e "$OUTPUT/manifest.json" ||
      -e "$OUTPUT/crossings.csv.gz" ||
      -e "$OUTPUT/crossings.csv.gz.incomplete" ]];then
  echo "SAFE STOP: existing result (or incomplete result) preserved: $OUTPUT" >&2
  exit 2
fi
echo "READ-ONLY p=4mm LUT20 analysis for 6 stations; field wires 0 / 75 um."
echo "Original 20-trial/cell maps, not the ongoing new p4 500-trial/cell maps."
echo "One physical muon gas crossing per observation; no parent pooling."
echo "NO new Garfield physics. Input generated J/psi parent slots: $((CHUNKS*50000))."
"$PY" "$SCRIPT" --chunks "$CHUNKS" --stations 0 1 2 3 4 5 \
  --pitches 4 --lut-trials 20 --max-crossings "$MAX" --out "$OUTPUT"
echo "PASS output: $OUTPUT"
if command -v wslpath >/dev/null 2>&1;then
  echo "Windows path: $(wslpath -w "$OUTPUT")"
fi

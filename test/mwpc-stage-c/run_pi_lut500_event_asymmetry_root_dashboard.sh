#!/usr/bin/env bash
# p=3mm LUT500 parent-event asymmetry for six MS stations and physical DETs.
# Never reruns Geant4 / Garfield, and never overwrites old LUT20 dashboards.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
if [[ "$#" -ne 1 || ! "$1" =~ ^(check|pilot|full|status|style)$ ]];then
 echo "Usage: bash $0 check|pilot|full|status|style" >&2
 exit 2
fi
MODE="$1"
BASE="$REPO/test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1"
PY="$REPO/test_runs/mwpc_stage_c/.venv_pi_slides/bin/python"
SCRIPT="$HERE/pi_lut500_event_asymmetry_root_dashboard.py"
TEST="$HERE/test_pi_lut500_event_asymmetry_root_dashboard.py"
STYLE="$HERE/pi_lut500_event_asymmetry_root_style_refresh.py"
if [[ "$MODE" == "status" ]];then
 for f in "$BASE"/pilot200/dashboard.html "$BASE"/dashboard.html;do
  if [[ -f "$f" ]];then
   echo "EXISTS: $f"
   if command -v wslpath >/dev/null 2>&1;then
    echo "Windows: $(wslpath -w "$f")"
   fi
  else
   echo "MISSING: $f"
  fi
 done
 exit
fi
if [[ ! -x "$PY" ]];then
 echo "Missing existing Python venv: $PY" >&2
 exit 2
fi
"$PY" -m py_compile "$SCRIPT" "$TEST" "$STYLE"
"$PY" "$TEST"
if [[ "$MODE" == "check" ]];then
 echo "PASS: LUT500 parent-asymmetry dashboard test and syntax"
 exit 0
fi
if [[ "$MODE" == "style" ]];then
 for folder in "$BASE/pilot200" "$BASE";do
  if [[ -f "$folder/manifest.json" ]];then
   "$PY" "$STYLE" --folder "$folder"
  else
   echo "SKIP no verified dashboard folder: $folder"
  fi
 done
 echo "PASS style-only: no Geant4, Garfield, Heed or LUT500 physics recalculation"
 exit 0
fi
if [[ "$MODE" == "pilot" ]];then
 OUTPUT="$BASE/pilot200"
 if [[ -f "$OUTPUT/dashboard.html" ]];then
  echo "SAFE SKIP: existing pilot preserved: $OUTPUT/dashboard.html"
  exit 0
 fi
 echo "PILOT: 200 accepted gas crossings/station, one chunk; NOT full parent exposure"
 "$PY" "$SCRIPT" --chunks 1 --stations 0 1 2 3 4 5 \
  --max-crossings 200 --out "$OUTPUT"
else
 PILOT="$BASE/pilot200/manifest.json"
 if [[ ! -f "$PILOT" ]];then
  echo "Run pilot first: bash $0 pilot" >&2
  exit 2
 fi
 OUTPUT="$BASE"
 if [[ -f "$OUTPUT/dashboard.html" ]];then
  echo "SAFE SKIP: existing full LUT500 dashboard preserved: $OUTPUT/dashboard.html"
  exit 0
 fi
 echo "FULL TWO CHUNKS: 100k generated Jpsi parent slots; p=3mm LUT500 only"
 "$PY" "$SCRIPT" --chunks 2 --stations 0 1 2 3 4 5 \
  --max-crossings 0 --out "$OUTPUT"
fi
echo "HTML: $OUTPUT/dashboard.html"
echo "No old LUT20 physics data were modified."
if command -v wslpath >/dev/null 2>&1;then
 echo "Windows path: $(wslpath -w "$OUTPUT/dashboard.html")"
fi

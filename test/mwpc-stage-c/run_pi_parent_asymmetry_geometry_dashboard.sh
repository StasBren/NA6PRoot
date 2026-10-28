#!/usr/bin/env bash
# 3mm LUT500 + 3/4mm LUT20 parent-level A and electron-counter dashboard.
# Existing parent-level and crossing-level physics are READ ONLY.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(check|build|replace-html|status)$ ]];then
  echo "Usage: bash $0 check|build|replace-html|status" >&2
  exit 2
fi
PY="$REPO/test_runs/mwpc_stage_c/.venv_pi_slides/bin/python"
SCRIPT="$HERE/pi_parent_asymmetry_geometry_dashboard.py"
TEST="$HERE/test_pi_parent_asymmetry_geometry_dashboard.py"
OLD="$REPO/test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1"
CROSS="$REPO/test_runs/mwpc_stage_c/pi_lut20_p4_muon_crossing_asymmetry_v1/full34"
OUTPUT="$OLD/dashboard_geometry_lut_compare.html"
if [[ "$MODE" == status ]];then
  for p in "$OLD/summary.json" "$OLD/manifest.json" \
           "$CROSS/summary.json" "$CROSS/manifest.json" \
           "$CROSS/crossings.csv.gz" "$OUTPUT";do
     if [[ -f "$p" ]];then echo "EXISTS $p";else echo "MISSING $p";fi
  done
  exit
fi
[[ -x "$PY" ]] || { echo "Missing scientific Python $PY" >&2;exit 2; }
"$PY" -m py_compile "$SCRIPT" "$TEST"
"$PY" "$TEST"
if [[ "$MODE" == check ]];then
  echo "PASS: parent-versus-crossing A, noFW/FW75 electron conservation regression"
  exit 0
fi
[[ -f "$OLD/manifest.json" && -f "$CROSS/manifest.json" ]] || {
  echo "Missing complete original p3 LUT500 or full34 LUT20 source." >&2
  exit 2
}
args=(--old "$OLD" --full34 "$CROSS" --out "$OUTPUT")
if [[ "$MODE" == replace-html ]];then
  args+=(--replace-html)
fi
echo "Grouping saved 100k parent J/psi crossing-level LUT20 raw data."
echo "No Geant4/Heed/Garfield physics rerun. Original data are immutable."
nice -n 15 "$PY" "$SCRIPT" "${args[@]}"
if command -v wslpath >/dev/null 2>&1;then
  echo "Windows HTML: $(wslpath -w "$OUTPUT")"
fi

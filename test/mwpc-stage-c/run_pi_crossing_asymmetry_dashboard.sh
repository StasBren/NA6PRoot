#!/usr/bin/env bash
# Dedicated CROSSING-level NA60+ dashboard; never overwrite the old parent-event dashboard.
# Before rendering LUT500 p3, use the existing crossing-level fold of one muon gas hit.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(check|status|lut20-full34|lut500-p3|lut500-p34)$ ]];then
  echo "Usage: bash $0 check|status|lut20-full34|lut500-p3|lut500-p34" >&2
  exit 2
fi

PY="$REPO/test_runs/mwpc_stage_c/.venv_pi_slides/bin/python"
if [[ "$MODE" == status ]];then
  for d in \
    "$REPO/test_runs/mwpc_stage_c/pi_lut20_p4_muon_crossing_asymmetry_v1/full34" \
    "$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full3" \
    "$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full34";do
    [[ -d "$d" ]] || continue
    printf "%s: raw=%s summary=%s manifest=%s new-dashboard=%s\n" "$d" \
       "$([[ -f "$d/crossings.csv.gz" ]] && echo YES || echo NO)" \
       "$([[ -f "$d/summary.json" ]] && echo YES || echo NO)" \
       "$([[ -f "$d/manifest.json" ]] && echo YES || echo NO)" \
       "$([[ -f "$d/crossing_dashboard.html" ]] && echo YES || echo NO)"
  done
  exit 0
fi
[[ -x "$PY" ]] || { echo "Missing scientific Python: $PY" >&2; exit 2; }
"$PY" -m py_compile \
    "$HERE/pi_crossing_asymmetry_dashboard.py" \
    "$HERE/test_pi_crossing_asymmetry_dashboard.py" \
    "$HERE/pi_lut500_muon_crossing_asymmetry.py"
"$PY" "$HERE/test_pi_crossing_asymmetry_dashboard.py"
if [[ "$MODE" == check ]];then
  echo "PASS: dashboard syntax, 6-station crossing histogram conservation and SHA regression."
  exit 0
fi

case "$MODE" in
 lut20-full34)
   DIR="$REPO/test_runs/mwpc_stage_c/pi_lut20_p4_muon_crossing_asymmetry_v1/full34"
   ;;
 lut500-p3)
   DIR="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full3"
   if [[ ! -f "$DIR/manifest.json" ]];then
     echo "Recomputing CROSSING-level p3 LUT500 from original native Geant4/Heed inputs."
     echo "This does NOT change or overwrite the original parent-event LUT500 dashboard."
     PILOT="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/pilot3/manifest.json"
     if [[ ! -f "$PILOT" ]];then
       bash "$HERE/run_pi_lut500_muon_crossing_asymmetry.sh" pilot3
     fi
     bash "$HERE/run_pi_lut500_muon_crossing_asymmetry.sh" full3
   fi
   ;;
 lut500-p34)
   DIR="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full34"
   bash "$HERE/run_pi_pitch_fieldwire_500_p4.sh" verify
   if [[ ! -f "$DIR/manifest.json" ]];then
     echo "Recomputing p3+p4 LUT500 CROSSING-level response."
     PILOT="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/pilot34/manifest.json"
     P3PILOT="$REPO/test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/pilot3/manifest.json"
     if [[ ! -f "$P3PILOT" ]];then
       bash "$HERE/run_pi_lut500_muon_crossing_asymmetry.sh" pilot3
     fi
     if [[ ! -f "$PILOT" ]];then
       bash "$HERE/run_pi_lut500_muon_crossing_asymmetry.sh" pilot34
     fi
     bash "$HERE/run_pi_lut500_muon_crossing_asymmetry.sh" full34
   fi
   ;;
esac

[[ -f "$DIR/manifest.json" ]] || {
  echo "No complete crossing-level scientific run: $DIR" >&2
  echo "No dashboard was published; original old parent-level dashboard is unchanged." >&2
  exit 2
}
echo "Rendering SHA-validated CROSSING-level dashboard from: $DIR"
"$PY" "$HERE/pi_crossing_asymmetry_dashboard.py" --input-dir "$DIR"
echo "PASS: $DIR/crossing_dashboard.html"
if command -v wslpath >/dev/null 2>&1;then
  echo "Windows file path: $(wslpath -w "$DIR/crossing_dashboard.html")"
fi

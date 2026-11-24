#!/usr/bin/env bash
# SYNTHETIC end-to-end test of the geometry-study chain (toy transport from synthetic_paths.py; NOT NA60+ physics):
# baseline and alternative geometry on the SAME toy paths: restrip -> S0 (frozen replay under the overlay) -> Green
# table from ONE shared grid -> G closure -> kernels -> K2A -> accuracy compare, geometry-effect compare, identities.
#   bash tests/test_pipeline_geometry.sh WORKDIR NUMW_DIR BASE_LUT ALT_LUT ALT_GEOM.json
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; SA="$HERE/../spatial-ablation"; TA="$HERE/../timing-ablation"
W="${1:?}"; NUMW="${2:?}"; L0="${3:?}"; L1="${4:?}"; G1="${5:?}"; G0="$HERE/geometries/baseline.json"; T=t
[[ ! -e "$W" ]] || { echo "$W exists"; exit 2; }
mkdir -p "$W/base/muon_0001" "$W/alt/muon_0001"; cp "$G1" "$W/alt/geometry.json"
python3 "$SA/synthetic_paths.py" muon "$W/base/muon_0001/T1_export" --seed 4 > /dev/null
echo '{"patch_shift_cm": 0.23, "u_in_mm": 7.0, "w_in_mm": 0.0, "base_U": 0, "base_V": 0, "stereo_tan": 0.1, "strip_pitch_mm": 1.7}' > "$W/frame.json"
python3 "$HERE/geom_restrip_export.py" "$W/base/muon_0001/T1_export" "$W/frame.json" "$G1" "$W/alt/muon_0001/T1_export"
# identity restrip of the toy export must be byte-identical
python3 "$HERE/geom_restrip_export.py" "$W/base/muon_0001/T1_export" "$W/frame.json" "$G0" "$W/identity" > /dev/null
cmp "$W/identity/vertices.f64" "$W/base/muon_0001/T1_export/vertices.f64" && cmp "$W/identity/ions.csv" "$W/base/muon_0001/T1_export/ions.csv"
echo "PASS identity restrip byte-identical"
python3 "$SA/s_grid_spec.py" "$W/base/muon_0001/T1_export" --nr 8 > "$W/base/grid_spec_$T.txt"
python3 "$SA/s_grid_spec.py" "$W/alt/muon_0001/T1_export" --nr 8 --strip-pitch-mm "$(python3 -c "import json;print(json.load(open('$G1'))['strip_pitch_mm'])")" > "$W/alt/grid_spec_$T.txt"
RU=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$W/base/grid_spec_$T.txt"); WI=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$W/base/grid_spec_$T.txt")
[[ "$(sed -n 's/^--wires//p' "$W/base/grid_spec_$T.txt")" == "$(sed -n 's/^--wires//p' "$W/alt/grid_spec_$T.txt")" ]] && echo "PASS grid spec (wires, radii) unchanged by restriping"
python3 "$SA/synthetic_paths.py" grid "$W/base/green_paths_$T" --wires="$WI" --r-um="$RU" --ntheta 16 > /dev/null
for g in base alt; do
  GJ=$G0; L=$L0; [[ $g == alt ]] && { GJ=$G1; L=$L1; }
  R="$W/$g/muon_0001"; DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$W/$g/grid_spec_$T.txt")
  python3 "$HERE/geom_run.py" "$GJ" "$NUMW" replay "$L" "$R/T1_export" "$R/S0/replay" --half-strips 5 > "$W/$g.s0.log"
  python3 "$HERE/geom_run.py" "$GJ" "$NUMW" prompt "$L" "$R/T1_export" "$R/S0/prompt" --replay-dir "$NUMW" --half-strips 5 > /dev/null
  python3 "$HERE/geom_run.py" "$GJ" "$NUMW" green "$L" "$W/base/green_paths_$T" "$W/$g/green_$T" --replay-dir "$NUMW" --d-max-mm "$DM" --ntau 120 | tail -1
  python3 "$SA/s_synthesize.py" "$W/$g/green_$T" "$R/T1_export" "$R/S0/replay/waveforms.csv" "$R" --variants G --suffix "_$T" > "$W/$g.synth.log"
  python3 "$SA/s_kernel.py" build "$W/$g/green_$T" "$W/$g/kernel_$T" > /dev/null
  python3 "$SA/s_kernel.py" moments "$R/T1_export" "$R/cloud_moments.csv" > /dev/null
  python3 "$SA/s_kernel.py" apply "$W/$g/kernel_$T" "$R/cloud_moments.csv" "$R/S0/replay/waveforms.csv" "$R" --model K2A --suffix "_$T" --repeat 1 > /dev/null
done
# negative control: the alternative geometry must refuse the baseline LUT (frozen loader would accept it silently)
if python3 "$HERE/geom_run.py" "$G1" "$NUMW" replay "$L0" "$W/alt/muon_0001/T1_export" "$W/neg" > "$W/neg.log" 2>&1; then
  echo "FAIL: alternative geometry accepted the baseline LUT"; exit 1; fi
grep -q "meta stereo_tan" "$W/neg.log" && echo "PASS LUT/geometry mismatch refused"
A1=$(python3 -c "import json;print(json.load(open('$G1'))['stereo_tan'])")
python3 "$TA/t_shaped_compare.py" "$W/alt" --variants "S0,G_$T,K2Akern_$T" --tp 100 --no-prompt --stereo-tan "$A1" --bootstrap 50 --out "$W/alt/acc" > /dev/null
python3 "$TA/t_shaped_compare.py" "$W/base" --variants "S0,G_$T,K2Akern_$T" --tp 100 --no-prompt --bootstrap 50 --out "$W/base/acc" > /dev/null
python3 - "$W" <<'PY'
import json, sys
W = sys.argv[1]; L = {}
for g in ("base", "alt"):
    for l in open(f"{W}/{g}/acc/per_muon_Tp100.jsonl"):
        r = json.loads(l)
        if r["variant"] != "S0":
            L[g, r["variant"]] = max(r["U_L2_rel"], r["V_L2_rel"])
            print(f"  {g:4s} {r['variant']:10s} shaped rel L2 vs S0 (max U,V): {L[g, r['variant']]:.4g}  (toy grid nr 8 x ntheta 16)")
for v in ("G_t", "K2Akern_t"):
    assert L["alt", v] <= 2 * L["base", v] + 1e-3, f"{v}: closure under the alternative geometry much worse than baseline"
print("PASS closure under the alternative geometry comparable to the baseline closure on the same toy grid")
PY
python3 "$HERE/geom_physics_compare.py" --base "$W/base" --alt "$W/alt" --base-geom "$G0" --alt-geom "$G1" --variants "S0,K2Akern_$T" --tp 100 --no-prompt --out "$W/effect" > /dev/null
grep -E "U/V charge ratio|u - truth|V mult" "$W/effect/summary.md"
# response engine: moments tagged with the baseline strip frame must be refused by the alternative kernel
EC="$SA/engine/engine_cli.py"; KID=$(python3 "$EC" info "$W/alt/kernel_$T" | python3 -c "import json,sys; print(json.load(sys.stdin)['kernel_id'])")
mkdir -p "$W/em_bad" "$W/em_ok"; cp "$W/alt/muon_0001/cloud_moments.csv" "$W/em_bad/m.moments.csv"; cp "$W/alt/muon_0001/cloud_moments.csv" "$W/em_ok/m.moments.csv"
echo '{"stereo_tan": 0.1, "strip_pitch_mm": 1.7}' > "$W/em_bad/m.moments.csv.frame.json"
echo "{\"stereo_tan\": $A1, \"strip_pitch_mm\": 1.7}" > "$W/em_ok/m.moments.csv.frame.json"
if python3 "$EC" run --kernel "$W/alt/kernel_$T" --kernel-id "$KID" --moments-dir "$W/em_bad" --out "$W/eo_bad" --impl vector > "$W/eo_bad.log" 2>&1; then
  echo "FAIL: engine accepted moments of another strip geometry"; exit 1; fi
python3 "$EC" run --kernel "$W/alt/kernel_$T" --kernel-id "$KID" --moments-dir "$W/em_ok" --out "$W/eo_ok" --impl vector --format replay > /dev/null
python3 "$SA/s_compare_waveforms.py" "$W/alt/muon_0001/K2Akern_$T" "$W/eo_ok/m" --report "$W/eo_ok.md" > /dev/null
echo "PASS engine: geometry-mismatched moments refused; matched moments reproduce s_kernel K2A"
python3 "$HERE/geom_manifest.py" "$W/base" "$W/alt" "$T" --out "$W/identities.json" | tail -1
echo "PASS synthetic geometry pipeline; work dir $W"

#!/usr/bin/env bash
# SYNTHETIC tests of the compiled absolute-strip K2A evaluator (toy kernel, toy moments: NOT physics).
#   bash tests/test_strip_eval_synthetic.sh WORKDIR        (WORKDIR must not exist)
# A  C++ adapter == offline Python recentering (k2a_recentered_k2a.py with impl loop, vector, compiled) on 4 events
#    (12-primary multiwire event with shifts m = -8..+7 in both planes, near-entry event, shift tie, late cloud);
#    labels and time grid identical, ion and electron L2 / strip charge / total charge / U-V ratio <= 1e-12
# B  imposed shift table (groups.json form): equal to the rule -> identical files; one shift changed -> equal to Python
#    run with the same table, and reported as differing from the rule
# C  displacement symmetry: the whole event displaced by (+3, -2) and (-4, +5) strips gives the same response on labels
#    shifted by the displacement (<= 1e-12)
# D  statuses: wire not in kernel, r_eff below / above the grid, local D outside the table (spread primary), negative
#    time, zero charge, NaN time, NaN D (undefined shift) -> BLOCKED (exit 4), no waveform file, blocked charge reported; the offline Python engine
#    refuses the same wire / r_eff / D cases; T_AFTER_WINDOW and SHIFT_NEAR_TIE are informational
# E  in-memory chain: births -> MomentAccumulator -> Event -> Evaluator == the same on the moments CSV (bitwise) and
#    == Python on that CSV (<= 1e-12); blocked moment event -> INPUT_BLOCKED
# F  refusals: wrong kernel id, existing output, unknown model; comparator negatives; repeats bit-identical
# G  run_m2_adapter.sh explicit coverage: E2b 0560 + 0997 (one grouping equal to the rule, one not), accepted T1 0000,
#    one validated multiwire + displaced regression -> PASS. Missing 0997, no T1 run, no regression run, a regression
#    that is neither multiwire nor displaced, or a both-refuse regression -> INCONCLUSIVE. A regression whose T1 verdict
#    is not PASS, a T1 item that is not slot 0000, a failing preload benchmark, a corrupted offline response -> FAIL
# H  k2a_select_regression.py: N = last included ORIGINAL primary index + 1 when primaries have zero ions (absent from
#    the CSV); budget over the same prefix; census checked against event.json (inconsistent -> skipped); sealed slot
#    skipped, sealed directory refused
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; CIR="$(cd "$HERE/.." && pwd)"
ENG="$CIR/spatial-ablation/engine"; REC="$CIR/k2a-campaign/recenter/k2a_recentered_k2a.py"
W="${1:?usage: test_strip_eval_synthetic.sh WORKDIR}"; [[ ! -e "$W" ]] || { echo "WORKDIR exists" >&2; exit 2; }
mkdir -p "$W"; W="$(cd "$W" && pwd)"; cd "$W"
die(){ echo "FAIL $*" >&2; exit 1; }
ok(){ echo "PASS $*"; }
CXX="${CXX:-g++}"
$CXX -O3 -std=c++17 -Wall -Wextra -o k2a_strip_eval "$HERE/k2a_strip_eval.cc" "$ENG/kernel_eval.cc"
$CXX -O3 -std=c++17 -Wall -Wextra -o inmem "$HERE/tests/strip_eval_inmemory_test.cc" "$ENG/kernel_eval.cc"
$CXX -O3 -std=c++17 -shared -fPIC -o libkernel_eval.so "$ENG/kernel_eval.cc"
export NA60_KERNEL_EVAL_LIB="$W/libkernel_eval.so"
python3 "$HERE/tests/make_toy_kernel.py" "$W/toy" > /dev/null
KD="$W/toy/kern"
KID=$(python3 -c "import sys; sys.path.insert(0, '$ENG'); import cathode_engine as ce; print(ce.KernelTable('$KD').kernel_id)")
CMP="$HERE/k2a_strip_compare.py"

# A
for ev in multi near tie late; do
  M="$W/toy/moments_$ev.csv"
  ./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "cc_$ev" --repeat 3 > "cc_$ev.log" || die "A $ev C++ ($(cat cc_$ev.log))"
  for impl in loop vector compiled; do
    python3 "$REC" "$KD" "$M" "py_${ev}_$impl" --impl "$impl" --kernel-id "$KID" > /dev/null
    python3 "$CMP" "py_${ev}_$impl" "cc_$ev" "cmp_${ev}_$impl.json" --name "$ev vs $impl" > "cmp_${ev}_$impl.log" || die "A $(cat cmp_${ev}_$impl.log)"
  done
  python3 - "cc_$ev/result.json" <<'PY' || die "A $ev result.json"
import json, sys
r = json.load(open(sys.argv[1])); assert r["usable"] and r["timing"]["repeats_bit_identical"] and r["kernel_verified"], r
PY
done
grep -q SHIFT_NEAR_TIE cc_tie/result.json || die "A tie not flagged"
grep -q T_AFTER_WINDOW cc_late/result.json || die "A late cloud not flagged"
python3 - <<'PY' || die "A multi shift range"
import csv
s = list(csv.DictReader(open("cc_multi/shifts.csv")))
mu = [int(r["m_U"]) for r in s]; mv = [int(r["m_V"]) for r in s]
assert min(mu) <= -6 and max(mu) >= 6 and min(mv) <= -4 and max(mv) >= 4, (mu, mv)
g = list(csv.DictReader(open("cc_multi/groups.csv"))); assert len(g) >= 10
c = list(csv.DictReader(open("cc_multi/clouds.csv"))); assert len({r["wire"] for r in c}) >= 4
PY
ok "A C++ adapter == Python offline recentering (loop, vector, compiled) on 4 events; worst ion L2 $(python3 -c "
import json, glob; print('%.1e' % max(max(json.load(open(f))['metrics']['ion'][k] for k in ('U_L2', 'V_L2')) for f in glob.glob('cmp_*.json')))")"

# B
python3 - <<'PY'
import csv, json
s = {r["primary"]: [int(r["m_U"]), int(r["m_V"])] for r in csv.DictReader(open("cc_multi/shifts.csv"))}
json.dump({"primary_shift": s}, open("same.json", "w"))
s2 = dict(s); k = sorted(s2, key=int)[3]; s2[k] = [s2[k][0] + 1, s2[k][1] - 1]
json.dump({"primary_shift": s2}, open("changed.json", "w"))
PY
./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_multi.csv" --out cc_same --shifts same.json > /dev/null
for f in replay/waveforms.csv prompt/prompt_waveforms.csv groups.csv; do cmp -s "cc_multi/$f" "cc_same/$f" || die "B imposed = rule differs in $f"; done
./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_multi.csv" --out cc_changed --shifts changed.json > /dev/null
python3 "$REC" "$KD" "$W/toy/moments_multi.csv" py_changed --groups changed.json --impl vector --kernel-id "$KID" > /dev/null
python3 "$CMP" py_changed cc_changed cmp_changed.json --name "imposed table" > /dev/null || die "B changed table vs Python"
grep -q '"shifts_differing_from_computed": 1' cc_changed/result.json || die "B differing shift not reported"
ok "B imposed shift table: equal to the rule -> identical files; changed table == Python with the same table, reported"

# C
for d in "3 -2" "-4 5"; do
  set -- $d
  ./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_multi.csv" --out "cc_disp_$1_$2" --displace-event "$1" "$2" > /dev/null
  python3 "$CMP" cc_multi "cc_disp_$1_$2" "cmp_disp_$1_$2.json" --label-shift "$1" "$2" --name "displaced $1 $2" > "cmp_disp.log" || die "C $(cat cmp_disp.log)"
done
ok "C displacement symmetry: events displaced by (+3,-2) and (-4,+5) strips reproduce the response on shifted labels"

# D
declare -A EXP=([bad_wire]=WIRE_NOT_IN_KERNEL [bad_r_low]=R_EFF_BELOW_GRID [bad_r_high]=R_EFF_ABOVE_GRID [bad_spread]=D_U_OUTSIDE_TABLE
                [bad_t]=T_NEGATIVE [bad_q]=NONPOSITIVE_CHARGE [bad_nan]=NONFINITE_MOMENT [bad_dnan]=SHIFT_UNDEFINED)
for b in "${!EXP[@]}"; do
  set +e; ./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_$b.csv" --out "cc_$b" > "cc_$b.log" 2>&1; rc=$?; set -e
  [[ $rc -eq 4 ]] || die "D $b: exit $rc ($(cat cc_$b.log))"
  grep -q "${EXP[$b]}" "cc_$b/clouds.csv" || die "D $b: status ${EXP[$b]} missing"
  [[ ! -e "cc_$b/replay" ]] || die "D $b: a waveform was written for a blocked event"
  python3 - "cc_$b/result.json" <<'PY' || die "D $b blocked charge"
import json, sys
r = json.load(open(sys.argv[1])); assert not r["usable"] and (r["blocked_charge_fC"] > 0 or "NONPOSITIVE" in r["status_names"]), r
PY
done
for b in bad_wire bad_r_low bad_r_high bad_spread; do
  if python3 "$REC" "$KD" "$W/toy/moments_$b.csv" "py_$b" --kernel-id "$KID" > "py_$b.log" 2>&1; then die "D Python accepted $b"; fi
done
ok "D statuses: 8 blocking cases -> exit 4, named status per cloud, no waveform, blocked charge; Python refuses the wire/r_eff/D cases too"

# E
./inmem "$KD" "$KID" im > inmem.log || { cat inmem.log; die "E in-memory"; }
python3 "$REC" "$KD" im/inmem.moments.csv py_inmem --impl vector --kernel-id "$KID" > /dev/null
python3 "$CMP" py_inmem im/inmem_result cmp_inmem.json --name "in-memory chain" > cmp_inmem.log || die "E $(cat cmp_inmem.log)"
ok "E in-memory chain == CSV path (bitwise) == Python ($(cat cmp_inmem.log | cut -c1-60))"

# F
set +e
./k2a_strip_eval --kernel "$KD" --kernel-id "$(printf '0%.0s' {1..64})" --moments "$W/toy/moments_near.csv" --out f1 > f1.log 2>&1; f1=$?
./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_near.csv" --out cc_near > f2.log 2>&1; f2=$?
./k2a_strip_eval --kernel "$KD" --kernel-id "$KID" --moments "$W/toy/moments_near.csv" --out f3 --model K3 > f3.log 2>&1; f3=$?
./k2a_strip_eval --kernel "$KD" --moments "$W/toy/moments_near.csv" --out f4 > f4.log 2>&1; f4=$?
cp -r cc_near bad_cmp; python3 - <<'PY'
import numpy as np
p = "bad_cmp/replay/waveforms.csv"; h = open(p).readline().strip()
a = np.loadtxt(p, delimiter=",", skiprows=1); a[100, 3] *= 1 + 1e-9
np.savetxt(p, a, delimiter=",", header=h, comments="")
PY
python3 "$CMP" cc_near bad_cmp c1.json > /dev/null; c1=$?
python3 "$CMP" cc_near cc_multi c2.json > /dev/null; c2=$?
python3 "$CMP" cc_near nowhere c3.json > /dev/null; c3=$?
set -e
[[ $f1 -eq 2 && $f2 -eq 2 && $f3 -eq 2 && $f4 -eq 2 ]] && grep -q "kernel id mismatch" f1.log || die "F refusals ($f1 $f2 $f3 $f4)"
[[ $c1 -eq 1 && $c2 -eq 1 && $c3 -eq 3 ]] || die "F comparator negatives ($c1 $c2 $c3)"
ok "F refusals: wrong kernel id, existing output, unknown model, missing id; comparator: 1e-9 perturbation FAIL, other labels FAIL, missing INCONCLUSIVE"

# G
mk_e2b(){ # slot shifts.json
  local d="$W/E2B/recenter_experiment_x/muon_$1"; mkdir -p "$d/prepared" "$d/groups"
  cp "$W/toy/moments_multi.csv" "$d/prepared/cloud_moments.csv"; cp "$2" "$d/groups/groups.json"
  python3 "$REC" "$KD" "$d/prepared/cloud_moments.csv" "$d/K2A_local" --groups "$d/groups/groups.json" --kernel-id "$KID" > /dev/null
}
mk_t1(){ # name slot moments verdict
  local d="$W/T1/$1"; mkdir -p "$d/R1"; cp "$3" "$d/R1/inline.moments.csv"
  printf 'label SYNTHETIC\nslot %s\nmax_heed_primaries 3\n' "$2" > "$d/provenance.txt"
  printf '{"label": "SYNTHETIC", "verdict": "%s"}\n' "$4" > "$d/T1_verdict.json"
}
mk_e2b 0560 same.json; mk_e2b 0997 changed.json
python3 "$HERE/tests/make_toy_kernel.py" --single-cloud "$W/toy/moments_near.csv" "$W/single.csv"
mk_t1 t1_0000 0000 "$W/toy/moments_near.csv" PASS
mk_t1 reg_0042 0042 im/inmem.moments.csv PASS
mk_t1 reg_0043 0043 "$W/single.csv" PASS
mk_t1 reg_0044 0044 "$W/toy/moments_bad_wire.csv" PASS
mk_t1 reg_0045 0045 im/inmem.moments.csv FAIL
mk_t1 t1_0007 0007 "$W/toy/moments_near.csv" PASS
mkdir -p "$W/failpre"
printf '#!/usr/bin/env bash\nfor a in "$@"; do [[ "$a" == --preload ]] && exit 2; done\nexec "%s/k2a_strip_eval" "$@"\n' "$W" > "$W/failpre/k2a_strip_eval"
chmod +x "$W/failpre/k2a_strip_eval"
E5="$W/E2B/recenter_experiment_x/muon_0560"; E9="$W/E2B/recenter_experiment_x/muon_0997"
m2(){ # outdir build e2b t1 regressions -> prints the exit code
  set +e
  MWPC_M2_SYNTHETIC=1 MWPC_K2A_INLINE_BUILD="$2" MWPC_K2A_KERNEL_DIR="$KD" MWPC_K2A_KERNEL_ID="$KID" MWPC_M2_REPEAT=2 \
    MWPC_M2_E2B_DIRS="$3" MWPC_M2_T1_RUN="$4" MWPC_M2_REGRESSION_RUNS="$5" bash "$HERE/run_m2_adapter.sh" "$W/$1" > "$W/$1.log" 2>&1
  local rc=$?; set -e; echo "$rc"
}
g=$(m2 m2_full "$W" "$E5 $E9" "$W/T1/t1_0000" "$W/T1/reg_0042")
[[ $g -eq 0 ]] || { cat "$W/m2_full.log"; die "G full coverage verdict ($g)"; }
grep -q "1 of 12 primaries: rule shift" "$W/m2_full/steps.tsv" && grep -q "0 of 12 primaries.*files identical" "$W/m2_full/steps.tsv" || die "G production-rule reporting"
python3 "$HERE/tests/check_m2_verdict.py" "$W/m2_full/m2_verdict.json" || die "G coverage record"
for k in "m2_no0997|$E5|$W/T1/t1_0000|$W/T1/reg_0042" "m2_noreg|$E5 $E9|$W/T1/t1_0000|" \
         "m2_nocov|$E5 $E9|$W/T1/t1_0000|$W/T1/reg_0043" "m2_refuse|$E5 $E9|$W/T1/t1_0000|$W/T1/reg_0044" \
         "m2_not1|$E5 $E9||$W/T1/reg_0042"; do
  IFS='|' read -r n e t r <<< "$k"
  g=$(m2 "$n" "$W" "$e" "$t" "$r"); [[ $g -eq 3 ]] || { tail -5 "$W/$n.log"; die "G $n: expected INCONCLUSIVE, exit $g"; }
done
for k in "m2_t1fail|$E5 $E9|$W/T1/t1_0000|$W/T1/reg_0045|$W" "m2_t1slot|$E5 $E9|$W/T1/t1_0007|$W/T1/reg_0042|$W" \
         "m2_prefail|$E5 $E9|$W/T1/t1_0000|$W/T1/reg_0042|$W/failpre"; do
  IFS='|' read -r n e t r bld <<< "$k"
  g=$(m2 "$n" "$bld" "$e" "$t" "$r"); [[ $g -eq 1 ]] || { tail -5 "$W/$n.log"; die "G $n: expected FAIL, exit $g"; }
done
grep -q "bench_preload.FAIL" <(tr '\t' '.' < "$W/m2_prefail/steps.tsv") || die "G preload failure not attributed"
cp -r "$W/E2B/recenter_experiment_x" "$W/E2B/bad"
python3 "$HERE/tests/make_toy_kernel.py" --corrupt "$W/E2B/bad/muon_0560/K2A_local/replay/waveforms.csv"
g=$(m2 m2_badoff "$W" "$W/E2B/bad/muon_0560 $W/E2B/bad/muon_0997" "$W/T1/t1_0000" "$W/T1/reg_0042")
[[ $g -eq 1 ]] && grep -q "reproduce_offline.FAIL" <(tr '\t' '.' < "$W/m2_badoff/steps.tsv") || die "G corrupted offline response ($g)"
ok "G run_m2_adapter.sh coverage: complete -> PASS; missing 0997 / no regression / no T1 / regression without coverage / both-refuse -> INCONCLUSIVE; T1 verdict FAIL / T1 not slot 0000 / preload benchmark failure / corrupted offline response -> FAIL"

# H
SD="$W/sel"; SN="$SD/moments_nonsealed/chunk_000"
mkdir -p "$SN" "$SD/moments_sealed" "$SD/inventory"
cp "$W/toy/moments_multi.csv" "$SN/muon_0001.moments.csv"
cp "$W/toy/moments_near.csv" "$SN/muon_0002.moments.csv"
cp "$W/toy/moments_multi.csv" "$SN/muon_0003.moments.csv"
# muon_0004: only ORIGINAL primaries 2 and 5 have ions (0, 1, 3, 4 and 6 are zero-ion, absent from the CSV)
python3 "$HERE/tests/make_toy_kernel.py" --gappy "$W/toy/moments_multi.csv" "$SN/muon_0004.moments.csv"
printf '{"n_primaries": 7, "n_primaries_zero_ions": 5}\n' > "$SN/muon_0004.moments.csv.event.json"
cp "$SN/muon_0004.moments.csv" "$SN/muon_0005.moments.csv"
printf '{"n_primaries": 7, "n_primaries_zero_ions": 3}\n' > "$SN/muon_0005.moments.csv.event.json"     # inconsistent census
printf 'slot,role\n0001,development\n0002,development\n0003,sealed_validation\n0004,development\n0005,development\n' > "$SD/inventory/slots.csv"
python3 "$HERE/k2a_select_regression.py" --moments-dir "$SD/moments_nonsealed" --slots "$SD/inventory/slots.csv" --kernel "$KD" \
  --out "$SD/cand.csv" --max-births 1e9 --top 50 > "$W/sel.log" || { cat "$W/sel.log"; die "H selector"; }
grep -q "1 files skipped: sealed" "$W/sel.log" && grep -q "1 skipped: census" "$W/sel.log" && ! grep -q "^0003," "$SD/cand.csv" || die "H skip counts ($(head -1 "$W/sel.log"))"
python3 "$HERE/tests/check_selector.py" "$SD/cand.csv" "$W/toy/moments_multi.csv" || die "H prefix N"
set +e; python3 "$HERE/k2a_select_regression.py" --moments-dir "$SD/moments_sealed" --slots "$SD/inventory/slots.csv" --kernel "$KD" --out "$SD/x.csv" > "$W/sel2.log" 2>&1; h2=$?; set -e
[[ $h2 -ne 0 ]] && grep -q "sealed moments directory" "$W/sel2.log" || die "H sealed directory accepted"
ok "H selector: zero-ion leading/intermediate primaries -> N = last index + 1 (3 and 6, not 1 and 2), budget over the same prefix; inconsistent census skipped; sealed slot skipped; sealed directory refused"
echo "PASS strip-evaluator synthetic tests (toy kernel, NOT physics); work dir $W"

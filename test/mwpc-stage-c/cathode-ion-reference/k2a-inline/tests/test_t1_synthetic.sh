#!/usr/bin/env bash
# SYNTHETIC end-to-end test of milestone 1 (stub Garfield++/ROOT: NOT physics, see tests/stubs/Garfield/StubCore.hh).
#   bash tests/test_t1_synthetic.sh WORKDIR      (WORKDIR must not exist)
# Builds the FROZEN and the K2A-inline generator sources, the archive dump and p3_birth_moments against the stubs, then:
#   A  T0 unit test (real header, no stubs needed)
#   B  allowed-diff checker: PASS on the inline copy; FAIL on four mutations (frozen line changed, //@frozen line
#      changed, broken block marker, K2A code outside a block)
#   C  frozen vs inline, no K2A options: every deterministic output byte-identical; with --k2a-moments-out the archive
#      and per-seed CSVs are still byte-identical to the frozen ones (the sink composition does not perturb anything)
#   D  the complete T1 runner in synthetic mode (R0-R4, comparisons, T2, frozen cross-check) -> PASS, MATCH
#   E  negatives: inside-anode birth -> exit 4, no moments file; outside-patch birth -> frozen abort kept (exit 2,
#      no K2A output); zero-ion primary -> informational, T1 comparison still PASS; corrupted dump / perturbed moment /
#      truncated batch -> FAIL; missing input -> INCONCLUSIVE; changed rerun -> identity FAIL; changed reservoir row
#      -> subset FAIL and frozen cross-check MISMATCH; bad options refused; existing outputs refused; sealed slot refused
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CIR="$(cd "$HERE/.." && pwd)"; AV="$(cd "$HERE/../../../mwpc-avalanche" && pwd)"
W="${1:?usage: test_t1_synthetic.sh WORKDIR}"; [[ ! -e "$W" ]] || { echo "WORKDIR exists" >&2; exit 2; }
mkdir -p "$W/build"; W="$(cd "$W" && pwd)"; BLD="$W/build"; ST="$HERE/tests/stubs"
ok(){ echo "PASS $*"; }
die(){ echo "FAIL $*" >&2; exit 1; }
CXX="${CXX:-g++}"; F="-std=c++17 -O2"
$CXX $F -I"$ST" -o "$BLD/frozen_stub" "$AV/phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx" 2>/dev/null
$CXX $F -I"$ST" -o "$BLD/mwpc_p3_noFW_penning057_k2a_inline_v1" "$AV/phaseB_native_muon_p3_noFW_penning057_k2a_inline_v1.cxx" 2>/dev/null
$CXX $F -I"$ST" -o "$BLD/k2a_archive_dump" "$HERE/k2a_archive_dump.cc"
$CXX $F -I"$ST" -I"$CIR" -I"$CIR/spatial-ablation" -o "$BLD/p3_birth_moments" "$CIR/spatial-ablation/p3_birth_moments.cc"
$CXX $F -Wall -Wextra -o "$BLD/moment_accumulator_t0" "$HERE/tests/moment_accumulator_t0.cc"
echo "SYNTHETIC stub build (not Garfield++/ROOT)" > "$BLD/k2a_inline_build_info.txt"
INL="$BLD/mwpc_p3_noFW_penning057_k2a_inline_v1"; FRZ="$BLD/frozen_stub"
awk 'NR==1||NR==2' "$CIR/../p3_wiretruth_1000_selection.csv" > "$W/input.csv"

# A
"$BLD/moment_accumulator_t0" "$W/t0" > "$W/t0.log" || { tail -20 "$W/t0.log"; die "T0"; }
ok "A T0: $(tail -1 "$W/t0.log")"

# B
python3 "$HERE/check_inline_generator_diff.py" > "$W/diff.log" || die "diff checker on the inline copy"
G="$AV/phaseB_native_muon_p3_noFW_penning057_k2a_inline_v1.cxx"
i=0
for m in 's/constexpr double startEnergy=0.1;/constexpr double startEnergy=0.2;/' \
         's|^//@frozen                   hitOut,nCl,discarded,failed,nullptr,|//@frozen                   hitOut,nCl,discarded,failed,k2aSubset,|' \
         's|// K2A-INLINE END exit|// K2A-INLINE END exitx|' '1035a\    k2a_extra();'; do
  i=$((i+1)); sed "$m" "$G" > "$W/mut$i.cxx"; ! cmp -s "$W/mut$i.cxx" "$G" || die "mutation $i did not apply"
  if python3 "$HERE/check_inline_generator_diff.py" --inline "$W/mut$i.cxx" > "$W/mut$i.log"; then die "diff checker accepted mutation $i"; fi
done
ok "B allowed-diff checker: inline copy PASS, 4/4 mutations FAIL"

# C
mkdir "$W/C" && cd "$W/C"
mkdir F I M
"$FRZ" --cathode-ion-reference "$W/input.csv" F/m 2100 0 full 0 F/ion.root > F/log 2>&1
"$INL" --cathode-ion-reference "$W/input.csv" I/m 2100 0 full 0 I/ion.root > I/log 2>&1
"$INL" --cathode-ion-reference "$W/input.csv" M/m 2100 0 full 0 M/ion.root --k2a-moments-out M/k > M/log 2>&1
for f in ion.root m_seed.csv m_crossing.csv m_wire_seed.csv m_wire_muon.csv m_ion_birth_sample.csv m_heed_primaries.csv m_heed_clusters.csv; do
  cmp -s F/$f I/$f || die "C: $f differs frozen vs inline (no options)"
  cmp -s F/$f M/$f || die "C: $f differs frozen vs inline (with moments)"
done
[[ -s M/k.moments.csv ]] || die "C: no moments"
cd "$W"
ok "C frozen vs inline: archive and 7 per-seed/per-muon CSVs byte-identical without and with --k2a-moments-out"

# frozen 'production' K256 archive of the same muon (all primaries), made with the frozen source
mkdir "$W/prod"; "$FRZ" --cathode-ion-reference "$W/input.csv" "$W/prod/m" 2100 0 weighted 256 "$W/prod/ion_births.root" > "$W/prod/log" 2>&1

# D
set +e
MWPC_T1_SYNTHETIC=1 MWPC_T1_INPUT="$W/input.csv" MWPC_T1_FROZEN_ARCHIVE="$W/prod/ion_births.root" MWPC_K2A_INLINE_BUILD="$BLD" \
  bash "$HERE/run_t1_inline.sh" SYN "$W/D" 3 > "$W/D.log" 2>&1
rc=$?; set -e
[[ $rc -eq 0 ]] || { cat "$W/D.log"; die "D runner verdict (exit $rc)"; }
python3 - "$W/D" <<'PY' || die "D verdict content"
import json, sys
v = json.load(open(sys.argv[1] + "/T1_verdict.json"))
assert v["label"] == "SYNTHETIC" and v["verdict"] == "PASS", v
assert all(s == "PASS" for s in v["required"].values()), v["required"]
assert v["diagnostic_frozen_production_crosscheck"] == "MATCH", v
assert v["diagnostic_T2"] == "DIAGNOSTIC", v
c = json.load(open(sys.argv[1] + "/T1_compare/compare.json"))
assert c["n_primaries"] == 3 and c["status"] == "PASS"
for f in ("N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A", "min_clearance_um", "t_raw_rms_ns"):
    b = c["bitwise_inline_vs_batch"][f]; assert b["equal"] == b["compared"] > 0, (f, b)
PY
ok "D runner (synthetic): verdict PASS, all required items PASS, frozen cross-check MATCH, T2 diagnostic; spatial moments bitwise = batch"

# E
mkdir "$W/E" && cd "$W/E"
mkdir ia op zp
set +e
K2A_STUB_INSIDE_ANODE=1 "$INL" --cathode-ion-reference "$W/input.csv" ia/m 2100 0 full 0 ia/ion.root --k2a-moments-out ia/k --max-heed-primaries 3 > ia/log 2>&1; rc=$?
set -e
[[ $rc -eq 4 && ! -e ia/k.moments.csv && -s ia/k.clouds.csv ]] && grep -q INSIDE_ANODE ia/k.inline.json || die "E inside anode (rc $rc)"
ok "E inside-anode birth: exit 4, INSIDE_ANODE, diagnostics only, no engine-facing moments file"
set +e
K2A_STUB_OUTSIDE_PATCH=1 "$INL" --cathode-ion-reference "$W/input.csv" op/m 2100 0 full 0 op/ion.root --k2a-moments-out op/k --max-heed-primaries 3 > op/log 2>&1; rc=$?
K2A_STUB_OUTSIDE_PATCH=1 "$FRZ" --cathode-ion-reference "$W/input.csv" op/f 2100 0 full 0 op/f.root > op/flog 2>&1; rcf=$?
K2A_STUB_OUTSIDE_PATCH=1 "$INL" --wiretruth-hv-p3-noFW "$W/input.csv" op/w 2100 0 --k2a-moments-out op/kw --max-heed-primaries 3 > op/wlog 2>&1; rcw=$?
set -e
[[ $rc -eq 2 && $rcf -eq 2 && $rcw -eq 2 && ! -e op/k.inline.json && ! -e op/kw.inline.json ]] || die "E outside patch (rc $rc frozen $rcf moments-only $rcw)"
grep -q "outside nine-wire patch" op/log && grep -q "outside nine-wire patch" op/flog && grep -q "Ion birth outside 9-wire analytic patch" op/wlog || die "E outside-patch messages"
ok "E outside-patch birth: frozen abort preserved in archive mode and moments-only mode (exit 2, same message), no K2A output"
K2A_STUB_ZERO_PRIMARY=1 "$INL" --cathode-ion-reference "$W/input.csv" zp/m 2100 0 full 0 zp/ion.root --k2a-moments-out zp/k --max-heed-primaries 3 > zp/log 2>&1
"$BLD/k2a_archive_dump" zp/ion.root zp/dump > /dev/null && "$BLD/p3_birth_moments" zp/ion.root zp/b.csv > /dev/null
python3 "$HERE/k2a_inline_compare.py" --inline zp/k --dump zp/dump --batch zp/b.csv --input "$W/input.csv" --max-heed-primaries 3 --out zp/cmp > zp/cmp.log \
  || { tail -20 zp/cmp.log; die "E zero-ion primary compare"; }
grep -q ZERO_ION_PRIMARY zp/k.inline.json || die "E zero-ion flag"
ok "E zero-ion primary: informational flag, moments usable, T1 comparison PASS"
D="$W/D"
cp -r "$D/dump_R1" bad_dump; sed -i '5d' bad_dump/births.csv
set +e
python3 "$HERE/k2a_inline_compare.py" --inline "$D/R1/inline" --dump bad_dump --batch "$D/batch_R1.moments.csv" --input "$W/input.csv" --max-heed-primaries 3 --out c1 > c1.log; r1=$?
cp "$D/R1/inline.moments.csv" pm.moments.csv; for s in inline.json clouds.csv primaries.csv; do cp "$D/R1/inline.$s" pm.$s; done
python3 - <<'PY'
import csv
r = list(csv.DictReader(open("pm.moments.csv"))); r[0]["D_U_mm"] = repr(float(r[0]["D_U_mm"]) * (1 + 1e-9))
w = csv.DictWriter(open("pm.moments.csv", "w", newline=""), fieldnames=list(r[0])); w.writeheader(); w.writerows(r)
PY
python3 "$HERE/k2a_inline_compare.py" --inline pm --dump "$D/dump_R1" --batch "$D/batch_R1.moments.csv" --input "$W/input.csv" --max-heed-primaries 3 --out c2 > c2.log; r2=$?
head -n 2 "$D/batch_R1.moments.csv" > tb.csv
python3 "$HERE/k2a_inline_compare.py" --inline "$D/R1/inline" --dump "$D/dump_R1" --batch tb.csv --input "$W/input.csv" --max-heed-primaries 3 --out c3 > c3.log; r3=$?
python3 "$HERE/k2a_inline_compare.py" --inline "$D/R1/inline" --dump "$D/dump_R1" --batch nothere.csv --input "$W/input.csv" --max-heed-primaries 3 --out c4 > c4.log; r4=$?
set -e
[[ $r1 -eq 1 && $r2 -eq 1 && $r3 -eq 1 && $r4 -eq 3 ]] || die "E comparator negatives ($r1 $r2 $r3 $r4)"
grep -q "CENSUS" c1/summary.md && grep -q "D_U_mm" c2/summary.md || die "E comparator failure attribution"
ok "E comparator: dropped birth -> FAIL (census), D_U x(1+1e-9) -> FAIL (moments), truncated batch -> FAIL, missing input -> INCONCLUSIVE"
cp "$D/R2/inline.moments.csv" r2m.moments.csv; for s in inline.json clouds.csv primaries.csv; do cp "$D/R2/inline.$s" r2m.$s; done
sed -i '2s/^0,/9,/' r2m.moments.csv; ! cmp -s r2m.moments.csv "$D/R2/inline.moments.csv" || die "E: rerun mutation did not apply"
cp -r "$D/dump_R4" bad_w; awk -F, 'BEGIN{OFS=","} NR==2{$3=$3"1"} {print}' "$D/dump_R4/births.csv" > bad_w/births.csv
set +e
python3 "$HERE/k2a_inline_identity.py" moments r2m "$D/R1/inline" i1.json > i1.log; r1=$?
python3 "$HERE/k2a_inline_identity.py" subset "$D/dump_R1" bad_w i2.json > i2.log; r2=$?
python3 "$HERE/k2a_inline_identity.py" subset "$D/dump_R1" bad_w i3.json --primaries 3 --label frozen-production > i3.log; r3=$?
set -e
[[ $r1 -eq 1 && $r2 -eq 1 && $r3 -eq 1 ]] && grep -q MISMATCH i3.log || die "E identity negatives ($r1 $r2 $r3)"
ok "E identity: changed rerun -> FAIL; changed reservoir row -> subset FAIL, frozen cross-check MISMATCH"
set +e
"$INL" --cathode-ion-reference "$W/input.csv" x 2100 0 full 0 x.root --max-heed-primaries 0 > o1.log 2>&1; o1=$?
"$INL" --wiretruth-p3-noFW "$W/input.csv" x2 --k2a-moments-out k > o2.log 2>&1; o2=$?
"$INL" --cathode-ion-reference "$W/input.csv" x3 2100 0 full 0 x3.root --k2a-moments-out "$D/R1/inline" > o3.log 2>&1; o3=$?
"$INL" --cathode-ion-reference "$W/input.csv" x4 2100 0 full 0 x4.root --max-heed-primaries 3x > o4.log 2>&1; o4=$?
set -e
[[ $o1 -eq 2 && $o2 -eq 2 && $o3 -eq 2 && $o4 -eq 2 ]] && grep -q "refusing to overwrite" o3.log || die "E options ($o1 $o2 $o3 $o4)"
ok "E options: N=0, non-integer N, unsupported mode and an existing moments prefix are refused"
mkdir -p pf/inventory; printf 'slot,role\n0007,sealed_validation\n0008,development\n' > pf/inventory/slots.csv
set +e
MWPC_K2A_PF="$PWD/pf" MWPC_K2A_INLINE_BUILD="$BLD" bash "$HERE/run_t1_inline.sh" 0007 "$PWD/sealed_out" 3 > s1.log 2>&1; s1=$?
MWPC_K2A_INLINE_BUILD="$BLD" bash "$HERE/run_t1_inline.sh" 0008 "$D" 3 > s2.log 2>&1; s2=$?
set -e
[[ $s1 -eq 2 && $s2 -eq 2 ]] && grep -q "sealed" s1.log && grep -q "exists" s2.log || die "E runner refusals ($s1 $s2)"
ok "E runner: sealed slot refused, existing OUTDIR refused"
echo "PASS milestone-1 synthetic test (stub Garfield/ROOT, NOT physics); work dir $W"

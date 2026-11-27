#!/usr/bin/env bash
# SYNTHETIC end-to-end test of the K2A-vs-S0 validation chain on a TOY run (toy transport, coarse toy kernel, fake
# exporter and a stand-in for the ROOT-dependent p3_birth_moments): NOT NA60+ physics, only the plumbing and the
# explicit-outcome logic.
#   bash tests/test_validation_synthetic.sh TOY_RUN TAG NUMW_DIR LUT WORK
# TOY_RUN as for k2a-campaign/tests/test_recenter_synthetic.sh (muon_0001 T1_export, S0, cloud_moments.csv, K2Akern_TAG,
# kernel_TAG, green_TAG). Nothing outside WORK is written; TOY_RUN is not modified.
#   A selector on the synthetic preflight fixture: no sealed / kernel / development / excluded slot, pilot strata named
#   B refusals before anything runs: sealed slot, kernel-reference slot, existing OUTDIR
#   C full chain on one far toy muon (half of the primaries 7 strips away; per-ion birth-time jitter with zero weighted
#     mean per primary, so T1 collapse recovers the toy export and S0_T0 differs from S0_T1 by a real time spread):
#     execution VALID; the PRIMARY K2A runs on its own production grouping (no S0 shifts), the assignment AGREEs with
#     the S0 grouping and the primary and S0-imposed responses are bit-identical; S0_T0 grouping == S0_T1 grouping;
#     compiled == Python K2A (1e-12); moment path (1e-7); six comparisons; four shaped comparisons; accuracy (v1 dry
#     run) reported separately from execution; runtime split with cathode-response evaluation separated from archive
#     I/O and moments; report COMPLETE, exit 0; archive unchanged; with MWPC_K2A_VAL_PRUNE_VERTICES=1 only the 5
#     derived path copies inside OUTDIR are removed (sha256 recorded)
#   D outcome negatives on copies (tests/validation_negatives.py): identity FAIL, domain changed by the moment path,
#     consistent refusal -> OUT_OF_DOMAIN, production rule blocked with S0 grouping usable -> OUT_OF_DOMAIN (primary
#     decides), agreeing assignments with different responses -> FAIL, assignment not compared / malformed, frozen
#     replay refusal -> S0_FAILED, archive changed, archive-vs-export moments FAIL, T0 failure -> VALID_T1_ONLY, missing
#     shaped / grouping-effect comparison -> INCONCLUSIVE, label-union failure -> FAIL; report: an intended-OOD slot is
#     listed separately and never counted in the accuracy set; missing outcome -> exit 3; combined report refuses
#     duplicates
#   F unit regression (tests/assignment_units.py) of k2a_shift_assignment.py (agree, disagree, K2A-only primary, imposed
#     run, missing) and k2a_label_union.py (identical -> symlinks, different -> bit-exact zero padding, inconsistent
#     variant refused)
#   G disagreement end-to-end (tests/disagreement_case.py) on a copy whose S0 grouping of one primary differs from the
#     production rule: reported with primary and charge, grouping effect measured, zero-padded labels, primary K2A on
#     its own grouping, equivalence tests PASS
#   H time grid (tests/grid_units.py): different dt, first-bin offset, too-short variant, one moved bin, non-uniform
#     spacing, replay != prompt axis, non-integer window, missing file -> explicit FAIL, nothing written for the shaped
#     comparison; longer S0 clock + K2A window -> one grid (S0 truncated bit-exactly); real run passes, shifted copy fails
#   E2 criteria with missing / None / NaN / non-integer lead or hit-flip metrics, missing peaks / threshold or an
#     inconsistent flip count -> FAIL, never AMBIGUOUS / MARGINAL / PASS (tests/criteria_invalid.py)
#   E criteria dry run unit cases: identical responses PASS; lead change with S0 top-2 within 1 % -> AMBIGUOUS (not a
#     failure); a real lead change and a 20 um w shift -> FAIL
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; CAMP="$(cd "$HERE/.." && pwd)"; CIR="$(cd "$CAMP/.." && pwd)"
RC="$CAMP/recenter"; TA="$CIR/timing-ablation"; SA="$CIR/spatial-ablation"; INL="$CIR/k2a-inline"
SRC="${1:?}"; TAG="${2:?}"; NUMW="${3:?}"; LUT="${4:?}"; W="${5:?}"
[[ ! -e "$W" ]] || { echo "$W exists"; exit 2; }; mkdir -p "$W"; W="$(cd "$W" && pwd)"
die(){ echo "FAIL $*" >&2; exit 1; }

# ---- A. selector
python3 "$HERE/tests/make_validation_fixture.py" "$W/fx" > /dev/null
python3 "$HERE/k2a_validation_select.py" --pf "$W/fx/pf" --recenter "$W/fx/rec" --out "$W/sel" --exclude-run "$W/fx/run" > "$W/sel.log"
python3 - "$W/sel" <<'PY'
import csv, sys
rows = [r for r in csv.DictReader(open(f"{sys.argv[1]}/selection.csv")) if r["slot"]]
sl = [r["slot"] for r in rows]
bad = [s for s in sl if int(s) < 5 or 50 <= int(s) < 60 or s in ("0033", "0000")]
assert not bad, f"forbidden slots selected: {bad}"
assert {r["role"] for r in rows} <= {"pool", "development"}, {r["role"] for r in rows}   # never previously_inspected / sealed
assert len(set(sl)) == len(sl)
pilot = [r for r in rows if r["set"] == "pilot"]; reasons = {r["reason"] for r in pilot}
assert 3 <= len(pilot) <= 5 and "typical" in reasons and "out_of_domain" in reasons, reasons
assert any(r["domain_after"] == "OUT_OF_DOMAIN" for r in pilot), "the pilot must keep an out-of-domain event"
pf = open(f"{sys.argv[1]}/pilot_slots.txt").read().split(); assert pf == [r["slot"] for r in pilot]
print(f"PASS selector (SYNTHETIC fixture): {len(rows)} slots, pilot {pf} ({', '.join(sorted(reasons))}); no sealed, kernel-reference, "
      "previously-inspected or excluded slot")
PY

# ---- reference run (symlinks to the toy run + provenance), far toy event, fake exporter / births / archive / preflight
RUN="$W/refrun"; mkdir -p "$RUN/muon_0001"
for f in kernel_$TAG green_$TAG; do ln -s "$SRC/$f" "$RUN/$f"; done
for f in T1_export S0 cloud_moments.csv K2Akern_$TAG; do ln -s "$SRC/muon_0001/$f" "$RUN/muon_0001/$f"; done
F="$W/far"; mkdir -p "$F"
python3 - "$RUN/muon_0001/T1_export" "$F" <<'PY'
import csv, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
src, F = Path(sys.argv[1]), Path(sys.argv[2])
seal = dict(l.split("=", 1) for l in (src / "complete.txt").read_text().split() if "=" in l)
V = np.fromfile(src / "vertices.f64", dtype="<f8").reshape(int(seal["vertices"]), 5)
rows = list(csv.DictReader(open(src / "ions.csv")))
prim = sorted({int(r["primary"]) for r in rows}); far = set(prim[::2])
for r in rows:                                     # half of the primaries 7 strips along w (w-invariant toy field)
    if int(r["primary"]) in far:
        s, n = int(r["first_vertex"]), int(r["vertices"]); V[s:s + n, 3] += 7 * 1.7; V[s:s + n, 4] += 7 * 1.7
tot = defaultdict(float)                           # integer census per primary (as real archives)
for r in rows:
    tot[int(r["primary"])] += float(r["weight"])
target = {p: max(1, round(t)) for p, t in tot.items()}
target[prim[0]] += int(round(sum(tot.values()))) - sum(target.values())
for r in rows:
    r["weight"] = repr(float(r["weight"]) * target[int(r["primary"])] / tot[int(r["primary"])])
ws = sum(float(r["weight"]) for r in rows)
# "frozen" export with true per-ion birth times: zero weighted-mean jitter per primary (T1 collapse recovers the toy times)
rng = np.random.default_rng(7)
jit = np.zeros(len(rows)); idx = defaultdict(list)
for j, r in enumerate(rows):
    idx[int(r["primary"])].append(j)
for p, js in idx.items():
    tb = min(float(rows[j]["birth_ns"]) for j in js)
    x = rng.uniform(0, min(40.0, 0.9 * tb) if tb > 0 else 0.0, len(js)); w = np.array([float(rows[j]["weight"]) for j in js])
    jit[js] = x - np.sum(w * x) / np.sum(w)
d = F / "frozen_like_export"; d.mkdir()
for j, r in enumerate(rows):
    s, n = int(r["first_vertex"]), int(r["vertices"]); V[s:s + n, 0] += jit[j]
    r["birth_ns"] = repr(float(r["birth_ns"]) + float(jit[j])); r["end_ns"] = repr(float(r["end_ns"]) + float(jit[j]))
with open(d / "ions.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
V.tofile(d / "vertices.f64")
(d / "complete.txt").write_text((src / "complete.txt").read_text().replace(f"weight_sum={seal['weight_sum']}", f"weight_sum={ws!r}"))
print(f"toy far event: {len(far)} of {len(prim)} primaries +7 strips; weighted rms birth jitter {np.sqrt(np.mean(jit ** 2)):.2f} ns")
PY
python3 "$TA/t_time_collapse_export.py" "$F/frozen_like_export" "$F/T1_ref" --mode primary > /dev/null
python3 "$SA/s_kernel.py" moments "$F/T1_ref" "$F/birth_moments_src.csv" > /dev/null
FX="$W/fake_exporter.sh"
cat > "$FX" <<FAKE
#!/usr/bin/env bash
[[ \$# -eq 3 ]] || { echo "usage: PROGRAM frozen.root MOBILITY NEW_DIR"; exit 2; }
cp -r "$F/frozen_like_export" "\$3" && echo "PASS archive census (fake exporter for the synthetic test)"
FAKE
FB="$W/fake_births.sh"
cat > "$FB" <<FAKE
#!/usr/bin/env bash
# SYNTHETIC stand-in for p3_birth_moments (ROOT archive reader): moments of the same toy event
[[ \$# -ge 2 ]] || { echo "usage: p3_birth_moments ion_births.root OUT_MOMENTS.csv [--repeat N]"; exit 2; }
[[ ! -e "\$2" ]] || { echo "output exists"; exit 2; }
cp "$F/birth_moments_src.csv" "\$2"
printf '{"tool": "p3_birth_moments SYNTHETIC stand-in", "archive": "%s", "archive_bytes": 9, "rows": 960, "clouds": 40, "repeat": 5, "read_validate_s_min": 0.004, "moments_s_min": 0.003, "peak_rss_kib": 9000, "moments_bytes": 1}\n' "\$1" > "\$2.json"
echo "PASS birth moments rows=960 clouds=40 (SYNTHETIC stand-in)"
FAKE
chmod +x "$FX" "$FB"; echo "toy mobility" > "$W/mobility.txt"
{ echo "lut_sha256=$(sha256sum "$LUT" | cut -d' ' -f1)"; echo "replay_sha256=$(sha256sum "$NUMW/p3_replay_full_spectral_muon.py" | cut -d' ' -f1)"
  echo "export_bin_sha256=$(sha256sum "$FX" | cut -d' ' -f1)"; echo "mobility_sha256=$(sha256sum "$W/mobility.txt" | cut -d' ' -f1)"; } > "$RUN/provenance_prepare_0001.txt"
mkdir -p "$W/pf/inventory"
for s in 0001 0002 0003; do mkdir -p "$W/archive/muon_$s"; echo "synthetic $s" > "$W/archive/muon_$s/ion_births.root"; (cd "$W/archive/muon_$s" && sha256sum ion_births.root > SHA256SUMS); done
printf 'slot,role,status\n0001,dev_reference,USABLE\n0002,pool,USABLE\n0003,sealed_validation,USABLE\n' > "$W/pf/inventory/slots.csv"
mkdir -p "$W/build"; g++ -O3 -std=c++17 -o "$W/build/k2a_strip_eval" "$INL/k2a_strip_eval.cc" "$SA/engine/kernel_eval.cc"
KID=$(python3 -c "import sys; sys.path.insert(0, '$SA/engine'); import cathode_engine as ce; print(ce.KernelTable('$RUN/kernel_$TAG').kernel_id)")
printf 'slot,reason,station,channel,abs_theta_u,abs_theta_w,log_gain,mean_abs_A,n_wires_5pct,max_abs_m,clouds\n0002,displaced,3,phi,4.1,2.2,9.3,0.2,1,7,40\n' > "$W/selection.csv"
ENVV=(MWPC_SPECTRAL_LUT="$LUT" MWPC_NUMW_DIR="$NUMW" MWPC_EXPORT_BIN="$FX" MWPC_ION_MOBILITY="$W/mobility.txt" MWPC_K2A_PF="$W/pf"
      MWPC_ARCHIVE_2100_K256="$W/archive" MWPC_K2A_KERNEL_ID="$KID" MWPC_BIRTHS_BIN="$FB" MWPC_K2A_INLINE_BUILD="$W/build"
      MWPC_K2A_FAR_MIN_FREE_GIB=0 MWPC_K2A_VAL_MIN_FREE_GIB=0 MWPC_K2A_VAL_SYNTHETIC=1 MWPC_K2A_VAL_REPEAT=5 MWPC_K2A_VAL_SELECTION="$W/selection.csv")
VR="$HERE/run_k2a_validation.sh"

# ---- B. refusals before anything runs
if env "${ENVV[@]}" bash "$VR" "$RUN" "$TAG" "$W/o_sealed" 0002 0003 > "$W/o_sealed.log" 2>&1; then die "sealed slot accepted"; fi
grep -q "is SEALED: nothing was run" "$W/o_sealed.log" && [[ ! -d "$W/o_sealed/muons" ]] || die "sealed refusal: $(tail -2 "$W/o_sealed.log")"
if env "${ENVV[@]}" bash "$VR" "$RUN" "$TAG" "$W/o_dev" 0001 > "$W/o_dev.log" 2>&1; then die "kernel-reference slot accepted"; fi
grep -q "is not independent" "$W/o_dev.log" && [[ ! -e "$W/o_dev" ]] || die "dev refusal: $(tail -2 "$W/o_dev.log")"
mkdir "$W/o_exists"
if env "${ENVV[@]}" bash "$VR" "$RUN" "$TAG" "$W/o_exists" 0002 > "$W/o_dev2.log" 2>&1; then die "existing OUTDIR accepted"; fi
grep -q "exists" "$W/o_dev2.log" || die "existing OUTDIR: $(tail -2 "$W/o_dev2.log")"
echo "PASS refusals: a sealed slot aborts the whole run before any muon; a kernel-reference slot is refused as not independent; an existing OUTDIR is refused"

# ---- C. full chain
rc=0; env "${ENVV[@]}" MWPC_K2A_VAL_PRUNE_VERTICES=1 bash "$VR" "$RUN" "$TAG" "$W/val" 0002 > "$W/val.log" 2>&1 || rc=$?
[[ $rc -eq 0 ]] || { tail -30 "$W/val.log"; die "validation run exit $rc"; }
python3 - "$W/val" <<'PY'
import json, sys
from pathlib import Path
O = Path(sys.argv[1]); o = json.load(open(O / "muons/muon_0002/outcome.json")); r = json.load(open(O / "validation_report.json"))
assert o["label"] == "SYNTHETIC" and r["label"] == "SYNTHETIC"
assert o["outcome"] == "VALID", (o["outcome"], o["failures"], o["inconclusive"])
v = o["steps"]["validation"]
for k in ("T0_groups", "T0_S0_local_census", "births_moments", "K2A_prod", "K2A_imposed", "K2A_cc_export", "label_union"):
    assert v[k]["status"] == "OK", (k, v[k])
for k in ("identity_cc_vs_python", "moment_path_archive_vs_export", "moments_archive_vs_export", "grouping_identity"):
    assert v[k]["status"] == "PASS", (k, v[k])
# the PRIMARY K2A ran WITHOUT the S0 grouping; agreement of assignments is checked, and then responses are bit-identical
rp = json.load(open(O / "muons/muon_0002/K2A_prod/result.json")); ri = json.load(open(O / "muons/muon_0002/K2A_imposed/result.json"))
assert rp["shifts_from"].startswith("production rule") and ri["shifts_from"].endswith("groups.json"), (rp["shifts_from"], ri["shifts_from"])
assert v["shift_assignment"]["status"] == "AGREE" and o["shift_assignment"]["n_disagree"] == 0
assert v["time_grid"]["status"] == "PASS", v["time_grid"]
assert o["k2a_domain"]["grouping"].startswith("production rule")
c = o["comparisons"]
for n in ("S0T1_vs_K2A", "S0T0_vs_K2A", "S0T1_vs_K2Aimp", "S0T0_vs_K2Aimp", "S0T0_vs_S0T1", "K2Aimp_vs_K2A"):
    assert c[n] and c[n]["status"] in ("OK", "WARN"), (n, c[n] and c[n]["status"])
assert c["K2Aimp_vs_K2A"]["worst_L2_rel"] == 0.0
t = c["S0T0_vs_S0T1"]["ion_only.U.unshaped_L2_rel"]
assert 1e-6 < t < 0.5, f"time jitter must change the waveform measurably: {t}"
exp = {"refS0T1_ionelec": {"K2A", "K2Aimp", "S0T0"}, "refS0T0_ionelec": {"K2A", "K2Aimp", "S0T1"},
       "refS0T1_ion": {"K2A", "K2Aimp", "S0T0"}, "refK2Aimp_ionelec": {"K2A"}}
for n, e in exp.items():
    assert set(o["shaped"][n]["metrics"]) == e, (n, set(o["shaped"][n]["metrics"]))
d = o["criteria_v1_dry_run"]["primary_vs_S0T1"]
assert d["decisive"] and d["candidate"] == "K2A" and d["n_checks"] == 2 * (2 * 7 + 3), d["n_checks"]
assert set(o["criteria_v1_dry_run"]) == {"primary_vs_S0T1", "primary_vs_S0T0", "imposed_vs_S0T1"}
assert o["accuracy"]["status"] == d["status"] and d["status"] in ("PASS", "FAIL")
assert "APPROXIMATE" in o["references"]["S0_T0"] and "DIAGNOSTIC" in o["references"]["K2Aimp"]
rt = o["runtime"]
assert "excluded" in rt["scope"]
for k in ("total_T0_single_pass", "response_T0_excl_archive_read_approx", "rkf_export_archive_read_plus_ion_drift"):
    assert rt["S0"][k] and rt["S0"][k] > 0, k
for k in ("response_eval_only", "moments_plus_response_eval", "archive_io_plus_moments_plus_eval", "single_event_cold_incl_kernel_open"):
    assert rt["K2A"][k] and rt["K2A"][k] > 0, k
sp = rt["speedup_measured"]
assert sp["S0_T0/cathode_response_eval"] > sp["S0_T0/response_from_births_incl_moments"] > sp["S0_T0/end_to_end_from_archive_incl_io"] > 1
assert r["status"] == "COMPLETE" and r["execution_counts"] == {"VALID": 1}, r["execution_counts"]
assert r["accuracy_set"] == ["0002"] and r["accuracy_v1_primary_counts"] == {d["status"]: 1}
assert (O / "per_muon.csv").read_text().count("\n") == 2
md = (O / "validation_report.md").read_text()
for x in ("SYNTHETIC test data", "execution is not accuracy", "DIAGNOSTIC dry run", "Avalanche generation is excluded", "Shift assignment"):
    assert x in md, x
b, a = (O / "muons/muon_0002/archive_sha256_before.txt").read_text(), (O / "muons/muon_0002/archive_sha256_after.txt").read_text()
assert a == b
pr = (O / "muons/muon_0002/pruned_vertices.sha256").read_text().splitlines()
assert len(pr) == 5 and all(l.split()[1].startswith(str(O)) for l in pr), pr      # frozen export + 2 T1 groups + 2 T0 groups
assert not list(O.rglob("frozen_export/vertices.f64")) and (Path(o["far_dir"]) / "prepared/T1_export/vertices.f64").is_file()
print(f"PASS full chain (SYNTHETIC toy): execution VALID, accuracy (v1 dry run, primary K2A vs S0_T1) {d['status']} ({d['n_checks']} checks); "
      f"primary K2A ran on its own grouping, assignments AGREE and K2A == K2Aimp bit-identically; compiled == Python, moment path; "
      f"model worst L2 {c['S0T1_vs_K2A']['worst_L2_rel']:.3g}, time collapse U L2 {t:.3g} (coarse toy kernel); "
      f"measured toy cathode-response-eval speedup {sp['S0_T0/cathode_response_eval']:.3g} (toy sizes: not a physics number)")
PY

# ---- D. outcome negatives on copies
python3 "$HERE/tests/validation_negatives.py" "$W/val" "$HERE" "$W/neg" "$CAMP/sealed/sealed_acceptance_criteria_v1.json"

# ---- F. shift assignment and label union: unit regression (agreement, disagreement, malformed, missing)
python3 "$HERE/tests/assignment_units.py" "$HERE" "$W/unit"

# ---- G. disagreement end-to-end: the S0 grouping of one primary differs from the production rule (copies; the S0
#      grouping file is edited and the Python K2A_local regenerated with it, so the equivalence tests still compare like
#      with like; the S0 waveforms are those of the original grouping: plumbing test of the disagreement path only)
python3 "$HERE/tests/disagreement_case.py" prepare "$W/val" "$W/disagree" "$RC/k2a_recentered_k2a.py" "$INL/tests/make_toy_kernel.py" "$SA/engine"
env "${ENVV[@]}" K2AV_EV="$W/build/k2a_strip_eval" K2AV_REPEAT=3 bash "$HERE/k2a_validation_k2a.sh" "$W/disagree/muon_0002" "$W/disagree/E" \
    "$(cat "$W/disagree/kernel_dir.txt")" "$(cat "$W/disagree/kernel_id.txt")" "$W/archive/muon_0002/ion_births.root" > "$W/disagree.log" 2>&1 || { tail -20 "$W/disagree.log"; die "K2A stage (disagreement)"; }
python3 "$HERE/k2a_validation_muon.py" "$W/disagree/muon_0002" --muon 0002 --label SYNTHETIC --independent yes --t0 1 \
    --criteria "$CAMP/sealed/sealed_acceptance_criteria_v1.json" > "$W/disagree_outcome.log" 2>&1
python3 "$HERE/tests/disagreement_case.py" check "$W/disagree"

# ---- E. criteria dry-run unit cases
python3 - "$HERE" "$CAMP/sealed/sealed_acceptance_criteria_v1.json" <<'PY'
import copy, json, sys
sys.path.insert(0, sys.argv[1]); import k2a_validation_muon as vm
crit = json.load(open(sys.argv[2]))
def lines(mod):
    out = {}
    for tp in ("100", "200"):
        ref = {"variant": "S0T1", "U_peaks": [10.0, 5.0, 1.0], "V_peaks": [8.0, 4.0, 0.5],
               "thresholds_fCeq": {s: {f"f{f:g}": f * p for f in (0.02, 0.05, 0.1, 0.2, 0.5)} for s, p in (("U", 10.0), ("V", 8.0))}}
        var = {"variant": "K2A", "U_peaks": [10.0, 5.0, 1.0], "V_peaks": [8.0, 4.0, 0.5], "w_mm": 0.0, "u_mm": 0.0, "UV_charge_ratio_rel": 0.0}
        for s in ("U", "V"):
            var.update({f"{s}_lead_index_changed": 0, f"{s}_hitflips_f0.05": 0, f"{s}_hitflips_f0.1": 0, f"{s}_t_f0.1_ns": 0.0,
                        f"{s}_charge_rel": 0.0, f"{s}_second_over_lead": 0.0, f"{s}_L2_rel": 0.0})
        mod(ref, var); out[tp] = {"S0T1": ref, "K2A": var}
    return out
r = vm.criteria_dry_run(crit, lines(lambda r, v: None), "S0T1", "K2A"); assert r["status"] == "PASS" and r["n_checks"] == 34, r
def amb(r, v):
    r["U_peaks"] = [10.0, 9.95, 1.0]; v["U_lead_index_changed"] = 1
r = vm.criteria_dry_run(crit, lines(amb), "S0T1", "K2A"); assert r["status"] == "PASS" and any("AMBIGUOUS" in e for e in r["exemptions"]), r
def bad(r, v):
    v["V_lead_index_changed"] = 1; v["w_mm"] = 0.020
r = vm.criteria_dry_run(crit, lines(bad), "S0T1", "K2A")
assert r["status"] == "FAIL" and any(f.startswith("V_lead_index_changed") for f in r["failed"]) and any(f.startswith("w_mm") for f in r["failed"]), r["failed"]
def marg(r, v):
    v["U_hitflips_f0.1"] = 1; v["U_peaks"] = [10.0, 5.0, 1.01]; r["U_peaks"] = [10.0, 5.0, 0.995]
r = vm.criteria_dry_run(crit, lines(marg), "S0T1", "K2A"); assert r["status"] == "PASS" and any("MARGINAL" in e for e in r["exemptions"]), r
r = vm.criteria_dry_run(crit, {"100": lines(lambda r, v: None)["100"]}, "S0T1", "K2A"); assert r["status"] == "MISSING"
print("PASS criteria v1 dry run: identical -> PASS (34 checks); ambiguous lead and marginal hit flip are exemptions; real lead "
      "change and 20 um w shift -> FAIL; a missing Tp -> MISSING")
PY
# ---- E2. criteria dry run: missing / invalid lead and hit-flip metrics are failures, never exemptions or passes
python3 "$HERE/tests/criteria_invalid.py" "$HERE" "$CAMP/sealed/sealed_acceptance_criteria_v1.json"

# ---- H. time grid: explicit FAIL on any cross-variant mismatch; the real synthetic variants share one grid
python3 "$HERE/tests/grid_units.py" "$HERE" "$W/grid" "$W/val"
echo "PASS K2A validation synthetic test; work dir $W"

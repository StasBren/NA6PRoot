#!/usr/bin/env bash
# Milestone 2 validation of the compiled absolute-strip K2A adapter on REAL moments with the REAL frozen kernel.
#
#   bash run_m2_adapter.sh OUTDIR
#
# env (all explicit, nothing guessed; every input is read-only):
#   MWPC_K2A_INLINE_BUILD     build directory containing k2a_strip_eval (default <this dir>/build)
#   MWPC_K2A_KERNEL_DIR       the validated kernel directory (RUN/kernel_TAG); MWPC_K2A_KERNEL_ID its validated id
#   MWPC_M2_E2B_DIRS          E2b muon directories of k2a_recenter_experiment.sh far mode; muon_0560 AND muon_0997 are
#                             REQUIRED (prepared/cloud_moments.csv, groups/groups.json, K2A_local/ = offline response)
#   MWPC_M2_T1_RUN            REQUIRED: the accepted run_t1_inline.sh directory of slot 0000 (T1_verdict.json PASS)
#   MWPC_M2_REGRESSION_RUNS   REQUIRED, >= 1: further run_t1_inline.sh directories (other non-sealed slots) whose T1
#                             verdict is PASS; at least one must show >= 2 wires (each >= 1 % of the charge) or a
#                             primary with m_U or m_V != 0 in its full-birth moments
#   NA60_KERNEL_EVAL_LIB      optional (Python 'compiled' reference); MWPC_M2_REPEAT (default 20)
#   MWPC_M2_SYNTHETIC=1       TEST MODE: accept T1 verdicts labelled SYNTHETIC; the verdict is labelled SYNTHETIC
#
# Required steps per E2b item (e2b_muon_0560, e2b_muon_0997):
#   evaluate, reproduce_offline (C++ with groups.json and the offline model / half strips / window == archived
#   K2A_local), reproduce_loop (== fresh k2a_recentered_k2a.py --impl loop), production_rule (INFO: rule shifts vs the
#   E2b grouping; identical files required when no shift differs), displace_{1_0,-1_0,0_1,0_-1,2_-2,-2_2} (displaced
#   event == response on shifted labels), bench_verified, bench_mapped, bench_preload, inputs_unchanged
# Required steps per full-birth item (t1_0000 and every regression_<slot>):
#   t1_validated (T1 verdict PASS, label REAL), slot (0000 for the T1 item; another slot for a regression item),
#   evaluate (C++ and Python both produce a response; both refusing = outside the kernel domain = INCONCLUSIVE, one
#   refusing alone = FAIL), full_birth_vs_python_vector, full_birth_vs_python_loop, bench_verified, bench_mapped,
#   bench_preload, inputs_unchanged; regression items also coverage (multiwire or displaced, else INCONCLUSIVE)
# Benchmarks (run only after a usable evaluation): each mode must exit 0 with a usable result and bit-identical repeats; 'verified' must have checked the
# kernel id, 'mapped' must have mapped K.npy, 'preload' must have copied it.
# Verdict (m2_verdict.json, m2_summary.md): FAIL if any step FAILs; INCONCLUSIVE if any required item or step is missing
# or not PASS, or no regression item fulfils the coverage; PASS only with the complete required coverage.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CIR="$(cd "$HERE/.." && pwd)"
REC="$CIR/k2a-campaign/recenter/k2a_recentered_k2a.py"; CMP="$HERE/k2a_strip_compare.py"
fail(){ echo "ERROR: $*" >&2; exit 2; }
OUT="${1:?usage: run_m2_adapter.sh OUTDIR}"; [[ ! -e "$OUT" ]] || fail "$OUT exists"
B="${MWPC_K2A_INLINE_BUILD:-$HERE/build}"; EV="$B/k2a_strip_eval"; [[ -x "$EV" ]] || fail "missing $EV"
KD="${MWPC_K2A_KERNEL_DIR:?set MWPC_K2A_KERNEL_DIR}"; KID="${MWPC_K2A_KERNEL_ID:?set MWPC_K2A_KERNEL_ID}"
[[ -s "$KD/kernel.json" && -s "$KD/K.npy" ]] || fail "not a kernel directory: $KD"
REP="${MWPC_M2_REPEAT:-20}"; SYN="${MWPC_M2_SYNTHETIC:-0}"
LABEL=REAL; [[ "$SYN" == 1 ]] && LABEL=SYNTHETIC
case "$(cd "$(dirname "$OUT")" 2>/dev/null && pwd)/" in
  "$HOME/na6p/wt-validation/"*|*test_runs/mwpc_stage_c/s_ladder_2100/*) fail "OUTDIR inside a frozen reference tree";;
esac
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
STEPS="$OUT/steps.tsv"; printf 'item\tstep\tstatus\tdetail\n' > "$STEPS"
step(){ printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$STEPS"; echo "[$3] $1 $2 ${4:-}"; }
{ echo "label $LABEL"; echo "date_utc $(date -u +%FT%TZ)"; echo "git_commit $(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "binary_sha256 $(sha256sum "$EV")"; echo "kernel_dir $KD"; echo "kernel_id_expected $KID"
  echo "e2b_dirs ${MWPC_M2_E2B_DIRS:-}"; echo "t1_run ${MWPC_M2_T1_RUN:-}"; echo "regression_runs ${MWPC_M2_REGRESSION_RUNS:-}"
  if [[ -s "$B/k2a_inline_build_info.txt" ]]; then sed 's/^/build_info /' "$B/k2a_inline_build_info.txt"; fi; } > "$OUT/provenance.txt"
# declared required coverage (the verdict checks it; anything absent is INCONCLUSIVE)
REQ="$OUT/required.tsv"; printf 'item\tkind\n' > "$REQ"
printf 'e2b_muon_0560\te2b\ne2b_muon_0997\te2b\nt1_0000\tt1\n' >> "$REQ"
hashes(){ for d in "$@"; do find -L "$d" -type f \( -name '*.csv' -o -name '*.json' \) 2>/dev/null | sort | xargs -r sha256sum; done; }
cmpstep(){ # item step ref test [compare args...]
  local it="$1" st="$2" r="$3" t="$4"; shift 4
  set +e; python3 "$CMP" "$r" "$t" "$OUT/$it/$st.json" --name "$it $st" "$@" > "$OUT/$it/$st.log" 2>&1; local rc=$?; set -e
  case $rc in 0) step "$it" "$st" PASS "$(tail -1 "$OUT/$it/$st.log")";; 3) step "$it" "$st" INCONCLUSIVE "$(tail -1 "$OUT/$it/$st.log")";;
              *) step "$it" "$st" FAIL "$(tail -1 "$OUT/$it/$st.log")";; esac
}
bench(){ # item moments [extra evaluator args]; three modes, each checked
  local it="$1" m="$2"; shift 2
  local mode rc
  for mode in verified mapped preload; do
    local args=(--no-verify); [[ $mode == verified ]] && args=(--kernel-id "$KID"); [[ $mode == preload ]] && args+=(--preload)
    set +e; "$EV" --kernel "$KD" "${args[@]}" --moments "$m" --out "$OUT/$it/bench_$mode" --repeat "$REP" "$@" > "$OUT/$it/bench_$mode.log" 2>&1; rc=$?; set -e
    if python3 - "$OUT/$it/bench_$mode/result.json" "$mode" "$KID" "$rc" > "$OUT/$it/bench_$mode.check" 2>&1 <<'PY'
import json, sys
p, mode, kid, rc = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
r = json.load(open(p))
bad = []
if rc != 0: bad.append(f"exit {rc}")
if not r["usable"]: bad.append("not usable: " + r["status_names"])
if not r["timing"]["repeats_bit_identical"]: bad.append("repeats differ")
if mode == "verified" and not (r["kernel_verified"] and r["kernel_id"] == kid): bad.append("kernel id not verified")
if mode == "mapped" and not r["memory"]["kernel_K_mapped"]: bad.append("K.npy not mapped")
if mode == "preload" and r["memory"]["kernel_K_mapped"]: bad.append("K.npy not preloaded")
print("; ".join(bad) if bad else f"cold {1e3 * r['timing']['evaluate_cold_s']:.3g} ms, warm {1e3 * (r['timing']['evaluate_warm_median_s'] or 0):.3g} ms, "
      f"kernel verify {r['timing']['kernel_verify_s']:.3g} s, peak RSS {r['memory']['peak_rss_kib_end']} KiB")
sys.exit(1 if bad else 0)
PY
    then step "$it" "bench_$mode" PASS "$(cat "$OUT/$it/bench_$mode.check")"
    else step "$it" "bench_$mode" FAIL "$(tail -1 "$OUT/$it/bench_$mode.check") (exit $rc: $(tail -1 "$OUT/$it/bench_$mode.log"))"; fi
  done
}

# ---------------------------------------------------------------- E2b muons: reproduction of the offline response
for E in ${MWPC_M2_E2B_DIRS:-}; do
  it="e2b_$(basename "$E")"; mkdir "$OUT/$it"
  M="$E/prepared/cloud_moments.csv"; G="$E/groups/groups.json"; OFF="$E/K2A_local"
  if [[ ! -s "$M" || ! -s "$G" || ! -s "$OFF/replay/waveforms.csv" ]]; then step "$it" inputs INCONCLUSIVE "missing $M / $G / $OFF"; continue; fi
  hashes "$M" "$G" "$OFF" > "$OUT/$it/inputs_sha256_before.txt"
  # evaluator settings of the archived offline run (recentered.json): model, half strips, window, kernel id
  read -r MODEL HSx WIN OKID < <(python3 -c "
import json, sys
j = json.load(open(sys.argv[1])); print(j.get('model', 'K2A'), j.get('half_strips', 5), j.get('window_ns', 7000.0), j.get('kernel_id') or '-')
" "$OFF/recentered.json" 2>/dev/null || echo "K2A 5 7000 -")
  if [[ "$OKID" != "-" && "$OKID" != "$KID" ]]; then step "$it" kernel_id FAIL "offline run used kernel $OKID"; continue; fi
  X=(--model "$MODEL" --half-strips "$HSx" --window-ns "$WIN")
  set +e; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --shifts "$G" --out "$OUT/$it/cc_e2b_groups" "${X[@]}" > "$OUT/$it/cc_e2b_groups.log" 2>&1; rc=$?; set -e
  if [[ $rc -ne 0 ]]; then step "$it" evaluate FAIL "exit $rc: $(tail -1 "$OUT/$it/cc_e2b_groups.log")"; continue; fi
  step "$it" evaluate PASS "$(tail -1 "$OUT/$it/cc_e2b_groups.log")"
  cmpstep "$it" reproduce_offline "$OFF" "$OUT/$it/cc_e2b_groups"
  if python3 "$REC" "$KD" "$M" "$OUT/$it/py_loop" --groups "$G" --impl loop --kernel-id "$KID" --model "$MODEL" --half-strips "$HSx" \
       --window-ns "$WIN" > "$OUT/$it/py_loop.log" 2>&1; then cmpstep "$it" reproduce_loop "$OUT/$it/py_loop" "$OUT/$it/cc_e2b_groups"
  else step "$it" reproduce_loop FAIL "python loop failed: $(tail -1 "$OUT/$it/py_loop.log")"; fi
  "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "$OUT/$it/cc_rule" "${X[@]}" > "$OUT/$it/cc_rule.log" 2>&1 || true
  if python3 - "$OUT/$it" "$G" <<'PY'
import csv, filecmp, json, sys
d, g = sys.argv[1], json.load(open(sys.argv[2]))["primary_shift"]
try:
    s = {r["primary"]: (int(r["m_U"]), int(r["m_V"])) for r in csv.DictReader(open(f"{d}/cc_rule/shifts.csv"))}
except FileNotFoundError:
    open(f"{d}/rule_vs_e2b.txt", "w").write("production-rule evaluation produced no shifts"); sys.exit(1)
diff = sorted(p for p in s if tuple(g[p]) != s[p])
msg = f"{len(diff)} of {len(s)} primaries: rule shift (cloud moments) != E2b export shift" + (f" {diff[:10]}" if diff else "")
same = not diff and all(filecmp.cmp(f"{d}/cc_rule/{f}", f"{d}/cc_e2b_groups/{f}", shallow=False) for f in ("replay/waveforms.csv", "prompt/prompt_waveforms.csv"))
if not diff and not same:
    msg += "; identical shifts but different files"; open(f"{d}/rule_vs_e2b.txt", "w").write(msg); sys.exit(1)
open(f"{d}/rule_vs_e2b.txt", "w").write(msg + ("; files identical" if same else ""))
m = sorted({tuple(v) for v in g.values()})
open(f"{d}/shift_range.txt", "w").write(f"m_U {min(x[0] for x in m)}..{max(x[0] for x in m)}, m_V {min(x[1] for x in m)}..{max(x[1] for x in m)}, {len(m)} groups")
PY
  then step "$it" production_rule INFO "$(cat "$OUT/$it/rule_vs_e2b.txt")"
  else step "$it" production_rule FAIL "$(cat "$OUT/$it/rule_vs_e2b.txt" 2>/dev/null)"; fi
  step "$it" shift_range INFO "$(cat "$OUT/$it/shift_range.txt" 2>/dev/null)"
  for d in "1 0" "-1 0" "0 1" "0 -1" "2 -2" "-2 2"; do
    set -- $d; nm="displace_${1}_${2}"
    if "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "$OUT/$it/cc_$nm" --displace-event "$1" "$2" "${X[@]}" > "$OUT/$it/cc_$nm.log" 2>&1; then
      cmpstep "$it" "$nm" "$OUT/$it/cc_rule" "$OUT/$it/cc_$nm" --label-shift "$1" "$2"
    else step "$it" "$nm" FAIL "$(tail -1 "$OUT/$it/cc_$nm.log")"; fi
  done
  bench "$it" "$M" --shifts "$G" "${X[@]}"
  hashes "$M" "$G" "$OFF" > "$OUT/$it/inputs_sha256_after.txt"
  if cmp -s "$OUT/$it/inputs_sha256_before.txt" "$OUT/$it/inputs_sha256_after.txt"; then step "$it" inputs_unchanged PASS; else step "$it" inputs_unchanged FAIL; fi
done

# ---------------------------------------------------------------- full-birth runs: accepted T1 0000 and regressions
full_birth(){ # run_dir kind(t1|regression)
  local R="$1" kind="$2" slot it M rcc rpv rpl
  slot=$(awk '$1=="slot"{print $2}' "$R/provenance.txt" 2>/dev/null | head -1); slot="${slot:-unknown_$(basename "$R")}"
  if [[ $kind == t1 ]]; then it="t1_$slot"; else it="regression_$slot"; printf '%s\tregression\n' "$it" >> "$REQ"; fi
  if [[ -e "$OUT/$it" ]]; then step "$it" duplicate FAIL "item given twice: $R"; return; fi
  mkdir "$OUT/$it"; echo "$R" > "$OUT/$it/source_run.txt"
  if [[ $kind == t1 && "$slot" != 0000 ]]; then step "$it" slot FAIL "MWPC_M2_T1_RUN is slot $slot, not the accepted 0000"
  elif [[ $kind == regression && "$slot" == 0000 ]]; then step "$it" slot FAIL "a regression run must be another slot than 0000"
  elif [[ "$slot" == unknown_* ]]; then step "$it" slot INCONCLUSIVE "no 'slot' line in $R/provenance.txt"
  else step "$it" slot PASS "$slot"; fi
  if [[ ! -s "$R/T1_verdict.json" ]]; then step "$it" t1_validated INCONCLUSIVE "no T1_verdict.json in $R"
  elif python3 - "$R/T1_verdict.json" "$LABEL" <<'PY'
import json, sys
v = json.load(open(sys.argv[1])); sys.exit(0 if v.get("verdict") == "PASS" and v.get("label") == sys.argv[2] else 1)
PY
  then step "$it" t1_validated PASS "T1 verdict PASS ($LABEL)"
  else step "$it" t1_validated FAIL "T1 verdict is not PASS/$LABEL: $(python3 -c "import json,sys; v=json.load(open(sys.argv[1])); print(v.get('verdict'), v.get('label'))" "$R/T1_verdict.json")"; fi
  M="$R/R1/inline.moments.csv"
  if [[ ! -s "$M" ]]; then step "$it" inputs INCONCLUSIVE "missing $M (blocked or absent full-birth moments)"; return; fi
  hashes "$M" > "$OUT/$it/inputs_sha256_before.txt"
  set +e; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "$OUT/$it/cc" > "$OUT/$it/cc.log" 2>&1; rcc=$?
  python3 "$REC" "$KD" "$M" "$OUT/$it/py_vector" --impl vector --kernel-id "$KID" > "$OUT/$it/py_vector.log" 2>&1; rpv=$?
  python3 "$REC" "$KD" "$M" "$OUT/$it/py_loop" --impl loop --kernel-id "$KID" > "$OUT/$it/py_loop.log" 2>&1; rpl=$?
  set -e
  if [[ $rcc -eq 0 && $rpv -eq 0 && $rpl -eq 0 ]]; then
    step "$it" evaluate PASS "$(tail -1 "$OUT/$it/cc.log")"
    cmpstep "$it" full_birth_vs_python_vector "$OUT/$it/py_vector" "$OUT/$it/cc"
    cmpstep "$it" full_birth_vs_python_loop "$OUT/$it/py_loop" "$OUT/$it/cc"
  elif [[ $rcc -eq 4 && $rpv -ne 0 && $rpl -ne 0 ]]; then
    step "$it" evaluate INCONCLUSIVE "both refuse (consistent, but no waveform comparison): C++ $(grep -o '"status_names": "[^"]*"' "$OUT/$it/cc/result.json"); Python: $(tail -1 "$OUT/$it/py_vector.log")"
  else
    step "$it" evaluate FAIL "inconsistent: C++ exit $rcc, Python vector $rpv, loop $rpl"
  fi
  if [[ $kind == regression && $rcc -eq 0 ]]; then
    if python3 - "$OUT/$it/cc" > "$OUT/$it/coverage.txt" <<'PY'
import csv, sys
d = sys.argv[1]
c = list(csv.DictReader(open(f"{d}/clouds.csv"))); s = list(csv.DictReader(open(f"{d}/shifts.csv")))
Q = sum(float(r["N_fC"]) for r in c); wq = {}
for r in c: wq[r["wire"]] = wq.get(r["wire"], 0.0) + float(r["N_fC"])
wires = sorted((w for w, q in wq.items() if q >= 0.01 * Q), key=int)
disp = [r["primary"] for r in s if int(r["m_U"]) or int(r["m_V"])]
print(f"wires >= 1 % of charge: {wires}; displaced primaries: {disp[:20]}; {len(c)} clouds, {len(s)} primaries")
sys.exit(0 if len(wires) >= 2 or disp else 1)
PY
    then step "$it" coverage PASS "$(cat "$OUT/$it/coverage.txt")"
    else step "$it" coverage INCONCLUSIVE "neither multiwire nor displaced: $(cat "$OUT/$it/coverage.txt")"; fi
  fi
  # benchmarks only for a usable evaluation (a refused event has no response to time; its steps stay missing)
  if [[ $rcc -eq 0 ]]; then bench "$it" "$M"; fi
  hashes "$M" > "$OUT/$it/inputs_sha256_after.txt"
  if cmp -s "$OUT/$it/inputs_sha256_before.txt" "$OUT/$it/inputs_sha256_after.txt"; then step "$it" inputs_unchanged PASS; else step "$it" inputs_unchanged FAIL; fi
}
if [[ -n "${MWPC_M2_T1_RUN:-}" ]]; then full_birth "$MWPC_M2_T1_RUN" t1; fi
for R in ${MWPC_M2_REGRESSION_RUNS:-}; do full_birth "$R" regression; done

# ---------------------------------------------------------------- verdict (explicit coverage)
python3 - "$OUT" "$LABEL" <<'PY'
import csv, json, sys
from pathlib import Path
O, label = Path(sys.argv[1]), sys.argv[2]
rows = list(csv.DictReader(open(O / "steps.tsv"), delimiter="\t"))
req = list(csv.DictReader(open(O / "required.tsv"), delimiter="\t"))
BENCH = ["bench_verified", "bench_mapped", "bench_preload"]
NEED = {"e2b": ["evaluate", "reproduce_offline", "reproduce_loop", "production_rule"] +
               [f"displace_{d}" for d in ("1_0", "-1_0", "0_1", "0_-1", "2_-2", "-2_2")] + BENCH + ["inputs_unchanged"],
        "t1": ["slot", "t1_validated", "evaluate", "full_birth_vs_python_vector", "full_birth_vs_python_loop"] + BENCH + ["inputs_unchanged"]}
NEED["regression"] = NEED["t1"] + ["coverage"]
INFO_OK = {"production_rule"}
st = {(r["item"], r["step"]): r["status"] for r in rows}
missing, notpass = [], []
item_ok = {}
for q in req:
    it, kind = q["item"], q["kind"]
    good = True
    for s in NEED[kind]:
        v = st.get((it, s))
        if v is None:
            missing.append(f"{it}/{s}"); good = False
        elif not (v == "PASS" or (s in INFO_OK and v == "INFO")):
            notpass.append(f"{it}/{s}={v}"); good = False
    item_ok[it] = good
regs = [q["item"] for q in req if q["kind"] == "regression"]
coverage = {"e2b_muon_0560": item_ok.get("e2b_muon_0560", False), "e2b_muon_0997": item_ok.get("e2b_muon_0997", False),
            "t1_0000": item_ok.get("t1_0000", False), "regression_validated": [r for r in regs if item_ok.get(r)]}
if not regs:
    missing.append("regression item (MWPC_M2_REGRESSION_RUNS)")
fails = [f"{r['item']}/{r['step']}" for r in rows if r["status"] == "FAIL"]
complete = coverage["e2b_muon_0560"] and coverage["e2b_muon_0997"] and coverage["t1_0000"] and bool(coverage["regression_validated"])
verdict = "FAIL" if fails else ("PASS" if complete and not missing and not notpass else "INCONCLUSIVE")
bench = {}
metrics = {}
for it in sorted({r["item"] for r in rows}):
    for kind in BENCH:
        p = O / it / kind / "result.json"
        if p.exists():
            r = json.loads(p.read_text()); bench.setdefault(it, {})[kind] = {"clouds": r["clouds"], "groups": r["groups"],
                "strips": [r["strips_U"], r["strips_V"]], **r["timing"], **r["memory"]}
    for f in (O / it).glob("*.json"):
        try:
            j = json.loads(f.read_text())
        except Exception:
            continue
        if isinstance(j, dict) and j.get("tool") == "k2a_strip_compare.py":
            mi = j["metrics"].get("ion", {}); me = j["metrics"].get("electron", {})
            metrics[f"{it}/{f.stem}"] = {"status": j["status"], "ion_L2": [mi.get("U_L2"), mi.get("V_L2")],
                                        "ion_strip_charge_max_rel": [mi.get("U_strip_charge_max_rel"), mi.get("V_strip_charge_max_rel")],
                                        "ion_total_charge_fC": [mi.get("U_total_charge_ref_fC"), mi.get("V_total_charge_ref_fC")],
                                        "electron_L2": [me.get("U_L2"), me.get("V_L2")], "strips": [mi.get("U_strips"), mi.get("V_strips")]}
v = {"label": label, "verdict": verdict, "coverage": coverage, "failed_steps": fails, "missing": missing, "not_pass": notpass,
     "required": req, "steps": rows, "comparisons": metrics, "benchmark": bench}
(O / "m2_verdict.json").write_text(json.dumps(v, indent=1) + "\n")
L = [f"# Milestone 2 adapter validation ({label}): **{verdict}**", "",
     f"Coverage: 0560 {coverage['e2b_muon_0560']}, 0997 {coverage['e2b_muon_0997']}, T1 0000 {coverage['t1_0000']}, "
     f"validated regression items {coverage['regression_validated'] or 'none'}"]
if fails: L.append(f"Failed: {', '.join(fails)}")
if missing: L.append(f"Missing: {', '.join(missing)}")
if notpass: L.append(f"Not PASS: {', '.join(notpass)}")
L += ["", "| item | step | status | detail |", "|---|---|---|---|"]
L += [f"| {r['item']} | {r['step']} | {r['status']} | {r['detail'][:150]} |" for r in rows]
L += ["", "| comparison | ion L2 U / V | max strip charge rel U / V | electron L2 U / V |", "|---|---|---|---|"]
f = lambda x: "-" if x is None else f"{x:.1e}"
L += [f"| {k} | {f(m['ion_L2'][0])} / {f(m['ion_L2'][1])} | {f(m['ion_strip_charge_max_rel'][0])} / {f(m['ion_strip_charge_max_rel'][1])} | "
      f"{f(m['electron_L2'][0])} / {f(m['electron_L2'][1])} |" for k, m in sorted(metrics.items())]
L += ["", "| item | mode | clouds | groups | kernel verify s | kernel load s | cold ms | warm ms | peak RSS after kernel / end KiB |", "|---|---|---|---|---|---|---|---|---|"]
for it, b in bench.items():
    for kind, x in b.items():
        w = x["evaluate_warm_median_s"]
        L.append(f"| {it} | {kind[6:]} | {x['clouds']} | {x['groups']} | {x['kernel_verify_s']:.3g} | {x['kernel_load_s']:.3g} | "
                 f"{1e3 * x['evaluate_cold_s']:.3g} | {'-' if w is None else format(1e3 * w, '.3g')} | {x['peak_rss_kib_after_kernel']} / {x['peak_rss_kib_end']} |")
(O / "m2_summary.md").write_text("\n".join(L) + "\n")
print("\n".join(L[:6]))
sys.exit({"PASS": 0, "FAIL": 1}.get(verdict, 3))
PY

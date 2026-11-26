#!/usr/bin/env bash
# Milestone 2 validation of the compiled absolute-strip K2A adapter on REAL moments with the REAL frozen kernel.
#
#   bash run_m2_adapter.sh OUTDIR
#
# env (all explicit, nothing guessed):
#   MWPC_K2A_INLINE_BUILD    build directory containing k2a_strip_eval (default <this dir>/build)
#   MWPC_K2A_KERNEL_DIR      the validated kernel directory (RUN/kernel_TAG); MWPC_K2A_KERNEL_ID its validated id
#   MWPC_M2_E2B_DIRS         space-separated E2b muon directories of k2a_recenter_experiment.sh far mode
#                            (.../recenter_experiment_*/muon_0560 .../muon_0997): prepared/cloud_moments.csv,
#                            groups/groups.json, K2A_local/{replay,prompt} (the offline K2A response), read-only
#   MWPC_M2_INLINE_PREFIXES  space-separated full-birth moment prefixes of the inline generator (PREFIX.moments.csv),
#                            e.g. the accepted T1 run R1/inline and the multiwire regression run R1/inline
#   NA60_KERNEL_EVAL_LIB     optional (Python 'compiled' reference); MWPC_M2_REPEAT (default 20)
# Per E2b muon (REQUIRED unless marked INFO):
#   reproduce_offline   C++ with the E2b grouping (--shifts groups.json) and the offline run's model / half strips /
#                       window (K2A_local/recentered.json, kernel id checked) == the archived offline K2A_local response
#   reproduce_loop      the same == a fresh run of the validated loop algorithm (k2a_recentered_k2a.py --impl loop)
#   production_rule     C++ with the production rule (shifts from the cloud moments): INFO - number of primaries whose
#                       rule shift differs from the E2b export-based shift; if none, the files must be identical
#   displace_<dU>_<dV>  the same event displaced by +-1 / +-2 strips per plane: same response on shifted labels
#   bench               --repeat N, kernel verified; and --no-verify, mapped vs --preload: timing and peak RSS
# Per inline prefix:
#   full_birth_vs_python C++ (production rule) == k2a_recentered_k2a.py (vector and loop) on the same full-birth
#                       moments; if the event is outside the kernel domain both must refuse (C++ exit 4 with the
#                       named status, Python an error) - a one-sided refusal is a FAIL
#   bench
# Verdict: m2_verdict.json / m2_summary.md, PASS / FAIL / INCONCLUSIVE (missing input). Read-only on every input.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CIR="$(cd "$HERE/.." && pwd)"
REC="$CIR/k2a-campaign/recenter/k2a_recentered_k2a.py"; CMP="$HERE/k2a_strip_compare.py"
fail(){ echo "ERROR: $*" >&2; exit 2; }
OUT="${1:?usage: run_m2_adapter.sh OUTDIR}"; [[ ! -e "$OUT" ]] || fail "$OUT exists"
B="${MWPC_K2A_INLINE_BUILD:-$HERE/build}"; EV="$B/k2a_strip_eval"; [[ -x "$EV" ]] || fail "missing $EV"
KD="${MWPC_K2A_KERNEL_DIR:?set MWPC_K2A_KERNEL_DIR}"; KID="${MWPC_K2A_KERNEL_ID:?set MWPC_K2A_KERNEL_ID}"
[[ -s "$KD/kernel.json" && -s "$KD/K.npy" ]] || fail "not a kernel directory: $KD"
REP="${MWPC_M2_REPEAT:-20}"
case "$(cd "$(dirname "$OUT")" 2>/dev/null && pwd)/" in
  "$HOME/na6p/wt-validation/"*|*test_runs/mwpc_stage_c/s_ladder_2100/*) fail "OUTDIR inside a frozen reference tree";;
esac
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
STEPS="$OUT/steps.tsv"; printf 'item\tstep\tstatus\tdetail\n' > "$STEPS"
step(){ printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$STEPS"; echo "[$3] $1 $2 ${4:-}"; }
{ echo "date_utc $(date -u +%FT%TZ)"; echo "git_commit $(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "binary_sha256 $(sha256sum "$EV")"; echo "kernel_dir $KD"; echo "kernel_id_expected $KID"
  if [[ -s "$B/k2a_inline_build_info.txt" ]]; then sed 's/^/build_info /' "$B/k2a_inline_build_info.txt"; fi; } > "$OUT/provenance.txt"
hashes(){ for d in "$@"; do find -L "$d" -type f \( -name '*.csv' -o -name '*.json' \) 2>/dev/null | sort | xargs -r sha256sum; done; }
cmpstep(){ # item step ref test [compare args...]
  local it="$1" st="$2" r="$3" t="$4"; shift 4
  set +e; python3 "$CMP" "$r" "$t" "$OUT/$it/$st.json" --name "$it $st" "$@" > "$OUT/$it/$st.log" 2>&1; local rc=$?; set -e
  case $rc in 0) step "$it" "$st" PASS "$(tail -1 "$OUT/$it/$st.log")";; 3) step "$it" "$st" INCONCLUSIVE "$(tail -1 "$OUT/$it/$st.log")";;
              *) step "$it" "$st" FAIL "$(tail -1 "$OUT/$it/$st.log")";; esac
}
bench(){ # item moments [extra]
  local it="$1" m="$2"; shift 2
  "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$m" --out "$OUT/$it/bench_verified" --repeat "$REP" "$@" > "$OUT/$it/bench_verified.log" 2>&1 || true
  "$EV" --kernel "$KD" --no-verify --moments "$m" --out "$OUT/$it/bench_mapped" --repeat "$REP" "$@" > "$OUT/$it/bench_mapped.log" 2>&1 || true
  "$EV" --kernel "$KD" --no-verify --preload --moments "$m" --out "$OUT/$it/bench_preload" --repeat "$REP" "$@" > "$OUT/$it/bench_preload.log" 2>&1 || true
  if [[ -s "$OUT/$it/bench_mapped/result.json" ]]; then step "$it" bench INFO "$(tail -1 "$OUT/$it/bench_mapped.log")"; else step "$it" bench FAIL "no timing"; fi
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
  python3 "$REC" "$KD" "$M" "$OUT/$it/py_loop" --groups "$G" --impl loop --kernel-id "$KID" --model "$MODEL" --half-strips "$HSx" \
    --window-ns "$WIN" > "$OUT/$it/py_loop.log" 2>&1 \
    && cmpstep "$it" reproduce_loop "$OUT/$it/py_loop" "$OUT/$it/cc_e2b_groups" || step "$it" reproduce_loop FAIL "python loop failed"
  "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "$OUT/$it/cc_rule" "${X[@]}" > "$OUT/$it/cc_rule.log" 2>&1 || true
  python3 - "$OUT/$it" "$G" <<'PY' && step "$it" production_rule INFO "$(cat "$OUT/$it/rule_vs_e2b.txt")" || step "$it" production_rule FAIL "$(cat "$OUT/$it/rule_vs_e2b.txt" 2>/dev/null)"
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

# ---------------------------------------------------------------- full-birth inline moments: C++ vs Python engine
for PFX in ${MWPC_M2_INLINE_PREFIXES:-}; do
  it="inline_$(basename "$(dirname "$(dirname "$PFX")")")_$(basename "$(dirname "$PFX")")"; mkdir "$OUT/$it"
  M="$PFX.moments.csv"
  if [[ ! -s "$M" ]]; then step "$it" inputs INCONCLUSIVE "missing $M"; continue; fi
  set +e; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$M" --out "$OUT/$it/cc" > "$OUT/$it/cc.log" 2>&1; rcc=$?
  python3 "$REC" "$KD" "$M" "$OUT/$it/py_vector" --impl vector --kernel-id "$KID" > "$OUT/$it/py_vector.log" 2>&1; rpv=$?
  python3 "$REC" "$KD" "$M" "$OUT/$it/py_loop" --impl loop --kernel-id "$KID" > "$OUT/$it/py_loop.log" 2>&1; rpl=$?
  set -e
  if [[ $rcc -eq 0 && $rpv -eq 0 && $rpl -eq 0 ]]; then
    step "$it" evaluate PASS "$(tail -1 "$OUT/$it/cc.log")"
    cmpstep "$it" full_birth_vs_python_vector "$OUT/$it/py_vector" "$OUT/$it/cc"
    cmpstep "$it" full_birth_vs_python_loop "$OUT/$it/py_loop" "$OUT/$it/cc"
  elif [[ $rcc -eq 4 && $rpv -ne 0 && $rpl -ne 0 ]]; then
    step "$it" evaluate PASS "both refuse: C++ $(grep -o '"status_names": "[^"]*"' "$OUT/$it/cc/result.json"); Python: $(tail -1 "$OUT/$it/py_vector.log")"
  else
    step "$it" evaluate FAIL "inconsistent: C++ exit $rcc, Python vector $rpv, loop $rpl"
  fi
  bench "$it" "$M"
done

# ---------------------------------------------------------------- verdict
python3 - "$OUT" <<'PY'
import csv, json, sys
from pathlib import Path
O = Path(sys.argv[1])
rows = list(csv.DictReader(open(O / "steps.tsv"), delimiter="\t"))
items = sorted({r["item"] for r in rows})
st = [r["status"] for r in rows if r["status"] != "INFO"]
verdict = "FAIL" if "FAIL" in st else ("INCONCLUSIVE" if "INCONCLUSIVE" in st or not items else "PASS")
bench = {}
for it in items:
    for kind in ("bench_verified", "bench_mapped", "bench_preload"):
        p = O / it / kind / "result.json"
        if p.exists():
            r = json.loads(p.read_text()); bench.setdefault(it, {})[kind] = {"clouds": r["clouds"], "groups": r["groups"],
                "strips": [r["strips_U"], r["strips_V"]], **r["timing"], **r["memory"]}
metrics = {}
for it in items:
    for f in (O / it).glob("*.json"):
        try:
            j = json.loads(f.read_text())
        except Exception:
            continue
        if j.get("tool") == "k2a_strip_compare.py":
            mi = j["metrics"].get("ion", {}); me = j["metrics"].get("electron", {})
            metrics[f"{it}/{f.stem}"] = {"status": j["status"], "ion_L2": [mi.get("U_L2"), mi.get("V_L2")],
                                        "ion_strip_charge_max_rel": [mi.get("U_strip_charge_max_rel"), mi.get("V_strip_charge_max_rel")],
                                        "ion_total_charge_fC": [mi.get("U_total_charge_ref_fC"), mi.get("V_total_charge_ref_fC")],
                                        "electron_L2": [me.get("U_L2"), me.get("V_L2")], "strips": [mi.get("U_strips"), mi.get("V_strips")]}
v = {"verdict": verdict, "items": items, "steps": rows, "comparisons": metrics, "benchmark": bench}
(O / "m2_verdict.json").write_text(json.dumps(v, indent=1) + "\n")
L = [f"# Milestone 2 adapter validation: **{verdict}**", "", "| item | step | status | detail |", "|---|---|---|---|"]
L += [f"| {r['item']} | {r['step']} | {r['status']} | {r['detail'][:150]} |" for r in rows]
L += ["", "| comparison | ion L2 U / V | max strip charge rel U / V | electron L2 U / V |", "|---|---|---|---|"]
f = lambda x: "-" if x is None else f"{x:.1e}"
L += [f"| {k} | {f(m['ion_L2'][0])} / {f(m['ion_L2'][1])} | {f(m['ion_strip_charge_max_rel'][0])} / {f(m['ion_strip_charge_max_rel'][1])} | "
      f"{f(m['electron_L2'][0])} / {f(m['electron_L2'][1])} |" for k, m in metrics.items()]
L += ["", "| item | mode | clouds | groups | kernel verify s | kernel load s | cold ms | warm ms | peak RSS after kernel / end KiB |", "|---|---|---|---|---|---|---|---|---|"]
for it, b in bench.items():
    for kind, x in b.items():
        L.append(f"| {it} | {kind[6:]} | {x['clouds']} | {x['groups']} | {x['kernel_verify_s']:.3g} | {x['kernel_load_s']:.3g} | "
                 f"{1e3 * x['evaluate_cold_s']:.3g} | {'-' if x['evaluate_warm_median_s'] is None else format(1e3 * x['evaluate_warm_median_s'], '.3g')} | {x['peak_rss_kib_after_kernel']} / {x['peak_rss_kib_end']} |")
(O / "m2_summary.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
sys.exit({"PASS": 0, "FAIL": 1}.get(verdict, 3))
PY

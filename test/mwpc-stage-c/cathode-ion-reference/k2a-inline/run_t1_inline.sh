#!/usr/bin/env bash
# T1: one small deterministic 2100 V avalanche, in-memory full-birth K2A moments vs independent full-birth references.
#
#   bash run_t1_inline.sh SLOT OUTDIR [N_HEED_PRIMARIES]        (N default 3)
#
# SLOT     4-digit NON-sealed slot of the frozen 997-muon K256 production archive; its input.csv (the exact frozen
#          one-muon Geant4 input) is copied, never modified. Sealed slots are refused (slots.csv role).
# OUTDIR   NEW directory; refused inside the frozen archive tree, inside ~/na6p/wt-validation, or inside any
#          test_runs/mwpc_stage_c/s_ladder_2100 reference tree.
# env      MWPC_K2A_INLINE_BUILD   k2a-inline CMake build directory (default: <this dir>/build)
#          MWPC_ARCHIVE_2100_K256  frozen production K256 archive root (default as the recentering runner)
#          MWPC_K2A_PF             preflight directory with inventory/slots.csv (required: sealed-slot check)
#          MWPC_T1_TIMEOUT_SEC     timeout per generator run (default 7200)
#          MWPC_T1_PREVIEW=1       only dump the frozen archive census of primaries 0..N-1 (cost preview), run nothing
#          MWPC_T1_MAX_BIRTHS      refuse if the frozen archive's census of primaries 0..N-1 exceeds this (default 2e6)
#          MWPC_T1_SYNTHETIC=1     TEST MODE ONLY (stub Garfield/ROOT binaries): MWPC_T1_INPUT gives the input CSV,
#                                  MWPC_T1_FROZEN_ARCHIVE an optional "frozen" archive; slot checks are skipped and
#                                  every output is labelled SYNTHETIC
# Runs (each in OUTDIR/Rx, /usr/bin/time -v (or rss_run.py if GNU time is absent), timeout; the Heed primaries 0..N-1 keep their original indices/seeds):
#   R1  full archive (every birth, weight 1) + in-memory moments, ONE run: both sinks see the same birth stream
#   R2  R1 again (reproducibility)
#   R3  production path: moments only, no archive, reservoir 0
#   R0  R3 without moments (memory baseline; frozen behaviour)
#   R4  weighted K=256 archive + moments (T2 input; the reservoir rows must be a bitwise subset of R1's births)
# Verdict (T1_verdict.json): PASS only if every REQUIRED item passes; any failure -> FAIL; anything missing ->
# INCONCLUSIVE. REQUIRED: allowed-diff check; R0-R4 exit 0; T1 comparison (census, seeds, identity, frame, statuses,
# moments vs both references, charge); R1 = R2 (moments and archive content); R3 = R1 and R4 = R1 (moments);
# R4 reservoir rows subset of R1; R0 and R3 per-seed physics CSVs identical. Reported separately, never part of the
# verdict: T2 (K256 vs full, diagnostic) and the cross-check against the frozen production K256 archive (MATCH /
# MISMATCH, environment-dependent). The frozen archive's sha256 is recorded before and after (must be unchanged).
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
fail(){ echo "ERROR: $*" >&2; exit 2; }
[[ $# -ge 2 ]] || fail "usage: run_t1_inline.sh SLOT OUTDIR [N_HEED_PRIMARIES]"
SLOT="$1"; OUT="$2"; N="${3:-3}"
[[ "$N" =~ ^[1-9][0-9]*$ ]] || fail "N must be a positive integer"
SYN="${MWPC_T1_SYNTHETIC:-0}"
B="${MWPC_K2A_INLINE_BUILD:-$HERE/build}"
EXE="$B/mwpc_p3_noFW_penning057_k2a_inline_v1"; DUMP="$B/k2a_archive_dump"; BM="$B/p3_birth_moments"
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
TMO="${MWPC_T1_TIMEOUT_SEC:-7200}"
for x in "$EXE" "$DUMP" "$BM"; do [[ -x "$x" ]] || fail "missing executable $x (build $HERE first)"; done
for c in python3 sha256sum timeout; do command -v "$c" >/dev/null || fail "missing $c"; done
if [[ -x /usr/bin/time ]]; then TIMEV=(/usr/bin/time -v -o); else TIMEV=(python3 "$HERE/rss_run.py"); fi   # peak RSS
python3 -c "import numpy" || fail "numpy missing"
[[ ! -e "$OUT" ]] || fail "$OUT exists (no overwrite)"
mkdir -p "$(dirname "$OUT")"; OUTR="$(cd "$(dirname "$OUT")" && pwd)/$(basename "$OUT")"
case "$OUTR/" in
  "$(cd "$ARCHIVE" 2>/dev/null && pwd || echo /nonexistent)/"*) fail "OUTDIR inside the frozen archive tree";;
  "$HOME/na6p/wt-validation/"*) fail "OUTDIR inside the frozen validation worktree";;
  *test_runs/mwpc_stage_c/s_ladder_2100/*) fail "OUTDIR inside the frozen S-ladder reference tree";;
esac
LABEL=REAL; [[ "$SYN" == 1 ]] && LABEL=SYNTHETIC

# ---- input, slot role, frozen archive integrity (all checked BEFORE anything is created) -----------------------------
FROZEN_ROOT=""; SRC_INPUT=""; ROLE=""
if [[ "$SYN" == 1 ]]; then
  [[ -s "${MWPC_T1_INPUT:-}" ]] || fail "MWPC_T1_INPUT required in synthetic mode"
  SRC_INPUT="$MWPC_T1_INPUT"; FROZEN_ROOT="${MWPC_T1_FROZEN_ARCHIVE:-}"; ROLE="synthetic (not checked)"
else
  [[ "$SLOT" =~ ^[0-9]{4}$ ]] || fail "SLOT must be 4 digits"
  SL="${MWPC_K2A_PF:-}/inventory/slots.csv"; [[ -s "$SL" ]] || fail "MWPC_K2A_PF/inventory/slots.csv required (sealed check)"
  ROLE=$(awk -F, -v s="$SLOT" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["slot"]==s{print $h["role"]}' "$SL")
  [[ -n "$ROLE" && "$ROLE" != sealed_validation ]] || fail "slot $SLOT is sealed or unknown (role '$ROLE')"
  D="$ARCHIVE/muon_$SLOT"
  [[ -s "$D/input.csv" && -s "$D/ion_births.root" && -s "$D/SHA256SUMS" ]] || fail "frozen archive slot $D incomplete"
  (cd "$D" && sha256sum -c --quiet SHA256SUMS) || fail "frozen archive slot fails its SHA256SUMS"
  SRC_INPUT="$D/input.csv"; FROZEN_ROOT="$D/ion_births.root"
fi
[[ "$(wc -l < "$SRC_INPUT")" -eq 2 ]] || fail "input must hold exactly one muon"
mkdir "$OUTR"; OUT="$OUTR"
STEPS="$OUT/steps.tsv"; printf 'step\tstatus\tdetail\n' > "$STEPS"
step(){ printf '%s\t%s\t%s\n' "$1" "$2" "${3:-}" >> "$STEPS"; echo "[$2] $1 ${3:-}"; }
echo "label=$LABEL slot=$SLOT N=$N out=$OUT"
step slot_role "$([[ "$SYN" == 1 ]] && echo SKIPPED || echo PASS)" "$ROLE"
[[ "$SYN" == 1 ]] || step frozen_archive_sha256sums PASS "$(dirname "$FROZEN_ROOT")"
cp "$SRC_INPUT" "$OUT/input.csv"; cmp -s "$SRC_INPUT" "$OUT/input.csv" || fail "input copy differs"
if [[ -n "$FROZEN_ROOT" ]]; then sha256sum "$FROZEN_ROOT" > "$OUT/frozen_archive_sha256_before.txt"; fi

# ---- cost preview from the frozen archive census (expected births of primaries 0..N-1; same seeds) ---------------
if [[ -n "$FROZEN_ROOT" && -s "$FROZEN_ROOT" ]]; then
  "$DUMP" "$FROZEN_ROOT" "$OUT/dump_frozen" > "$OUT/dump_frozen.log" 2>&1 || fail "cannot dump the frozen archive"
  EXPB=$(awk -F, -v n="$N" 'NR>1 && $1<n {s+=$5} END{print s+0}' "$OUT/dump_frozen/primaries.csv")
  MAXB="${MWPC_T1_MAX_BIRTHS:-2000000}"
  step expected_births INFO "primaries 0..$((N-1)) of the frozen archive: $EXPB ion births (limit $MAXB)"
  (( EXPB <= MAXB )) || fail "expected $EXPB births > MWPC_T1_MAX_BIRTHS=$MAXB; choose a smaller N or another slot"
  awk -F, -v n="$N" 'NR==1 || $1<n' "$OUT/dump_frozen/primaries.csv" | cut -d, -f1,3,5
fi
if [[ "${MWPC_T1_PREVIEW:-0}" == 1 ]]; then echo "PREVIEW only (MWPC_T1_PREVIEW=1): nothing was simulated; $OUT holds the preview"; exit 0; fi

# ---- provenance ---------------------------------------------------------------------------------------------------
{
  echo "label $LABEL"; echo "date_utc $(date -u +%FT%TZ)"; echo "host $(hostname)"
  echo "git_commit $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_branch $(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
  echo "git_dirty_paths $(git -C "$REPO" status --porcelain -- test/mwpc-stage-c/cathode-ion-reference test/mwpc-avalanche 2>/dev/null | wc -l)"
  for x in "$EXE" "$DUMP" "$BM"; do echo "binary_sha256 $(sha256sum "$x")"; done
  echo "input_sha256 $(sha256sum "$OUT/input.csv" | awk '{print $1}')"
  echo "slot $SLOT"; echo "max_heed_primaries $N"
  if [[ -s "$B/k2a_inline_build_info.txt" ]]; then sed 's/^/build_info /' "$B/k2a_inline_build_info.txt"; fi
  grep -E '^CMAKE_BUILD_TYPE|^CMAKE_CXX_FLAGS:' "$B/CMakeCache.txt" 2>/dev/null | sed 's/^/cmake_cache /' || true
} > "$OUT/provenance.txt"
if python3 "$HERE/check_inline_generator_diff.py" --report "$OUT/diff_check.json" > "$OUT/diff_check.log" 2>&1; then
  step diff_check PASS "$(tail -1 "$OUT/diff_check.log")"; else step diff_check FAIL "$(tail -1 "$OUT/diff_check.log")"; fi

# ---- generator runs -----------------------------------------------------------------------------------------------
run(){  # name, args...
  local r="$1"; shift; mkdir "$OUT/$r"
  set +e
  "${TIMEV[@]}" "$OUT/$r/resource_time.txt" timeout --signal=TERM --kill-after=30s "$TMO" "$EXE" "$@" > "$OUT/$r/physics.log" 2>&1
  local rc=$?; set -e
  echo "$rc" > "$OUT/$r/exit_code.txt"
  if [[ $rc -eq 0 ]]; then step "run_$r" PASS "rss_kib=$(awk -F': ' '/Maximum resident/{print $2}' "$OUT/$r/resource_time.txt")"
  else step "run_$r" FAIL "exit $rc: $(tail -1 "$OUT/$r/physics.log")"; fi
}
I="$OUT/input.csv"
run R1 --cathode-ion-reference "$I" "$OUT/R1/microscopic" 2100 0 full 0 "$OUT/R1/ion_births.root" --k2a-moments-out "$OUT/R1/inline" --max-heed-primaries "$N"
run R2 --cathode-ion-reference "$I" "$OUT/R2/microscopic" 2100 0 full 0 "$OUT/R2/ion_births.root" --k2a-moments-out "$OUT/R2/inline" --max-heed-primaries "$N"
run R3 --wiretruth-hv-p3-noFW "$I" "$OUT/R3/microscopic" 2100 0 --k2a-moments-out "$OUT/R3/inline" --max-heed-primaries "$N"
run R0 --wiretruth-hv-p3-noFW "$I" "$OUT/R0/microscopic" 2100 0 --max-heed-primaries "$N"
run R4 --cathode-ion-reference "$I" "$OUT/R4/microscopic" 2100 0 weighted 256 "$OUT/R4/ion_births.root" --k2a-moments-out "$OUT/R4/inline" --max-heed-primaries "$N"

# ---- references and comparisons -----------------------------------------------------------------------------------
tool(){  # step name, command...
  local s="$1"; shift
  set +e; "$@" > "$OUT/$s.log" 2>&1; local rc=$?; set -e
  case $rc in 0) step "$s" PASS "$(tail -1 "$OUT/$s.log")";; 3) step "$s" INCONCLUSIVE "$(tail -1 "$OUT/$s.log")";;
              *) step "$s" FAIL "exit $rc: $(tail -1 "$OUT/$s.log")";; esac
}
for r in R1 R2 R4; do
  if [[ -s "$OUT/$r/ion_births.root" ]]; then tool "dump_$r" "$DUMP" "$OUT/$r/ion_births.root" "$OUT/dump_$r"; fi
done
if [[ -s "$OUT/R1/ion_births.root" ]]; then tool batch_R1 "$BM" "$OUT/R1/ion_births.root" "$OUT/batch_R1.moments.csv"; fi
if [[ -s "$OUT/R4/ion_births.root" ]]; then tool batch_R4 "$BM" "$OUT/R4/ion_births.root" "$OUT/batch_R4.moments.csv"; fi
tool T1_compare python3 "$HERE/k2a_inline_compare.py" --inline "$OUT/R1/inline" --dump "$OUT/dump_R1" --batch "$OUT/batch_R1.moments.csv" \
     --input "$I" --max-heed-primaries "$N" --out "$OUT/T1_compare"
ID="$HERE/k2a_inline_identity.py"
tool identity_moments_R1_R2 python3 "$ID" moments "$OUT/R1/inline" "$OUT/R2/inline" "$OUT/identity_moments_R1_R2.json"
tool identity_archive_R1_R2 python3 "$ID" dumps "$OUT/dump_R1" "$OUT/dump_R2" "$OUT/identity_archive_R1_R2.json"
tool identity_moments_R3_R1 python3 "$ID" moments "$OUT/R3/inline" "$OUT/R1/inline" "$OUT/identity_moments_R3_R1.json"
tool identity_moments_R4_R1 python3 "$ID" moments "$OUT/R4/inline" "$OUT/R1/inline" "$OUT/identity_moments_R4_R1.json"
tool subset_R4_in_R1 python3 "$ID" subset "$OUT/dump_R1" "$OUT/dump_R4" "$OUT/subset_R4_in_R1.json"
seedcmp(){ local ok=1; for f in microscopic_seed.csv microscopic_wire_seed.csv microscopic_ion_birth_sample.csv; do cmp -s "$OUT/R0/$f" "$OUT/R3/$f" || ok=0; done; [[ $ok == 1 ]]; }
if [[ -s "$OUT/R0/microscopic_seed.csv" && -s "$OUT/R3/microscopic_seed.csv" ]]; then
  if seedcmp; then step physics_csv_R0_R3 PASS "per-seed CSVs identical with and without moments"; else step physics_csv_R0_R3 FAIL "per-seed CSVs differ"; fi
else step physics_csv_R0_R3 INCONCLUSIVE "missing per-seed CSVs"; fi
# diagnostics (never part of the verdict)
if [[ -s "$OUT/batch_R4.moments.csv" ]]; then
  python3 "$HERE/k2a_inline_t2.py" "$OUT/R1/inline" "$OUT/batch_R4.moments.csv" "$OUT/T2" > "$OUT/T2.log" 2>&1 || true
fi
if [[ -n "$FROZEN_ROOT" && -s "$FROZEN_ROOT" ]]; then
  python3 "$ID" subset "$OUT/dump_R1" "$OUT/dump_frozen" "$OUT/frozen_crosscheck.json" --primaries "$N" --label frozen-production > "$OUT/frozen_crosscheck.log" 2>&1 || true
  sha256sum "$FROZEN_ROOT" > "$OUT/frozen_archive_sha256_after.txt"
  if cmp -s "$OUT/frozen_archive_sha256_before.txt" "$OUT/frozen_archive_sha256_after.txt"; then step frozen_archive_unchanged PASS
  else step frozen_archive_unchanged FAIL "frozen archive sha256 changed"; fi
fi

# ---- verdict ------------------------------------------------------------------------------------------------------
python3 - "$OUT" "$LABEL" <<'PY'
import csv, json, re, sys
from pathlib import Path
O, label = Path(sys.argv[1]), sys.argv[2]
st = {r["step"]: r for r in csv.DictReader(open(O / "steps.tsv"), delimiter="\t")}
req = ["diff_check", "run_R1", "run_R2", "run_R3", "run_R0", "run_R4", "dump_R1", "dump_R2", "dump_R4", "batch_R1", "batch_R4",
       "T1_compare", "identity_moments_R1_R2", "identity_archive_R1_R2", "identity_moments_R3_R1", "identity_moments_R4_R1",
       "subset_R4_in_R1", "physics_csv_R0_R3"] + (["frozen_archive_unchanged"] if "frozen_archive_unchanged" in st else [])
s = {k: st[k]["status"] if k in st else "MISSING" for k in req}
verdict = "FAIL" if any(v == "FAIL" for v in s.values()) else ("INCONCLUSIVE" if any(v != "PASS" for v in s.values()) else "PASS")
def rss(r):
    p = O / r / "resource_time.txt"
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", p.read_text()) if p.exists() else None
    return int(m.group(1)) if m else None
def jl(p):
    try: return json.loads(Path(p).read_text())
    except Exception: return None
mem = {r: rss(r) for r in ("R0", "R1", "R2", "R3", "R4")}
inl = jl(O / "R3" / "inline.inline.json") or {}
cmpj = jl(O / "T1_compare" / "compare.json") or {}
out = {"label": label, "verdict": verdict, "required": s,
       "memory": {"peak_rss_kib": mem, "rss_overhead_kib_R3_minus_R0": (mem["R3"] - mem["R0"]) if mem["R3"] and mem["R0"] else None,
                  "accumulator_state_bytes_peak": inl.get("accumulator_state_bytes_peak")},
       "t1": {k: cmpj.get(k) for k in ("n_primaries", "births", "n_clouds", "bitwise_inline_vs_batch", "worst_diff_over_tolerance")},
       "diagnostic_T2": (jl(O / "T2" / "t2.json") or {}).get("status", "not run"),
       "diagnostic_frozen_production_crosscheck": (jl(O / "frozen_crosscheck.json") or {}).get("status", "not run")}
(O / "T1_verdict.json").write_text(json.dumps(out, indent=1) + "\n")
L = [f"# T1 inline full-birth moments ({label}): **{verdict}**", "", "| required item | status |", "|---|---|"]
L += [f"| {k} | {v} |" for k, v in s.items()]
L += ["", f"Peak RSS (kB): {mem}; R3 - R0 = {out['memory']['rss_overhead_kib_R3_minus_R0']} kB; "
      f"accumulator state {out['memory']['accumulator_state_bytes_peak']} bytes",
      f"Diagnostic T2 (K256 vs full): {out['diagnostic_T2']}; frozen production cross-check: {out['diagnostic_frozen_production_crosscheck']}",
      "", "See T1_compare/summary.md (cloud-by-cloud: T1_compare/clouds_comparison.csv, census: T1_compare/primary_census.csv)."]
(O / "T1_summary.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
sys.exit({"PASS": 0, "FAIL": 1}.get(verdict, 3))
PY

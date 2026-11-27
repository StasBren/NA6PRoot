#!/usr/bin/env bash
# K2A-versus-S0 validation on independent NON-SEALED muons of the frozen 2100 V K256 archive (no new avalanche).
#
#   bash run_k2a_validation.sh RUN TAG OUTDIR SLOT [SLOT ...]
#
# RUN/TAG   the validated reference run and kernel tag (e.g. .../s_ladder_2100/run_20261124_124740 r14t72edge5)
# OUTDIR    NEW directory (refused if it exists or lies inside RUN, the archive, ~/na6p/wt-validation or any
#           s_ladder_2100 reference tree). Nothing outside OUTDIR is written.
# SLOT      4-digit preflight slots (k2a_validation_select.py pilot_slots.txt / full_slots.txt). ALL slots are checked
#           before anything runs: a sealed_validation slot aborts the whole run; a slot used to build or tune the
#           kernel/recentering (RUN/provenance_prepare_<slot>.txt, or MWPC_K2A_VAL_DEV_SLOTS) is refused unless
#           MWPC_K2A_VAL_ALLOW_DEV=1, and is then labelled NOT independent.
#
# Per muon, from the SAME archived weighted ion births (no regeneration, frozen tools unchanged):
#   S0_T1   k2a_recenter_experiment.sh far (frozen RKF exporter, T1 collapse, per-primary groups, frozen replay with the
#           late-step ladder + approximate electron step per group, k2a_local_sum census) = the established E2b
#           reference; also the Python K2A_local, truncation check and far analysis. Run per muon in OUTDIR/far/.
#   S0_T0   the same groups built from the UNCOLLAPSED frozen export (true per-ion birth times; group membership must
#           equal S0_T1's), same ladder replays, electron step and census sum. Disable with MWPC_K2A_VAL_T0=0.
#   K2A     production path: p3_birth_moments (archive -> compact moments, the frozen archive reader) -> compiled
#           k2a_strip_eval (validated kernel, id verified) with the S0 grouping imposed (--shifts groups.json); the
#           production rule is also run and any differing primary shift reported.
#   checks  K2A(C++, export moments) == Python K2A_local (1e-12: implementation identity);
#           K2A(C++, archive moments) == K2A(C++, export moments) (1e-7: moment path; s_compare_moments at 1e-9 first)
#   compare k2a_recenter_compare (k2a_vs_s0): S0_T1 vs K2A, S0_T0 vs K2A, S0_T0 vs S0_T1 (time-collapse part);
#           t_shaped_compare (provisional CR-RC^3, Tp 100/200 ns): ion+electron with reference S0_T1 and S0_T0, and
#           ion-only with reference S0_T1.
#   outcome k2a_validation_muon.py: VALID / OUT_OF_DOMAIN (K2A blocking statuses named, S0 still computed) /
#           S0_FAILED (reference unavailable, cause named) / FAIL (numerical or integrity failure) / INCONCLUSIVE
#           (something missing); runtime split; DIAGNOSTIC dry run of the frozen sealed criteria v1. Nothing is
#           clamped, excluded or re-run selectively.
# Then k2a_validation_report.py OUTDIR -> validation_report.{json,md}, per_muon.csv.
#
# env: as k2a_recenter_experiment.sh far (MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR, MWPC_EXPORT_BIN, MWPC_ION_MOBILITY,
#      MWPC_K2A_PF, optional MWPC_K2A_KERNEL_ID, MWPC_ARCHIVE_2100_K256, MWPC_K2A_LATE_STEPS), plus
#      MWPC_BIRTHS_BIN        built spatial-ablation p3_birth_moments (required)
#      MWPC_K2A_INLINE_BUILD  directory with the built k2a_strip_eval (default ../../k2a-inline/build)
#      MWPC_K2A_VAL_T0        1 (default) build S0_T0; 0 = S0_T1 only (outcome then VALID_T1_ONLY)
#      MWPC_K2A_VAL_REPEAT    k2a_strip_eval repetitions (cold + warm median; default 20)
#      MWPC_K2A_VAL_BIRTHS_REPEAT  p3_birth_moments repetitions (min time; default 5)
#      MWPC_K2A_VAL_DEV_SLOTS default 0000,0858,0560,0997,0133 (used for development decisions)
#      MWPC_K2A_VAL_SELECTION selection.csv of k2a_validation_select.py (features in the report; optional)
#      MWPC_K2A_VAL_MIN_FREE_GIB  free space required before each muon (default 24)
#      MWPC_K2A_VAL_PRUNE_VERTICES=1  after a muon's outcome is written, delete ONLY the derived ion-path copies made by
#                             this run (far/.../prepared/frozen_export, groups/g_*, T0/groups/g_* vertices.f64; ~1 export
#                             size each), recording their sha256 in muons/muon_M/pruned_vertices.sha256. They are rebuilt
#                             bit-identically from the frozen archive by rerunning; T1_export, every waveform, comparison
#                             and log are kept. Default 0 (keep everything).
#      MWPC_K2A_VAL_SYNTHETIC=1   TEST MODE: every output labelled SYNTHETIC
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CAMP="$(cd "$HERE/.." && pwd)"; CIR="$(cd "$CAMP/.." && pwd)"
RC="$CAMP/recenter"; TA="$CIR/timing-ablation"; SA="$CIR/spatial-ablation"; INL="$CIR/k2a-inline"
set -u
umask 077
fail(){ echo "K2A VALIDATION ABORTED: $*" >&2; exit 4; }
[[ $# -ge 4 ]] || fail "usage: run_k2a_validation.sh RUN TAG OUTDIR SLOT [SLOT ...]"
RUN="$(cd "$1" && pwd)"; TAG="$2"; OUT="$3"; shift 3
KD="$RUN/kernel_$TAG"; [[ -s "$KD/kernel.json" && -s "$KD/K.npy" ]] || fail "missing kernel $KD"
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
HS="${MWPC_HALF_STRIPS:-5}"; LATE_STEPS="${MWPC_K2A_LATE_STEPS:-100,50,20,10}"
T0="${MWPC_K2A_VAL_T0:-1}"; REP="${MWPC_K2A_VAL_REPEAT:-20}"; BREP="${MWPC_K2A_VAL_BIRTHS_REPEAT:-5}"
LABEL=REAL; [[ "${MWPC_K2A_VAL_SYNTHETIC:-0}" == 1 ]] && LABEL=SYNTHETIC
B="${MWPC_K2A_INLINE_BUILD:-$INL/build}"; EV="$B/k2a_strip_eval"; [[ -x "$EV" ]] || fail "missing $EV (set MWPC_K2A_INLINE_BUILD)"
[[ -n "${MWPC_BIRTHS_BIN:-}" && -x "${MWPC_BIRTHS_BIN}" ]] || fail "set MWPC_BIRTHS_BIN to the built p3_birth_moments"
[[ -n "${MWPC_K2A_PF:-}" && -s "$MWPC_K2A_PF/inventory/slots.csv" ]] || fail "set MWPC_K2A_PF (preflight with inventory/slots.csv)"
for v in MWPC_SPECTRAL_LUT MWPC_NUMW_DIR MWPC_EXPORT_BIN MWPC_ION_MOBILITY; do [[ -n "${!v:-}" && -e "${!v}" ]] || fail "set $v"; done
[[ "$T0" == 0 || "$T0" == 1 ]] || fail "MWPC_K2A_VAL_T0 must be 0 or 1"

# ---- output location
[[ ! -e "$OUT" ]] || fail "$OUT exists"
mkdir -p "$(dirname "$OUT")"; OUTP="$(cd "$(dirname "$OUT")" && pwd)/$(basename "$OUT")"
case "$OUTP/" in
  "$RUN/"*|"$ARCHIVE/"*|"$HOME/na6p/wt-validation/"*|*test_runs/mwpc_stage_c/s_ladder_2100/*) fail "OUTDIR inside a frozen reference tree";;
esac

# ---- every slot checked before anything runs
DEV=",${MWPC_K2A_VAL_DEV_SLOTS:-0000,0858,0560,0997,0133},"
SLOTS=(); declare -A INDEP=()
for I in "$@"; do
  [[ "$I" =~ ^[0-9]{1,4}$ ]] || fail "slot '$I' is not a number"
  M=$(printf %04d "$((10#$I))")
  role=$(awk -F, -v s="$M" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["slot"]==s{print $h["role"]}' "$MWPC_K2A_PF/inventory/slots.csv")
  [[ -n "$role" ]] || fail "slot $M is not in the preflight inventory"
  [[ "$role" != sealed_validation ]] && [[ "$role" != *sealed* ]] || fail "slot $M is SEALED: nothing was run"
  why=""
  [[ -s "$RUN/provenance_prepare_$M.txt" ]] && why="kernel/reference muon of RUN (provenance_prepare_$M.txt)"
  [[ "$DEV" == *",$M,"* ]] && why="${why:+$why; }used for development decisions (MWPC_K2A_VAL_DEV_SLOTS)"
  if [[ -n "$why" ]]; then
    [[ "${MWPC_K2A_VAL_ALLOW_DEV:-0}" == 1 ]] || fail "slot $M is not independent: $why (set MWPC_K2A_VAL_ALLOW_DEV=1 to run it labelled NOT independent)"
    INDEP[$M]="no: $why"
  else
    INDEP[$M]="yes"
  fi
  [[ -s "$ARCHIVE/muon_$M/ion_births.root" ]] || fail "archive muon_$M missing in $ARCHIVE"
  [[ " ${SLOTS[*]-} " == *" $M "* ]] && fail "slot $M listed twice"
  SLOTS+=("$M")
done

mkdir -p "$OUTP"; OUT="$OUTP"
CRIT="$CAMP/sealed/sealed_acceptance_criteria_v1.json"
python3 "$CAMP/sealed/check_criteria_freeze.py" > "$OUT/criteria_freeze_check.txt" 2>&1 || fail "sealed criteria v1 freeze check failed (see $OUT/criteria_freeze_check.txt)"
{ echo "label=$LABEL"; echo "date_utc=$(date -u +%FT%TZ)"; echo "git_commit=$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_dirty=$(git -C "$HERE" status --porcelain 2>/dev/null | wc -l)"
  echo "run=$RUN"; echo "tag=$TAG"; echo "kernel=$KD"; echo "archive=$ARCHIVE"; echo "half_strips=$HS"; echo "late_steps=$LATE_STEPS"
  echo "t0_reference=$T0"; echo "repeat=$REP"; echo "births_repeat=$BREP"; echo "prune_vertices=${MWPC_K2A_VAL_PRUNE_VERTICES:-0}"
  echo "births_bin=$MWPC_BIRTHS_BIN"; echo "births_bin_sha256=$(sha256sum "$MWPC_BIRTHS_BIN" | cut -d' ' -f1)"
  echo "strip_eval=$EV"; echo "strip_eval_sha256=$(sha256sum "$EV" | cut -d' ' -f1)"
  [[ -s "$B/k2a_inline_build_info.txt" ]] && sed 's/^/strip_eval_build_info: /' "$B/k2a_inline_build_info.txt"
  echo "criteria_v1_sha256=$(sha256sum "$CRIT" | cut -d' ' -f1)"
  echo "selection=${MWPC_K2A_VAL_SELECTION:-}"
  for M in "${SLOTS[@]}"; do echo "slot $M independent=${INDEP[$M]}"; done
  for f in "$HERE"/*.py "$HERE"/*.sh "$RC"/*.py "$RC"/*.sh "$TA/t_shaped_compare.py" "$TA/t_prompt_electron_term.py" \
           "$TA/t_time_collapse_export.py" "$SA/s_compare_moments.py" "$SA/s_kernel.py" "$INL/k2a_strip_compare.py"; do
    echo "tool_sha256 $(sha256sum "$f")"; done; } > "$OUT/provenance.txt"
printf '%s\n' "${SLOTS[@]}" > "$OUT/slots.txt"

# ---- helpers (each muon: V = OUT/muons/muon_M)
V=/dev/null
vstep(){ printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$V/steps.tsv"; echo "  [$2] $1: $3"; }
TF=/dev/null
tm(){ # tm LABEL CMD...: CMD unchanged; appends label, wall_s, peak_rss_kib (GNU time only), exit to $TF
  local label="$1"; shift
  local rc=0 t0 t1 rss=""
  if [[ $(type -t "$1") != function && -x /usr/bin/time ]]; then
    local tf; tf=$(mktemp)
    /usr/bin/time -f '%e %M' -o "$tf" "$@" || rc=$?
    read -r t1 rss < <(tail -1 "$tf"); rm -f "$tf"
  else
    t0=$(date +%s.%N); "$@" || rc=$?; t1=$(date +%s.%N)
    t1=$(python3 -c "print(f'{$t1 - $t0:.3f}')")
  fi
  printf '%s\t%s\t%s\t%s\n' "$label" "$t1" "$rss" "$rc" >> "$TF"
  return $rc
}
stripcmp(){ # NAME REF TEST TOL
  local rc=0
  python3 "$INL/k2a_strip_compare.py" "$2" "$3" "$V/$1.json" --tol "$4" --name "$1" > "$V/$1.log" 2>&1 || rc=$?
  case $rc in 0) vstep "$1" PASS "$(tail -1 "$V/$1.log")" "$V/$1.json";; 3) vstep "$1" INCONCLUSIVE "$(tail -1 "$V/$1.log")" "$V/$1.json";;
              *) vstep "$1" FAIL "$(tail -1 "$V/$1.log")" "$V/$1.json";; esac
}
physcmp(){ # NAME REF CAND
  local rc=0
  python3 "$RC/k2a_recenter_compare.py" "$2" "$3" --out "$V/cmp_$1.json" --kind k2a_vs_s0 --half-strips "$HS" --label "$1" > "$V/cmp_$1.txt" 2>&1 || rc=$?
  local st; st=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['status'])" "$V/cmp_$1.json" 2>/dev/null || echo FAIL)
  [[ $rc -le 3 ]] || st=FAIL
  vstep "cmp_$1" "$st" "$(head -c 300 "$V/cmp_$1.txt" | tr '\t\n' '  ')" "$V/cmp_$1.json"
}
ladder_groups(){ # GROUPS_DIR S0GROUPS_DIR PREFIX -> 0 all OK
  local GD="$1" SG="$2" P="$3" bad=0 rc g G
  mkdir -p "$SG"
  for g in "$GD"/g_*; do
    G="$SG/${g##*/}"; mkdir -p "$G"; rc=0
    tm "${P}replay_${g##*/}" python3 "$RC/k2a_timing_ladder.py" "$g" "$G" --lut "$MWPC_SPECTRAL_LUT" --numw "$MWPC_NUMW_DIR" \
        --steps "$LATE_STEPS" --half-strips "$HS" --link-replay > "$G.ladder.log" 2>&1 || rc=$?
    if [[ $rc -eq 0 ]]; then
      tm "${P}electron_step_${g##*/}" python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$g" "$G/prompt" \
          --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$G.prompt.log" 2>&1 || rc=4
    fi
    case $rc in
      0) vstep "${P}replay_${g##*/}" OK "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
      3) bad=1; vstep "${P}replay_${g##*/}" INCONCLUSIVE "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
      4) bad=1; vstep "${P}replay_${g##*/}" FAIL "electron step failed ($G.prompt.log)" "$G/timing_attempts.json";;
      *) bad=1; vstep "${P}replay_${g##*/}" FAIL "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
    esac
  done
  return $bad
}
shaped(){ # NAME VARIANTS [--no-prompt]  (layout OUT/muons/muon_M/shaped_layout/muon_M/<variant>)
  local n="$1" vars="$2"; shift 2
  local rc=0
  python3 "$TA/t_shaped_compare.py" "$V/shaped_layout" --variants "$vars" --out "$V/shaped_$n" --bootstrap 0 "$@" > "$V/shaped_$n.log" 2>&1 || rc=$?
  if [[ $rc -eq 0 ]]; then vstep "shaped_$n" OK "variants $vars ${*:-ion+electron}" "$V/shaped_$n/per_muon_Tp100.jsonl"
  else vstep "shaped_$n" FAIL "t_shaped_compare exit $rc: $(tail -1 "$V/shaped_$n.log" | head -c 250)" "$V/shaped_$n.log"; fi
}

# ---- muons
for M in "${SLOTS[@]}"; do
  V="$OUT/muons/muon_$M"; mkdir -p "$V"; printf 'step\tstatus\tdetail\tfile\n' > "$V/steps.tsv"
  TF="$V/timing.tsv"; printf 'label\twall_s\tpeak_rss_kib\texit\n' > "$TF"
  echo "muon $M (independent: ${INDEP[$M]}):"
  avail=$(( $(df -Pk "$OUT" | awk 'NR==2{print $4}') >> 20 )); minf="${MWPC_K2A_VAL_MIN_FREE_GIB:-24}"
  if [[ $avail -lt $minf ]]; then
    vstep disk INCONCLUSIVE "only $avail GiB free in $OUT (need >= $minf GiB; MWPC_K2A_VAL_MIN_FREE_GIB); muon not run"
  else
    sha256sum "$ARCHIVE/muon_$M/ion_births.root" > "$V/archive_sha256_before.txt"
    # 1. S0_T1 reference + Python K2A_local (frozen far-mode experiment, instrumented, uncollapsed export kept)
    rc=0
    env MWPC_K2A_OUT="$OUT/far/muon_$M" MWPC_K2A_KEEP_FROZEN_VERTICES="$T0" MWPC_ARCHIVE_2100_K256="$ARCHIVE" \
        bash "$RC/k2a_recenter_experiment.sh" far "$RUN" "$TAG" "$M" > "$V/far.log" 2>&1 || rc=$?
    X=$( (ls -d "$OUT/far/muon_$M"/recenter_experiment_* 2>/dev/null || true) | tail -1); E="$X/muon_$M"
    echo "$E" > "$V/far_dir.txt"
    if [[ -z "$X" || ! -s "$E/steps.tsv" ]]; then
      vstep far_experiment FAIL "far experiment aborted before muon $M (exit $rc): $(grep -m1 ABORTED "$V/far.log" | head -c 250)" "$V/far.log"
    else
      vstep far_experiment "$(case $rc in 0) echo PASS;; 1) echo FAIL;; 2) echo WARN;; 3) echo INCONCLUSIVE;; *) echo ABORTED;; esac)" "aggregate exit $rc" "$X/results.json"
      R="$E/prepared"
      KID=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['kernel_id'])" "$X/env_check.json")
      HAVE_T1=0; [[ -s "$E/S0_local/replay/waveforms.csv" && -s "$E/S0_local/prompt/prompt_waveforms.csv" ]] && HAVE_T1=1
      # 2. S0_T0 reference (uncollapsed frozen export, same groups)
      HAVE_T0=0
      if [[ $T0 == 1 ]]; then
        if [[ -s "$R/frozen_export/vertices.f64" && -s "$E/groups/groups.json" ]]; then
          rc=0; tm T0_group_export python3 "$RC/k2a_group_export.py" "$R/frozen_export" "$V/T0/groups" --mode primary > "$V/T0_groups.log" 2>&1 || rc=$?
          if [[ $rc -ne 0 ]]; then vstep T0_groups FAIL "group export of the uncollapsed export failed ($V/T0_groups.log)"
          elif python3 -c "
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:3])
sys.exit(0 if a['primary_shift'] == b['primary_shift'] and sorted(a['groups']) == sorted(b['groups']) else 1)" "$E/groups/groups.json" "$V/T0/groups/groups.json"; then
            vstep T0_groups OK "same primaries -> (m_U, m_V) as S0_T1: $(tail -1 "$V/T0_groups.log" | head -c 200)" "$V/T0/groups/groups.json"
            if ladder_groups "$V/T0/groups" "$V/T0/S0groups" T0_; then
              rc=0; tm T0_local_sum python3 "$RC/k2a_local_sum.py" "$V/T0/groups" "$V/T0/S0groups" "$V/T0/S0_local" > "$V/T0_local_sum.log" 2>&1 || rc=$?
              vstep T0_S0_local_census "$([[ $rc -eq 0 ]] && echo OK || echo FAIL)" "$(tail -1 "$V/T0_local_sum.log")" "$V/T0/S0_local/local_sum.json"
              [[ $rc -eq 0 ]] && HAVE_T0=1
            fi
          else
            vstep T0_groups FAIL "uncollapsed export grouped differently from S0_T1 (primary_shift differs)" "$V/T0/groups/groups.json"
          fi
        else
          vstep T0_groups FAIL "no uncollapsed frozen export or no S0_T1 groups (far experiment incomplete)"
        fi
      fi
      # 3. production K2A: archive -> compact moments -> compiled kernel evaluation
      mkdir -p "$V/prod"; rc=0
      tm births_moments "$MWPC_BIRTHS_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$V/prod/birth_moments.csv" --repeat "$BREP" > "$V/prod/births.log" 2>&1 || rc=$?
      if [[ $rc -eq 0 ]] && grep -q '^PASS birth moments' "$V/prod/births.log"; then
        vstep births_moments OK "$(grep -m1 '^PASS birth moments' "$V/prod/births.log")" "$V/prod/birth_moments.csv.json"
        if [[ -s "$R/cloud_moments.csv" ]]; then
          rc=0; python3 "$SA/s_compare_moments.py" "$R/cloud_moments.csv" "$V/prod/birth_moments.csv" --report "$V/prod/compare_moments.md" > "$V/prod/compare_moments.log" 2>&1 || rc=$?
          vstep moments_archive_vs_export "$([[ $rc -eq 0 ]] && echo PASS || echo FAIL)" "$(tail -1 "$V/prod/compare_moments.log" | head -c 250)" "$V/prod/compare_moments.md"
        else
          vstep moments_archive_vs_export INCONCLUSIVE "no export-path moments (far experiment incomplete)"
        fi
        if [[ -s "$E/groups/groups.json" ]]; then
          rc=0; tm k2a_eval "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$V/prod/birth_moments.csv" --shifts "$E/groups/groups.json" \
                  --half-strips "$HS" --out "$V/K2A_prod" --repeat "$REP" > "$V/K2A_prod.log" 2>&1 || rc=$?
          case $rc in
            0) vstep K2A_prod OK "$(tail -1 "$V/K2A_prod.log" | head -c 300)" "$V/K2A_prod/result.json";;
            4) vstep K2A_prod BLOCKED "$(tail -1 "$V/K2A_prod.log" | head -c 300)" "$V/K2A_prod/result.json";;
            *) vstep K2A_prod FAIL "k2a_strip_eval exit $rc: $(tail -1 "$V/K2A_prod.log" | head -c 250)" "$V/K2A_prod.log";;
          esac
          rc=0; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$V/prod/birth_moments.csv" --half-strips "$HS" --out "$V/K2A_rule" > "$V/K2A_rule.log" 2>&1 || rc=$?
          # production rule (no imposed grouping): status, groups; how many primaries the S0 grouping assigns differently
          if [[ -s "$V/K2A_rule/result.json" ]]; then
            vstep K2A_rule INFO "$(python3 -c "
import json, sys
r = json.load(open(sys.argv[1])); p = sys.argv[2]
d = json.load(open(p))['shifts_differing_from_computed'] if p and __import__('os').path.isfile(p) else 'n/a'
print(f\"production rule: status {r['status_names']}, {r['groups']} groups, usable {r['usable']}; S0 grouping differs from the rule for {d} primaries\")" \
              "$V/K2A_rule/result.json" "$V/K2A_prod/result.json")" "$V/K2A_rule/result.json"
          else vstep K2A_rule INFO "production rule run exit $rc, no result ($V/K2A_rule.log)"; fi
          rc=0; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$R/cloud_moments.csv" --shifts "$E/groups/groups.json" --half-strips "$HS" \
                  --out "$V/K2A_cc_export" > "$V/K2A_cc_export.log" 2>&1 || rc=$?
          vstep K2A_cc_export "$(case $rc in 0) echo OK;; 4) echo BLOCKED;; *) echo FAIL;; esac)" "$(tail -1 "$V/K2A_cc_export.log" | head -c 300)" "$V/K2A_cc_export/result.json"
          if [[ -s "$E/K2A_local/replay/waveforms.csv" && -s "$V/K2A_cc_export/replay/waveforms.csv" ]]; then
            stripcmp identity_cc_vs_python "$E/K2A_local" "$V/K2A_cc_export" 1e-12
          fi
          if [[ -s "$V/K2A_cc_export/replay/waveforms.csv" && -s "$V/K2A_prod/replay/waveforms.csv" ]]; then
            stripcmp moment_path_archive_vs_export "$V/K2A_cc_export" "$V/K2A_prod" 1e-7
          fi
        else
          vstep K2A_prod FAIL "no S0 grouping (groups.json) to evaluate with: far experiment incomplete"
        fi
      else
        vstep births_moments FAIL "p3_birth_moments exit $rc: $(tail -1 "$V/prod/births.log" | head -c 250)" "$V/prod/births.log"
      fi
      # 4. physics comparisons (K2A production response vs the detailed references)
      HAVE_K=0; [[ -s "$V/K2A_prod/replay/waveforms.csv" ]] && HAVE_K=1
      [[ $HAVE_T1 == 1 && $HAVE_K == 1 ]] && physcmp S0T1_vs_K2A "$E/S0_local" "$V/K2A_prod"
      [[ $HAVE_T0 == 1 && $HAVE_K == 1 ]] && physcmp S0T0_vs_K2A "$V/T0/S0_local" "$V/K2A_prod"
      [[ $HAVE_T0 == 1 && $HAVE_T1 == 1 ]] && physcmp S0T0_vs_S0T1 "$V/T0/S0_local" "$E/S0_local"
      # 5. shaped observables (layout of symlinks; nothing copied)
      L="$V/shaped_layout/muon_$M"; mkdir -p "$L"
      [[ $HAVE_T1 == 1 ]] && ln -s "$E/S0_local" "$L/S0T1"
      [[ $HAVE_T0 == 1 ]] && ln -s "$V/T0/S0_local" "$L/S0T0"
      [[ $HAVE_K == 1 ]] && ln -s "$V/K2A_prod" "$L/K2A"
      if [[ $HAVE_T1 == 1 && $HAVE_K == 1 ]]; then
        if [[ $HAVE_T0 == 1 ]]; then
          shaped refS0T1_ionelec S0T1,K2A,S0T0; shaped refS0T0_ionelec S0T0,K2A,S0T1; shaped refS0T1_ion S0T1,K2A,S0T0 --no-prompt
        else
          shaped refS0T1_ionelec S0T1,K2A; shaped refS0T1_ion S0T1,K2A --no-prompt
        fi
      fi
    fi
    sha256sum "$ARCHIVE/muon_$M/ion_births.root" > "$V/archive_sha256_after.txt"
  fi
  # 6. explicit outcome, runtime split, criteria dry run
  python3 "$HERE/k2a_validation_muon.py" "$V" --muon "$M" --label "$LABEL" --independent "${INDEP[$M]}" --t0 "$T0" \
      --criteria "$CRIT" > "$V/outcome.log" 2>&1 || echo "  (k2a_validation_muon.py exit $?; see $V/outcome.log)"
  tail -1 "$V/outcome.log"
  if [[ "${MWPC_K2A_VAL_PRUNE_VERTICES:-0}" == 1 && -s "$V/outcome.json" ]]; then
    E=$(cat "$V/far_dir.txt" 2>/dev/null || true)
    for f in ${E:+"$E/prepared/frozen_export/vertices.f64"} ${E:+"$E"/groups/g_*/vertices.f64} "$V"/T0/groups/g_*/vertices.f64; do
      case "$f" in "$OUT"/*) ;; *) continue;; esac                       # never anything outside this run's OUTDIR
      [[ -f "$f" ]] || continue
      sha256sum "$f" >> "$V/pruned_vertices.sha256" && rm -f "$f"
    done
  fi
done
rc=0; python3 "$HERE/k2a_validation_report.py" "$OUT" ${MWPC_K2A_VAL_SELECTION:+--selection "$MWPC_K2A_VAL_SELECTION"} || rc=$?
echo "return for inspection: $OUT/validation_report.md $OUT/validation_report.json $OUT/per_muon.csv $OUT/muons/muon_*/outcome.json"
exit $rc

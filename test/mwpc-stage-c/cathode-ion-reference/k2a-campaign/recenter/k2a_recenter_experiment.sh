#!/usr/bin/env bash
# Reference experiment: does per-primary strip recentering + absolute-strip accumulation reproduce the FROZEN S0
# weighting response? Design: ../RECENTERING_DESIGN_REVIEW.md section 5. READ-ONLY on the reference run, archive, LUT:
# all outputs go to a NEW directory $MWPC_K2A_OUT/recenter_experiment_<timestamp>/; reference files are hashed before
# and after and the aggregate fails if any changed.
#
#   bash k2a_recenter_experiment.sh envcheck RUN TAG MUON ...     environment vs the validated run's own records
#   bash k2a_recenter_experiment.sh refs     RUN TAG MUON ...     E1 + E2a on development reference muons
#   bash k2a_recenter_experiment.sh far      RUN TAG INDEX ...    E2b on NON-SEALED far-D muons (frozen export + T1 into
#                                                                 the new directory; needs MWPC_EXPORT_BIN,
#                                                                 MWPC_ION_MOBILITY, MWPC_K2A_PF)
# E1  symmetry: whole T1 export relabelled by m = +1 / -1 strips, frozen replay + electron step, compared with the
#     frozen S0 on the overlapping strips: PASS / FAIL / INCONCLUSIVE (replay refusal = INCONCLUSIVE, with its reason)
# E2  per-primary groups -> frozen replay + electron step per group -> S0_local (absolute strips, census checked);
#     K2A_local with the unchanged engine and the validated kernel (kernel id checked); comparisons with documented
#     warning thresholds (thresholds.json): OK / WARN / FAIL / INCONCLUSIVE
# Every step is recorded in muon_*/steps.tsv; k2a_recenter_aggregate.py writes results.json / results.md: overall PASS
# only if every required step and report of the mode is present and none FAILs or is INCONCLUSIVE (exit 0/1/3).
# Memory: exports are streamed (k2a_group_export.py); disk: ~3 x T1_export per reference muon (checked before start).
# Environment (from the validated setup, e.g. `source ~/t_ladder_deps.env`; verified, not assumed):
#   MWPC_SPECTRAL_LUT  sha256 must equal RUN/green_TAG/axes.json lut_sha256 and RUN/provenance_prepare_<M>.txt
#   MWPC_NUMW_DIR      p3_replay_full_spectral_muon.py sha256 must equal replay_sha256 in provenance_prepare_<M>.txt
#   far mode also: MWPC_EXPORT_BIN (export_bin_sha256), MWPC_ION_MOBILITY (mobility_sha256), ROOT/Garfield runtime
#   MWPC_K2A_LATE_STEPS (far mode group replays: frozen --late-step-ns ladder, default 100,50,20,10; only the frozen
#   time-compression rejection is retried, see k2a_timing_ladder.py)
#   MWPC_K2A_OUT (output root), MWPC_HALF_STRIPS (5), MWPC_K2A_KERNEL_ID (default: latest RUN/engine_TAG_*/kernel_id.txt)
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; TA="$HERE/../../timing-ablation"; SA="$HERE/../../spatial-ablation"
REPO="$(cd "$HERE/../../../../.." && pwd)"
set -u
umask 077
HS="${MWPC_HALF_STRIPS:-5}"; LATE_STEPS="${MWPC_K2A_LATE_STEPS:-100,50,20,10}"; OUTROOT="${MWPC_K2A_OUT:-$REPO/test_runs/mwpc_stage_c/k2a_campaign_2100}"
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
fail(){ echo "RECENTER EXPERIMENT ABORTED: $*" >&2; exit 4; }
need(){ [[ -n "${!1:-}" && -e "${!1}" ]] || fail "set $1 to an existing path (validated setup, e.g. source ~/t_ladder_deps.env)"; }
sha(){ sha256sum "$1" | cut -d' ' -f1; }
prov(){ sed -n "s/^$2=//p" "$1" | head -1; }

MODE="${1:-}"; RUN="${2:?RUN}"; TAG="${3:?TAG}"; shift 3
KD="$RUN/kernel_$TAG"; [[ -s "$KD/kernel.json" ]] || fail "missing $KD"
[[ $# -ge 1 ]] || fail "list the muons explicitly"
need MWPC_SPECTRAL_LUT; need MWPC_NUMW_DIR
[[ $MODE != far ]] || { need MWPC_EXPORT_BIN; need MWPC_ION_MOBILITY; need MWPC_K2A_PF; }
(cd "$MWPC_NUMW_DIR" && python3 -c "import p3_replay_full_spectral_muon") > /dev/null 2>&1 || fail "replay module does not import from MWPC_NUMW_DIR"

envcheck(){ # writes $1 (json) ; muons in remaining args; aborts on mismatch
  local OUTJ="$1"; shift
  local L; L=$(sha "$MWPC_SPECTRAL_LUT"); local RP; RP=$(sha "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py")
  local AX; AX=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['lut_sha256'])" "$RUN/green_$TAG/axes.json")
  [[ "$L" == "$AX" ]] || fail "MWPC_SPECTRAL_LUT sha256 $L != green_$TAG lut_sha256 $AX"
  local KID; KID="${MWPC_K2A_KERNEL_ID:-}"
  if [[ -z "$KID" ]]; then
    local V; V=$( (ls -d "$RUN"/engine_"${TAG}"_* 2>/dev/null || true) | sort | tail -1)
    [[ -n "$V" && -s "$V/kernel_id.txt" ]] || fail "no validated engine_${TAG}_*/kernel_id.txt; set MWPC_K2A_KERNEL_ID"
    KID=$(sed -n 's/^kernel_id=//p' "$V/kernel_id.txt")
  fi
  local rows=""; EXTRA=""
  for M in "$@"; do
    M=$(printf %04d "$((10#$M))"); local P="$RUN/provenance_prepare_$M.txt"
    if [[ -s "$P" ]]; then
      [[ "$(prov "$P" lut_sha256)" == "$L" ]] || fail "muon $M: LUT differs from the one its S0 was made with ($P)"
      [[ "$(prov "$P" replay_sha256)" == "$RP" ]] || fail "muon $M: replay code differs from the one its S0 was made with ($P)"
      if [[ $MODE == far ]]; then
        [[ "$(prov "$P" export_bin_sha256)" == "$(sha "$MWPC_EXPORT_BIN")" ]] || fail "exporter differs from $P"
        [[ "$(prov "$P" mobility_sha256)" == "$(sha "$MWPC_ION_MOBILITY")" ]] || fail "mobility file differs from $P"
      fi
      rows+="\"$M\": \"$P\","
    elif [[ $MODE == far ]]; then
      rows+="\"$M\": \"(new muon: checked against the reference provenance)\","
    else
      fail "muon $M: no $P (not a prepared reference of this run)"
    fi
  done
  if [[ $MODE == far ]]; then
    local XB; XB=$(sha "$MWPC_EXPORT_BIN"); local MB; MB=$(sha "$MWPC_ION_MOBILITY"); local nref=0
    for REFP in "$RUN"/provenance_prepare_*.txt; do                 # every reference must agree with the environment
      [[ -s "$REFP" ]] || continue; nref=$((nref + 1))
      [[ "$(prov "$REFP" export_bin_sha256)" == "$XB" ]] || fail "MWPC_EXPORT_BIN sha256 differs from the exporter recorded in $REFP"
      [[ "$(prov "$REFP" mobility_sha256)" == "$MB" ]] || fail "MWPC_ION_MOBILITY sha256 differs from the mobility recorded in $REFP"
      [[ "$(prov "$REFP" lut_sha256)" == "$L" ]] || fail "LUT differs from $REFP"
    done
    [[ $nref -ge 1 ]] || fail "no provenance_prepare_*.txt in $RUN to verify the exporter against"
    EXTRA="\"far_provenance\": {\"export_bin_sha256\": \"$XB\", \"mobility_sha256\": \"$MB\", \"reference_files_checked\": $nref}, "
    local usage; usage=$("$MWPC_EXPORT_BIN" 2>&1 || true)
    [[ $usage == *"usage: PROGRAM frozen.root MOBILITY NEW_DIR"* ]] || fail "MWPC_EXPORT_BIN is not the full-spectral exporter or does not start (ROOT/Garfield runtime?): ${usage:0:200}"
  fi
  local K; K=$(python3 - "$KD" <<'PY'
import hashlib, sys
from pathlib import Path
d = Path(sys.argv[1]); h = hashlib.sha256()
for f in ("K.npy", "KPHI0.npy", "kernel.json"):
    s = hashlib.sha256()
    with open(d / f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            s.update(b)
    h.update(s.hexdigest().encode())
print(h.hexdigest())
PY
)
  [[ "$K" == "$KID" ]] || fail "kernel $KD id $K != validated $KID"
  python3 -c "import numpy, sys; print('numpy', numpy.__version__, 'python', sys.version.split()[0])" > /dev/null || fail "python/numpy"
  printf '{"mode": "%s", "run": "%s", "tag": "%s", "lut": "%s", "lut_sha256": "%s", "numw_dir": "%s", "replay_sha256": "%s", "kernel": "%s", "kernel_id": "%s", "half_strips": %s, %s"provenance": {%s}, "status": "PASS"}\n' \
    "$MODE" "$RUN" "$TAG" "$MWPC_SPECTRAL_LUT" "$L" "$MWPC_NUMW_DIR" "$RP" "$KD" "$KID" "$HS" "$EXTRA" "${rows%,}" > "$OUTJ"
  echo "PASS environment: LUT $L, replay $RP, kernel $KID match the validated run"
  KERNEL_ID="$KID"
}

step(){ printf '%s\t%s\t%s\t%s\n' "$2" "$3" "$4" "${5:-}" >> "$1/steps.tsv"; echo "  [$3] $2: $4"; }

replay(){ # EXPORT OUTDIR -> 0 ok, 3 refused (window/phase range), 1 other failure
  local rc=0
  (cd "$MWPC_NUMW_DIR" && python3 p3_replay_full_spectral_muon.py "$MWPC_SPECTRAL_LUT" "$1" "$2/replay" --half-strips "$HS") > "$2.replay.log" 2>&1 || rc=$?
  if [[ $rc -ne 0 ]]; then
    grep -q -E "selected channels omit the landing strip|query outside declared strip phase range|query outside validated cross-section" "$2.replay.log" && return 3
    return 1
  fi
  python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$1" "$2/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$2.prompt.log" 2>&1 || {
    grep -q -E "selected channels omit|query outside declared strip phase range|query outside validated cross-section" "$2.prompt.log" && return 3; return 1; }
  return 0
}
why(){ grep -h -m1 -o -E "selected channels omit the landing strip[^\"]*|query outside declared strip phase range|query outside validated cross-section|[A-Za-z]*Error: .*" "$@" 2>/dev/null | head -1; }

compare(){ # E NAME REF CAND KIND [extra]
  local E="$1" N="$2" R="$3" C="$4" K="$5"; shift 5
  local rc=0
  python3 "$HERE/k2a_recenter_compare.py" "$R" "$C" --out "$E/cmp_$N.json" --kind "$K" --half-strips "$HS" --label "$N" "$@" > "$E/cmp_$N.txt" 2>&1 || rc=$?
  local st; st=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['status'])" "$E/cmp_$N.json" 2>/dev/null || echo FAIL)
  [[ $rc -le 3 ]] || st=FAIL
  step "$E" "$N" "$st" "$(head -c 300 "$E/cmp_$N.txt" | tr '\t\n' '  ')" "$E/cmp_$N.json"
}

e1(){ # R E
  local R="$1" E="$2"
  for m in 1 -1; do
    python3 "$HERE/k2a_group_export.py" "$R/T1_export" "$E/shift_$m" --mode shift --m-U "$m" --m-V "$m" > "$E/shift_$m.log" 2>&1 \
      || { step "$E" "E1_shift_$m" FAIL "relabelled export not written ($E/shift_$m.log)"; continue; }
    local rc=0; replay "$E/shift_$m/g_${m}_${m}" "$E/S0shift_$m" || rc=$?
    if [[ $rc -eq 3 ]]; then step "$E" "E1_shift_$m" INCONCLUSIVE "frozen replay refused: $(why "$E/S0shift_$m.replay.log" "$E/S0shift_$m.prompt.log")"; continue
    elif [[ $rc -ne 0 ]]; then step "$E" "E1_shift_$m" FAIL "frozen replay failed: $(why "$E/S0shift_$m.replay.log")"; continue; fi
    compare "$E" "E1_shift_$m" "$R/S0" "$E/S0shift_$m" e1 --shift "$m"
  done
}

e2(){ # R E HAVE_BASE_S0
  local R="$1" E="$2" HAVE="$3" rc
  python3 "$HERE/k2a_group_export.py" "$R/T1_export" "$E/groups" --mode primary > "$E/groups.log" 2>&1 \
    || { step "$E" E2_groups FAIL "group export failed ($E/groups.log)"; return; }
  step "$E" E2_groups OK "$(tail -1 "$E/groups.log")" "$E/groups/groups.json"
  mkdir -p "$E/S0groups"; local bad=0
  for g in "$E"/groups/g_*; do
    if [[ $MODE == far ]]; then
      # far mode: frozen replay with the late-step convergence ladder (retries ONLY the frozen time-compression
      # rejection, each attempt in a new directory; replay -> chosen attempt is a symlink), then the electron step
      local G="$E/S0groups/${g##*/}"; mkdir -p "$G"; rc=0
      python3 "$HERE/k2a_timing_ladder.py" "$g" "$G" --lut "$MWPC_SPECTRAL_LUT" --numw "$MWPC_NUMW_DIR" \
          --steps "$LATE_STEPS" --half-strips "$HS" --link-replay > "$G.ladder.log" 2>&1 || rc=$?
      if [[ $rc -eq 0 ]]; then
        python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$g" "$G/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$G.prompt.log" 2>&1 || rc=4
      fi
      case $rc in
        0) step "$E" "E2_replay_${g##*/}" OK "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
        3) bad=1; step "$E" "E2_replay_${g##*/}" INCONCLUSIVE "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
        4) bad=1; step "$E" "E2_replay_${g##*/}" FAIL "electron step failed ($G.prompt.log)" "$G/timing_attempts.json";;
        *) bad=1; step "$E" "E2_replay_${g##*/}" FAIL "$(tail -1 "$G.ladder.log" | head -c 300)" "$G/timing_attempts.json";;
      esac
      continue
    fi
    rc=0; replay "$g" "$E/S0groups/${g##*/}" || rc=$?
    if [[ $rc -ne 0 ]]; then bad=1; step "$E" "E2_replay_${g##*/}" "$([[ $rc -eq 3 ]] && echo INCONCLUSIVE || echo FAIL)" "$(why "$E/S0groups/${g##*/}.replay.log")"; fi
  done
  [[ $bad -eq 0 ]] || return
  rc=0; python3 "$HERE/k2a_local_sum.py" "$E/groups" "$E/S0groups" "$E/S0_local" > "$E/local_sum.log" 2>&1 || rc=$?
  step "$E" E2_S0_local_census "$([[ $rc -eq 0 ]] && echo OK || echo FAIL)" "$(tail -1 "$E/local_sum.log")" "$E/S0_local/local_sum.json"
  [[ $rc -eq 0 ]] || return
  local PG=(); [[ $MODE == far ]] && PG=(--per-group-out "$E/K2A_groups")
  rc=0; python3 "$HERE/k2a_recentered_k2a.py" "$KD" "$R/cloud_moments.csv" "$E/K2A_local" --groups "$E/groups/groups.json" \
        --half-strips "$HS" --kernel-id "$KERNEL_ID" "${PG[@]}" > "$E/k2a_local.log" 2>&1 || rc=$?
  step "$E" E2_K2A_local "$([[ $rc -eq 0 ]] && echo OK || echo FAIL)" "$(tail -1 "$E/k2a_local.log")" "$E/K2A_local/recentered.json"
  [[ $rc -eq 0 ]] || return
  local single=0; [[ $(ls -d "$E"/groups/g_* | wc -l) -eq 1 && -d "$E/groups/g_0_0" ]] && single=1
  compare "$E" E2_K2Alocal_vs_S0local "$E/S0_local" "$E/K2A_local" k2a_vs_s0
  if [[ $MODE == far ]]; then
    rc=0; python3 "$HERE/k2a_truncation_check.py" "$E/groups" "$MWPC_SPECTRAL_LUT" "$MWPC_NUMW_DIR" "$E/truncation.json" --half-strips "$HS" > "$E/truncation.log" 2>&1 || rc=$?
    step "$E" E2_truncation "$([[ $rc -eq 0 ]] && echo OK || echo FAIL)" "$(tr '\n' ' ' < "$E/truncation.log" | head -c 300)" "$E/truncation.json"
    rc=0; python3 "$HERE/k2a_far_analysis.py" "$E" --kernel "$KD" --half-strips "$HS" > "$E/far_analysis.log" 2>&1 || rc=$?
    step "$E" E2_far_analysis "$(case $rc in 0) echo OK;; 2) echo WARN;; *) echo FAIL;; esac)" "$(tail -1 "$E/far_analysis.log" | head -c 300)" "$E/far_analysis.json"
  fi
  if [[ $HAVE == 1 ]]; then
    compare "$E" E2_S0local_vs_S0 "$R/S0" "$E/S0_local" "$([[ $single == 1 ]] && echo identity || echo window)"
    compare "$E" E2_K2Alocal_vs_S0 "$R/S0" "$E/K2A_local" k2a_vs_s0
    if [[ $MODE == far ]]; then :                                                         # no base-centred K2A for new muons
    elif [[ -d "$R/K2Akern_$TAG" ]]; then
      compare "$E" E2_K2Abase_vs_S0 "$R/S0" "$R/K2Akern_$TAG" k2a_vs_s0                 # the validated baseline, for context
      compare "$E" E2_K2Alocal_vs_K2Abase "$R/K2Akern_$TAG" "$E/K2A_local" "$([[ $single == 1 ]] && echo identity || echo window)"
    else
      step "$E" E2_K2Alocal_vs_K2Abase INCONCLUSIVE "no $R/K2Akern_$TAG"
    fi
  fi
}

refhash(){ # R OUT  -> sha256 of every reference file read (never written)
  ( cd "$1" && find -L T1_export S0 cloud_moments.csv "K2Akern_$TAG" -type f 2>/dev/null | sort | xargs -r sha256sum ) > "$2"
}

X="$OUTROOT/recenter_experiment_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$X"
case "$X" in "$RUN"*) fail "output inside the reference run";; esac
git -C "$REPO" rev-parse HEAD > "$X/SOURCE_GIT_SHA.txt" 2>/dev/null || true
cp "$HERE/thresholds.json" "$X/thresholds.json"
KERNEL_ID=""
envcheck "$X/env_check.json" "$@"
[[ $MODE == envcheck ]] && { echo "environment OK: $X/env_check.json"; exit 0; }
case "$MODE" in
  refs)
    for M in "$@"; do
      M=$(printf %04d "$((10#$M))"); R="$RUN/muon_$M"; E="$X/muon_$M"; mkdir -p "$E"; : > "$E/steps.tsv"
      [[ -s "$R/T1_export/vertices.f64" && -s "$R/S0/replay/waveforms.csv" && -s "$R/S0/prompt/prompt_waveforms.csv" && -s "$R/cloud_moments.csv" ]] \
        || fail "$R is not a prepared reference (T1_export, S0 replay+prompt, cloud_moments.csv)"
      need_b=$(( 3 * $(stat -c %s "$R/T1_export/vertices.f64") + (512 << 20) ))         # 2 relabelled copies + groups + outputs
      avail=$(( $(df -Pk "$X" | awk 'NR==2{print $4}') * 1024 ))
      [[ $avail -gt $need_b ]] || fail "muon $M needs ~$(( need_b >> 30 )) GiB free in $X, only $(( avail >> 30 )) GiB available"
      refhash "$R" "$E/reference_sha256_before.txt"
      echo "muon $M:"; e1 "$R" "$E"; e2 "$R" "$E" 1
      refhash "$R" "$E/reference_sha256_after.txt"
    done;;
  far)
    for I in "$@"; do
      M=$(printf %04d "$((10#$I))")
      role=$(awk -F, -v s="$M" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["slot"]==s{print $h["role"]}' "$MWPC_K2A_PF/inventory/slots.csv")
      [[ -n "$role" && "$role" != sealed_validation ]] || fail "$M is sealed or unknown (role '$role')"
      E="$X/muon_$M"; R="$E/prepared"; mkdir -p "$R"; : > "$E/steps.tsv"; echo "muon $M ($role):"
      avail=$(( $(df -Pk "$X" | awk 'NR==2{print $4}') >> 20 )); minf="${MWPC_K2A_FAR_MIN_FREE_GIB:-8}"
      [[ $avail -ge $minf ]] || fail "muon $M: only $avail GiB free in $X (need >= $minf GiB; set MWPC_K2A_FAR_MIN_FREE_GIB)"
      [[ -s "$ARCHIVE/muon_$M/ion_births.root" && -s "$ARCHIVE/muon_$M/SHA256SUMS" ]] || fail "archive muon_$M missing"
      (cd "$ARCHIVE/muon_$M" && sha256sum -c --quiet SHA256SUMS) || fail "archive muon_$M fails its SHA256SUMS"
      sha256sum "$ARCHIVE/muon_$M/ion_births.root" > "$E/archive_sha256_before.txt"
      "$MWPC_EXPORT_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$MWPC_ION_MOBILITY" "$R/frozen_export" > "$R/export.log" 2>&1 \
        && grep -q '^PASS archive census' "$R/export.log" || { step "$E" export FAIL "exporter ($R/export.log)"; continue; }
      python3 "$TA/t_time_collapse_export.py" "$R/frozen_export" "$R/T1_export" --mode primary > "$R/T1.log" 2>&1 || { step "$E" T1 FAIL "$R/T1.log"; continue; }
      rm -f "$R/frozen_export/vertices.f64"          # our own new copy; the T1 export keeps the paths
      python3 "$SA/s_kernel.py" moments "$R/T1_export" "$R/cloud_moments.csv" > "$R/moments.log" 2>&1 || { step "$E" moments FAIL "$R/moments.log"; continue; }
      rc=0; replay "$R/T1_export" "$R/S0" || rc=$?
      case $rc in
        0) HAVE=1; step "$E" base_centred_S0 OK "frozen base-centred S0 ran";;
        3) HAVE=0; step "$E" base_centred_S0 EXPECTED_REFUSAL "frozen base-centred S0 refused, as expected for far D (not a failure): $(why "$R/S0.replay.log" "$R/S0.prompt.log")";;
        *) HAVE=0; step "$E" base_centred_S0 FAIL "frozen base-centred S0 failed: $(why "$R/S0.replay.log")";;
      esac
      e2 "$R" "$E" "$HAVE"
      sha256sum "$ARCHIVE/muon_$M/ion_births.root" > "$E/archive_sha256_after.txt"
    done;;
  *) fail "mode $MODE";;
esac
rc=0; python3 "$HERE/k2a_recenter_aggregate.py" "$X" || rc=$?
echo "return for inspection: $X/results.md $X/results.json $X/muon_*/steps.tsv $X/muon_*/cmp_*.json"
exit $rc

#!/usr/bin/env bash
# OPTIONAL diagnostic reference: rebuild one muon's S0 at the frozen replay's NATIVE 10 ns step (no late-step time
# compression), with the same frozen code, LUT, groups and electron step. It tells whether a K2A-vs-S0 residual is S0's
# own time-compression error or K2A's. Weighting replay only: no avalanche, no ion transport rerun (the RKF paths of the
# existing groups are reused). Expensive (the native step evaluates every 10 ns bin); run on the worst muon only, with
# approval.
#
#   bash k2a_s0_native_step.sh MUON_VALIDATION_DIR NEW_OUT [T1|T0]
#
# MUON_VALIDATION_DIR  OUT/muons/muon_M of run_k2a_validation.sh (far_dir.txt; T0/groups for T0)
# NEW_OUT              new directory: S0groups/g_*/ (ladder with --steps 10 only), S0_local/, timing.tsv, provenance.txt
# The group exports must still hold their vertices.f64 (not pruned). Use the result with
#   k2a_residual_diagnostics.py DIAG --from RUN --native S0T1native=M:NEW_OUT/S0_local
# env: MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR (as the run), MWPC_HALF_STRIPS (5)
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CAMP="$(cd "$HERE/.." && pwd)"; CIR="$(cd "$CAMP/.." && pwd)"
RC="$CAMP/recenter"; TA="$CIR/timing-ablation"; INL="$CIR/k2a-inline"
set -u
V="${1:?usage: k2a_s0_native_step.sh MUON_VALIDATION_DIR NEW_OUT [T1|T0]}"; OUT="${2:?NEW_OUT}"; REF="${3:-T1}"
[[ ! -e "$OUT" ]] || { echo "refusing: $OUT exists" >&2; exit 2; }
for v in MWPC_SPECTRAL_LUT MWPC_NUMW_DIR; do [[ -n "${!v:-}" && -e "${!v}" ]] || { echo "set $v" >&2; exit 2; }; done
HS="${MWPC_HALF_STRIPS:-5}"
E=$(cat "$V/far_dir.txt")
case "$REF" in T1) GD="$E/groups";; T0) GD="$V/T0/groups";; *) echo "REF must be T1 or T0" >&2; exit 2;; esac
[[ -s "$GD/groups.json" ]] || { echo "no groups in $GD" >&2; exit 2; }
for g in "$GD"/g_*; do [[ -s "$g/vertices.f64" ]] || { echo "$g/vertices.f64 missing (pruned?): rerun the muon without pruning" >&2; exit 2; }; done
umask 077; mkdir -p "$OUT/S0groups"; OUT="$(cd "$OUT" && pwd)"
V=/dev/null; TF="$OUT/timing.tsv"; printf 'label\twall_s\tpeak_rss_kib\texit\n' > "$TF"
source "$HERE/k2a_validation_lib.sh"
V="$OUT"; printf 'step\tstatus\tdetail\tfile\n' > "$V/steps.tsv"
{ echo "tool=k2a_s0_native_step.sh"; echo "date_utc=$(date -u +%FT%TZ)"; echo "git_commit=$(git -C "$HERE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "groups=$GD"; echo "reference=S0_$REF native 10 ns late step"; echo "lut_sha256=$(sha256sum "$MWPC_SPECTRAL_LUT" | cut -d' ' -f1)"
  echo "replay_sha256=$(sha256sum "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py" | cut -d' ' -f1)"; } > "$OUT/provenance.txt"
bad=0
for g in "$GD"/g_*; do
  G="$OUT/S0groups/${g##*/}"; mkdir -p "$G"; rc=0
  tm "replay10_${g##*/}" python3 "$RC/k2a_timing_ladder.py" "$g" "$G" --lut "$MWPC_SPECTRAL_LUT" --numw "$MWPC_NUMW_DIR" --steps 10 \
      --half-strips "$HS" --link-replay > "$G.ladder.log" 2>&1 || rc=$?
  [[ $rc -eq 0 ]] && { tm "electron_step_${g##*/}" python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$g" "$G/prompt" \
      --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$G.prompt.log" 2>&1 || rc=4; }
  if [[ $rc -eq 0 ]]; then vstep "replay10_${g##*/}" OK "$(tail -1 "$G.ladder.log" | head -c 200)"
  else bad=1; vstep "replay10_${g##*/}" FAIL "exit $rc ($G.ladder.log)"; fi
done
[[ $bad -eq 0 ]] || { echo "native-step S0 incomplete" >&2; exit 1; }
rc=0; tm local_sum python3 "$RC/k2a_local_sum.py" "$GD" "$OUT/S0groups" "$OUT/S0_local" > "$OUT/local_sum.log" 2>&1 || rc=$?
vstep S0_local_census "$([[ $rc -eq 0 ]] && echo OK || echo FAIL)" "$(tail -1 "$OUT/local_sum.log")"
exit $rc

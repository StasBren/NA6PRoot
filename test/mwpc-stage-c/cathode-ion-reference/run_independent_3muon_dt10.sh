#!/usr/bin/env bash
# Safe isolated 3-event validation; reads only SEALED original uncapped Garfield
# K256 muons and the frozen 256-site 10ns kernel trained on DIFFERENT muon0.
# Never writes into or changes the still-active 1000-muon WSL campaign.
set -Eeuo pipefail
umask 027
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
MODE="${1:-preflight}"
DATA="${MWPC_3MUON_DATA:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
TRAIN="${MWPC_3MUON_TRAIN:-$HOME/na6p/cathode-timing-2100/test_runs/mwpc_stage_c/plane_timing_dt10_v1/templates256_dt10_training.root}"
MOBILITY="${MWPC_3MUON_MOBILITY:-$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt}"
OUT="${MWPC_3MUON_OUT:-$REPO/test_runs/mwpc_stage_c/independent_muon_dt10_v1}"
BIN="$REPO/test/mwpc-avalanche/build-independent"
EXE="$BIN/mwpc_plane_response_from_ions_root"
FAST="$BIN/mwpc_plane_muon_kernel_compare"
AUDIT="$BIN/mwpc_audit_independent_muon_kernel_10ns"
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:$HOME/na6p/install/root/lib:${LD_LIBRARY_PATH:-}"
ids=(0301 0408 0446)
fail(){ echo "INDEPENDENT 3MUON STOP: $*" >&2; exit 2; }
sealed(){
  local dir="$1"
  [[ -s "$dir/complete.txt" && -s "$dir/SHA256SUMS" &&
     -s "$dir/ion_births.root" ]] || return 1
  grep -qx 'PASS CATHODE ROOT PHYSICS AUDIT' "$dir/root_audit.log" || return 1
  grep -qx 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$dir/physics.log" || return 1
  (cd "$dir" && sha256sum -c SHA256SUMS >/dev/null 2>&1) || return 1
}
preflight(){
  [[ -s "$TRAIN" && -s "$MOBILITY" ]] ||
    fail "missing frozen template ROOT or ion mobility"
  [[ -x "$EXE" && -x "$FAST" && -x "$AUDIT" ]] ||
    fail "build the three dedicated tools in $BIN before running"
  local tag dir
  for tag in "${ids[@]}";do
    dir="$DATA/muon_$tag"
    sealed "$dir" || fail "muon_$tag is absent, unsealed, invalid, or has bad checksum"
    echo "SEALED muon_$tag $(du -h "$dir/ion_births.root"|cut -f1)"
  done
  [[ "$(df -PB1 "$REPO"|awk 'NR==2{print $4}')" -ge 2147483648 ]] ||
    fail "less than 2 GiB free space; preserve active campaign disk reserve"
  echo "PASS PRECHECK 3 completed independent muons; active microsimulation untouched"
}
run(){
  preflight
  mkdir -p "$OUT"
  exec 9>"$OUT/.analysis.lock"
  flock -n 9 || fail "another analysis process holds the analysis lock"
  [[ ! -e "$OUT/RESULTS_COMPLETE.txt" ]] || fail "results already complete"
  local tag in prefix
  for tag in "${ids[@]}";do
    in="$DATA/muon_$tag/ion_births.root"
    prefix="$OUT/muon_$tag"
    for name in "${prefix}_direct_dt10_summary.csv" \
      "${prefix}_fast_dt10_summary.csv" \
      "${prefix}_fast_dt10_waveforms.csv";do
      [[ ! -e "$name" ]] || fail "already exists; refusing overwrite $name"
    done
    echo "BEGIN independent muon_$tag $(date -Is)"
    "$EXE" --ion-file "$in" \
      --out-prefix "${prefix}_direct_dt10" \
      --mobility "$MOBILITY" --mode direct \
      --dt-ns 10 --window-us 100 --max-step-mm 0.05 \
      --max-direct-ions 20000
    "$FAST" --ion-file "$in" \
      --training-root "$TRAIN" \
      --reference-waveforms "${prefix}_direct_dt10_waveforms.csv" \
      --reference-summary "${prefix}_direct_dt10_summary.csv" \
      --out-prefix "${prefix}_fast_dt10" \
      --neighbors 4 --dt-ns 10 --window-us 100
    echo "DONE independent muon_$tag $(date -Is)"
  done
  "$AUDIT" --data "$DATA" --results "$OUT" \
    --training-root "$TRAIN" \
    --out-prefix "$OUT/independent3"
  {
   printf 'git_commit=%s\n' "$(git -C "$REPO" rev-parse HEAD)"
   printf 'validation=3_independent_muons_fixed_training_256_10ns\n'
   sha256sum "$TRAIN" "$MOBILITY" "$EXE" "$FAST" "$AUDIT"
   for tag in "${ids[@]}";do
     printf 'selected_muon_%s=' "$tag"
     sha256sum "$DATA/muon_$tag/ion_births.root"
   done
  } > "$OUT/PROVENANCE.txt"
  echo "PASS 3 independent muon timing tests, plot/csv written." > "$OUT/RESULTS_COMPLETE.txt"
  cat "$OUT/independent3_events.csv"
  echo "RESULTS $OUT/independent3_windows.csv"
}
case "$MODE" in
 preflight)preflight;;
 run)run;;
 status)
  for tag in "${ids[@]}";do
    if [[ -s "$OUT/muon_${tag}_fast_dt10_summary.csv" ]];then
      echo "OUTPUT_EXISTS muon_$tag"
    else echo "PENDING muon_$tag";fi
  done
  [[ ! -f "$OUT/independent3_events.csv" ]] ||
    cat "$OUT/independent3_events.csv"
  ;;
 *)fail "usage $0 preflight|run|status";;
esac

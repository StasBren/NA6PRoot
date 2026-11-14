#!/usr/bin/env bash
# Explicit gated PBS submission only. Never called by the physics script.
# This launcher prints and records the REAL IDs returned by /opt/pbs/bin/qsub.
set -Eeuo pipefail
STAGE=${1:?Usage: bash submit_wiretruth30k_stage_v1.sh build|pilot20|pilot1000|full}
REPO=${WT_REPO:-"$HOME/NA6PRoot-wiretruth30k"}
ROOT=${WT_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}
QSUB=/opt/pbs/bin/qsub
PROGRAM="$REPO/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
CHUNK_SCRIPT="$REPO/test/mwpc-stage-c/pbs_wiretruth30k_chunk_v1.sh"
BUILD_SCRIPT="$REPO/test/mwpc-stage-c/pbs_wiretruth30k_build_v1.sh"
[[ -x "$QSUB" ]] || { echo "PBS qsub missing: $QSUB" >&2;exit 2; }
[[ -f "$PROGRAM" ]] || { echo "Expected wiretruth production worktree missing" >&2;exit 2; }
[[ "$ROOT" == /storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1 ]] || {
  echo "Refusing unreviewed campaign path; edit intentionally if needed" >&2;exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/submissions"
LEDGER="$ROOT/submissions/qsub_ids.tsv"
if [[ -f "$LEDGER" ]] && grep -q $'\t'"$STAGE"$'\t' "$LEDGER"; then
  if [[ "${WT_RESUBMIT:-0}" != 1 ]]; then
    echo "Stage $STAGE was submitted before. Check qstat and completed chunks." >&2
    echo "Set WT_RESUBMIT=1 ONLY for an explicitly reviewed failed/subset rerun." >&2
    exit 3
  fi
fi
if [[ "$STAGE" == "pilot20" || "$STAGE" == "pilot1000" || "$STAGE" == "full" ]]; then
  [[ -f "$ROOT/locked/binary.json" ]] || { echo "Build lock absent" >&2;exit 4; }
fi
if [[ "$STAGE" == "pilot1000" ]]; then
  [[ -f "$ROOT/gates/gate_20.json" ]] || { echo "20-muon gate absent" >&2;exit 4; }
elif [[ "$STAGE" == full ]]; then
  [[ -f "$ROOT/gates/gate_1000.json" ]] || { echo "1000-muon gate absent" >&2;exit 4; }
fi
submit_one(){
  local name=$1 channel=$2 range=$3 stage=$4 script=$5
  local vars="WT_REPO=$REPO,WT_ROOT=$ROOT"
  if [[ "$stage" != build ]]; then
    vars="$vars,CHANNEL=$channel,WT_STAGE=$stage"
  fi
  local args=(-N "$name" -v "$vars" -o "$ROOT/logs" -e "$ROOT/logs")
  if [[ -n "$range" ]]; then
    args+=(-J "$range")
  elif [[ "$stage" == pilot20 ]]; then
    args+=(-v "$vars,CHUNK_INDEX=0")
  fi
  local jobid
  jobid=$("$QSUB" "${args[@]}" "$script")
  [[ -n "$jobid" ]] || { echo "qsub returned no job ID" >&2;exit 5; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u +%FT%TZ)" "$stage" "$channel" "$range" "$name" "$jobid" >> "$LEDGER"
  printf 'PBS_SUBMITTED stage=%s channel=%s range=%s name=%s id=%s\n' \
    "$stage" "$channel" "${range:-single}" "$name" "$jobid"
}
case "$STAGE" in
  build)
    submit_one w30build all "" build "$BUILD_SCRIPT";;
  pilot20)
    submit_one w30pilot20 Jpsi "" pilot20 "$CHUNK_SCRIPT";;
  pilot1000)
    submit_one w30pjpsi Jpsi "1-19%2" pilot1000 "$CHUNK_SCRIPT"
    submit_one w30pphi Phi "0-14%2" pilot1000 "$CHUNK_SCRIPT"
    submit_one w30pomega Omega "0-14%2" pilot1000 "$CHUNK_SCRIPT";;
  full)
    submit_one w30fjpsi Jpsi "20-499%2" full "$CHUNK_SCRIPT"
    submit_one w30fphi Phi "15-499%2" full "$CHUNK_SCRIPT"
    submit_one w30fomega Omega "15-499%2" full "$CHUNK_SCRIPT";;
  *)
    echo "Unknown stage. Choose build|pilot20|pilot1000|full." >&2;exit 2;;
esac
echo "Job IDs ledger: $LEDGER"
echo "Monitor with /usr/wipp/bin/qtop or qstat -t -u $USER"

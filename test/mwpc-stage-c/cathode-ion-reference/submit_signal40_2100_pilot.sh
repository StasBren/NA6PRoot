#!/usr/bin/env bash
# Prepare separate 2100V ROOT-ion archive first with:
# MWPC_DETAIL_OUTPUT_DIR=/storage/cgrp/stanis/mwpc_penning057/phase_space_2100_v1 \
#   bash test/mwpc-stage-c/cathode-ion-reference/manage_cathode_ion_reference.sh prepare
# This script submits just TWO *selected* Jpsi K512 original muons, NO 1000 array.
set -Eeuo pipefail
umask 077
: "${MWPC_SIGNAL40_SELECTION:?Copy selected40.csv from local verified 997 K256 survey to WIPP}"
: "${MWPC_DETAIL_QUEUE:?Inspect qstat -Qf and provide a compatible PBS queue}"
: "${MWPC_DETAIL_MEM:?Explicitly choose queue-approved amount e.g. 16gb}"
: "${MWPC_DETAIL_WALL:?Explicitly choose queue-approved walltime}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
ROOT="${MWPC_DETAIL_OUTPUT_DIR:-/storage/cgrp/stanis/mwpc_penning057/phase_space_2100_v1}"
[[ "$ROOT" != *hv2100_uncapped_1000_v1* &&
   "$ROOT" != *wiretruth30k_v1* ]] || {
  echo "Refusing to overwrite original 2100/1800 production" >&2;exit 2;
}
[[ -f "$MWPC_SIGNAL40_SELECTION" &&
   "$(head -n1 "$MWPC_SIGNAL40_SELECTION")" == role,index,* ]] || {
  echo "Missing/malformed frozen selected40.csv" >&2;exit 2;
}
[[ "$MWPC_DETAIL_MEM" =~ ^[0-9]+(gb|mb)$ &&
   "$MWPC_DETAIL_WALL" =~ ^[0-9][0-9]:[0-5][0-9]:[0-5][0-9]$ ]] || {
  echo "Invalid MEM/WALL format" >&2;exit 2;
}
ORIG="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
EXE="$REPO/test/mwpc-avalanche/build-cathode-2100/mwpc_p3_noFW_penning057_cathode_ion_reference_v1"
AUDITOR="$REPO/test/mwpc-avalanche/build-cathode-2100/verify_cathode_ion_root"
SCRIPT="$REPO/test/mwpc-stage-c/cathode-ion-reference/pbs_signal40_selected.sh"
[[ -s "$ORIG" && -x "$EXE" && -x "$AUDITOR" && -s "$ROOT/config.lock" ]] || {
  echo "Independent C++/ROOT target not prepared; run manager prepare first" >&2;exit 3;
}
GIT="$(git -C "$REPO" rev-parse HEAD)"
SOURCESHA="$(sha256sum "$ORIG" | cut -d' ' -f1)"
BINHASH="$(sha256sum "$EXE" | cut -d' ' -f1)"
SELECTSHA="$(sha256sum "$MWPC_SIGNAL40_SELECTION" | cut -d' ' -f1)"
[[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ]] || {
  echo "Refusing PBS submission with uncommitted physics code changes" >&2;exit 3;
}
grep -Fxq "git_sha=$GIT" "$ROOT/config.lock"
grep -Fxq "selection_sha256=$SOURCESHA" "$ROOT/config.lock"
grep -Fxq "binary_sha256=$BINHASH" "$ROOT/config.lock"
[[ "$(awk -F, '$1=="development"&&$3=="jpsi"{n++}END{print n+0}' \
    "$MWPC_SIGNAL40_SELECTION")" == 20 ]] || {
  echo "Expected exactly 20 original Jpsi development muons" >&2;exit 3;
}
mkdir -p "$ROOT/pbs-logs" "$ROOT/weighted_k512"
LEDGER="$ROOT/signal40_submissions.tsv"
[[ ! -e "$LEDGER" ]] || {
  echo "Pilot already submitted or recorded; no duplicate; inspect $LEDGER" >&2;exit 3;
}
VARS="MWPC_SIGNAL40_SELECTION=$MWPC_SIGNAL40_SELECTION,MWPC_SIGNAL40_SELECTION_SHA=$SELECTSHA,MWPC_DETAIL_REPO=$REPO,MWPC_DETAIL_OUT=$ROOT/weighted_k512,MWPC_DETAIL_MODE=weighted,MWPC_DETAIL_K=512,MWPC_DETAIL_GIT_SHA=$GIT,MWPC_DETAIL_SELECTION_SHA=$SOURCESHA,MWPC_DETAIL_BINARY_SHA=$BINHASH"
echo "Pilot only: 2 existing-original muons, K512, queue=$MWPC_DETAIL_QUEUE mem=$MWPC_DETAIL_MEM wall=$MWPC_DETAIL_WALL"
JOBID=$(/opt/pbs/bin/qsub -N cSig21 -J '0-1%1' -q "$MWPC_DETAIL_QUEUE" \
  -l "select=1:ncpus=1:mem=$MWPC_DETAIL_MEM" \
  -l "walltime=$MWPC_DETAIL_WALL" \
  -v "$VARS" -o "$ROOT/pbs-logs" -e "$ROOT/pbs-logs" "$SCRIPT")
[[ -n "$JOBID" ]] || { echo "qsub returned empty job ID" >&2;exit 5; }
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%FT%TZ)" \
  signal40-dev-K512 '0-1%1' "$JOBID" "$GIT" "$SELECTSHA" > "$LEDGER"
echo "PBS_SUBMITTED pilot=$JOBID (NO other jobs submitted)"
echo "Verify status with /opt/pbs/bin/qstat -x -f $JOBID"

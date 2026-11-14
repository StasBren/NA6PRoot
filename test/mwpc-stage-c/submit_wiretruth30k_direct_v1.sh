#!/usr/bin/env bash
# EXPLICIT direct 3x10,000 user-requested production, WITHOUT a fabricated gate1000.
# Separate Git-tracked control worktree; frozen physics worktree is unchanged.
# Usage: bash submit_wiretruth30k_direct_v1.sh Jpsi|Phi|Omega
set -Eeuo pipefail
CHANNEL=${1:?Usage: bash submit_wiretruth30k_direct_v1.sh Jpsi|Phi|Omega}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid channel $CHANNEL" >&2; exit 2;; esac
CONTROL=${WT_CONTROL_REPO:-"$HOME/NA6PRoot-wiretruth30k-direct-control"}
REPO=${WT_REPO:-"$HOME/NA6PRoot-wiretruth30k"}
ROOT=${WT_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}
QSUB=/opt/pbs/bin/qsub
WORKER="$CONTROL/test/mwpc-stage-c/pbs_wiretruth30k_direct_v1.sh"
PROGRAM="$REPO/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
EXPECTED_PHYSICS_COMMIT=599e29522edb6fcec9dfcca2904308cb3186edce

[[ -x "$QSUB" && -f "$WORKER" && -f "$PROGRAM" ]] || {
  echo "Missing qsub, Git-tracked direct worker, or frozen physics checker" >&2;exit 2;
}
[[ "$ROOT" == /storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1 ]] || {
  echo "Refusing to run on unreviewed production root" >&2;exit 2;
}
[[ "$(git -C "$REPO" rev-parse HEAD)" == "$EXPECTED_PHYSICS_COMMIT" ]] || {
  echo "Frozen physics repo HEAD changed; do not change existing locked campaign" >&2;exit 2;
}
CONTROL_SHA=$(git -C "$CONTROL" rev-parse HEAD)
[[ -z "$(git -C "$CONTROL" status --porcelain --untracked-files=no)" ]] || {
  echo "Direct-control branch has uncommitted tracked changes" >&2;exit 2;
}
[[ -f "$ROOT/gates/gate_20.json" && -f "$ROOT/locked/binary.json" ]] || {
  echo "Actual measured pilot20 gate and locked Garfield binary are required" >&2;exit 2;
}
# Explicit independent check; never write or invent gate_1000.json.
python3 - "$ROOT" "$EXPECTED_PHYSICS_COMMIT" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1]); expected=sys.argv[2]
gate=json.loads((root/"gates/gate_20.json").read_text())
binary=json.loads((root/"locked/binary.json").read_text())
selection=json.loads((root/"locked/selection.json").read_text())
assert gate["gate"]=="20" and gate["chunks"]==1
assert gate["counts"]["muons"]==20 and gate["counts"]["capped"]==0
assert gate["source_git_commit"]==binary["git_commit"]==selection["git_commit"]==expected
assert gate["binary_sha256"]==binary["binary_sha256"]
assert len(selection["input_sha256"])==1500
print("PASS direct 30k opt-in: genuine gate20 and 3x500 source chunks; 1000 pilot waived")
PY

mkdir -p "$ROOT/logs" "$ROOT/submissions"
LEDGER="$ROOT/submissions/qsub_ids.tsv"
# Avoid re-submitting a channel (even after a partial 3-channel submit).
# Any recovery must first review job IDs, claims and verified outputs.
if [[ -f "$LEDGER" ]] && grep -q $'\t'direct30k$'\t'"$CHANNEL"$'\t' "$LEDGER"; then
  echo "DIRECT30K $CHANNEL previously submitted. Review qstat/status; no duplicate array." >&2
  exit 3
fi
case "$CHANNEL" in
  Jpsi) NAME=w30djpsi;;
  Phi) NAME=w30dphi;;
  Omega) NAME=w30domega;;
esac
RANGE="0-499%2"
VARS="WT_REPO=$REPO,WT_ROOT=$ROOT,WT_CONTROL_REPO=$CONTROL,WT_CONTROL_SHA=$CONTROL_SHA,CHANNEL=$CHANNEL,WT_STAGE=direct30k"
JOBID=$("$QSUB" -N "$NAME" -J "$RANGE" -v "$VARS" \
  -o "$ROOT/logs" -e "$ROOT/logs" "$WORKER")
[[ -n "$JOBID" ]] || { echo "PBS returned no array job ID" >&2; exit 5; }
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$(date -u +%FT%TZ)" "direct30k" "$CHANNEL" "$RANGE" "$NAME" "$JOBID" "$CONTROL_SHA" >> "$LEDGER"
printf 'PBS_SUBMITTED stage=direct30k channel=%s range=%s job_id=%s control_git=%s physics_git=%s\n' \
  "$CHANNEL" "$RANGE" "$JOBID" "$CONTROL_SHA" "$EXPECTED_PHYSICS_COMMIT"
echo "Resources: each 20-muon worker 1 CPU / 4 GB / 8 h; at most 2 active chunks/channel."
echo "Scientific gate1000 intentionally not created, as explicitly requested."

#!/usr/bin/env bash
#PBS -q N
#PBS -N w18hs
#PBS -m n
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1
# Standard-library-only stage: DO NOT source CVMFS or run Garfield++/ROOT.
set -Eeo pipefail
set -u
: "${PBS_JOBID:?PBS job required}"
: "${PBS_ARRAY_INDEX:?PBS array index required}"
: "${TMPDIR:?PBS scratch required}"
: "${MWPC_1800_ENRICHED:?output required}"
: "${MWPC_1800_DASH_WORKTREE:?worktree required}"
LOG="$MWPC_1800_ENRICHED/pbs-logs/shard_array_${PBS_ARRAY_INDEX}.progress"
printf 'SHARD_WORKER_ENTRY utc=%s job=%s host=%s array=%s\n' \
  "$(date -u +%FT%TZ)" "$PBS_JOBID" "$(hostname)" "$PBS_ARRAY_INDEX" >> "$LOG"
exec >> "$LOG" 2>&1
trap 'rc=$?; echo "SHARD_WORKER_EXIT utc=$(date -u +%FT%TZ) code=$rc"' EXIT
python3 -u "$MWPC_1800_DASH_WORKTREE/test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.py" \
  stage --output "$MWPC_1800_ENRICHED" --array-index "$PBS_ARRAY_INDEX"

#!/usr/bin/env bash
#PBS -q N
#PBS -N w18hr
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1
# ROOT-only reduction from locally staged, SHA-audited normalized Heed records.
set -Eeo pipefail
set -u
: "${PBS_JOBID:?PBS job required}"
: "${TMPDIR:?PBS scratch required}"
: "${MWPC_1800_ENRICHED:?output required}"
: "${MWPC_1800_DASH_WORKTREE:?worktree required}"
LOG="$MWPC_1800_ENRICHED/pbs-logs/final_root.progress"
printf 'HEED_REDUCE_WORKER_ENTRY utc=%s job=%s host=%s\n' \
  "$(date -u +%FT%TZ)" "$PBS_JOBID" "$(hostname)" >> "$LOG"
exec >> "$LOG" 2>&1
trap 'rc=$?; echo "HEED_REDUCE_WORKER_EXIT utc=$(date -u +%FT%TZ) code=$rc"' EXIT
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
echo "HEED_REDUCE_ENV_READY utc=$(date -u +%FT%TZ)"
python3 -u "$MWPC_1800_DASH_WORKTREE/test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.py" \
  combine --output "$MWPC_1800_ENRICHED"

#!/usr/bin/env bash
#PBS -q N
#PBS -N w18all
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -l io=0.1
# No original Garfield, Heed, Geant4, detector or production outputs are modified.
set -Eeo pipefail
: "${MWPC_1800_FULL:?Missing frozen original 1800V analysis campaign}"
: "${MWPC_1800_WORKTREE:?Missing isolated worktree}"
: "${PBS_ARRAY_INDEX:?PBS shard index required}"
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
echo "W1800_PBS_JOB=${PBS_JOBID:-unknown} INDEX=$PBS_ARRAY_INDEX HOST=$(hostname)"
echo "FROZEN_CAMPAIGN=$MWPC_1800_FULL"
echo "SOURCE_BRANCH=$(git -C "$MWPC_1800_WORKTREE" rev-parse HEAD)"
echo "START_UTC=$(date -u +%FT%TZ)"
python3 -u "$MWPC_1800_WORKTREE/test/mwpc-stage-c/wiretruth1800/full_wiretruth1800_campaign.py" \
  run-shard --campaign "$MWPC_1800_FULL" --shard "$PBS_ARRAY_INDEX"
echo "PASS_PBS_SHARD_INDEX=$PBS_ARRAY_INDEX"
echo "FINISH_UTC=$(date -u +%FT%TZ)"

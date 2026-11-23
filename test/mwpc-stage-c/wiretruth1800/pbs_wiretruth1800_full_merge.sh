#!/usr/bin/env bash
#PBS -q N
#PBS -N w18merge
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -l io=0.1
# Exact TFileMerger + distinct crossing IDs, ROOT hist and TProfile merge,
# real original Geant4/Heed/anode counts, final plots once.
set -Eeo pipefail
: "${MWPC_1800_FULL:?Missing frozen original 1800V analysis campaign}"
: "${MWPC_1800_WORKTREE:?Missing isolated worktree}"
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
echo "W1800_MERGE_JOB=${PBS_JOBID:-unknown} HOST=$(hostname)"
echo "FROZEN_CAMPAIGN=$MWPC_1800_FULL"
echo "MERGE_GIT=$(git -C "$MWPC_1800_WORKTREE" rev-parse HEAD)"
echo "START_UTC=$(date -u +%FT%TZ)"
python3 -u "$MWPC_1800_WORKTREE/test/mwpc-stage-c/wiretruth1800/full_wiretruth1800_campaign.py" \
  merge --campaign "$MWPC_1800_FULL"
echo "PASS_FULL_ROOT_MERGE_COMPLETE=$MWPC_1800_FULL/merged"
echo "FINISH_UTC=$(date -u +%FT%TZ)"

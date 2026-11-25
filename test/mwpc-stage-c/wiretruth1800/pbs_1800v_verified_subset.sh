#!/usr/bin/env bash
# PBS worker for separate, frozen 5620-muon C++/ROOT merged plots.
set -Eeo pipefail
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
: "${PBS_JOBID:?PBS job required}"
: "${TMPDIR:?PBS scratch required}"
: "${MWPC_1800_SUBSET:?frozen subset required}"
: "${MWPC_1800_DASH_WORKTREE:?dashboard worktree required}"
python3 -u "$MWPC_1800_DASH_WORKTREE/test/mwpc-stage-c/wiretruth1800/merge_1800v_verified_subset.py" run --output "$MWPC_1800_SUBSET"

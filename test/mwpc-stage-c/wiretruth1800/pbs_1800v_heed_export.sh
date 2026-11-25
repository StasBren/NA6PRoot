#!/usr/bin/env bash
#PBS -q N
#PBS -N w18heed
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -l io=0.1
# No Garfield++ rerun: enrich SHA-sealed ROOT COPY using original Heed CSV.
set -Eeo pipefail
: "${PBS_JOBID:?PBS submission required}"
: "${TMPDIR:?PBS scratch required}"
: "${MWPC_1800_ENRICHED:?Prepared new output required}"
: "${MWPC_1800_DASH_WORKTREE:?Dashboard worktree required}"
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
echo "HEED_EXPORT_BEGIN_UTC=$(date -u +%FT%TZ) PBS=$PBS_JOBID HOST=$(hostname)"
python3 -u "$MWPC_1800_DASH_WORKTREE/test/mwpc-stage-c/wiretruth1800/export_1800v_heed_root.py" \
  run --output "$MWPC_1800_ENRICHED"
echo "HEED_EXPORT_DONE_UTC=$(date -u +%FT%TZ)"

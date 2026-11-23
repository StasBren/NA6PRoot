#!/usr/bin/env bash
#PBS -N w1800ana
#PBS -m n
#PBS -j oe
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
# PBS ROOT/C++ analysis only: never reruns microscopic Garfield++ avalanches.
set -Eeo pipefail
: "${MWPC_1800_ANALYSIS_REPO:?Missing isolated analysis worktree path}"
: "${MWPC_1800_ANALYSIS_OUT:?Missing unique analysis output}"
REPO="$MWPC_1800_ANALYSIS_REPO"
SCRIPT="$REPO/test/mwpc-stage-c/wiretruth1800/run_wiretruth1800_analysis.sh"
[[ -f "$SCRIPT" ]] || { echo "Missing ROOT/C++ analysis runner" >&2;exit 2; }
[[ -x "$REPO/test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_distributions" ]] || {
 echo "Build ROOT-only executable on login before submitting" >&2;exit 2;
}
echo "PBS_JOB=${PBS_JOBID:-manual} NODE=$(hostname)"
echo "SOURCE_1800_V=${MWPC_1800_ROOT:-unset}"
echo "ANALYSIS_DEST=$MWPC_1800_ANALYSIS_OUT"
echo "ANALYSIS_GIT=$(git -C "$REPO" rev-parse HEAD)"
echo "START_UTC=$(date -u +%FT%TZ)"
bash "$SCRIPT" run
echo "FINISH_UTC=$(date -u +%FT%TZ)"

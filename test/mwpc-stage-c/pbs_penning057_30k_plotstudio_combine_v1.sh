#!/bin/bash
#PBS -N plot30kAll
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=02:00:00
#PBS -l io=0.1

# Combine already SHA-verified channel caches; never touches raw Garfield.
set -eo pipefail
BASE=/storage/cgrp/stanis/mwpc_penning057
LOGDIR="$BASE/plotstudio_30k_v1/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/combine_${PBS_JOBID:-manual}.log"
exec >"$LOG" 2>&1
START=$(date +%s)
echo "START combined 30,000-muon cache on $(hostname) $(date -u)"
if python3 "$HOME/NA6PRoot-mwpc/test/mwpc-stage-c/pi_penning057_30k_plotstudio_cache_v1.py" combine; then
  rc=0
else
  rc=$?
fi
echo "FINISH combined 30k cache exit=$rc elapsed_seconds=$(( $(date +%s)-START )) $(date -u)"
exit "$rc"

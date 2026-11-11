#!/bin/bash
#PBS -N plot30k
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1

# Read-only, original-2M, Penning ON; produces a NEW compact plot cache.
set -eo pipefail
CHANNEL=${CHANNEL:?Set CHANNEL=Jpsi, Phi or Omega using qsub -v}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid CHANNEL=$CHANNEL"; exit 2;; esac
BASE=/storage/cgrp/stanis/mwpc_penning057
LOGDIR="$BASE/plotstudio_30k_v1/logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/prepare_${CHANNEL}_${PBS_JOBID:-manual}.log"
exec >"$LOG" 2>&1
START=$(date +%s)
echo "START Plot Studio 30k cache $CHANNEL on $(hostname) $(date -u)"
PYTHON="$BASE/plotstudio_venv/bin/python"
if ! "$PYTHON" -c 'import numpy; print("Worker NumPy", numpy.__version__)'; then
  echo "ERROR: NumPy unavailable in cluster plotstudio_venv" >&2
  exit 10
fi
if "$PYTHON" "$HOME/NA6PRoot-mwpc/test/mwpc-stage-c/pi_penning057_30k_plotstudio_cache_v1.py" prepare --channel "$CHANNEL"; then
  rc=0
else
  rc=$?
fi
echo "FINISH Plot Studio cache $CHANNEL exit=$rc elapsed_seconds=$(( $(date +%s)-START )) $(date -u)"
exit "$rc"

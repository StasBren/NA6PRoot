#!/bin/bash
#PBS -N penningSHA
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1

# Read-only full scientific SHA audit on a PBS worker, never login node.
set -eo pipefail
CHANNEL=${CHANNEL:?Submit using -v CHANNEL=Jpsi/Phi/Omega}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid CHANNEL=$CHANNEL" >&2; exit 2;; esac
ROOT=/storage/cgrp/stanis/mwpc_penning057
LOGDIR="$ROOT/audit_logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/full_sha_${CHANNEL}_${PBS_JOBID:-manual}.log"
exec >"$LOG" 2>&1
START=$(date +%s)
echo "START full SHA-256 audit channel=$CHANNEL job=$PBS_JOBID node=$(hostname) date=$(date -u)"
if python3 "$HOME/NA6PRoot-mwpc/test/mwpc-stage-c/pi_penning057_original2m_full_recovery_audit_v1.py" --channel "$CHANNEL" --report-every 25; then
    rc=0
else
    rc=$?
fi
echo "FINISH channel=$CHANNEL rc=$rc wall=$(( $(date +%s)-START )) seconds at $(date -u)"
exit "$rc"

#!/usr/bin/env bash
# PBS Pro one-job-per-chunk wrapper; submit ONLY after a 5-voltage pilot passed.
# 49 production chunks per HV x 5 voltages = 245 PBS array indices.
# Typical PBS submission (queue limits/site policy may differ):
# qsub -J 0-244%12 -v NA60_HV_OUT=/storage/cgrp/stanis/p3_hv_scan_wiretruth_1000_v1 \
#      test/mwpc-stage-c/hv_scan_v1/pbs_hv_scan_chunk.sh
#PBS -N na60_p3_hvscan
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=06:00:00
#PBS -j oe
set -Eeuo pipefail
# PBS often copies the submitted script into its spool directory; do not
# derive the Git worktree from BASH_SOURCE[0] inside a job.
REPO="${NA60_HV_REPO:-${PBS_O_WORKDIR:-$HOME/na6p/p3-hv-scan}}"
HERE="$REPO/test/mwpc-stage-c/hv_scan_v1"
[[ -f "$HERE/run_hv_scan_v1.sh" ]] || { echo "Worktree unavailable: $HERE" >&2;exit 2; }
[[ -n "${PBS_ARRAY_INDEX:-}" ]] || { echo 'PBS_ARRAY_INDEX is missing' >&2;exit 2; }
[[ "${PBS_ARRAY_INDEX}" =~ ^[0-9]+$ ]] || exit 2
idx=$((10#$PBS_ARRAY_INDEX))
((idx>=0 && idx<245)) || { echo 'Array index out of range' >&2;exit 2; }
volts=(1700 1750 1800 1850 1900)
hv="${volts[$((idx/49))]}"
chunk="$((1+idx%49))"
export NA60_HV_OUT="${NA60_HV_OUT:-/storage/cgrp/stanis/p3_hv_scan_wiretruth_1000_v1}"
# Do not build in parallel: login node builds binary before pilot/array launch.
exec bash "$HERE/run_hv_scan_v1.sh" chunk "$hv" "$chunk"

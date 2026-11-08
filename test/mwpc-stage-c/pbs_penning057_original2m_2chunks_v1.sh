#!/bin/bash
#PBS -N penning2x
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1

# 250 array subjobs/channel, instead of 499, each processing up to two
# sequential 20-muon chunks. This avoids the cluster's Q-state job cap.
# Index 1 -> chunk 1,2; ... index 249 -> 497,498; index 250 -> 499.
# Chunk 0000 is the completed pilot and is never re-run.
set -euo pipefail

CHANNEL=${CHANNEL:?Submit with qsub -v CHANNEL=Jpsi or Phi or Omega}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid channel: $CHANNEL"; exit 2 ;; esac
INDEX=${PBS_ARRAY_INDEX:?Must submit as a PBS array with -J 1-250}
case "$INDEX" in ''|*[!0-9]*) echo "Invalid PBS array index: $INDEX" >&2; exit 2 ;; esac
if (( INDEX < 1 || INDEX > 250 )); then
    echo "ERROR: PBS array index must be 1..250" >&2
    exit 2
fi

REPO="$HOME/NA6PRoot-mwpc"
WORKER="$REPO/test/mwpc-stage-c/pbs_penning057_original2m_10k_v1.sh"
ROOT=/storage/cgrp/stanis/mwpc_penning057
lc_channel=${CHANNEL,,}
RESULTS="$ROOT/original2m_${lc_channel}_10000_v1/results"
START=$((2*INDEX-1))
END=$((2*INDEX))

echo "PBS packed array: channel=$CHANNEL index=$INDEX chunks $START..$END node=$(hostname) $(date -u)"
for ((I=START; I<=END; I++)); do
    if ((I > 499)); then break; fi
    printf -v FOUR '%04d' "$I"
    OUT="$RESULTS/chunk_$FOUR"
    if [[ -d "$OUT" ]]; then
        if [[ -f "$OUT/exit_status.txt" &&
              "$(cat "$OUT/exit_status.txt")" == "0" &&
              -f "$OUT/verified_report.json" &&
              -f "$OUT/outputs.sha256" ]]; then
            # Existing completed job may have been produced by previous PBS
            # array. Verify checksums before skipping it.
            (cd "$OUT" && sha256sum -c outputs.sha256 >/dev/null)
            echo "SKIP SHA-verified existing complete chunk $FOUR"
            continue
        fi
        echo "ERROR: existing incomplete/failed chunk $FOUR; NEVER OVERWRITE" >&2
        exit 3
    fi
    echo "START original Geant4 20-muon chunk $FOUR $(date -u)"
    # Per-chunk worker expects PBS_ARRAY_INDEX by default; override with an
    # isolated CHUNK_INDEX and leave the actual parent array untouched.
    env -u PBS_ARRAY_INDEX CHUNK_INDEX="$I" /bin/bash "$WORKER"
    echo "PASS completed chunk $FOUR $(date -u)"
done
echo "PASS packed worker completed channel=$CHANNEL chunks $START..$END $(date -u)"

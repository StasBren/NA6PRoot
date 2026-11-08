#!/bin/bash
#PBS -N penning10k
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1

# A normal qsub executes pilot chunk 0.
# qsub -J 1-499%10 executes remaining 9,980 muons, max ten concurrent jobs.
set -eo pipefail
CAMPAIGN=/storage/cgrp/stanis/mwpc_penning057/jpsi10000_penning057_v1
REPO="$HOME/NA6PRoot-mwpc"
SCRIPT="$REPO/test/mwpc-stage-c/pi_p3_noFW_penning057_cluster10k_v1.py"
BIN="$REPO/test/mwpc-avalanche/build-cvmfs-gcc13/mwpc_p3_noFW_jpsi_penning057_highheed_v1"
IDX=${PBS_ARRAY_INDEX:-${CHUNK_INDEX:-0}}

case "$IDX" in ''|*[!0-9]*)
  echo "ERROR: invalid array index $IDX" >&2
  exit 2 ;;
esac
if ((IDX < 0 || IDX > 499)); then
  echo "ERROR: chunk index must be 0..499" >&2
  exit 2
fi
printf -v FOUR '%04d' "$IDX"
IN="$CAMPAIGN/chunks/chunk_$FOUR.csv"
OUT="$CAMPAIGN/results/chunk_$FOUR"

# Deliberately refuse any existing task folder (completed or partial).
mkdir "$OUT" || {
  echo "ERROR: result directory exists: $OUT ; inspect, never overwrite" >&2
  exit 3
}
exec > "$OUT/worker.log" 2>&1
echo "START PBS job $PBS_JOBID chunk $FOUR node $(hostname) date $(date -u)"
START_TIME=$(date +%s)
WORK=""

finalize() {
  rc=$?
  trap - EXIT
  set +e
  if [[ -n "$WORK" && -d "$WORK" ]]; then
    cp -a "$WORK/." "$OUT/"
    copy_rc=$?
    if ((copy_rc != 0)); then
      echo "ERROR: scratch -> Lustre copy failed: $copy_rc" >&2
      rc=90
    fi
  fi
  wall=$(( $(date +%s) - START_TIME ))
  printf '%s\n' "$wall" > "$OUT/wall_seconds.txt"
  printf '%s\n' "$rc" > "$OUT/exit_status.txt"
  echo "FINISH chunk $FOUR exit=$rc worker_wall_seconds=$wall at $(date -u)"
  exit "$rc"
}
trap finalize EXIT

# The same CVMFS environment was confirmed on wn024 with the 1-muon test.
export GARFIELD_INSTALL=/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt
source "$GARFIELD_INSTALL/Garfield++-env.sh"
set -u
echo "ROOT $(root-config --version), GCC $(g++ -dumpfullversion)"
echo "GARFIELD_INSTALL=$GARFIELD_INSTALL"
sha256sum -c "$CAMPAIGN/binary.sha256"
ldd "$BIN" > "$OUT/ldd.txt"
if grep -q 'not found' "$OUT/ldd.txt"; then
  echo "ERROR: unresolved runtime libraries" >&2
  exit 4
fi
python3 "$SCRIPT" check-input --campaign "$CAMPAIGN" --index "$IDX" --binary "$BIN"

# Run physics and produce large CSVs on per-worker local scratch only.
WORK="${TMPDIR:?}/penning_jpsi_$FOUR"
mkdir -p "$WORK"
cp "$IN" "$WORK/native_muons.csv"
"$BIN" --validate-input "$WORK/native_muons.csv"
echo "=== START 20 real muons: TrackHeed + Penning r=0.57 + full microscopic Garfield ==="
"$BIN" --muon-cache-physics-p3-noFW "$WORK/native_muons.csv" "$WORK/penning"
echo "=== Validate and seal the output set ==="
python3 "$SCRIPT" verify-chunk --campaign "$CAMPAIGN" --index "$IDX" \
  --binary "$BIN" --workdir "$WORK"
(cd "$WORK" && sha256sum penning_*.csv verified_report.json > outputs.sha256)
echo "PASS: 20 true J/psi muons complete; copying scratch data to persistent Lustre"

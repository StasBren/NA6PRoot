#!/bin/bash
#PBS -N penningRec
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=04:00:00
#PBS -l io=0.1

# Independent full-microscopic recovery, no original result overwritten.
set -euo pipefail
CHANNEL=${CHANNEL:?Set CHANNEL=Jpsi/Phi/Omega}
INDEX=${CHUNK_INDEX:?Set CHUNK_INDEX=0..499}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Bad CHANNEL" >&2; exit 2;; esac
case "$INDEX" in ''|*[!0-9]*) echo "Bad index" >&2; exit 2;; esac
(( INDEX>=0 && INDEX<=499 )) || exit 2
printf -v FOUR '%04d' "$INDEX"
lower=${CHANNEL,,}
REPO="$HOME/NA6PRoot-mwpc"
CAMPAIGN="/storage/cgrp/stanis/mwpc_penning057/original2m_${lower}_10000_v1"
PY="$REPO/test/mwpc-stage-c/pi_penning057_original2m_cluster10k_v1.py"
BIN="$REPO/test/mwpc-avalanche/build-cvmfs-gcc13/mwpc_p3_noFW_jpsi_penning057_highheed_v1"
ORIG="$CAMPAIGN/results/chunk_$FOUR"
IN="$CAMPAIGN/chunks/chunk_$FOUR.csv"
OUT="$CAMPAIGN/recovery_runtime_v1/chunk_$FOUR"
if [[ -f "$ORIG/exit_status.txt" &&
      "$(cat "$ORIG/exit_status.txt")" == "0" &&
      -f "$ORIG/verified_report.json" && -f "$ORIG/outputs.sha256" ]]; then
    echo "ERROR: original chunk $FOUR already complete; refusing duplicate" >&2
    exit 3
fi
mkdir -p "$CAMPAIGN/recovery_runtime_v1"
mkdir "$OUT" || { echo "ERROR: recovery already exists $OUT" >&2; exit 3; }
exec >"$OUT/worker.log" 2>&1
START=$(date +%s)
printf '%s\n' "$START" >"$OUT/started_epoch.txt"
WORK=""
finalize() {
    rc=$?
    trap - EXIT
    set +e
    if [[ -n "$WORK" && -d "$WORK" ]]; then
        cp -a "$WORK/." "$OUT/" || rc=90
    fi
    printf '%s\n' "$(( $(date +%s) - START ))" >"$OUT/wall_seconds.txt"
    printf '%s\n' "$rc" >"$OUT/exit_status.txt"
    echo "FINISH $CHANNEL chunk $FOUR exit=$rc $(date -u)"
    exit "$rc"
}
trap finalize EXIT
echo "START recovery $CHANNEL:$FOUR on $(hostname), job $PBS_JOBID, $(date -u)"
export GARFIELD_INSTALL=/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt
source "$GARFIELD_INSTALL/Garfield++-env.sh"
ROOT_LIB="$(root-config --libdir)"
export LD_LIBRARY_PATH="$ROOT_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
echo "ROOT $(root-config --version), ROOT_LIB=$ROOT_LIB"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
sha256sum -c "$CAMPAIGN/binary.sha256"
python3 "$PY" check-input --campaign "$CAMPAIGN" --index "$INDEX" --binary "$BIN"
WORK="${TMPDIR:?}/penning_rec_${lower}_$FOUR"
mkdir -p "$WORK"
cp "$IN" "$WORK/native_muons.csv"
ready=0
for attempt in 1 2 3; do
    echo "ROOT/CVMFS preflight attempt=$attempt $(date -u)"
    if [[ -r "$ROOT_LIB/libTree.so" ]]; then
        ldd "$BIN" > "$OUT/ldd.txt"
        if ! grep -q 'not found' "$OUT/ldd.txt" &&
           "$BIN" --validate-input "$WORK/native_muons.csv"; then
            ready=1
            break
        fi
    fi
    if ((attempt<3)); then sleep 15; fi
done
((ready==1)) || { echo "ROOT/CVMFS validation failed after three tries" >&2; exit 127; }
echo "START FULL microscopic TrackHeed + Garfield Penning r=0.57"
"$BIN" --muon-cache-physics-p3-noFW "$WORK/native_muons.csv" "$WORK/penning"
python3 "$PY" verify-chunk --campaign "$CAMPAIGN" --index "$INDEX" \
  --binary "$BIN" --workdir "$WORK"
(cd "$WORK" && sha256sum penning_*.csv verified_report.json > outputs.sha256)
echo "PASS recovery 20 genuine muons; original evidence untouched"

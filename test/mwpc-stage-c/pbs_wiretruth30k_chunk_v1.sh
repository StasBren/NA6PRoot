#!/usr/bin/env bash
#PBS -N wt30chunk
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -l io=0.1

# FULL native Garfield++ per-wire production, one immutable 20-muon chunk.
set -Eeuo pipefail
CHANNEL=${CHANNEL:?Set CHANNEL=Jpsi/Phi/Omega with qsub -v}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid CHANNEL" >&2; exit 2;; esac
IDX=${PBS_ARRAY_INDEX:-${CHUNK_INDEX:?PBS_ARRAY_INDEX or CHUNK_INDEX required}}
[[ "$IDX" =~ ^[0-9]+$ ]] && (( IDX >= 0 && IDX <= 499 )) || {
  echo "Invalid 20-muon chunk index: $IDX" >&2;exit 2;
}
printf -v FOUR '%04d' "$IDX"

STAGE=${WT_STAGE:?Set WT_STAGE=pilot20/pilot1000/full}
case "$STAGE" in
  pilot20)
    [[ "$CHANNEL" == Jpsi && "$IDX" == 0 ]] || {
      echo "pilot20 is exactly Jpsi chunk 0000" >&2;exit 2;};;
  pilot1000)
    case "$CHANNEL" in
      Jpsi) (( IDX <= 19 )) || exit 2;;
      Phi|Omega) (( IDX <= 14 )) || exit 2;;
    esac;;
  full) ;;
  *) echo "Invalid production stage $STAGE" >&2;exit 2;;
esac

REPO=${WT_REPO:-"$HOME/NA6PRoot-wiretruth30k"}
ROOT=${WT_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}
OLD=/storage/cgrp/stanis/mwpc_penning057
PROGRAM="$REPO/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
BIN="$REPO/test/mwpc-avalanche/build-wiretruth30k-cvmfs/mwpc_p3_noFW_penning057_wiretruth_v1"
INPUT="$ROOT/inputs/${CHANNEL,,}/chunk_$FOUR.csv"
FINAL="$ROOT/results/${CHANNEL,,}/chunk_$FOUR"
JOB=${PBS_JOBID:-manual_$$}

# Skip a SHA-verified success. Never recalculate completed physics.
if [[ -f "$FINAL/verified_report.json" ]]; then
  python3 "$PROGRAM" audit-one --root "$ROOT" --channel "$CHANNEL" --index "$IDX"
  echo "SKIP completed and SHA-verified $CHANNEL $FOUR"
  exit 0
fi
if [[ -e "$FINAL" ]]; then
  echo "INCOMPLETE canonical result; inspect manually: $FINAL" >&2;exit 3
fi

if [[ "$STAGE" == "pilot1000" ]]; then
  [[ -f "$ROOT/gates/gate_20.json" ]] || {
    echo "MISSING verified 20-muon gate" >&2;exit 3; }
elif [[ "$STAGE" == "full" ]]; then
  [[ -f "$ROOT/gates/gate_1000.json" ]] || {
    echo "MISSING verified 1000-muon gate" >&2;exit 3; }
fi
[[ -f "$INPUT" && -f "$BIN" ]] || {
  echo "Frozen input or binary missing" >&2;exit 4; }

export GARFIELD_INSTALL=/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt
source "$GARFIELD_INSTALL/Garfield++-env.sh"
python3 "$PROGRAM" runtime-check --root "$ROOT"

CLAIM="$ROOT/claims/${CHANNEL,,}/chunk_$FOUR"
if ! mkdir "$CLAIM" 2>/dev/null; then
  echo "Another PBS attempt owns $CLAIM. Do not run chunk twice." >&2;exit 5
fi
printf '%s\n' "$JOB" > "$CLAIM/owner_job_id"
cleanup(){
  rc=$?
  trap - EXIT INT TERM
  rm -f "$CLAIM/owner_job_id"
  rmdir "$CLAIM" 2>/dev/null || true
  exit "$rc"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

SAFE_JOB=${JOB//[^A-Za-z0-9._-]/_}
ATTEMPT="$ROOT/attempts/${CHANNEL,,}/chunk_${FOUR}_$SAFE_JOB"
mkdir "$ATTEMPT"
exec > >(tee -a "$ATTEMPT/worker.log") 2>&1
echo "START channel=$CHANNEL chunk=$FOUR job=$JOB node=$(hostname)"
date -u +%FT%TZ
echo "FROZEN_G4_INPUT=$INPUT"
echo "GIT_COMMIT=$(git -C "$REPO" rev-parse HEAD)"
sha256sum "$BIN" "$INPUT"

WORK="${TMPDIR:?PBS TMPDIR is required}/wiretruth_${CHANNEL}_${FOUR}_$SAFE_JOB"
mkdir "$WORK"
cp "$INPUT" "$WORK/native_muons.csv"
"$BIN" --validate-input "$WORK/native_muons.csv"
START=$(date +%s)
if [[ "$STAGE" == "pilot20" ]]; then
  "$BIN" --validate-heed "$WORK/native_muons.csv" \
    > "$ATTEMPT/preflight_heed.log" 2>&1
  echo "PASS original first-20 Geant4 TrackHeed preflight"
fi
echo "START Garfield++ full microscopic per-wire"
if [[ ! -x /usr/bin/time ]]; then
  echo "GNU /usr/bin/time required for hard 4GB/8h resource pilot gate" >&2
  exit 6
fi
/usr/bin/time -v -o "$ATTEMPT/resource_time.txt" \
  "$BIN" --wiretruth-p3-noFW "$WORK/native_muons.csv" "$WORK/microscopic"
echo "GARFIELD_WALL_SECONDS=$(( $(date +%s) - START ))"

# Inspect ALL independent original endpoints, ion birth censuses, wire map,
# same-Heed replay, and real reservoir weights before persistence.
python3 "$PROGRAM" verify --root "$ROOT" --old-root "$OLD" \
  --channel "$CHANNEL" --index "$IDX" --work "$WORK"

STAGING="$ATTEMPT/staged"
mkdir "$STAGING"
cp -a "$WORK/microscopic_"*.csv "$WORK/verified_report.json" "$STAGING/"
python3 "$PROGRAM" audit-work --root "$ROOT" --channel "$CHANNEL" \
  --index "$IDX" --work "$STAGING"
if [[ -f "$ATTEMPT/resource_time.txt" ]]; then
  cp "$ATTEMPT/resource_time.txt" "$STAGING/resource_time.txt"
fi
if [[ -e "$FINAL" ]]; then
  echo "Duplicate canonical result, refusing overwrite." >&2;exit 7
fi
mv "$STAGING" "$FINAL"
echo "PASS COMPLETE 20 muons: $FINAL"
date -u +%FT%TZ

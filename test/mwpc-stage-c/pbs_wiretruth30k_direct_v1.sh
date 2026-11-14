#!/usr/bin/env bash
#PBS -N wt30chunk
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=08:00:00
#PBS -l io=0.1

# Explicitly authorized DIRECT 30k production after the successful 20-muon gate.
# Separate CONTROL branch; physical C++ and Python verification remain pinned
# to the immutable 599e295 worktree with the exact locked binary SHA.
set -Eeo pipefail
CHANNEL=${CHANNEL:?Set CHANNEL=Jpsi/Phi/Omega with qsub -v}
case "$CHANNEL" in Jpsi|Phi|Omega) ;; *) echo "Invalid CHANNEL" >&2; exit 2;; esac
IDX=${PBS_ARRAY_INDEX:-${CHUNK_INDEX:?PBS_ARRAY_INDEX or CHUNK_INDEX required}}
[[ "$IDX" =~ ^[0-9]+$ ]] && (( IDX >= 0 && IDX <= 499 )) || {
  echo "Invalid 20-muon chunk index: $IDX" >&2;exit 2;
}
printf -v FOUR '%04d' "$IDX"

STAGE=${WT_STAGE:?Set WT_STAGE=direct30k}
[[ "$STAGE" == direct30k ]] || { echo "DIRECT worker requires direct30k" >&2;exit 2; }

REPO=${WT_REPO:-"$HOME/NA6PRoot-wiretruth30k"}
ROOT=${WT_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}
OLD=/storage/cgrp/stanis/mwpc_penning057
PROGRAM="$REPO/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
BIN="$REPO/test/mwpc-avalanche/build-wiretruth30k-cvmfs/mwpc_p3_noFW_penning057_wiretruth_v1"
# Point directly at the old, SHA-sealed Geant4 source selection.
# Use read-only original inputs, never resimulate Geant4 and never bulk copy.
INPUT="$OLD/original2m_${CHANNEL,,}_10000_v1/chunks/chunk_$FOUR.csv"
FINAL="$ROOT/results/${CHANNEL,,}/chunk_$FOUR"
JOB=${PBS_JOBID:-manual_$}

CONTROL=${WT_CONTROL_REPO:?Set WT_CONTROL_REPO to separate direct-control worktree}
CONTROL_SHA=${WT_CONTROL_SHA:?Set frozen control Git commit SHA}
[[ "$(git -C "$CONTROL" rev-parse HEAD)" == "$CONTROL_SHA" ]] || {
  echo "Direct-control worker Git commit changed since PBS submission" >&2; exit 3;
}
[[ -z "$(git -C "$CONTROL" status --porcelain --untracked-files=no)" ]] || {
  echo "Direct-control tracked code has uncommitted changes" >&2; exit 3;
}
# The original physics worktree must remain at the only already approved
# physics binary/source lock. The original Python runtime-check checks again.
[[ "$(git -C "$REPO" rev-parse HEAD)" == "599e29522edb6fcec9dfcca2904308cb3186edce" ]] || {
  echo "Original physics worktree is no longer at the verified frozen commit" >&2; exit 3;
}
[[ -f "$ROOT/gates/gate_20.json" && -f "$ROOT/locked/binary.json" ]] || {
  echo "Missing actual physical pilot20 gate / binary SHA lock" >&2; exit 3;
}
python3 - "$ROOT" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1])
g=json.loads((root/"gates/gate_20.json").read_text())
b=json.loads((root/"locked/binary.json").read_text())
s=json.loads((root/"locked/selection.json").read_text())
assert g["gate"]=="20" and g["chunks"]==1
assert g["counts"]["muons"]==20 and g["counts"]["capped"]==0
assert g["source_git_commit"]==s["git_commit"]==b["git_commit"]
assert g["binary_sha256"]==b["binary_sha256"]
assert len(s["input_sha256"])==1500
print("PASS DIRECT30K actual gate20 + frozen input/binary provenance",flush=True)
PY

# Skip a SHA-verified success. Never recalculate completed physics.
if [[ -f "$FINAL/verified_report.json" ]]; then
  python3 "$PROGRAM" audit-one --root "$ROOT" --channel "$CHANNEL" --index "$IDX"
  echo "SKIP completed and SHA-verified $CHANNEL $FOUR"
  exit 0
fi
if [[ -e "$FINAL" ]]; then
  echo "INCOMPLETE canonical result; inspect manually: $FINAL" >&2;exit 3
fi

[[ -f "$INPUT" && -f "$BIN" ]] || {
  echo "Frozen input or binary missing" >&2;exit 4; }

export GARFIELD_INSTALL=/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt
source "$GARFIELD_INSTALL/Garfield++-env.sh"
# CVMFS initialization may reference an unset ZSH_NAME; activate nounset only
# after sourcing the validated Garfield++ environment.
set -u
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
echo "START DIRECT30K channel=$CHANNEL chunk=$FOUR job=$JOB node=$(hostname)"
echo "DIRECT_CONTROL_GIT=$CONTROL_SHA"
date -u +%FT%TZ
echo "FROZEN_G4_INPUT=$INPUT"
echo "GIT_COMMIT=$(git -C "$REPO" rev-parse HEAD)"
sha256sum "$BIN" "$INPUT"

WORK="${TMPDIR:?PBS TMPDIR is required}/wiretruth_${CHANNEL}_${FOUR}_$SAFE_JOB"
mkdir "$WORK"
cp "$INPUT" "$WORK/native_muons.csv"
"$BIN" --validate-input "$WORK/native_muons.csv"
START=$(date +%s)
echo "START Garfield++ full microscopic per-wire"
if [[ ! -x /usr/bin/time ]]; then
  echo "GNU /usr/bin/time required for 4GB/8h per-chunk resource audit" >&2
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

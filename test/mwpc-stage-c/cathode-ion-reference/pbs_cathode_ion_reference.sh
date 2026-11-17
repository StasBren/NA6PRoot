#!/usr/bin/env bash
#PBS -N cathIon21
#PBS -m n
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -l walltime=12:00:00
#PBS -j oe
# Dedicated one-real-muon uncapped ROOT reference. All previous runs read-only.
set -Eeo pipefail
export GARFIELD_INSTALL=${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}
[[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]] || {
 echo "Pinned CVMFS Garfield++ missing" >&2; exit 2; }
source "$GARFIELD_INSTALL/Garfield++-env.sh"
set -u
umask 077
: "${MWPC_DETAIL_REPO:?}"
: "${MWPC_DETAIL_OUT:?}"
: "${MWPC_DETAIL_MODE:?}"
: "${MWPC_DETAIL_K:?}"
: "${MWPC_DETAIL_GIT_SHA:?}"
: "${MWPC_DETAIL_SELECTION_SHA:?}"
: "${MWPC_DETAIL_BINARY_SHA:?}"
: "${PBS_ARRAY_INDEX:?Requires PBS array}"
IDX=$PBS_ARRAY_INDEX
[[ "$IDX" =~ ^[0-9]+$ ]] && (( IDX >= 0 && IDX < 1000 )) || {
  echo "Invalid frozen muon index $IDX" >&2;exit 3; }
case "$MWPC_DETAIL_MODE:$MWPC_DETAIL_K" in
 weighted:256|weighted:512|weighted:1024|full:0) ;;
 *) echo "Invalid mode/reservoir" >&2;exit 3;;
esac
REPO="$MWPC_DETAIL_REPO"
OUT="$MWPC_DETAIL_OUT"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BUILD="$REPO/test/mwpc-avalanche/build-cathode-2100"
EXE="$BUILD/mwpc_p3_noFW_penning057_cathode_ion_reference_v1"
CHECK="$BUILD/verify_cathode_ion_root"
LOCK="$OUT/../config.lock"
[[ -s "$SELECTION" && -x "$EXE" && -x "$CHECK" && -s "$LOCK" ]] || {
  echo "Missing detail executable, ROOT auditor, selection or lock" >&2;exit 3; }
[[ "$(git -C "$REPO" rev-parse HEAD)" == "$MWPC_DETAIL_GIT_SHA" &&
   "$(sha256sum "$SELECTION" | cut -d' ' -f1)" == "$MWPC_DETAIL_SELECTION_SHA" &&
   "$(sha256sum "$EXE" | cut -d' ' -f1)" == "$MWPC_DETAIL_BINARY_SHA" ]] || {
  echo "Frozen source/input/binary changed after submission" >&2;exit 3; }
grep -Fxq "git_sha=$MWPC_DETAIL_GIT_SHA" "$LOCK"
grep -Fxq "selection_sha256=$MWPC_DETAIL_SELECTION_SHA" "$LOCK"
grep -Fxq "binary_sha256=$MWPC_DETAIL_BINARY_SHA" "$LOCK"
grep -Fxq "garfield_install_realpath=$(readlink -f "$GARFIELD_INSTALL")" "$LOCK"
GARFIELD_SO="$(ldd "$EXE" | awk '$1 ~ /Garfield/ && $3 ~ /^\\/ {print $3;exit}')"
[[ -n "$GARFIELD_SO" && -f "$GARFIELD_SO" ]] || {
 echo "Missing linked Garfield shared library" >&2;exit 3; }
grep -Fxq "garfield_library_realpath=$(readlink -f "$GARFIELD_SO")" "$LOCK"
grep -Fxq "garfield_library_sha256=$(sha256sum "$GARFIELD_SO" | cut -d' ' -f1)" "$LOCK"
grep -Fxq "root_version=$(root-config --version)" "$LOCK"
grep -Fxq "gpp_version=$(c++ -dumpfullversion)" "$LOCK"
printf -v TAG 'muon_%04d' "$IDX"
DEST="$OUT/$TAG"
if [[ -f "$DEST/complete.txt" ]]; then
 (cd "$DEST" && sha256sum -c SHA256SUMS >/dev/null) || {
  echo "CORRUPTED sealed output: $DEST" >&2;exit 3; }
 echo "SKIP SHA-verified $DEST";exit 0
fi
if ! mkdir "$DEST" 2>/dev/null;then
 echo "Refusing to overwrite incomplete or concurrent $DEST" >&2;exit 4
fi
echo "job_id=${PBS_JOBID:-unknown}" > "$DEST/metadata.txt"
echo "array_index=$IDX" >> "$DEST/metadata.txt"
echo "mode=$MWPC_DETAIL_MODE reservoir=$MWPC_DETAIL_K" >> "$DEST/metadata.txt"
echo "selection_sha256=$MWPC_DETAIL_SELECTION_SHA" >> "$DEST/metadata.txt"
echo "binary_sha256=$MWPC_DETAIL_BINARY_SHA" >> "$DEST/metadata.txt"
echo "git_sha=$MWPC_DETAIL_GIT_SHA" >> "$DEST/metadata.txt"
echo "start_utc=$(date -u +%FT%TZ)" >> "$DEST/metadata.txt"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/cathode2100.XXXXXXXX")
cleanup(){
 rc=$?
 trap - EXIT INT TERM
 if ((rc!=0)); then
   printf 'FAILED exit=%s at %s\n' "$rc" "$(date -u +%FT%TZ)" > "$DEST/FAILED.txt"
   for p in physics.log root_audit.log resource_time.txt;do
     [[ ! -f "$WORK/$p" ]] || cp "$WORK/$p" "$DEST/$p" || true
   done
 fi
 rm -rf "$WORK"
 exit "$rc"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
awk -v target="$((IDX+2))" 'NR==1||NR==target' "$SELECTION" > "$WORK/input.csv"
[[ "$(wc -l < "$WORK/input.csv")" == 2 ]] || exit 5
"$EXE" --validate-input "$WORK/input.csv" > "$WORK/input_check.log" 2>&1
cp "$WORK/input.csv" "$DEST/input.csv"
cp "$WORK/input_check.log" "$DEST/input_check.log"
echo "START true Geant4 muon index=$IDX HV=2100 no avalanche size cap mode=$MWPC_DETAIL_MODE K=$MWPC_DETAIL_K"
# Resource failure is NOT a valid finite gain or truncated ion population.
/usr/bin/time -v -o "$WORK/resource_time.txt" \
  timeout --signal=TERM --kill-after=30s 11h30m \
  "$EXE" --cathode-ion-reference "$WORK/input.csv" \
  "$WORK/microscopic" 2100 0 "$MWPC_DETAIL_MODE" "$MWPC_DETAIL_K" \
  "$WORK/ion_births.root" > "$WORK/physics.log" 2>&1 || {
    tail -n 60 "$WORK/physics.log" >&2;exit 6;
  }
grep -Fxq 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$WORK/physics.log" || {
 echo "Missing real uncapped Garfield initialization proof" >&2;exit 6; }
"$CHECK" "$WORK/ion_births.root" "$WORK/microscopic" \
  > "$WORK/root_audit.log" 2>&1 || {
    cat "$WORK/root_audit.log" >&2;exit 7;
  }
cp "$WORK/ion_births.root" "$WORK/microscopic_"*.csv \
 "$WORK/physics.log" "$WORK/root_audit.log" \
 "$WORK/resource_time.txt" "$DEST/"
(cd "$DEST" && sha256sum input.csv input_check.log ion_births.root \
 microscopic_*.csv physics.log root_audit.log resource_time.txt metadata.txt \
 > SHA256SUMS)
(cd "$DEST" && sha256sum -c SHA256SUMS >/dev/null)
printf 'PASS index=%s mode=%s K=%s\n' "$IDX" "$MWPC_DETAIL_MODE" "$MWPC_DETAIL_K" > "$DEST/complete.txt"
echo "PASS cathode-ion reference muon=$IDX"

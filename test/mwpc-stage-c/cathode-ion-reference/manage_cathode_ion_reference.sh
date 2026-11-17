#!/usr/bin/env bash
# 2100 V ion-birth archive: prepare | status | pilot-weighted | pilot-full.
# User selects PBS queue/resources based on qstat -Qf; no mass submissions.
set -Eeo pipefail
export GARFIELD_INSTALL=${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}
[[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]] || {
 echo "Missing pinned Garfield++ CVMFS install" >&2;exit 2; }
source "$GARFIELD_INSTALL/Garfield++-env.sh"
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BUILD="$REPO/test/mwpc-avalanche/build-cathode-2100"
INPUT_AUDITOR="$BUILD/audit_geant4_1000_root"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BIN="$BUILD/mwpc_p3_noFW_penning057_cathode_ion_reference_v1"
AUDITOR="$BUILD/verify_cathode_ion_root"
OUT="${MWPC_DETAIL_OUTPUT_DIR:-/storage/cgrp/stanis/mwpc_penning057/cathode_2100_uncapped_v1}"
LEDGER="$OUT/submissions.tsv"
MODE="${1:-status}"
# PBS request can be overridden by exporting following values explicitly.
QUEUE="${MWPC_DETAIL_QUEUE:-}"
MEM="${MWPC_DETAIL_MEM:-32gb}"
WALL="${MWPC_DETAIL_WALL:-12:00:00}"
if [[ "$OUT" == */hv2100_uncapped_1000_v1* ||
      "$OUT" == */mwpc_penning057_wiretruth30k_v1* ]];then
 echo "Refusing to write to historical/ongoing physics output" >&2;exit 3;
fi
need_input(){
 [[ -s "$SELECTION" ]] || { echo "Original 1000 Geant4 selection missing";exit 3; }
 [[ "$(wc -l < "$SELECTION")" == 1001 ]] || {
  echo "Expected 1000 frozen Geant4 crossings";exit 3; }
}
lock_data(){
 printf 'git_sha=%s\n' "$(git -C "$REPO" rev-parse HEAD)"
 printf 'selection_sha256=%s\n' "$(sha256sum "$SELECTION" | cut -d' ' -f1)"
 printf 'binary_sha256=%s\n' "$(sha256sum "$BIN" | cut -d' ' -f1)"
 printf 'physics=2100V;ArCO2_70_30;Penning0.57;p3_noFW;30um;3mm;cathodes_pm3mm;B0;UNCAPPED\n'
 printf 'root_archive=ion_callbacks_all_or_unbiased_per_primary_reservoir;no_species_inferred\n'
 printf 'garfield_install=%s\n' "$GARFIELD_INSTALL"
 printf 'garfield_install_realpath=%s\n' "$(readlink -f "$GARFIELD_INSTALL")"
 local garfield_so
 garfield_so="$(ldd "$BIN" | awk '$1 ~ /Garfield/ && index($3,"/")==1 {print $3;exit}')"
 [[ -n "$garfield_so" && -f "$garfield_so" ]] || {
   echo "Could not locate actual linked Garfield++ shared library" >&2; exit 4; }
 printf 'garfield_library_realpath=%s\n' "$(readlink -f "$garfield_so")"
 printf 'garfield_library_sha256=%s\n' "$(sha256sum "$garfield_so" | cut -d' ' -f1)"
 printf 'root_version=%s\n' "$(root-config --version)"
 printf 'gpp_version=%s\n' "$(c++ -dumpfullversion)"
}
check_lock(){
 [[ -x "$BIN" && -x "$AUDITOR" && -s "$OUT/config.lock" ]] || {
   echo "Missing sealed independent build or config.lock; prepare first" >&2;exit 4; }
 [[ "$(cat "$OUT/config.lock")" == "$(lock_data)" ]] || {
   echo "Source/selection/binary lock mismatch. Do not change campaign midrun." >&2;exit 4; }
}
prepare(){
 need_input
 mkdir -p "$OUT/pbs-logs"
 echo "Pinned compiler: $(c++ -dumpfullversion)"
 echo "Pinned ROOT: $(root-config --version)"
 echo "Pinned Garfield: $GARFIELD_INSTALL"
 cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD"
 cmake --build "$BUILD" --target \
    mwpc_p3_noFW_penning057_cathode_ion_reference_v1 \
    verify_cathode_ion_root audit_geant4_1000_root compare_sampling_to_full_root select_detailed_1000_root -j2
 [[ -x "$BIN" && -x "$AUDITOR" && -x "$INPUT_AUDITOR" ]] || {
   echo "ROOT/C++ detailed target build failed" >&2;exit 4; }
 TEMP=$(mktemp "$OUT/.inputcheck.XXXXXX.csv")
 { head -n 1 "$SELECTION"; sed -n '2p' "$SELECTION"; } > "$TEMP"
 "$BIN" --validate-input "$TEMP"
 rm -f "$TEMP"
 # This audit only reads immutable Geant4 geometry and draws ROOT-only plots.
 # Actual Heed counts, gain and multi-anode fractions require Garfield.
 if [[ ! -s "$OUT/geant4_input_audit/audit.root" ]];then
   "$INPUT_AUDITOR" "$SELECTION" "$OUT/geant4_input_audit/audit.root"
 fi
 if [[ ! -s "$OUT/geant4_input_audit/selection.root" ]];then
   "$BUILD/select_detailed_1000_root" \
    "$OUT/geant4_input_audit/audit.root" \
    "$OUT/geant4_input_audit/selection.root"
 fi
 expected=$(lock_data)
 if [[ -f "$OUT/config.lock" ]];then
   [[ "$(cat "$OUT/config.lock")" == "$expected" ]] || {
    echo "Existing detail archive lock cannot be overwritten" >&2;exit 4; }
 else
   printf '%s\n' "$expected" > "$OUT/config.lock"
 fi
 echo "PASS prepared isolated ROOT ion archive target, 1000 genuine Geant4 selections"
}
status(){
 echo "Output: $OUT"
 if [[ -f "$OUT/config.lock" ]];then cat "$OUT/config.lock";fi
 if [[ -f "$LEDGER" ]];then cat "$LEDGER";fi
 for spec in weighted_k256 weighted_k512 weighted_k1024 full_k0;do
    local p="$OUT/$spec"
    local n=0,bad=0,d
    if [[ -d "$p" ]];then
     for d in "$p"/muon_*;do
       [[ -d "$d" ]] || continue
       if [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" ]] &&
          (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1);then
          ((n+=1))
       else
          ((bad+=1))
          echo "INVALID/UNSEALED $d"
       fi
     done
    fi
    echo "$spec verified=$n incomplete_or_failed=$bad"
 done
}
submit(){
 check_lock
 command -v /opt/pbs/bin/qsub >/dev/null || {
   echo "qsub not available" >&2;exit 4; }
 local physics="$1" size="$2" range="$3" label="$4"
 local target="$OUT/${physics}_k$size"
 mkdir -p "$target" "$OUT/pbs-logs"
 if [[ -f "$LEDGER" ]] && awk -F'\t' -v label="$label" '$2==label {found=1} END{exit !found}' "$LEDGER";then
   echo "Already submitted $label. Review PBS; no duplicate array." >&2;exit 4
 fi
 local vars qargs=()
 vars="MWPC_DETAIL_REPO=$REPO,MWPC_DETAIL_OUT=$target,MWPC_DETAIL_MODE=$physics,MWPC_DETAIL_K=$size,MWPC_DETAIL_GIT_SHA=$(git -C "$REPO" rev-parse HEAD),MWPC_DETAIL_SELECTION_SHA=$(sha256sum "$SELECTION" | cut -d' ' -f1),MWPC_DETAIL_BINARY_SHA=$(sha256sum "$BIN" | cut -d' ' -f1)"
 if [[ -n "$QUEUE" ]];then qargs+=(-q "$QUEUE");fi
 echo "Submitting $label range=$range queue=${QUEUE:-PBS-default} mem=$MEM wall=$WALL"
 local jobid
 jobid=$(/opt/pbs/bin/qsub "${qargs[@]}" -N cathIon21 \
   -J "$range" -l "select=1:ncpus=1:mem=$MEM" \
   -l "walltime=$WALL" -v "$vars" \
   -o "$OUT/pbs-logs" "$HERE/pbs_cathode_ion_reference.sh")
 printf '%s\t%s\t%s\t%s\t%s\n' \
   "$(date -u +%FT%TZ)" "$label" "$range" "$jobid" \
   "$(git -C "$REPO" rev-parse HEAD)" >> "$LEDGER"
 echo "PBS_SUBMITTED $label $jobid"
}
case "$MODE" in
 prepare) prepare;;
 status) status;;
 pilot-weighted) submit weighted 256 0-1%1 pilot-weighted-256;;
 pilot-full) submit full 0 0-1%1 pilot-full-0;;
 weighted512) submit weighted 512 0-1%1 pilot-weighted-512;;
 weighted1024) submit weighted 1024 0-1%1 pilot-weighted-1024;;
 *) echo "Usage: bash manage_cathode_ion_reference.sh {prepare|status|pilot-weighted|pilot-full|weighted512|weighted1024}" >&2;exit 2;;
esac

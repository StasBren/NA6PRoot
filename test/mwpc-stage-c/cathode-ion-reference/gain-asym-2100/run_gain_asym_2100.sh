#!/usr/bin/env bash
# Existing 2100-V K256 source archives only. No Garfield++, RKF or 1800V runs.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
export PATH="$HOME/na6p/install/root/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:${LD_LIBRARY_PATH:-}"
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
BUILD="${MWPC_GAINASYM_BUILD:-$REPO/test/mwpc-avalanche/build-gain-asym-v0}"
BIN="$BUILD/mwpc_v0_archive_stats"
OUTBASE="${MWPC_GAINASYM_OUT:-$REPO/test_runs/mwpc_stage_c/2100_gain_asym}"
MODE="${1:-locate}"
fail(){ echo "NA60 2100 V analysis: $*" >&2; exit 2; }
sealed(){
 local d="$1"
 [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" &&
    -s "$d/root_audit.log" && -s "$d/ion_births.root" &&
    -s "$d/microscopic_seed.csv" && -s "$d/microscopic_wire_seed.csv" ]] &&
 grep -q "PASS CATHODE ROOT PHYSICS AUDIT" "$d/root_audit.log"
}
inventory(){
 [[ -d "$ARCHIVE" ]] || fail "cannot find archive: $ARCHIVE"
 local n=0 other=0 k name
 LIST_RECORDS=()
 for ((k=0;k<1000;k++));do
   printf -v name 'muon_%04d' "$k"
   if sealed "$ARCHIVE/$name";then
     LIST_RECORDS+=("$ARCHIVE/$name/ion_births.root"); ((n+=1))
   else ((other+=1));fi
 done
 printf "ARCHIVE=%s\nSEALED_AND_ROOT_AUDITED=%s\nOTHER_SLOTS=%s\n" "$ARCHIVE" "$n" "$other"
 echo "Inventory is not a fresh full-file SHA rehash (pilot does recheck five slots)."
}
locate_drift(){
 echo "Ion drift/RKF and U/V waveforms are separate from ion_births.root:"
 for dir in \
  "$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c" \
  "$HOME/na6p/cathode-plane-response-2100/test_runs/mwpc_stage_c" \
  "$HOME/na6p/cathode-timing-2100/test_runs/mwpc_stage_c"; do
  if [[ -d "$dir" ]];then
    echo "RESPONSE_DIR=$dir"
    find "$dir" -maxdepth 7 -type f \( -name '*waveforms.csv' -o -name '*waveform*.root' -o -name '*ion*response*.root' \) -print | head -n 25 || true
  fi
 done
}
case "$MODE" in
 locate)
   inventory
   for id in 0000 0301 0408 0446 0878;do
     d="$ARCHIVE/muon_$id"
     if [[ -d "$d" ]];then
       echo "MUON_$id:"
       find "$d" -maxdepth 1 -type f -printf '  %f\n' | head -n 20
     fi
   done
   locate_drift;;
 prepare)
   command -v root-config >/dev/null || fail "ROOT not available in PATH"
   command -v c++ >/dev/null || fail "C++ compiler not available"
   # Compile the ROOT-only archive reader directly. CMake for the avalanche
   # project globally calls find_package(Garfield), which is unnecessary here.
   # Include "V0Core.hh" resolves next to v0_archive_stats.cc.
   mkdir -p "$BUILD"
   SRC="$REPO/test/mwpc-stage-c/cathode-ion-reference/v0/v0_archive_stats.cc"
   [[ -s "$SRC" ]] || fail "missing existing V0 C++ reader: $SRC"
   c++ -std=c++17 -O2 -DV0_ROOT $(root-config --cflags) "$SRC" \
       -o "$BIN" $(root-config --libs)
   [[ -x "$BIN" ]] || fail "ROOT executable was not built"
   echo "PASS compiled existing ROOT/C++ v0 archive reader without Garfield/CMake (no physics rerun).";;
 pilot|population)
   [[ -x "$BIN" ]] || fail "run 'bash $0 prepare' first"
   command -v root >/dev/null || fail "ROOT not available"
   inventory
   (("${#LIST_RECORDS[@]}" >= 5)) || fail "not enough verified archived muons"
   if [[ "$MODE" == population && "${#LIST_RECORDS[@]}" -ne 997 ]];then
     fail "expected 997 seal/audit-pass slots, found ${#LIST_RECORDS[@]}; investigate before running"
   fi
   out="$OUTBASE/${MODE}_$(date +%Y%m%d_%H%M%S)_$$"
   [[ ! -e "$out" ]] || fail "output already exists"
   mkdir -p "$out"
   list="$out/rootfiles.txt"
   if [[ "$MODE" == pilot ]];then
     for ((i=0;i<5;i++));do
       path="${LIST_RECORDS[$i]}"; d="${path%/ion_births.root}"
       (cd "$d" && sha256sum -c SHA256SUMS >/dev/null) || fail "SHA mismatch at $d"
       printf '%s\n' "$path" >> "$list"
     done
   else
     printf '%s\n' "${LIST_RECORDS[@]}" > "$list"
   fi
   printf 'git_sha=%s\narchive=%s\nmode=%s\n' "$(git -C "$REPO" rev-parse HEAD)" "$ARCHIVE" "$MODE" > "$out/provenance.txt"
   echo "Read-only ROOT and exact wire-census analysis; NO Garfield, drift, FEM or 1800 V data."
   "$BIN" --list "$list" --no-holdout-exclusion --out-dir "$out/census" | tee "$out/v0_console.txt"
   root -l -b -q "$HERE/PlotGainAsym2100.C(\"$out/census\",\"$out/figures\")" | tee "$out/plots_console.txt"
   [[ -s "$out/figures/gain_asym_2100.root" ]] || fail "ROOT plot output not created"
   echo "OUTPUT=$out"
   echo "PNG/PDF=$out/figures/"
   echo "ROOT=$out/figures/gain_asym_2100.root"
   echo "SUMMARY=$out/figures/summary.txt"
   echo "NOTE: full 997 uses holdout events for population DESCRIPTION, not for model calibration.";;
 *)
   fail "Usage: bash $0 locate|prepare|pilot|population";;
esac

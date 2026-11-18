#!/usr/bin/env bash
# NA60+ 2100 V laptop pilot: same Garfield++ ROOT ion-birth archive as PBS.
# Usage: bash run_local_cathode_reference.sh prepare|weighted INDEX K|full INDEX|compare INDEX K|status
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BUILD="$REPO/test/mwpc-avalanche/build-cathode-2100-wsl"
OUT="${MWPC_LOCAL_OUTPUT_DIR:-$REPO/test_runs/mwpc_stage_c/cathode_2100_laptop_v1}"
SEL="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
EXE="$BUILD/mwpc_p3_noFW_penning057_cathode_ion_reference_v1"
CHECK="$BUILD/verify_cathode_ion_root"
COMPARE="$BUILD/compare_sampling_to_full_root"
AUDIT="$BUILD/audit_geant4_1000_root"
SELECT="$BUILD/select_detailed_1000_root"
if [[ -f "$HOME/na6p/env.sh" ]];then source "$HOME/na6p/env.sh";fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
set -u
umask 077
fail(){ echo "ERROR: $*" >&2;exit 2; }
lock_data(){
  printf 'git=%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'input=%s\n' "$(sha256sum "$SEL" | awk '{print $1}')"
  printf 'binary=%s\n' "$(sha256sum "$EXE" | awk '{print $1}')"
  printf 'checker=%s\n' "$(sha256sum "$CHECK" | awk '{print $1}')"
  printf 'root=%s\n' "$(root-config --version)"
  printf 'physics=2100V,P3_noFW,ArCO2_70_30,Penning0.57,UNCAPPED\n'
}
verify_lock(){
  [[ -s "$SEL" && "$(wc -l < "$SEL")" -eq 1001 ]] || fail 'Frozen input must have 1000 muons'
  [[ -s "$OUT/config.lock" && -x "$EXE" && -x "$CHECK" ]] || fail 'Run prepare first'
  [[ "$(cat "$OUT/config.lock")" == "$(lock_data)" ]] || fail 'Frozen build/input changed; choose NEW output folder'
}
valid(){
  local d="$1"
  [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" && -s "$d/ion_births.root" ]] || return 1
  (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1) || return 1
  grep -Fxq 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$d/physics.log"
}
prepare(){
  for c in cmake c++ root-config timeout sha256sum;do command -v "$c" >/dev/null || fail "Missing $c";done
  [[ -s "$SEL" && "$(wc -l < "$SEL")" -eq 1001 ]] || fail 'Frozen 1000-muon selection not present'
  mkdir -p "$OUT/geant4_input_audit"
  local -a cmargs=()
  if [[ -n "${Garfield_DIR:-}" ]];then cmargs+=("-DGarfield_DIR=$Garfield_DIR");
  elif [[ -f "$HOME/na6p/install/garfieldpp/lib/cmake/Garfield/GarfieldConfig.cmake" ]];then
    cmargs+=("-DGarfield_DIR=$HOME/na6p/install/garfieldpp/lib/cmake/Garfield")
  fi
  cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" "${cmargs[@]}"
  cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_cathode_ion_reference_v1 \
    verify_cathode_ion_root compare_sampling_to_full_root \
    audit_geant4_1000_root select_detailed_1000_root -j2
  local temp="$OUT/geant4_input_audit/first.csv"
  { head -n1 "$SEL";sed -n '2p' "$SEL"; } > "$temp"
  "$EXE" --validate-input "$temp"
  if [[ ! -s "$OUT/geant4_input_audit/audit.root" ]];then
    "$AUDIT" "$SEL" "$OUT/geant4_input_audit/audit.root"
  fi
  if [[ ! -s "$OUT/geant4_input_audit/selection.root" ]];then
    "$SELECT" "$OUT/geant4_input_audit/audit.root" "$OUT/geant4_input_audit/selection.root"
  fi
  local expected="$(lock_data)"
  if [[ -e "$OUT/config.lock" ]];then
    [[ "$(cat "$OUT/config.lock")" == "$expected" ]] || fail 'Existing frozen run differs'
  else printf '%s\n' "$expected" > "$OUT/config.lock";fi
  echo "PREPARED $OUT"
}
run_one(){
  local mode="$1" idx="$2" k="$3"
  [[ "$idx" =~ ^[0-9]+$ ]] || fail 'Index must be 0..999'
  idx="$((10#$idx))";((idx<1000)) || fail 'Index out of range'
  case "$mode:$k" in weighted:256|weighted:512|weighted:1024|full:0);;*)fail 'Invalid mode/K';;esac
  local timeout_s="${MWPC_LOCAL_TIMEOUT_SEC:-3600}"
  [[ "$timeout_s" =~ ^[0-9]+$ ]] && ((timeout_s>=60)) || fail 'Invalid timeout'
  verify_lock
  local tag dest prior
  printf -v tag 'muon_%04d' "$idx"
  dest="$OUT/${mode}_k${k}/$tag"
  if [[ -f "$dest/complete.txt" ]];then valid "$dest" || fail 'Checksum failure';echo "SKIP $dest";return;fi
  if [[ "$mode" == full ]];then
    prior="$OUT/weighted_k256/$tag"
    valid "$prior" || fail 'First finish weighted K=256 for this same muon'
  fi
  mkdir -p "$(dirname "$dest")";mkdir "$dest" 2>/dev/null || fail "Existing output $dest; no overwrite"
  trap 'rc=$?;if ((rc!=0));then echo "FAILED $rc" > "$dest/FAILED.txt";echo "INCOMPLETE $dest" >&2;fi' EXIT
  awk -v n="$((idx+2))" 'NR==1||NR==n' "$SEL" > "$dest/input.csv"
  [[ "$(wc -l < "$dest/input.csv")" -eq 2 ]] || fail 'Bad one-muon extraction'
  "$EXE" --validate-input "$dest/input.csv" > "$dest/input_check.log" 2>&1
  { echo "index=$idx mode=$mode K=$k start=$(date -u +%FT%TZ)";cat "$OUT/config.lock"; } > "$dest/metadata.txt"
  echo "START $mode K=$k muon=$idx; $dest"
  (
    if [[ -n "${MWPC_LOCAL_MAX_VMEM_KB:-}" ]];then ulimit -v "$MWPC_LOCAL_MAX_VMEM_KB";fi
    /usr/bin/time -v -o "$dest/resource_time.txt" \
      timeout --signal=TERM --kill-after=30s "$timeout_s" \
      "$EXE" --cathode-ion-reference "$dest/input.csv" "$dest/microscopic" \
      2100 0 "$mode" "$k" "$dest/ion_births.root" > "$dest/physics.log" 2>&1
  ) || { tail -n 30 "$dest/physics.log" >&2;fail 'Garfield failed / timed out. No valid result'; }
  grep -Fxq 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$dest/physics.log" || fail 'Uncapped proof missing'
  "$CHECK" "$dest/ion_births.root" "$dest/microscopic" > "$dest/root_audit.log" 2>&1 || {
    cat "$dest/root_audit.log" >&2;fail 'ROOT physics audit failed'; }
  (cd "$dest" && sha256sum input.csv input_check.log ion_births.root microscopic_*.csv \
    physics.log resource_time.txt root_audit.log metadata.txt > SHA256SUMS)
  (cd "$dest" && sha256sum -c SHA256SUMS >/dev/null)
  printf 'PASS HV=2100 mode=%s K=%s original_index=%d\n' "$mode" "$k" "$idx" > "$dest/complete.txt"
  trap - EXIT
  tail -n 2 "$dest/root_audit.log"
  echo "PASS $dest"
}
compare(){
  local idx="$1" k="$2" tag full weighted
  verify_lock
  [[ "$idx" =~ ^[0-9]+$ ]] || fail 'Bad index'
  idx="$((10#$idx))";((idx<1000)) || fail 'Bad index'
  case "$k" in 256|512|1024);;*)fail 'K=256|512|1024';;esac
  printf -v tag 'muon_%04d' "$idx"
  full="$OUT/full_k0/$tag";weighted="$OUT/weighted_k$k/$tag"
  valid "$full" && valid "$weighted" || fail 'Full/weighted result not both verified'
  mkdir -p "$OUT/comparisons"
  "$COMPARE" "$full/ion_births.root" "$weighted/ion_births.root" | tee "$OUT/comparisons/${tag}_k${k}.log"
}
status(){
  echo "Laptop output: $OUT"
  local kind d ok bad
  for kind in weighted_k256 weighted_k512 weighted_k1024 full_k0;do
    ok=0;bad=0
    for d in "$OUT/$kind"/muon_*;do
      [[ -d "$d" ]] || continue
      if valid "$d";then ((ok+=1));else ((bad+=1));fi
    done
    echo "$kind VERIFIED=$ok INCOMPLETE_OR_BAD=$bad"
  done
}
case "${1:-status}" in
  prepare) prepare;;
  weighted) [[ $# -eq 3 ]] || fail 'weighted INDEX K';run_one weighted "$2" "$3";;
  full) [[ $# -eq 2 ]] || fail 'full INDEX';run_one full "$2" 0;;
  compare) [[ $# -eq 3 ]] || fail 'compare INDEX K';compare "$2" "$3";;
  status) status;;
  *) fail 'Usage: prepare|weighted INDEX K|full INDEX|compare INDEX K|status';;
esac
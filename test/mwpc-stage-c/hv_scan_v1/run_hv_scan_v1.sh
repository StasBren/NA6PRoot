#!/usr/bin/env bash
# Complete Garfield++ 1000-muon / HV campaign, ROOT-only physics audits.
# Usage: bash run_hv_scan_v1.sh prepare|pilot|run|analyse|status|chunk VOLTAGE CHUNK_INDEX
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
META="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.json"
OUT="${NA60_HV_OUT:-$REPO/test_runs/mwpc_stage_c/hv_scan_wiretruth_1000_v1}"
BUILD="$REPO/test/mwpc-avalanche/build"
TARGET="mwpc_p3_noFW_penning057_hvscan_v1"
ROOT_ARGS="-l -b -q"
VOLTAGES="${NA60_HV_VOLTAGES:-1700 1750 1800 1850 1900}"
MODE="${1:-status}"

# WSL defaults plus explicit cluster overrides. A cluster may load the
# Garfield environment from CVMFS rather than from the local WSL installation.
if [[ -n "${NA60_HV_GARFIELD_ENV:-}" ]]; then
  [[ -f "$NA60_HV_GARFIELD_ENV" ]] || { echo "Garfield env script not found: $NA60_HV_GARFIELD_ENV" >&2;exit 2; }
  source "$NA60_HV_GARFIELD_ENV"
else
  if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
  if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]]; then
    source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
  fi
fi

require_root(){ command -v root >/dev/null || { echo "ROOT executable not found; source ROOT/Garfield environment" >&2;exit 2; }; }
require_cpp(){ command -v cmake >/dev/null || { echo "cmake missing" >&2;exit 2; }; }

build(){
  require_root;require_cpp
  if [[ ! -f "$BUILD/CMakeCache.txt" ]]; then
    local opts=()
    if [[ -n "${NA60_HV_GARFIELD_DIR:-}" ]]; then
      opts+=("-DGarfield_DIR=$NA60_HV_GARFIELD_DIR")
    elif [[ -d "$HOME/na6p/install/garfieldpp/lib/cmake/Garfield" ]]; then
      opts+=("-DGarfield_DIR=$HOME/na6p/install/garfieldpp/lib/cmake/Garfield")
    fi
    cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" "${opts[@]}"
  fi
  cmake --build "$BUILD" --target "$TARGET" -j2
  test -x "$BUILD/$TARGET"
}

check_voltage(){
  case "$1" in 1700|1750|1800|1850|1900) ;; *) echo "Unsupported HV: $1" >&2;exit 2;; esac
}
check_index(){ [[ "$1" =~ ^[0-9]+$ ]] && (( 10#$1 < 50 )) || { echo 'Chunk must be 0..49' >&2;exit 2; }; }
prepare(){
  mkdir -p "$OUT/inputs" "$OUT/voltages" "$OUT/reports"
  [[ -f "$SELECTION" && -f "$META" ]] || { echo 'Frozen 1000 selection missing' >&2;exit 2; }
  [[ "$(sha256sum "$SELECTION" | cut -d' ' -f1)" == 'e1de4ce0e58645cc833a9979bf6a2a21419adfd4d74a1fa11f78fe887ffb010c' ]] || {
     echo 'Frozen selection SHA256 mismatch' >&2;exit 2;
  }
  [[ "$(wc -l < "$SELECTION")" -eq 1001 ]] || { echo '1000 row count mismatch' >&2; exit 2; }
  cp -n "$META" "$OUT/selection_frozen.json"
  cp -n "$SELECTION" "$OUT/selection_frozen.csv"
  local h;h="$(head -n1 "$SELECTION")"
  for i in $(seq 0 49); do
    local id d tmp
    printf -v id '%04d' "$i"
    d="$OUT/inputs/chunk_${id}.csv";tmp="$d.tmp"
    printf '%s\n' "$h" > "$tmp"
    sed -n "$((2+20*i)),$((21+20*i))p" "$SELECTION" >> "$tmp"
    if [[ -e "$d" ]]; then
      cmp -s "$d" "$tmp" || { echo "Frozen input mismatch: $d" >&2;exit 2; }
      rm -f "$tmp"
    else mv "$tmp" "$d"; fi
    [[ "$(wc -l < "$d")" -eq 21 ]] || { echo "Bad chunk: $d" >&2;exit 2; }
  done
  {
    echo 'PHYSICS=P3_NOFW_ARCO2_7030_PENNING_057_B0'
    echo 'MUON_SAMPLE=FROZEN_1000_REUSED_AT_EACH_VOLTAGE'
    echo 'MICROSCOPIC_ENGINE=Garfield::AvalancheMicroscopic'
    echo 'GAIN_CAP=200000_WARNING_IF_ANY_CAPPED_STOP'
    echo 'HV_VALUES=1700,1750,1800,1850,1900'
    echo 'CONDITIONAL_ISSUE=HEED_SEEDS_SHARED_AT_ALL_VOLTAGES'
    sha256sum "$SELECTION"
  } > "$OUT/campaign_metadata.tmp"
  if [[ -f "$OUT/campaign_metadata.txt" ]]; then
    cmp -s "$OUT/campaign_metadata.tmp" "$OUT/campaign_metadata.txt" || { echo 'Campaign metadata changed: refusing overwrite' >&2;exit 2; }
    rm "$OUT/campaign_metadata.tmp"
  else mv "$OUT/campaign_metadata.tmp" "$OUT/campaign_metadata.txt"; fi
  echo "PREPARED 50 chunks x 20 real Geant4 muons in $OUT"
}

# Recompute only absent runs; never clobber a partial/failed physical output.
run_chunk(){
  local hv="$1" i="$2" id d inp pref
  check_voltage "$hv";check_index "$i"
  printf -v id '%04d' "$i"
  inp="$OUT/inputs/chunk_${id}.csv"
  d="$OUT/voltages/V${hv}/chunk_${id}"
  pref="$d/microscopic"
  if [[ -f "$d/ROOT_AUDIT_PASS.txt" ]]; then
    echo "SKIP sealed V${hv} chunk ${id}";return 0
  fi
  if [[ -e "$d" ]]; then echo "PARTIAL/FAILED OUTPUT: $d -- inspect, do not overwrite" >&2;exit 3;fi
  mkdir -p "$d"
  cp "$inp" "$d/input.csv"
  sha256sum "$inp" "$BUILD/$TARGET" > "$d/input_and_binary_SHA256.txt"
  echo "RUN V${hv} / chunk ${id} 20 crossings"
  date -Is > "$d/started_at.txt"
  "$BUILD/$TARGET" --wiretruth-hv-p3-noFW "$inp" "$pref" "$hv" > "$d/garfield.log" 2>&1 || {
    echo "Garfield failed. See $d/garfield.log" >&2;exit 3;
  }
  # TrackHeed ionisation must be identical at all HV values (fixed muon and seed).
  # This is a strict cross-voltage input equality test, not a physics observable.
  local common="$OUT/voltages/V1700/chunk_${id}/microscopic"
  if [[ "$hv" != 1700 && -f "${common}_heed_primaries.csv" ]]; then
    cmp -s "${common}_heed_primaries.csv" "${pref}_heed_primaries.csv" || {
      echo "FATAL: TrackHeed primaries changed across HV for chunk ${id}" >&2;exit 3;
    }
    cmp -s "${common}_heed_clusters.csv" "${pref}_heed_clusters.csv" || {
      echo "FATAL: TrackHeed clusters changed across HV for chunk ${id}" >&2;exit 3;
    }
  fi
  root $ROOT_ARGS "$HERE/VoltageChunkAudit_ROOT.C(\"$pref\",\"$inp\",20,\"$d/ROOT_AUDIT_PASS.txt\")" \
    > "$d/root_audit.log" 2>&1 || { echo "ROOT audit failed: $d/root_audit.log" >&2;exit 3; }
  test -s "$d/ROOT_AUDIT_PASS.txt" || { echo "Audit seal missing" >&2;exit 3; }
  sha256sum "${pref}"_*.csv > "$d/output_SHA256.txt"
  date -Is > "$d/completed_at.txt"
  echo "PASS V${hv} chunk ${id}"
}

run_pilot(){
  build;prepare
  for hv in $VOLTAGES;do run_chunk "$hv" 0;done
  echo 'PILOT PASSED: inspect V1700..V1900 results and no capped avalanches before full scan'
}
run_all(){
  build;prepare
  for hv in $VOLTAGES;do
    # No full V scan until its pilot has been validated.
    run_chunk "$hv" 0
  done
  for hv in $VOLTAGES;do
    for i in $(seq 1 49);do run_chunk "$hv" "$i";done
  done
  echo "FULL SCAN DONE: $OUT"
}
status(){
  for hv in $VOLTAGES;do
    local n=0
    for i in $(seq 0 49);do
      local id;printf -v id '%04d' "$i"
      [[ -f "$OUT/voltages/V${hv}/chunk_${id}/ROOT_AUDIT_PASS.txt" ]] && n=$((n+1)) || true
    done
    echo "HV ${hv}: sealed ${n}/50 chunks, muons=$((n*20))"
  done
}
case "$MODE" in
  prepare) build;prepare;;
  pilot) run_pilot;;
  run) run_all;;
  chunk)
    require_root
    [[ $# -eq 3 ]] || { echo 'chunk VOLTAGE CHUNK_INDEX' >&2;exit 2; }
    [[ -x "$BUILD/$TARGET" ]] || { echo 'Build once with: bash run_hv_scan_v1.sh prepare' >&2;exit 2; }
    [[ -f "$OUT/inputs/chunk_0000.csv" ]] || { echo 'Prepare the immutable chunks before PBS submission' >&2;exit 2; }
    run_chunk "$2" "$3";;
  status) status;;
  analyse) require_root;root $ROOT_ARGS "$HERE/AnalyzeVoltageScan_ROOT.C(\"$OUT\")";;
  *) echo 'Usage: run_hv_scan_v1.sh prepare|pilot|run|chunk V chunk|status|analyse' >&2;exit 2;;
esac

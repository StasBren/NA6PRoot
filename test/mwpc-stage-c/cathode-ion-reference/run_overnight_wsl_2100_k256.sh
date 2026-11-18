#!/usr/bin/env bash
# NA60+ Stage C — resumable 2100V K256 weighted ion birth WSL campaign.
# Genuine 1000 frozen Geant4 crossings, sequential FULL uncapped Garfield
# avalanche transport, one ROOT archive per muon, per-event physics audit.
# Timeout/memory/disk guards stop failed jobs: NEVER treat a cap as physics.
set -Eeuo pipefail
umask 027
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
SELECT="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BUILD="$REPO/test/mwpc-avalanche/build-overnight-2100-k256"
EXE="$BUILD/mwpc_p3_noFW_penning057_cathode_ion_reference_v1"
CHECK="$BUILD/verify_cathode_ion_root"
OUT="${MWPC_OVERNIGHT_OUT:-$REPO/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
MODE="${1:-status}"
MAX_EVENTS="${MWPC_OVERNIGHT_MAX_EVENTS:-1000}"
MAX_HOURS="${MWPC_OVERNIGHT_MAX_HOURS:-12}"
TIMEOUT_MIN="${MWPC_OVERNIGHT_EVENT_TIMEOUT_MIN:-90}"
KEEP_GB="${MWPC_OVERNIGHT_KEEP_FREE_GB:-8}"
VM_GB="${MWPC_OVERNIGHT_MAX_VM_GB:-12}"
GARFIELD_INSTALL="${GARFIELD_INSTALL:-$HOME/na6p/install/garfieldpp}"
export GARFIELD_INSTALL
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
export HEED_DATABASE="${HEED_DATABASE:-$GARFIELD_INSTALL/share/Heed/database}"
fail(){ echo "OVERNIGHT STOP: $*" >&2; exit 2; }
number(){ [[ "$1" =~ ^[0-9]+$ ]]; }
for x in "$MAX_EVENTS" "$MAX_HOURS" "$TIMEOUT_MIN" "$KEEP_GB" "$VM_GB"; do
  number "$x" || fail "invalid numeric limit $x"
done
(( MAX_EVENTS>=1 && MAX_EVENTS<=1000 && MAX_HOURS>=1 && MAX_HOURS<=168 && TIMEOUT_MIN>=1 && TIMEOUT_MIN<=720 && KEEP_GB>=2 && VM_GB>=4 )) || fail "unsafe limits"
[[ -s "$SELECT" && "$(wc -l < "$SELECT")" -eq 1001 ]] || fail "missing or incomplete 1000-entry selection"
[[ -f "$HEED_DATABASE/cbdel.dat" ]] || fail "HEED_DATABASE not found: $HEED_DATABASE/cbdel.dat"
[[ -d "$GARFIELD_INSTALL/lib" ]] || fail "Garfield++ libraries missing"
command -v sha256sum >/dev/null && command -v flock >/dev/null && command -v timeout >/dev/null || fail "required coreutils missing"
lock_text(){
  printf 'git_sha=%s\nselection_sha256=%s\nbinary_sha256=%s\nverifier_sha256=%s\n' \
    "$(git -C "$REPO" rev-parse HEAD)" "$(sha256sum "$SELECT"|cut -d' ' -f1)" \
    "$(sha256sum "$EXE"|cut -d' ' -f1)" "$(sha256sum "$CHECK"|cut -d' ' -f1)"
  printf 'configuration=2100V;noFW;3mm_pitch;pm3mm_gap;ArCO2_70_30;Penning0.57;B0;uncapped;K256\n'
}
check_lock(){
  [[ -s "$OUT/config.lock" && -x "$EXE" && -x "$CHECK" ]] || fail "run prepare first"
  [[ "$(cat "$OUT/config.lock")" == "$(lock_text)" ]] || fail "binary/Git/input changed after configuration freeze"
}
input_for(){
  local i="$1" dest="$2"
  awk -v row="$((i+2))" 'NR==1 || NR==row' "$SELECT" > "$dest"
  [[ "$(wc -l < "$dest")" -eq 2 ]] || fail "invalid input index $i"
}
finished(){
  local d="$1"
  [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" && -s "$d/ion_births.root" ]] || return 1
  (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1) || return 1
  grep -qx 'PASS CATHODE ROOT PHYSICS AUDIT' "$d/root_audit.log" || return 1
  grep -qx 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$d/physics.log" || return 1
  [[ "$(wc -l < "$d/microscopic_crossing.csv")" -eq 2 ]] || return 1
}
disk_ok(){
  local available
  available="$(df -PB1 "$OUT" | awk 'NR==2{print $4}')"
  [[ "$available" =~ ^[0-9]+$ ]] || fail "cannot inspect disk"
  (( available >= KEEP_GB*1024*1024*1024 )) || fail "free disk below ${KEEP_GB}GiB reserve"
}
prepare(){
  command -v cmake >/dev/null || fail "cmake not found"
  mkdir -p "$OUT"
  if [[ -e "$OUT/config.lock" ]];then
    check_lock
  else
    cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" "-DCMAKE_PREFIX_PATH=$GARFIELD_INSTALL"
    cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_cathode_ion_reference_v1 verify_cathode_ion_root -j2
    [[ -x "$EXE" && -x "$CHECK" ]] || fail "missing compiled binary"
    [[ ! -e "$OUT/.input_preflight.csv" ]] || fail "stale preflight input exists"
    input_for 0 "$OUT/.input_preflight.csv"
    "$EXE" --validate-input "$OUT/.input_preflight.csv"
    rm "$OUT/.input_preflight.csv"
    lock_text > "$OUT/config.lock"
  fi
  disk_ok
  echo "PREPARED $OUT; 1000 genuine muon selection, uncapped microscopic Garfield, per-primary K256."
}
status(){
  [[ -d "$OUT" ]] || fail "no output directory"
  local good=0 bad=0 missing=0 idx tag d
  for ((idx=0;idx<MAX_EVENTS;++idx));do
    printf -v tag 'muon_%04d' "$idx"
    d="$OUT/$tag"
    if finished "$d";then good=$((good+1))
    elif [[ -d "$d" ]];then bad=$((bad+1));if ((bad<=5));then echo "UNSEALED $d";fi
    else missing=$((missing+1));fi
  done
  echo "STATUS VERIFIED=$good/$MAX_EVENTS ABSENT=$missing UNSEALED=$bad DATA=$OUT"
}
run(){
  check_lock
  exec 9>"$OUT/.runner.lock"
  flock -n 9 || fail "another worker holds the run lock"
  local started idx tag dir rc completed=0
  started="$(date +%s)"
  echo "START $(date -Is) max_events=$MAX_EVENTS max_hours=$MAX_HOURS; PID=$$"
  for ((idx=0;idx<MAX_EVENTS;++idx));do
    if (( $(date +%s)-started >= MAX_HOURS*3600 ));then
      echo "STOP: time budget reached before event $idx";break
    fi
    check_lock
    printf -v tag 'muon_%04d' "$idx"
    dir="$OUT/$tag"
    if finished "$dir";then completed=$((completed+1));continue;fi
    [[ ! -e "$dir" ]] || fail "refusing to overwrite incomplete/invalid $dir"
    disk_ok
    mkdir "$dir"
    input_for "$idx" "$dir/input.csv"
    {
      printf 'index=%s\nstart_utc=%s\n' "$idx" "$(date -u -Is)"
      printf 'lock_sha256=%s\n' "$(sha256sum "$OUT/config.lock"|cut -d' ' -f1)"
      printf 'uncapped_Garfield=1\nreservoir_per_primary=256\n'
    } > "$dir/task_metadata.txt"
    echo "BEGIN $tag $(date -Is)"
    rc=0
    (
      ulimit -v "$((VM_GB*1024*1024))"
      timeout --signal=TERM --kill-after=30s "${TIMEOUT_MIN}m" \
        "$EXE" --cathode-ion-reference "$dir/input.csv" "$dir/microscopic" \
        2100 0 weighted 256 "$dir/ion_births.root" > "$dir/physics.log" 2>&1
    ) || rc=$?
    if ((rc!=0));then
      printf 'FAILED exit=%s\n' "$rc" > "$dir/FAILED"
      tail -n 25 "$dir/physics.log" || true
      fail "microscopic $tag failed (timeout, out-of-memory or Garfield); no truncated physics accepted"
    fi
    if ! "$CHECK" "$dir/ion_births.root" "$dir/microscopic" > "$dir/root_audit.log" 2>&1;then
      echo "FAILED root audit" > "$dir/FAILED"
      tail -n 30 "$dir/root_audit.log" || true
      fail "ROOT physics audit failed for $tag"
    fi
    grep -qx 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$dir/physics.log" || fail "$tag did not prove uncapped Garfield"
    [[ "$(wc -l < "$dir/microscopic_crossing.csv")" -eq 2 ]] || fail "$tag missing complete crossing"
    for item in input.csv physics.log root_audit.log ion_births.root task_metadata.txt microscopic_crossing.csv microscopic_seed.csv microscopic_wire_seed.csv microscopic_wire_muon.csv;do
      [[ -s "$dir/$item" ]] || fail "missing $tag/$item"
    done
    (cd "$dir" && sha256sum input.csv physics.log root_audit.log ion_births.root task_metadata.txt microscopic_*.csv > SHA256SUMS)
    printf 'VERIFIED INDEX=%s uncapped=1 weighted_K256\n' "$idx" > "$dir/complete.txt"
    completed=$((completed+1))
    echo "DONE $tag $(date -Is), count_this_run=$completed"
  done
  echo "STOP $(date -Is) completed_or_skipped=$completed"
  status
}
case "$MODE" in
 prepare) prepare;;
 run) run;;
 status) status;;
 *) fail "usage: $0 prepare|run|status";;
esac

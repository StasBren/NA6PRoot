#!/usr/bin/env bash
# NA60+ Stage C — resumable 2100V K256 weighted ion birth WSL campaign.
# Genuine 1000 frozen Geant4 crossings, sequential FULL uncapped Garfield
# avalanche transport, one ROOT archive per muon, per-event physics audit.
# Per-event physics, audit, and timeout failures are recorded and SKIPPED.
# Configuration/disk failure stops the campaign. NEVER cap avalanche size.
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
MAX_HOURS="${MWPC_OVERNIGHT_MAX_HOURS:-0}"
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
(( MAX_EVENTS>=1 && MAX_EVENTS<=1000 && MAX_HOURS>=0 && MAX_HOURS<=168 && TIMEOUT_MIN>=1 && TIMEOUT_MIN<=720 && KEEP_GB>=2 && VM_GB>=4 )) || fail "unsafe limits"
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
  echo "PREPARED $OUT; 1000 muons, continue-after-event-errors, uncapped microscopic Garfield, K256."
}
status(){
  [[ -d "$OUT" ]] || fail "no output directory"
  local good=0 failed=0 unsealed=0 missing=0 idx tag d
  for ((idx=0;idx<MAX_EVENTS;++idx));do
    printf -v tag 'muon_%04d' "$idx"
    d="$OUT/$tag"
    if finished "$d";then good=$((good+1))
    elif [[ -s "$d/FAILED" ]];then failed=$((failed+1))
    elif [[ -d "$d" ]];then unsealed=$((unsealed+1))
    else missing=$((missing+1));fi
  done
  echo "STATUS VERIFIED=$good/$MAX_EVENTS FAILED=$failed UNSEALED=$unsealed ABSENT=$missing DATA=$OUT"
}
mark_event_failed(){
  local dir="$1" reason="$2" tag
  tag="$(basename "$dir")"
  # NEVER erase original simulation or audit logs; failed != valid physics.
  printf 'classification=%s\nindex=%s\nrecorded_utc=%s\n' \
    "$reason" "$tag" "$(date -u -Is)" > "$dir/FAILED"
  echo "SKIP_FAILED $tag reason=$reason ; continuing to next input index"
}
run(){
  # Blocks safely if the ORIGINAL runner is processing an event now.
  # Takes over after the old worker exits; only one worker runs at a time.
  mkdir -p "$OUT"
  exec 9>"$OUT/.runner.lock"
  echo "WAIT_FOR_RUNNING_WORKER $(date -Is)"
  flock 9 || fail "cannot acquire campaign lock"
  check_lock
  local started idx tag dir rc audit_rc item missing_item
  local good=0 failed=0 existing_failed=0 existing_unsealed=0
  started="$(date +%s)"
  echo "RESILIENT_START $(date -Is) max_events=$MAX_EVENTS max_hours=$MAX_HOURS (0=no time limit); PID=$$"
  for ((idx=0;idx<MAX_EVENTS;++idx));do
    if (( MAX_HOURS>0 && $(date +%s)-started >= MAX_HOURS*3600 ));then
      echo "RESILIENT_STOP time budget reached before event $idx"
      break
    fi
    check_lock
    printf -v tag 'muon_%04d' "$idx"
    dir="$OUT/$tag"
    if finished "$dir";then good=$((good+1));continue;fi
    if [[ -e "$dir" ]];then
      if [[ -s "$dir/FAILED" ]];then
        existing_failed=$((existing_failed+1))
        echo "SKIP_PREVIOUS_FAILURE $tag (preserving evidence)"
      else
        existing_unsealed=$((existing_unsealed+1))
        echo "SKIP_PREVIOUS_UNSEALED $tag (preserving evidence)"
      fi
      continue
    fi
    # These are global safeguards: STOP on disk exhaustion or altered physics.
    disk_ok
    mkdir "$dir" || fail "cannot allocate event directory $dir"
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
      timeout --signal=TERM --kill-after=30s "$TIMEOUT_MIN"m \
        "$EXE" --cathode-ion-reference "$dir/input.csv" "$dir/microscopic" \
        2100 0 weighted 256 "$dir/ion_births.root" > "$dir/physics.log" 2>&1
    ) || rc=$?
    if ((rc!=0));then
      mark_event_failed "$dir" "microscopic_exit_$rc"
      tail -n 8 "$dir/physics.log" || true
      failed=$((failed+1))
      continue
    fi
    audit_rc=0
    "$CHECK" "$dir/ion_births.root" "$dir/microscopic" \
      > "$dir/root_audit.log" 2>&1 || audit_rc=$?
    if ((audit_rc!=0));then
      mark_event_failed "$dir" "root_audit_exit_$audit_rc"
      tail -n 8 "$dir/root_audit.log" || true
      failed=$((failed+1))
      continue
    fi
    if ! grep -qx 'PASS CATHODE ROOT PHYSICS AUDIT' "$dir/root_audit.log";then
      mark_event_failed "$dir" "root_audit_no_pass"
      failed=$((failed+1))
      continue
    fi
    if ! grep -qx 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$dir/physics.log";then
      mark_event_failed "$dir" "no_uncapped_proof"
      failed=$((failed+1))
      continue
    fi
    if [[ ! -s "$dir/microscopic_crossing.csv" ]] ||
       [[ "$(wc -l < "$dir/microscopic_crossing.csv")" -ne 2 ]];then
      mark_event_failed "$dir" "incomplete_crossing"
      failed=$((failed+1))
      continue
    fi
    missing_item=""
    for item in input.csv physics.log root_audit.log ion_births.root \
      task_metadata.txt microscopic_crossing.csv microscopic_seed.csv \
      microscopic_wire_seed.csv microscopic_wire_muon.csv;do
      if [[ ! -s "$dir/$item" ]];then missing_item="$item";break;fi
    done
    if [[ -n "$missing_item" ]];then
      mark_event_failed "$dir" "missing_$missing_item"
      failed=$((failed+1))
      continue
    fi
    if ! (cd "$dir" && sha256sum input.csv physics.log root_audit.log \
       ion_births.root task_metadata.txt microscopic_*.csv > SHA256SUMS \
       && sha256sum -c SHA256SUMS > /dev/null 2>&1);then
      mark_event_failed "$dir" "checksum_failure"
      failed=$((failed+1))
      continue
    fi
    printf 'VERIFIED INDEX=%s uncapped=1 weighted_K256\n' "$idx" > "$dir/complete.txt"
    good=$((good+1))
    echo "DONE $tag $(date -Is) newly_verified=$good"
  done
  echo "RESILIENT_END $(date -Is) verified_or_skipped=$good newly_failed=$failed previously_failed=$existing_failed previously_unsealed=$existing_unsealed"
  status
}
case "$MODE" in
 prepare) prepare;;
 run) run;;
 status) status;;
 *) fail "usage: $0 prepare|run|status";;
esac

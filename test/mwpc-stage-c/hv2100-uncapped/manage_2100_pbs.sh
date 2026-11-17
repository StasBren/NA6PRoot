#!/usr/bin/env bash
# Prepare, submit and audit 1000 true Geant4 muons at 2100 V, UNCAPPED Garfield++.
# Usage (on PBS login): bash manage_2100_pbs.sh prepare|pilot|status|full
# Never modifies 1800-V references, any other HV scan or prior PBS campaigns.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BUILD="$REPO/test/mwpc-avalanche/build-hv2100-uncapped"
EXE="$BUILD/mwpc_p3_noFW_penning057_2100_uncapped_v1"
OUT="${MWPC_2100_OUTPUT_DIR:-/storage/cgrp/stanis/mwpc_penning057/hv2100_uncapped_ion128_1000_v2}"
MODE="${1:-status}"

if [[ -f "$HOME/na6p/env.sh" ]];then source "$HOME/na6p/env.sh";fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
need_input() {
  [[ -s "$SELECTION" ]] || { echo "Missing frozen 1000 Geant4 selection: $SELECTION" >&2;exit 2; }
  local lines
  lines="$(wc -l < "$SELECTION")"
  [[ "$lines" -eq 1001 ]] || {
    echo "EXPECTED 1000 frozen muons (+header); actual lines=$lines" >&2;exit 2;
  }
}
lock_data() {
  local input_sha binary_sha git_sha
  input_sha="$(sha256sum "$SELECTION" | awk '{print $1}')"
  binary_sha="$(sha256sum "$EXE" | awk '{print $1}')"
  git_sha="$(git -C "$REPO" rev-parse HEAD)"
  printf 'selection_sha256=%s\nbinary_sha256=%s\ngit_revision=%s\n' \
      "$input_sha" "$binary_sha" "$git_sha"
  printf 'voltage_V=2100\navalanche_size_limit=DISABLED\nion_birth_reservoir_per_primary=128\nmuons=1000\n'
  printf 'configuration=P3_noFW;30um;3mm;gap_pm3mm;ArCO2_70_30;Penning0.57;B0\n'
}
check_lock() {
  [[ -x "$EXE" && -s "$OUT/config.lock" ]] || {
    echo "Missing compiled executable or config.lock; run 'prepare' first" >&2;exit 3;
  }
  [[ "$(cat "$OUT/config.lock")" == "$(lock_data)" ]] || {
    echo "Frozen config differs (source/input/binary). Do not mix versions; use a new output folder." >&2
    exit 3
  }
}
prepare() {
  need_input
  command -v cmake >/dev/null || { echo "CMake unavailable on this host" >&2;exit 3; }
  command -v c++ >/dev/null || { echo "C++ compiler unavailable on this host" >&2;exit 3; }
  mkdir -p "$OUT/pbs-logs"
  echo "Cluster host: $(hostname -f)"
  echo "Build directory: $BUILD"
  cmake_args=()
  if [[ -n "${Garfield_DIR:-}" ]];then
    cmake_args+=("-DGarfield_DIR=$Garfield_DIR")
  elif [[ -f "$HOME/na6p/install/garfieldpp/lib/cmake/Garfield/GarfieldConfig.cmake" ]];then
    cmake_args+=("-DGarfield_DIR=$HOME/na6p/install/garfieldpp/lib/cmake/Garfield")
  fi
  cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" "${cmake_args[@]}"
  cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_2100_uncapped_v1 -j2
  [[ -x "$EXE" ]] || { echo "Uncapped executable was not produced" >&2;exit 3; }
  # The legacy CLI accepts at most 24 rows. Validate the first frozen row.
  tmp="$(mktemp "$OUT/.preflight.XXXXXX.csv")"
  { head -n1 "$SELECTION"; sed -n '2p' "$SELECTION"; } > "$tmp"
  "$EXE" --validate-input "$tmp"
  rm -f "$tmp"
  local expected
  expected="$(lock_data)"
  if [[ -e "$OUT/config.lock" ]];then
    [[ "$(cat "$OUT/config.lock")" == "$expected" ]] || {
      echo "Refusing to change already-frozen run config.lock" >&2;exit 3;
    }
  else
    printf '%s\n' "$expected" > "$OUT/config.lock"
  fi
  echo "PREPARED: $OUT (1000 frozen muons, 2100 V, Garfield avalanche cap DISABLED)"
  echo "Selection $(sha256sum "$SELECTION")"
  echo "Binary $(sha256sum "$EXE")"
}
audit_one() {
  local i="$1" d
  printf -v d '%s/muon_%04d' "$OUT" "$i"
  [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" ]] || return 1
  (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1) || return 1
  [[ "$(wc -l < "$d/microscopic_crossing.csv")" -eq 2 ]] || return 1
  grep -q '^GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1$' "$d/physics.log" || return 1
  grep -q '^ION_BIRTH_RESERVOIR_PER_PRIMARY=128$' "$d/physics.log" || return 1
  [[ -s "$d/microscopic_ion_birth_sample.csv" ]] || return 1
  head -n1 "$d/microscopic_ion_birth_sample.csv" | grep -q 'nearest_wire_u_global_mm,reservoir_sample_size' || return 1
  return 0
}
status() {
  echo "Output root: $OUT"
  if command -v qstat >/dev/null;then qstat -u "${USER:-stanis}" || true;fi
  local ok=0 missing=0 bad=0 i d
  for ((i=0;i<1000;++i));do
    if audit_one "$i";then
      ((ok+=1))
    else
      printf -v d '%s/muon_%04d' "$OUT" "$i"
      if [[ -d "$d" ]];then
        ((bad+=1))
        if ((bad<=12));then echo "UNSEALED/CORRUPT $i $d";fi
      else
        ((missing+=1))
      fi
    fi
  done
  echo "VERIFIED_COMPLETE=$ok /1000; ABSENT=$missing; UNSEALED_OR_INVALID=$bad"
  if ((ok==1000));then echo "PASS: every muon SHA-verified";else echo "NOT COMPLETE";fi
}
submit() {
  need_input
  check_lock
  command -v qsub >/dev/null || { echo "OpenPBS qsub unavailable" >&2;exit 3; }
  local range="$1" label="$2"
  if [[ "$label" == "full" ]];then
    local i
    for ((i=0;i<5;++i));do
      audit_one "$i" || {
        echo "Pilot muon $i NOT SHA-verified; refusing to submit 995 full jobs" >&2
        exit 4
      }
    done
  fi
  local selection_sha binary_sha git_sha vars
  selection_sha="$(sha256sum "$SELECTION" | awk '{print $1}')"
  binary_sha="$(sha256sum "$EXE" | awk '{print $1}')"
  git_sha="$(git -C "$REPO" rev-parse HEAD)"
  vars="MWPC_2100_REPO=$REPO,MWPC_2100_OUT=$OUT,MWPC_2100_SELECTION_SHA=$selection_sha,MWPC_2100_BINARY_SHA=$binary_sha,MWPC_2100_GIT_SHA=$git_sha"
  echo "Submitting $label array [$range] to OpenPBS; per-task resources: 1 CPU, 32 GB, 12 h"
  qsub -J "$range" -v "$vars" -o "$OUT/pbs-logs" "$HERE/pbs_one_muon_2100.sh"
  echo "Check with: bash $HERE/manage_2100_pbs.sh status"
}
case "$MODE" in
  prepare) prepare;;
  pilot) submit '0-4%2' pilot;;
  full) submit '5-999%8' full;;
  status) status;;
  *) echo "Usage: bash $0 {prepare|pilot|status|full}" >&2;exit 2;;
esac

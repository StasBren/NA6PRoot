#!/usr/bin/env bash
# Survey 2100V K256 997/1000 archived Garfield events without re-simulation.
# Modes: prepare (ROOT/C++ compile only), survey (40-event frozen selection).
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
# WSL uses the existing local Garfield++/ROOT setup; WIPP can use CVMFS.
if [[ -f "$HOME/na6p/env.sh" ]];then source "$HOME/na6p/env.sh";fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
if [[ -n "${GARFIELD_INSTALL:-}" &&
      -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
umask 077
MODE="${1:-survey}"
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
OUTPUT="${MWPC_SIGNAL40_OUTPUT:-$REPO/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
BUILD="$REPO/test/mwpc-avalanche/build-signal40-root"
BIN="$BUILD/mwpc_survey_and_select_2100_signal"
fail(){ echo "SIGNAL40 SETUP FAILED: $*" >&2;exit 2; }
case "$MODE" in
  prepare)
    for p in cmake root-config c++;do command -v "$p" >/dev/null || fail "missing $p";done
    cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD"
    cmake --build "$BUILD" --target mwpc_survey_and_select_2100_signal -j2
    [[ -x "$BIN" ]] || fail "ROOT/C++ survey target not built"
    echo "PASS C++/ROOT 2100V phase-space survey compiled (no Garfield jobs)"
    ;;
  survey)
    [[ -d "$ARCHIVE" ]] || fail "missing 2100V real weighted K256 archive: $ARCHIVE"
    [[ -x "$BIN" ]] || fail "run 'bash $0 prepare' first"
    [[ ! -e "$OUTPUT" ]] || fail "output exists; use a new MWPC_SIGNAL40_OUTPUT, never overwrite"
    "$BIN" "$ARCHIVE" "$OUTPUT"
    [[ -s "$OUTPUT/selected40.csv" && -s "$OUTPUT/survey.root" ]] ||
      fail "ROOT survey output not produced"
    # Independently SHA-verify all 40 actual development/heldout ROOT archives.
    count=0
    while IFS=, read -r role idx rest;do
      [[ "$role" == role ]] && continue
      [[ "$idx" =~ ^[0-9]+$ ]] || fail "malformed original index $idx"
      printf -v tag 'muon_%04d' "$idx"
      d="$ARCHIVE/$tag"
      [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" ]] ||
        fail "unsealed $d (cannot include in 40 selected)"
      (cd "$d" && sha256sum -c SHA256SUMS >/dev/null) ||
        fail "SHA mismatch in actual 2100V signal40 selected $d"
      ((count+=1))
    done < "$OUTPUT/selected40.csv"
    [[ "$count" -eq 40 ]] || fail "expected exactly 40 SHA-verified selected muons"
    (
      cd "$OUTPUT"
      sha256sum selected40.csv development20.csv sealed_validation20.csv \
          survey.root archive_inventory.csv \
          channel_station_coverage.csv excluded_slots.txt \
          previously_inspected_not_sealed.txt selection_scope.txt \
          > FROZEN_SOURCE_SHA256.txt
      git -C "$REPO" rev-parse HEAD > SOURCE_GIT_SHA.txt
    )
    echo "PASS 40 unique chosen events SHA-verified, original K256 archive untouched"
    echo "FROZEN_SELECTION $OUTPUT/selected40.csv"
    echo "COVERAGE $OUTPUT/channel_station_coverage.csv"
    echo "Do NOT use sealed_validation for training or model parameter tuning."
    ;;
  *) fail "Usage: bash run_signal40_selection.sh prepare|survey";;
esac

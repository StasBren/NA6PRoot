#!/usr/bin/env bash
# V0 archive statistics for the frozen 2100 V K256 campaign (overnight_2100_weighted_k256_v1).
# READ-ONLY on archives. Writes only to a new, timestamped output directory. Never opens sealed validation slots.
#
#   bash run_v0_archive_stats.sh prepare      # build mwpc_v0_archive_stats with the local ROOT (no Garfield run)
#   bash run_v0_archive_stats.sh selftest     # synthetic estimator checks (not NA60+ data)
#   bash run_v0_archive_stats.sh schema       # --dump-schema on the first development archive
#   bash run_v0_archive_stats.sh pilot        # the 20 J/psi development archives of development20.csv
#   bash run_v0_archive_stats.sh population   # every sealed-file archive except sealed_validation20.csv
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
if [[ -f "$HOME/na6p/env.sh" ]];then source "$HOME/na6p/env.sh";fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
DEV="$SELECTION/development20.csv"
SEALED="$SELECTION/sealed_validation20.csv"
BUILD="${MWPC_V0_BUILD:-$REPO/test/mwpc-avalanche/build-v0-root}"
BIN="$BUILD/mwpc_v0_archive_stats"
OUTROOT="${MWPC_V0_OUT:-$REPO/test_runs/mwpc_stage_c/v0_archive_stats_2100}"
fail(){ echo "V0 RUNNER FAILED: $*" >&2; exit 2; }

column_index(){ # prints the 'index' column of a selection CSV, zero padded; never prints other columns
  awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i; if(!c)exit 3; next} NF{printf "%04d\n",$c}' "$1"
}
verified(){ local d="$ARCHIVE/muon_$1"
  [[ -s "$d/complete.txt" && -s "$d/SHA256SUMS" && -s "$d/ion_births.root" && -s "$d/microscopic_wire_seed.csv" ]] || return 1
  (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1)
}
new_out(){ local o="$OUTROOT/$1_$(date +%Y%m%d_%H%M%S)"; [[ ! -e "$o" ]] || fail "output exists"; mkdir -p "$o"; echo "$o"; }
provenance(){ { echo "git_sha=$(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "git_dirty=$(git -C "$REPO" status --porcelain 2>/dev/null | wc -l)"
  echo "binary_sha256=$(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "archive_root=$ARCHIVE"; echo "selection_dir=$SELECTION"
  [[ -s "$SEALED" ]] && echo "sealed_list_sha256=$(sha256sum "$SEALED" | cut -d' ' -f1)"
  [[ -s "$DEV" ]] && echo "development_list_sha256=$(sha256sum "$DEV" | cut -d' ' -f1)"
  echo "root_version=$(root-config --version 2>/dev/null || echo unknown)"; } > "$1/provenance.txt"; }

case "${1:-}" in
  prepare)
    for p in cmake root-config c++;do command -v "$p" >/dev/null || fail "missing $p";done
    cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD"
    cmake --build "$BUILD" --target mwpc_v0_archive_stats -j2
    [[ -x "$BIN" ]] || fail "mwpc_v0_archive_stats not built"
    echo "PASS built $BIN (ROOT $(root-config --version))";;
  selftest)
    [[ -x "$BIN" ]] || fail "run prepare first"
    OUT="$(new_out selftest)"
    "$BIN" --synthetic 300 --seed 1 --out-dir "$OUT/independent" > "$OUT/independent.txt"
    "$BIN" --synthetic 300 --seed 2 --rho-side 0.25 --out-dir "$OUT/correlated" > "$OUT/correlated.txt"
    grep -q '^PASS census-CSV reader round trip' "$OUT/independent.txt" || fail "census CSV round trip"
    echo "PASS synthetic self-test (estimator checks only, not NA60+ data) -> $OUT";;
  schema)
    [[ -x "$BIN" && -s "$DEV" ]] || fail "need built binary and $DEV"
    first="$(column_index "$DEV" | head -1)"; verified "$first" || fail "development muon_$first not SHA-verified"
    OUT="$(new_out schema)"; "$BIN" --dump-schema "$ARCHIVE/muon_$first/ion_births.root" | tee "$OUT/schema_muon_$first.txt"
    head -1 "$ARCHIVE/muon_$first/microscopic_wire_seed.csv" > "$OUT/wire_seed_header.txt"
    grep -q '^PASS V0 required-branch schema check' "$OUT/schema_muon_$first.txt" || fail "schema mismatch, see $OUT"
    echo "OUTPUT $OUT";;
  pilot|population)
    [[ -x "$BIN" ]] || fail "run prepare first"
    [[ -s "$SEALED" && -s "$DEV" ]] || fail "selection files missing in $SELECTION (run run_signal40_selection.sh survey first)"
    mapfile -t SEAL < <(column_index "$SEALED"); [[ ${#SEAL[@]} -eq 20 ]] || fail "expected 20 sealed indices"
    declare -A IS_SEALED=(); for i in "${SEAL[@]}"; do IS_SEALED[$i]=1; done
    OUT="$(new_out "$1")"; LIST="$OUT/input_list.txt"; : > "$LIST"; n_unverified=0
    if [[ "$1" == pilot ]]; then mapfile -t IDX < <(column_index "$DEV")
    else mapfile -t IDX < <(cd "$ARCHIVE" && ls -d muon_[0-9][0-9][0-9][0-9] | sed 's/muon_//'); fi
    for i in "${IDX[@]}"; do
      [[ -n "${IS_SEALED[$i]:-}" ]] && continue               # sealed slots are skipped BEFORE any file access
      if verified "$i"; then echo "$ARCHIVE/muon_$i/ion_births.root" >> "$LIST"; else n_unverified=$((n_unverified+1)); fi
    done
    echo "unverified_or_incomplete_skipped=$n_unverified" > "$OUT/selection_note.txt"
    [[ -s "$LIST" ]] || fail "no verified archives"
    provenance "$OUT"
    "$BIN" --list "$LIST" --exclude-index-csv "$SEALED" --out-dir "$OUT" | tee "$OUT/console.txt"
    echo "OUTPUT $OUT  (return v0_summary.txt, v0_manifest.txt, provenance.txt, selection_note.txt)";;
  *) echo "usage: $0 prepare|selftest|schema|pilot|population"; exit 1;;
esac

#!/usr/bin/env bash
# 2100V 997-muon original source data -> Claude compact ZIP.
# No Garfield, no ROOT/Garfield rebuild, no source-file modification.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
VENV="${NA60_PLOTSTUDIO_VENV:-$HOME/na6p/venvs/na60-plotstudio-2100}"
OUTBASE="${MWPC_GAINASYM_OUT:-$REPO/test_runs/mwpc_stage_c/2100_gain_asym}"
DEST="${MWPC_CLAUDE_2100_ZIP:-$HOME/NA60_2100V_997_for_Claude.zip}"
MODE="${1:-locate}"
fail(){ echo "NA60+ CLAUDE HANDOFF: $*" >&2;exit 2; }
pick(){
  shopt -s nullglob
  local d
  for d in $(ls -dt "$OUTBASE"/population_* 2>/dev/null); do
    if [[ -s "$d/census/v0_muons.csv" && -s "$d/census/v0_fragments.csv" &&
          -s "$d/rootfiles.txt" ]]; then
      printf '%s\n' "$d"
      return 0
    fi
  done
  fail "No full population export. First run: bash test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_gain_asym_2100.sh population"
}
case "$MODE" in
  locate)
    echo "ORIGINAL=$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1"
    echo "V0_ANALYSIS_BASE=$OUTBASE"
    echo "DESTINATION=$DEST"
    pick;;
  check)
    [[ -x "$VENV/bin/python" ]] || fail "Plot Studio venv missing: $VENV"
    "$VENV/bin/python" -m py_compile "$HERE/claude_handoff_2100.py"
    "$VENV/bin/python" -c 'import pandas,numpy; print("PASS pandas and numpy")'
    echo "SOURCE=$(pick)"
    echo "PASS Python syntax and compact handoff dependencies";;
  export)
    [[ -x "$VENV/bin/python" ]] || fail "Plot Studio venv missing; run Plot Studio install"
    [[ ! -e "$DEST" ]] || fail "Output exists (no overwrite): $DEST"
    SRC="$(pick)"
    "$VENV/bin/python" "$HERE/claude_handoff_2100.py" --run "$SRC" --out "$DEST"
    echo "Use in Windows: explorer.exe \"$(wslpath -w "$(dirname "$DEST")")\"";;
  *)
    fail "Usage: bash $0 locate|check|export";;
esac

#!/usr/bin/env bash
# 1800V 5620 muons: build/prepare/submit/check one-file Heed+Garfield ROOT.
set -Eeo pipefail
export GARFIELD_INSTALL="\${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
HERE="$(dirname "$(readlink -f "\${BASH_SOURCE[0]}")")"
REPO="$(cd "$HERE/../../.." && pwd)"
PROG="$HERE/export_1800v_heed_root.py"
BINARY="$REPO/test/mwpc-stage-c/build-wiretruth1800-heed-enrichment/enrich_heed_root"
QUEUE="\${MWPC_1800_PBS_QUEUE:-N}"
OUT="\${MWPC_1800_ENRICHED:-}"
SUBSET="\${MWPC_1800_SUBSET:-}"
check_out(){
  [[ -n "$OUT" ]] || { echo "Export MWPC_1800_ENRICHED to a NEW test_runs path" >&2;exit 2; }
}
case "\${1:-}" in
  build)
    command -v root-config >/dev/null || { echo "ROOT environment is not loaded" >&2;exit 2; }
    mkdir -p "$(dirname "$BINARY")"
    c++ -std=c++17 -O2 -Wall -Wextra \
      $(root-config --cflags) "$HERE/enrich_heed_root.cc" \
      $(root-config --libs) -o "$BINARY"
    test -x "$BINARY"
    python3 -m py_compile "$PROG"
    echo "HEED_ROOT_BUILT=$BINARY"
    sha256sum "$BINARY"
    ;;
  prepare)
    check_out
    [[ -n "$SUBSET" ]] || { echo "Set MWPC_1800_SUBSET" >&2;exit 2; }
    python3 "$PROG" prepare --output "$OUT" --subset "$SUBSET"
    ;;
  submit)
    check_out
    [[ -s "$OUT/locked/export.json" && ! -e "$OUT/merged" && ! -e "$OUT/enriched" ]] ||
      { echo "prepare first; do not overwrite existing result" >&2;exit 2; }
    [[ ! -e "$OUT/pbs-logs/submission_ledger.tsv" ]] ||
      { echo "Already submitted: check logs/status instead" >&2;exit 2; }
    [[ "$QUEUE" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] ||
      { echo "Invalid PBS queue" >&2;exit 2; }
    job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18heed \
      -l select=1:ncpus=1:mem=4gb -l walltime=08:00:00 -l io=0.1 \
      -v "MWPC_1800_ENRICHED=$OUT,MWPC_1800_DASH_WORKTREE=$REPO" \
      -o "$OUT/pbs-logs" -e "$OUT/pbs-logs" \
      "$HERE/pbs_1800v_heed_export.sh")
    [[ -n "$job" ]] || { echo "Empty job ID from PBS" >&2;exit 2; }
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$job" \
      > "$OUT/pbs-logs/submission_ledger.tsv"
    echo "HEED_EXPORT_PBS_SUBMITTED=$job"
    ;;
  status)
    check_out
    python3 "$PROG" status --output "$OUT"
    ;;
  *)
    echo "Usage: bash export_1800v_heed.sh {build|prepare|submit|status}" >&2
    exit 2
    ;;
esac

#!/usr/bin/env bash
# Separate read-only 48-shard ROOT merge; NO Garfield++ and NO original edits.
set -Eeo pipefail
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
HERE="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
PY="$HERE/merge_1800v_verified_subset.py"
OUT="${MWPC_1800_SUBSET:-}"
QUEUE="${MWPC_1800_PBS_QUEUE:-N}"
case "${1:-}" in
  prepare)
    [[ -n "$OUT" ]] || { echo "Set MWPC_1800_SUBSET to a NEW test_runs path" >&2; exit 2; }
    python3 "$PY" prepare --output "$OUT"
    ;;
  submit)
    [[ -n "$OUT" && -s "$OUT/locked/subset.json" ]] ||
      { echo "Prepare the frozen independent subset first" >&2; exit 2; }
    [[ ! -e "$OUT/pbs-logs/submission_ledger.tsv" ]] ||
      { echo "Refusing to submit twice; inspect existing PBS job first" >&2; exit 2; }
    [[ "$QUEUE" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] ||
      { echo "Invalid PBS queue" >&2; exit 2; }
    job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18dash \
      -l select=1:ncpus=1:mem=4gb -l walltime=08:00:00 -l io=0.1 \
      -v "MWPC_1800_SUBSET=$OUT" \
      -o "$OUT/pbs-logs" -e "$OUT/pbs-logs" \
      "$HERE/pbs_1800v_verified_subset.sh")
    [[ -n "$job" ]] || { echo "Empty PBS job ID" >&2; exit 2; }
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$job" > "$OUT/pbs-logs/submission_ledger.tsv"
    printf 'SUBSET_ROOT_MERGE_SUBMITTED=%s\n' "$job"
    ;;
  status)
    [[ -n "$OUT" ]] || { echo "Set MWPC_1800_SUBSET" >&2; exit 2; }
    python3 "$PY" status --output "$OUT"
    ;;
  *)
    echo "Usage: bash subset_1800v.sh {prepare|submit|status}" >&2
    exit 2
    ;;
esac

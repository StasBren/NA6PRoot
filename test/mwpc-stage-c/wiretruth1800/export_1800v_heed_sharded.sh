#!/usr/bin/env bash
# Resumable 1800 V 5620 Geant4 muons -> 48 Heed shard tables -> one ROOT.
# No Garfield++ physics is run and no original result is changed.
set -Eeo pipefail
set -u
HERE="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO="$(cd "$HERE/../../.." && pwd)"
S="$HERE/export_1800v_heed_sharded.py"
P="$HERE/export_1800v_heed_root.py"
OUT="${MWPC_1800_ENRICHED:-}"
SUBSET="${MWPC_1800_SUBSET:-}"
QUEUE="${MWPC_1800_PBS_QUEUE:-N}"
[[ "$QUEUE" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] || { echo "invalid queue" >&2;exit 2; }
check_out(){
  [[ -n "$OUT" ]] || { echo "Set MWPC_1800_ENRICHED to a new test_runs folder" >&2;exit 2; }
}
case "${1:-}" in
  prepare)
    check_out
    [[ -n "$SUBSET" ]] || { echo "Set MWPC_1800_SUBSET to verified 5620 ROOT folder" >&2;exit 2; }
    python3 "$P" prepare --output "$OUT" --subset "$SUBSET"
    ;;
  stage-submit)
    check_out
    [[ -s "$OUT/locked/export.json" && ! -e "$OUT/enriched" ]] ||
      { echo "prepare first; no overwrites" >&2;exit 2; }
    [[ ! -e "$OUT/pbs-logs/shard_submission.tsv" ]] ||
      { echo "shard array already submitted; do not duplicate" >&2;exit 2; }
    job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18hs \
      -J '0-47%2' \
      -l select=1:ncpus=1:mem=2gb -l walltime=04:00:00 -l io=0.1 \
      -v "MWPC_1800_ENRICHED=$OUT,MWPC_1800_DASH_WORKTREE=$REPO" \
      -o "$OUT/pbs-logs" -e "$OUT/pbs-logs" \
      "$HERE/pbs_1800v_heed_shard.sh")
    [[ -n "$job" ]] || { echo "empty PBS job ID" >&2;exit 2; }
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$job" \
      > "$OUT/pbs-logs/shard_submission.tsv"
    echo "HEED_SHARDED_PBS_ARRAY=$job"
    ;;
  status)
    check_out
    python3 "$S" status --output "$OUT"
    ;;
  final-submit)
    check_out
    [[ -s "$OUT/locked/export.json" && ! -e "$OUT/enriched" ]] ||
      { echo "invalid/finished export area" >&2;exit 2; }
    [[ ! -e "$OUT/pbs-logs/final_submission.tsv" ]] ||
      { echo "final ROOT job already submitted" >&2;exit 2; }
    python3 "$S" gate --output "$OUT"
    job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18hr \
      -l select=1:ncpus=1:mem=4gb -l walltime=04:00:00 -l io=0.1 \
      -v "MWPC_1800_ENRICHED=$OUT,MWPC_1800_DASH_WORKTREE=$REPO" \
      -o "$OUT/pbs-logs" -e "$OUT/pbs-logs" \
      "$HERE/pbs_1800v_heed_reduce.sh")
    [[ -n "$job" ]] || { echo "empty PBS job ID" >&2;exit 2; }
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$job" \
      > "$OUT/pbs-logs/final_submission.tsv"
    echo "HEED_FINAL_PBS_SUBMITTED=$job"
    ;;
  *)
    echo "Usage: bash export_1800v_heed_sharded.sh {prepare|stage-submit|status|final-submit}" >&2
    exit 2
    ;;
esac

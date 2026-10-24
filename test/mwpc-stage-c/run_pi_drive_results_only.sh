#!/usr/bin/env bash
# Upload NEW Stage-C physics result shards/reports only, not full ROOT/archive.
# source local WSL; destination same connected NA60plus Google Drive folder.
# Non-deleting; copy/checksum; safe to rerun after every analysis batch.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
BASE="$REPO/test_runs/mwpc_stage_c"
REMOTE="${NA60_DRIVE_REMOTE:-na60drive:}"
TARGET="${REMOTE}test_runs/mwpc_stage_c"
MODE="${1:-check}"
LOG="$HOME/na60plus_new_results.log"
[[ "$#" -eq 1 && "$MODE" =~ ^(check|copy|verify)$ ]] || {
  echo "Usage: bash $0 check|copy|verify" >&2;exit 2;
}
command -v rclone >/dev/null || {
  echo "ERROR: Install and configure rclone" >&2;exit 3;
}
[[ "$REMOTE" == *: ]] || {
  echo "ERROR: remote must be na60drive: style" >&2;exit 3;
}
echo "New-results-only cloud destination: $TARGET"
echo "No original ROOT or native crossings CSV is transferred here."
CATEGORIES=(
  "pi_six_geometry_2m_timing_v1/full"
  "pi_weekend_2m_operations_v1"
  "pi_geometric_vs_lut_benchmark_v1"
)
for group in "${CATEGORIES[@]}"; do
  if [[ -d "$BASE/$group" ]];then du -sh "$BASE/$group"
  else echo "Not yet present: $group";fi
done
if [[ "$MODE" == check ]];then
  rclone lsf "$REMOTE" --max-depth 1 >/dev/null
  echo "PASS: research-only Drive storage reachable";exit 0
fi
FILTER=(
  --exclude '*:Zone.Identifier'
  --exclude '*.incomplete'
  --exclude '*.tmp'
  --exclude '*.lock'
  --exclude '*.pid'
  --exclude '*.part'
)
COMMON=(
  --checksum --checkers 8 --transfers 4
  --drive-chunk-size 32M --tpslimit 7
  --low-level-retries 10 --retries 5
  --log-level INFO --log-file "$LOG"
)
for group in "${CATEGORIES[@]}"; do
  [[ -d "$BASE/$group" ]] || continue
  origin="$BASE/$group"
  target="$TARGET/$group"
  if [[ "$MODE" == copy ]];then
    echo "COPY: $group (skip files younger than 2 min to avoid writers)"
    rclone copy "$origin" "$target" "${FILTER[@]}" "${COMMON[@]}" \
      --min-age 2m --stats 60s
  else
    echo "VERIFY: $group (run after production and copy complete)"
    rclone check "$origin" "$target" "${FILTER[@]}" \
      --checksum --checkers 8 --tpslimit 7 \
      --log-level INFO --log-file "$LOG"
  fi
done
echo "PASS: $MODE finished for research results only."
echo "Log: $LOG"

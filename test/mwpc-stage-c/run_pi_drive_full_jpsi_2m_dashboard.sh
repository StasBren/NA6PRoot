#!/usr/bin/env bash
# Curated read-only 2M Jpsi dashboard publication; no 240 gzip folds,
# original ROOT, crossing CSV or 6M calculation is uploaded or changed.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
SRC="$REPO/test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1"
DEST="na60drive:test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1"
MODE="$1"
if [[ "$#" -ne 1 || ! "$MODE" =~ ^(check|copy|verify)$ ]];then
  echo "Usage: bash $0 check|copy|verify" >&2;exit 2
fi
command -v rclone >/dev/null || {
  echo "ERROR: rclone unavailable" >&2;exit 127;
}
python3 "$HERE/pi_full_jpsi_2m_dashboard.py" verify
echo "Verified FULL two-million-parent Jpsi dashboard."
echo "Drive destination: $DEST"
echo "Original native ROOT/CSV and frozen Garfield files remain in WSL."
if [[ "$MODE" == check ]];then
  du -sh "$SRC/pages" "$SRC/station_summaries" "$SRC/summary.csv"
  rclone lsf na60drive: --max-depth 1
  exit 0
fi
# Only publish interactive self-contained HTML and small summary/manifests.
# All large raw merged station JSONs and source folds remain local.
if [[ "$MODE" == copy ]];then
  rclone copy "$SRC" "$DEST" \
    --include '/pages/**' \
    --include '/station_summaries/**' \
    --include '/complete_dashboard_manifest.json' \
    --include '/summary.csv' \
    --exclude '**' \
    --checksum --transfers 2 --checkers 4 \
    --drive-chunk-size 16M --tpslimit 3 --retries 5 \
    --low-level-retries 8 --stats 45s -P
  echo "PASS: complete Jpsi 2M dashboard copied; WSL originals kept."
else
  rclone check "$SRC" "$DEST" \
    --include '/pages/**' \
    --include '/station_summaries/**' \
    --include '/complete_dashboard_manifest.json' \
    --include '/summary.csv' \
    --exclude '**' \
    --checksum --checkers 4 --tpslimit 3 --retries 4
  echo "PASS: cloud files check out against local published dashboard."
fi

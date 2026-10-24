#!/usr/bin/env bash
# Safe Stage-C archive to the user's connected Google Drive folder.
# Set rclone OAuth remote na60drive root_folder_id to project folder ID.
# NEVER delete local data; copy/check only, no sync/move.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
MODE="${1:-check}"
REMOTE="${NA60_DRIVE_REMOTE:-na60drive:}"
SRC="$REPO/test_runs/mwpc_stage_c"
DEST="${REMOTE}test_runs/mwpc_stage_c"
LOG="$HOME/na60plus_drive_archive.log"
[[ "$#" -eq 1 && "$MODE" =~ ^(check|copy|verify)$ ]] || {
  echo "Usage: bash $0 check|copy|verify" >&2;exit 2;
}
if ! command -v rclone >/dev/null;then
  echo "ERROR: rclone missing. In WSL: sudo apt update && sudo apt install rclone" >&2
  exit 3
fi
if [[ ! -d "$SRC" ]];then echo "ERROR: no results at $SRC" >&2;exit 4;fi
if [[ "$REMOTE" != *: ]];then
  echo "ERROR: NA60_DRIVE_REMOTE must be an rclone root such as na60drive:" >&2
  exit 4
fi
echo "STORAGE: local experiment files $SRC"
du -sh "$SRC"
echo "Google Drive archive destination: $DEST"
echo "Expected project Drive folder: https://drive.google.com/drive/folders/1Svbp3E230YmVk53D0M1AHwTGo6cBxy38"
echo "Set rclone remote root_folder_id=1Svbp3E230YmVk53D0M1AHwTGo6cBxy38"
echo "Credentials stay in local rclone OAuth configuration."
if ! rclone lsf "$REMOTE" --max-depth 1 >/dev/null;then
  echo "ERROR: Drive OAuth remote unavailable. Run rclone config."
  exit 5
fi
echo "REMOTE ACCESS PASS"
rclone about "$REMOTE" 2>/dev/null || echo "Quota information unavailable."
if [[ "$MODE" == check ]];then
  echo "PASS: no files uploaded or deleted. Next use copy mode."
  exit 0
fi
SKIP=(
  --exclude '*.incomplete'
  --exclude '*.tmp'
  --exclude '*.lock'
  --exclude '*.pid'
  --exclude '*.part'
)
COMMON=(
  --checksum
  --checkers 8
  --transfers 2
  --drive-chunk-size 32M
  --tpslimit 7
  --retries 6
  --low-level-retries 10
)
if [[ "$MODE" == copy ]];then
  echo "COPY: adds/updates cloud files, NO deletion locally or remotely."
  rclone copy "$SRC" "$DEST" "${SKIP[@]}" "${COMMON[@]}" \
      --create-empty-src-dirs --stats 60s --stats-one-line \
      --log-level INFO --log-file "$LOG"
  echo "COPY PASS: completed results archived. Log: $LOG"
  echo "Repeat copy after more folds; run verify at campaign completion."
  exit 0
fi
echo "VERIFY: compare complete source file checksums with stored cloud bytes."
rclone check "$SRC" "$DEST" "${SKIP[@]}" \
    --checksum --checkers 8 --tpslimit 7 \
    --stats 60s --log-level INFO --log-file "$LOG"
echo "VERIFY PASS: cloud completed dataset matches local WSL source."

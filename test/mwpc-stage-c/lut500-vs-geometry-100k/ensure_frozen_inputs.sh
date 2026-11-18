#!/usr/bin/env bash
# Retrieve exactly the frozen Fast4/Fast5 LUT500 and J/psi donor inputs.
# Reuses the historical verified Drive archive and refuses hash mismatches.
set -euo pipefail

kind="${1:-}"
case "$kind" in
  lut)
    filename="NA60_LUT500.zip"
    varname="LUT_ZIP"
    expected="6058c83d2d9a49fcab0867ba635f1318d40b6f8aad9d3e058ed6646561e40ee5"
    ;;
  calibration)
    filename="model_calibration_jpsi_v1.npz"
    varname="CALIBRATION_NPZ"
    expected="7eb519f3c62e495b750eff744751865c9563e2d5b0078b05f1b8d1c0f6fdc74e"
    ;;
  *)
    echo "Usage: bash ensure_frozen_inputs.sh lut|calibration" >&2
    exit 2
    ;;
esac

verify() {
  local file="$1" actual
  [[ -f "$file" ]] || { echo "Missing frozen input: $file" >&2; return 1; }
  actual="$(sha256sum "$file" | cut -d' ' -f1)"
  if [[ "$actual" != "$expected" ]]; then
    echo "SHA256 MISMATCH for $file" >&2
    echo "Expected: $expected" >&2
    echo "Actual:   $actual" >&2
    return 1
  fi
  echo "[SHA256 OK] $file" >&2
}

if [[ -n "${!varname:-}" ]]; then
  path="${!varname}"
  verify "$path"
  printf '%s\n' "$path"
  exit 0
fi

cache="${NA60_FROZEN_INPUT_DIR:-$HOME/na6p/studies/lut500_vs_geometry_100k/frozen_inputs}"
mkdir -p "$cache"
path="$cache/$filename"
if [[ -f "$path" ]]; then
  verify "$path"
  printf '%s\n' "$path"
  exit 0
fi

remote="${NA60_RCLONE_REMOTE:-na60drive}"
folder_id="1hrU0Puv1I8ADC1udhNsLdBXZoQbuok5B"
if ! command -v rclone >/dev/null 2>&1; then
  echo "rclone is not installed. Set $varname to the exact frozen input path, or install/configure rclone." >&2
  exit 2
fi
if ! rclone listremotes | grep -Fxq "${remote}:"; then
  echo "rclone remote ${remote}: not configured. Check: rclone listremotes" >&2
  echo "Alternatively set $varname=/absolute/path/$filename" >&2
  exit 2
fi

tmp="$(mktemp "$cache/.$filename.download.XXXXXX")"
trap 'rm -f "$tmp"' EXIT
echo "[rclone] Fetching frozen $filename from verified Fast4 validation Drive folder" >&2
rclone copyto "${remote}:$filename" "$tmp" \
  --drive-root-folder-id "$folder_id" \
  --transfers 1 --retries 5 --low-level-retries 5 --timeout 3m -P >&2
verify "$tmp"
mv "$tmp" "$path"
trap - EXIT
printf '%s\n' "$path"

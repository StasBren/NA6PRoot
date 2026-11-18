#!/usr/bin/env bash
# NA60+/DiCE: reproducible LUT500 vs midline, from verified frozen WSL inputs.
# Does not reset, overwrite, or re-run any historical Geant4/Garfield campaign.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
case "${1:-}" in
  test)
    "$PYTHON_BIN" -m unittest discover -s "$HERE" -p 'test_*.py' -v
    ;;
  export)
    "$PYTHON_BIN" "$REPO/test/mwpc-stage-c/pi_2m_export_completed_chunks.py" \
      --max-chunks-per-source 40 --root "${ROOT_BIN:-$HOME/na6p/install/root/bin/root}"
    ;;
  select)
    "$PYTHON_BIN" "$HERE/run_comparison.py" select \
      --repo "$REPO" --chunks 40 --n 100000
    ;;
  pilot|full)
    if [[ "$1" == pilot ]]; then n=10000; else n=100000; fi
    find_unique() {
      local filename="$1" varname="$2" path="" found=0 hit
      if [[ -n "${!varname:-}" ]]; then
        path="${!varname}"
        [[ -f "$path" ]] || { printf 'Missing %s: %s\n' "$varname" "$path" >&2; return 1; }
      else
        while IFS= read -r -d '' hit; do
          ((found+=1))
          path="$hit"
        done < <(find "$HOME/na6p" -type f -name "$filename" -print0 2>/dev/null)
        if (( found != 1 )); then
          printf 'Expected exactly one %s under ~/na6p, found %s.\n' "$filename" "$found" >&2
          printf 'Set %s=/absolute/path/to/%s and rerun.\n' "$varname" "$filename" >&2
          return 1
        fi
      fi
      printf '%s\n' "$path"
    }
    lut="$(find_unique 'NA60_LUT500.zip' LUT_ZIP)"
    calibration="$(find_unique 'model_calibration_jpsi_v1.npz' CALIBRATION_NPZ)"
    "$PYTHON_BIN" "$HERE/run_comparison.py" compare \
       --repo "$REPO" --selected-n 100000 --limit "$n" \
       --lut "$lut" --calibration "$calibration" --weights donor
    ;;
  *)
    cat <<'USAGE'
Usage: bash test/mwpc-stage-c/lut500-vs-geometry-100k/run.sh COMMAND
  test    -- seven synthetic tests (no physics data)
  export  -- verified/restartable existing 2M Geant4 ROOT-to-CSV export
  select  -- nine hash-ranked 100k crossing samples (3 channels x 3 stations)
  pilot   -- run first 10k in each sample with frozen LUT500 + J/psi donors
  full    -- run all 100k in each sample with frozen LUT500 + J/psi donors
Optional env vars: ROOT_BIN, LUT_ZIP, CALIBRATION_NPZ, PYTHON_BIN.
Requires the pre-existing 2M Geant4 parent chunks and Heed calibration files in WSL.
USAGE
    exit 2
    ;;
esac
#!/usr/bin/env bash
# VALIDATED symmetric Garfield 3+3mm + Heed on ACTUAL first-chunk
# Geant4 native gas-hit segments, all Jpsi/Phi/Omega; no ROOT or transport.
# Existing legacy 2+4mm outputs stay untouched.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/pi_centered_electron_sampler.py" \
 "$HERE/pi_native_centered_timing_pilot.py" \
 "$HERE/test_pi_native_centered_timing_pilot.py"
python3 "$HERE/test_pi_native_centered_timing_pilot.py"
echo "Reading pre-existing true Geant4 first-chunk native gas hits, Heed library,"
echo "and the finished MANIFEST-VERIFIED symmetric 3+3mm Garfield map."
echo "This may take several minutes; no new Garfield or Geant4 transport."
python3 "$HERE/pi_native_centered_timing_pilot.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_timing_pilot.html"
echo "PASS: Native Geant4/Heed/NEW centered Garfield timing pilot:"
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
  echo "Windows browser path: $(wslpath -w "$OUT")"
  if command -v explorer.exe >/dev/null 2>&1; then
    explorer.exe "$(wslpath -w "$OUT")" || true
  fi
fi

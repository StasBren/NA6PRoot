#!/usr/bin/env bash
# PI: true 2-anode Geant4 event tracks in local u-v plane and separate
# collected-primary-electron arrival-time histograms, no waveforms.
# Reuse completed 3+3mm viewer electron arrival times and native Geant4
# gas-hit TGeo-local entry/exit coordinates, absolutely NO transports.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/build_pi_two_anode_electron_explorer.py" \
 "$HERE/test_pi_two_anode_electron_explorer.py"
python3 "$HERE/test_pi_two_anode_electron_explorer.py"
python3 "$HERE/build_pi_two_anode_electron_explorer.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_two_anode_electron_explorer.html"
python3 - "$REPO" <<'PY'
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import sys
d=Path(sys.argv[1])/"test_runs/mwpc_stage_c/pi_fixed_wire_maps"
archive=d/"pi_two_anode_electron_PI_review.zip"
with ZipFile(archive,"w",ZIP_DEFLATED,compresslevel=6) as z:
    for name in ("pi_two_anode_electron_explorer.html",
                 "pi_two_anode_electron_explorer_manifest.json"):
        p=d/name
        if not p.is_file():raise FileNotFoundError(p)
        z.write(p,arcname=name)
print("PI electron-only event review ZIP:",archive)
PY
echo "PASS: all 4,632 curated TWO-ANODE Geant4 gas events, electron-count-only time histograms and local u-v tracks."
echo "  $OUT"
if command -v wslpath >/dev/null 2>&1; then
 echo "Windows file path: $(wslpath -w "$OUT")"
 if command -v explorer.exe >/dev/null 2>&1; then
  explorer.exe "$(wslpath -w "$OUT")" || true
 fi
fi

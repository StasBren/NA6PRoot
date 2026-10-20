#!/usr/bin/env bash
# Native Geant4 true muon angle -> independent Heed primary electrons ->
# NEW microscopic Garfield true drift-line points -> normalized 1D
# P(u(t)), P(v(t)), P(w(t)) for electron cloud STILL DRIFTING at t.
#
# No re-run of the frozen 600k Geant4 campaign.
# Unlike endpoint-map visualisations this runs a SMALL NEW Garfield
# microscopic trajectory PILOT (~selected thousands of electrons).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/pi_centered_native_drift_births.py" \
 "$HERE/pi_centered_true_1d_cloud.py" \
 "$HERE/test_pi_centered_true_1d_cloud.py"
python3 "$HERE/test_pi_centered_true_1d_cloud.py"
python3 "$HERE/pi_centered_native_drift_births.py" "$@"

echo "Compiling existing Garfield Stage-C transport executable (new trace-mode only)."
cmake -S "$HERE" -B "$HERE/build" \
 -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$HERE/build" -j4 \
 --target mwpc_stage_c_transport_box_map

python3 "$HERE/refresh_pi_dashboard_navigation.py"
OUT="$REPO/test_runs/mwpc_stage_c/pi_centered_native_drift_trajectories"
TRACE="$OUT/centered_native_microscopic_drift_points.csv"
echo "New Garfield microscopic drifts from TRUE native-angle Heed origins..."
"$HERE/build/mwpc_stage_c_transport_box_map" \
 --gap-minus-mm 3 --gap-plus-mm 3 \
 --wire-pitch-mm 4 --hv 1800 --b-tesla 0 \
 --wire-diam-um 30 \
 --trace-starts "$OUT/native_heed_starts.csv" \
 --trace-output "$TRACE" \
 --trace-collision-steps 20

python3 "$HERE/pi_centered_true_1d_cloud.py"
VIEW="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_true_1d_cloud.html"
python3 - "$REPO" <<'PY'
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import sys
root=Path(sys.argv[1])/"test_runs/mwpc_stage_c"
d=root/"pi_centered_native_drift_trajectories"
v=root/"pi_fixed_wire_maps"
archive=d/"pi_centered_true_1d_PI_review.zip"
with ZipFile(archive,"w",ZIP_DEFLATED,compresslevel=6) as z:
    for path in (
      v/"pi_centered_true_1d_cloud.html",
      v/"pi_centered_true_1d_cloud_manifest.json",
      d/"native_heed_trace_manifest.json",
      d/"native_heed_start_metadata.csv",
      d/"centered_native_microscopic_drift_points.csv",
    ):
        if not path.is_file():raise FileNotFoundError(path)
        z.write(path,arcname=path.name)
print("PI true microscopic 1D review ZIP:",archive)
PY
echo "PASS: actual microscopic intermediate-coordinate pilot vs native muon angles"
echo "  $VIEW"
if command -v wslpath >/dev/null 2>&1; then
 echo "Windows browser: $(wslpath -w "$VIEW")"
 if command -v explorer.exe >/dev/null 2>&1; then
  explorer.exe "$(wslpath -w "$VIEW")" || true
 fi
fi

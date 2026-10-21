#!/usr/bin/env bash
# Review already completed true native first-chunk six-geometry output.
# READ ONLY on original 496k-group full-fold JSON and the frozen Garfield
# and native Geant4/Heed inputs. No ROOT, Garfield, or Geant4 executables.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"

MODE="pages"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 ]] || [[ "$MODE" != "pages" && "$MODE" != "bootstrap" && "$MODE" != "all" ]]; then
  echo "Usage: bash $0 pages|bootstrap|all" >&2
  exit 2
fi
python3 -m py_compile \
 "$HERE/pi_six_geometry_native_dashboard_pages.py" \
 "$HERE/test_pi_six_geometry_native_dashboard_pages.py" \
 "$HERE/pi_six_geometry_event_bootstrap.py" \
 "$HERE/pi_six_geometry_parent_bootstrap_html.py" \
 "$HERE/test_pi_six_geometry_event_bootstrap.py"
if [[ "$MODE" == "pages" || "$MODE" == "all" ]]; then
  echo "SMOKE: streaming JSON and 18 offline HTML sharded dashboards."
  python3 "$HERE/test_pi_six_geometry_native_dashboard_pages.py"
  echo "READ ONLY: split completed native first chunk into browser-sized pages."
  python3 "$HERE/pi_six_geometry_native_dashboard_pages.py"
fi
if [[ "$MODE" == "bootstrap" || "$MODE" == "all" ]]; then
  echo "SMOKE: parent-event bootstrap and exact full-fold conservation."
  python3 "$HERE/test_pi_six_geometry_event_bootstrap.py"
  REFERENCE="$REPO/test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/pages/pi_six_geometry_station_reference.json"
  if [[ ! -f "$REFERENCE" ]]; then
    echo "Need station reference from pages mode first; building it..."
    python3 "$HERE/pi_six_geometry_native_dashboard_pages.py"
  fi
  echo "REFOLD first 3000 generated parents/source with original shared-Heed seeds."
  echo "No new Garfield or Geant4, but Python event refolding can take time."
  python3 "$HERE/pi_six_geometry_event_bootstrap.py"
  python3 "$HERE/pi_six_geometry_parent_bootstrap_html.py"
fi
echo "PASS: six-geometry first-chunk postprocessing $MODE finished."

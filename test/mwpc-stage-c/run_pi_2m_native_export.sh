#!/usr/bin/env bash
# After first 100k/channel Geant4 benchmark: verify/read frozen ROOT.
# preflight does not run ROOT, Garfield or Geant4. benchmark launches
# ROOT ONLY to export native 298-detector gas hit CSV in six original chunks.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="preflight"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 || ! "$MODE" =~ ^(preflight|benchmark)$ ]]; then
  echo "Usage: bash $0 preflight|benchmark" >&2
  exit 2
fi
python3 -m py_compile \
  "$HERE/pi_2m_export_completed_chunks.py" \
  "$HERE/test_pi_2m_export_completed_chunks.py"
python3 "$HERE/test_pi_2m_export_completed_chunks.py"
echo "Checking frozen 100k/Geant4 benchmark checkpoint ROOT SHA-256..."
if [[ "$MODE" == "preflight" ]]; then
  python3 "$HERE/pi_2m_export_completed_chunks.py" \
    --max-chunks-per-source 2 --dry-run
  echo "PASS: six native Geant4 chunk export preflight (no ROOT process)."
  exit 0
fi
NA6P_STAGEC_REQUIRE_ROOT_HADD=1
source "$HERE/activate_pi_stagec_env.sh"
python3 "$HERE/pi_2m_export_completed_chunks.py" \
  --root "$(command -v root)" --max-chunks-per-source 2
echo "PASS: six verified native CSV, each source+chunk separately identified."
echo "  $REPO/test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1"

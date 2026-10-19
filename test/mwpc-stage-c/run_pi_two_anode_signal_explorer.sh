#!/usr/bin/env bash
# Curated TRUE native-Geant4 two physical anode signal gallery:
# realistic first-chunk MCKine muon gas-hit entrance/exit geometry,
# EXACT same Heed sample + fresh verified centered Garfield 3+3mm response,
# INDEPENDENT uncalibrated anode-only electronics for two wire channels.
# A+B is a hypothetical arithmetic comparison, no real summed ASIC channel.
#
# No Geant4 transport, ROOT re-export, Garfield drift generation or
# cathode-strip simulation. Uses the already completed 74,322-hit pilot CSV.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile \
 "$HERE/pi_native_two_anode_signal_explorer.py" \
 "$HERE/test_pi_native_two_anode_signal_explorer.py" \
 "$HERE/pi_anode_digitizer.py"
python3 "$HERE/test_pi_anode_digitizer.py"
python3 "$HERE/test_pi_native_two_anode_signal_explorer.py"
echo "Replaying representative N_all=2 native Geant4 muon hits by source/chamber/angle."
echo "NO new Geant4, ROOT or microscopic Garfield simulation."
python3 "$HERE/pi_native_two_anode_signal_explorer.py" "$@"
OUTPUT="$REPO/test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_two_anode_signal_explorer.html"
echo "PASS: two separate physical anode readouts + illustrative sum, pulse times, XY/angle/event selector"
echo "OPEN: $OUTPUT"
if command -v wslpath >/dev/null 2>&1; then
 echo "Windows file: $(wslpath -w "$OUTPUT")"
 if command -v explorer.exe >/dev/null 2>&1; then
  explorer.exe "$(wslpath -w "$OUTPUT")" || true
 fi
fi

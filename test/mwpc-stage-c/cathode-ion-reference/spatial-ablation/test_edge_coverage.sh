#!/usr/bin/env bash
# SYNTHETIC regression test: ions born 0.001 um above the 15 um anode must be inside the generated radial grid,
# no grid start point may be on/inside the anode, and the strict no-extrapolation check in s_synthesize.py must pass
# (G closure is then computed). Also checks that the strict check still REJECTS a grid that misses those ions.
#   bash test_edge_coverage.sh LUT.npz NUMW_DIR [WORKDIR]
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LUT="${1:?spectral LUT .npz}"; NUMW="${2:?dir with p3_replay_full_spectral_muon.py}"
W="${3:-$(mktemp -d)}"; mkdir -p "$W"; [[ -z "$(ls -A "$W")" ]] || { echo "workdir $W not empty" >&2; exit 2; }
python3 "$HERE/synthetic_paths.py" muon "$W/muon" --seed 3 --edge-ions 6 > /dev/null
python3 "$HERE/s_grid_spec.py" "$W/muon" --nr 8 | tee "$W/spec.txt"
grep -q '^PASS radial coverage' "$W/spec.txt"
RU=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$W/spec.txt"); WI=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$W/spec.txt")
DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$W/spec.txt")
python3 - "$RU" <<'PY'
import sys
r = [float(x) for x in sys.argv[1].split(",")]
assert r[0] > 15.0 and r[0] <= 15.001 + 1e-12, f"innermost node {r[0]!r} not in (15, 15.001]"
print(f"PASS innermost node {r[0]!r} um: above the wire surface, below the 15.001 um edge ions")
PY
python3 "$HERE/synthetic_paths.py" grid "$W/grid" --wires="$WI" --r-um="$RU" --ntheta 16 > /dev/null
python3 "$HERE/s_green_table.py" "$LUT" "$W/grid" "$W/green" --replay-dir "$NUMW" --d-max-mm "$DM" --ntau 120 | tail -1
(cd "$NUMW" && python3 p3_replay_full_spectral_muon.py "$LUT" "$W/muon" "$W/S0" --half-strips 5 > "$W/s0.log")
python3 "$HERE/s_synthesize.py" "$W/green" "$W/muon" "$W/S0/waveforms.csv" "$W/out" --variants G | tee "$W/synth.log"
grep -q "ions outside the grid: 0" "$W/synth.log"
# negative control: the old lower bound (15.05 um) must still be rejected by the strict check
OLD=$(python3 -c "import sys; r=sys.argv[1].split(','); r[0]='15.05'; print(','.join(r))" "$RU")
python3 "$HERE/synthetic_paths.py" grid "$W/grid_old" --wires="$WI" --r-um="$OLD" --ntheta 16 > /dev/null
python3 "$HERE/s_green_table.py" "$LUT" "$W/grid_old" "$W/green_old" --replay-dir "$NUMW" --d-max-mm "$DM" --ntau 120 | tail -1
if python3 "$HERE/s_synthesize.py" "$W/green_old" "$W/muon" "$W/S0/waveforms.csv" "$W/out_old" --variants G > "$W/synth_old.log" 2>&1; then
  echo "FAIL: strict no-extrapolation check accepted a grid that misses the edge ions"; exit 1
fi
grep -q "outside the radial grid" "$W/synth_old.log"
echo "PASS edge-coverage regression (SYNTHETIC): edge ions covered, old grid still rejected; work dir $W"

#!/usr/bin/env bash
# Small Gmsh geometry gate for NA60+ P3 stereo-strip FEM.
# Isolated 3-wire subset: geometry and Jacobian feasibility ONLY.
# No Elmer, Garfield, microscopic avalanche, or existing data touched.
# Never overwrites an existing output directory.
set -Eeuo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="$REPO/test_runs/mwpc_stage_c/p3_2100_fem_3wire_fastcurve_v4"
if [[ "$#" -ge 1 ]]; then OUT="$1"; fi
[[ "$OUT" = /* ]] || { echo "Provide an absolute new output directory" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "OUTPUT EXISTS: $OUT (refusing overwrite)" >&2; exit 2; }
for cmd in python3 gmsh; do
  command -v "$cmd" >/dev/null 2>&1 ||
    { echo "Missing command: $cmd" >&2; exit 2; }
done
[[ -s "$HERE/p3_2100_stereo_width_along_w_1p7mm.json" &&
   -s "$HERE/generate_weighting_model.py" ]] ||
  { echo "Missing source FEM generator/config" >&2; exit 2; }
if pgrep -x gmsh >/dev/null || pgrep -x ElmerGrid >/dev/null ||
   pgrep -x ElmerSolver >/dev/null; then
  echo "FEM process already running: do not compete for 7.4-GiB RAM" >&2
  exit 2
fi
mkdir -p "$OUT"
cp -- "$HERE/p3_2100_stereo_width_along_w_1p7mm.json" \
  "$OUT/input_config.json"
python3 - "$OUT/input_config.json" <<'PY'
import json
import pathlib
import sys
p = pathlib.Path(sys.argv[1])
c = json.loads(p.read_text())
assert c["geometry"]["half_wires"] == 4
assert c["geometry"]["wire_diameter_um"] == 30
assert c["mesh"]["wire_size_mm"] == 0.012
g = c["geometry"]
g["half_wires"] = 1  # exactly three wires: -3, 0, +3 mm
g["u_half_span_mm"] = 5.0
g["w_half_span_mm"] = 4.0
c["mesh"]["curved_wire_surfaces"] = True
p.write_text(json.dumps(c, indent=2) + "\n")
print("PHYSICS TEST: three 30-um wires at 3-mm pitch, +/-3-mm cathodes")
print("Finite box now 10 x 6 x 8 mm; NOT the nine-wire detector field.")
PY
python3 "$HERE/generate_weighting_model.py" \
  --config "$OUT/input_config.json" \
  --output "$OUT/map" > "$OUT/generate.log" 2>&1
python3 - "$OUT/map/chamber.geo" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
s = p.read_text()
old = "Mesh.HighOrderOptimize = 1;"
assert s.count(old) == 1, "Unexpected Gmsh high-order option"
assert "Mesh.SecondOrderLinear = 0;" in s
s = s.replace(old, "Mesh.HighOrderOptimize = 4;")
p.write_text(s)
PY
cat > "$OUT/SCOPE.txt" <<'EOF'
Exploratory Gmsh-only three-wire curved P2 geometry test.
Same 30um anode wire diameter, 3mm wire pitch, +/-3mm cathodes,
12um near-wire mesh size and U/V strip boundaries; reduced box and
3 wires instead of 9. Mesh.HighOrderOptimize=4 (fast curving).
This is NOT a full nine-wire FEM convergence test or a physical
field calculation. This directory must not be used as a validated
Garfield/Elmer weighting map.
No muon or avalanche processing. Existing archives unchanged.
EOF
echo "Gmsh-only fast-curving mesh test: $OUT"
set +e
(
  cd "$OUT/map" &&
  gmsh chamber.geo -3 -order 2 -format msh2 \
    -o chamber.msh > gmsh.log 2>&1
)
status=$?
set -e
echo "Gmsh exit code: $status"
tail -n 22 "$OUT/map/gmsh.log"
if [[ $status -ne 0 ]]; then
  echo "FAIL: gmsh exited nonzero; preserve logs and investigate" >&2
  exit 3
fi
if grep -Eq '^Error[[:space:]]*:' "$OUT/map/gmsh.log"; then
  echo "FAIL: Gmsh log contains error(s), even though mesh was written" >&2
  exit 3
fi
[[ -s "$OUT/map/chamber.msh" ]] ||
  { echo "FAIL: No resulting .msh file" >&2; exit 3; }
if ! python3 -c 'import gmsh' >/dev/null 2>&1; then
  echo "PENDING: Gmsh Python API unavailable in this Python environment."
  echo "Mesh generated but NOT accepted. Preserve output; do not run Elmer."
  echo "We need independent curved-Tet10 Jacobian validation."
  exit 4
fi
echo "Auditing curved-element Jacobians..."
python3 "$HERE/audit_p3_curved_tet10_jacobians.py" \
  "$OUT/map/chamber.msh" | tee "$OUT/jacobian_audit.log"
echo "GMSH-ONLY ACCEPTANCE: no reported Gmsh errors or inverted Tet10"
echo "Not a field/Elmer/Garfield validation."

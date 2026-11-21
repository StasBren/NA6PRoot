#!/usr/bin/env bash
# NA60+ P3 geometry-only three-wire finer straight-sided P2 mesh pilot.
# New directory only; no Elmer, Garfield, microscopic data or prior result changes.
set -Eeuo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="$REPO/test_runs/mwpc_stage_c/p3_2100_fem_3wire_fine_straight_p2_v5"
[[ $# -le 1 ]] || { echo "Usage: $0 [absolute_new_output_dir]" >&2; exit 2; }
if [[ $# -eq 1 ]]; then OUT="$1"; fi
[[ "$OUT" = /* && ! -e "$OUT" ]] ||
  { echo "STOP: output must be new absolute path: $OUT" >&2; exit 2; }

for cmd in python3 gmsh; do command -v "$cmd" >/dev/null ||
  { echo "Missing $cmd" >&2; exit 2; }; done

if pgrep -x gmsh >/dev/null ||
   pgrep -x ElmerGrid >/dev/null ||
   pgrep -x ElmerSolver >/dev/null; then
  echo "STOP: FEM job already running; do not compete for RAM" >&2
  exit 2
fi

ENV_PY="$HOME/na6p/.venv-gmsh4121/bin/python"
if [[ -n "${MWPC_GMSH_PYTHON:-}" ]]; then ENV_PY="$MWPC_GMSH_PYTHON"; fi
[[ -x "$ENV_PY" ]] ||
  { echo "STOP: missing Gmsh Python API environment: $ENV_PY" >&2; exit 2; }
"$ENV_PY" -c 'import gmsh; gmsh.initialize(); gmsh.finalize()' ||
  { echo "STOP: Gmsh Python API is unavailable" >&2; exit 2; }

mkdir -p "$OUT"
cp "$HERE/p3_2100_stereo_width_along_w_1p7mm.json" \
  "$OUT/input_config.json"

python3 - "$OUT/input_config.json" <<'PY'
import json
from pathlib import Path
import sys
p = Path(sys.argv[1])
cfg = json.loads(p.read_text())
g, m = cfg["geometry"], cfg["mesh"]
assert g["half_wires"] == 4
assert g["wire_diameter_um"] == 30
assert g["wire_pitch_mm"] == 3
assert m["wire_size_mm"] == 0.012
g["half_wires"] = 1
g["u_half_span_mm"] = 5.0
g["w_half_span_mm"] = 4.0
m["wire_size_mm"] = 0.004
m["wire_curvature_points"] = 32
m["curved_wire_surfaces"] = False
p.write_text(json.dumps(cfg, indent=2) + "\n")
print("Three-wire reduced box: 10 x 6 x 8 mm")
print("Straight-sided P2. Target wire mesh 4 um, curvature points 32")
print("No full-detector field-accuracy claims.")
PY

python3 "$HERE/generate_weighting_model.py" \
  --config "$OUT/input_config.json" \
  --output "$OUT/map" > "$OUT/generate.log" 2>&1

grep -Fq 'Mesh.SecondOrderLinear = 1;' "$OUT/map/chamber.geo" ||
  { echo "STOP: not straight-sided P2" >&2; exit 2; }

cat > "$OUT/SCOPE.txt" <<'EOF'
Exploratory three-wire straight-sided quadratic tetrahedral Gmsh-only trial.
Anode cylinder 30um diameter, pitch 3mm, cathodes +/-3mm; finite reduced
10 x 6 x 8mm box. Local target 4um, curvature 32 points per full circle.
Unlike the prior failed curved-P2 cases, mid-edge nodes remain linear.
No Elmer, Garfield, microscopic avalanches, frozen muon data, or electronics.
This reduced geometry is NOT a physically validated full-chamber U/V map.
EOF

echo "Gmsh-only geometry test. New output: $OUT"
set +e
(
  cd "$OUT/map"
  gmsh chamber.geo -3 -order 2 -format msh2 \
    -o chamber.msh > gmsh.log 2>&1
)
status=$?
set -e
echo "Gmsh status: $status"
tail -n 20 "$OUT/map/gmsh.log"
if [[ $status -ne 0 ]]; then
  echo "FAIL: Gmsh nonzero exit; preserve all outputs" >&2
  exit 3
fi
if grep -Eq '^Error[[:space:]]*:' "$OUT/map/gmsh.log"; then
  echo "FAIL: Gmsh errors in log; mesh rejected" >&2
  exit 3
fi
[[ -s "$OUT/map/chamber.msh" ]] ||
  { echo "FAIL: no mesh" >&2; exit 3; }

echo "Independent full Tet10 minDetJac / minSJ scan:"
set -o pipefail
"$ENV_PY" "$HERE/audit_p3_curved_tet10_jacobians.py" \
  "$OUT/map/chamber.msh" | tee "$OUT/jacobian_audit.log"

echo "PASS: small straight-sided P2 geometry-only Jacobian gate."
echo "NOT a Garfield/Elmer or nine-wire geometry certification."

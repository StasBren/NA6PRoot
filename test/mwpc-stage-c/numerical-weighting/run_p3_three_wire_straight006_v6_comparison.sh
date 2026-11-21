#!/usr/bin/env bash
# NA60+ P3 geometry-only three-wire 6um comparison, FEM solutions and Garfield field survey.
# New directory only; no microscopic data or prior result changes.
set -Eeuo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="$REPO/test_runs/mwpc_stage_c/p3_2100_fem_3wire_straight006_v6"
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

ELMER="$HOME/na6p/install/elmer/bin"
BUILD="$REPO/test/mwpc-stage-c/build-fem-strip-probe"
BIN="$BUILD/mwpc_stage_c_fem_stereo_near_wire_probe"
REF="$REPO/test_runs/mwpc_stage_c/p3_2100_fem_3wire_v5_elmer_import_v1/map/v5_near_wire_probe_retry1.csv"
[[ -s "$REF" && -x "$ELMER/ElmerGrid" && -x "$ELMER/ElmerSolver" && -x "$BIN" ]] ||
  { echo "STOP: missing V5 probe or Elmer/Garfield executables" >&2; exit 2; }
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
if ldd "$BIN" | grep -q 'not found'; then
  echo "STOP: Garfield probe has unresolved libraries" >&2
  exit 2
fi
echo "Disk space and RAM:"
df -h "$REPO"
free -h
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
m["wire_size_mm"] = 0.006
m["wire_curvature_points"] = 32
m["curved_wire_surfaces"] = False
p.write_text(json.dumps(cfg, indent=2) + "\n")
print("Three-wire reduced box: 10 x 6 x 8 mm")
print("Straight-sided P2. Target wire mesh 6 um, curvature points 32")
print("No full-detector field-accuracy claims.")
PY

python3 "$HERE/generate_weighting_model.py" \
  --config "$OUT/input_config.json" \
  --output "$OUT/map" > "$OUT/generate.log" 2>&1

grep -Fq 'Mesh.SecondOrderLinear = 1;' "$OUT/map/chamber.geo" ||
  { echo "STOP: not straight-sided P2" >&2; exit 2; }

cat > "$OUT/SCOPE.txt" <<'EOF'
Exploratory three-wire straight-sided quadratic tetrahedral comparison trial.
Anode cylinder 30um diameter, pitch 3mm, cathodes +/-3mm; finite reduced
10 x 6 x 8mm box. Local target 6um, curvature 32 points per full circle.
Mid-edge nodes remain linear. All differences vs V5 are exploratory.
No microscopic avalanches, frozen muon data, or electronics.
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

echo "Auditing wire surface:"
"$ENV_PY" "$HERE/audit_p3_straight_wire_surface.py" \
  "$OUT/map/chamber.msh" | tee "$OUT/wire_surface_audit.log"

# Change logging only, not the solver or approximation space.
python3 - "$OUT/map" <<'PY'
import sys
from pathlib import Path
for p in Path(sys.argv[1]).glob("*.sif"):
    s=p.read_text()
    old="  Linear System Max Iterations = 3000"
    assert s.count(old)==1
    assert "Linear System Residual Output" not in s
    p.write_text(s.replace(old,old+"\n  Linear System Residual Output = 1"))
PY

(
  cd "$OUT/map"
  nice -n 10 "$ELMER/ElmerGrid" 14 2 chamber.msh -autoclean \
    > elmergrid.log 2>&1
  cat chamber/mesh.header
  for name in drift weight_plus_0 weight_minus_0; do
    echo "Running $name with residual logging"
    if ! nice -n 10 /usr/bin/time -v "$ELMER/ElmerSolver" "$name.sif" \
      > "$name.log" 2>&1; then
      tail -n 60 "$name.log"
      exit 3
    fi
    [[ -s "chamber/$name.result" ]] || { echo "Missing $name.result"; exit 3; }
    echo "=== $name convergence lines (inspect actual residuals) ==="
    grep -Ei 'IterSolve|residu|converg|fatal|error' "$name.log" | tail -n 18 || true
    head -n 9 "chamber/$name.result"
    tail -n 4 "$name.log"
  done
)
# The same points as V5; no particle tracking in this probe.
nice -n 10 "$BIN" --map-dir "$OUT/map" \
  --out-csv "$OUT/near_wire_probe.csv" \
  --tag v6_three_wire > "$OUT/near_wire_probe.log" 2>&1 || {
    tail -n 45 "$OUT/near_wire_probe.log"; exit 3;
  }
tail -n 22 "$OUT/near_wire_probe.log"
python3 - "$REF" "$OUT/near_wire_probe.csv" <<'PY' | tee "$OUT/v5_v6_comparison.log"
import csv, math, sys
def read(p):
    with open(p) as f:
        rows=list(csv.DictReader(f))
    keys=("side","wire_index","r_over_R","theta_index","w_mm")
    d={tuple(r[k] for k in keys):r for r in rows}
    assert len(rows)==len(d)==480
    return d
ref,new=map(read,sys.argv[1:])
assert ref.keys()==new.keys()
coords=("u_mm","v_mm","w_mm","h_mm")
for key in ref:
    assert all(abs(float(ref[key][c])-float(new[key][c]))<1e-12 for c in coords)
def l2(rows,cols):
    den=sum(float(a[c])**2 for a,b in rows for c in cols)
    num=sum((float(a[c])-float(b[c]))**2 for a,b in rows for c in cols)
    return math.sqrt(num/den) if den>0 else float("nan")
fields=("Eu_halfstep_per_mm","Ev_halfstep_per_mm","Ew_halfstep_per_mm")
print("V6 vs V5: relative L2 differences (V5 reference, not truth)")
for side in ("plus","minus"):
    for shell in (1.1,1.5,2.,5.,20.):
        rows=[(ref[k],new[k]) for k in ref if k[0]==side and float(k[2])==shell]
        assert len(rows)==48
        mx=max(abs(float(a["phi"])-float(b["phi"])) for a,b in rows)
        print(f"{side:5s} r/R={shell:4g} n={len(rows)} "
              f"phi_L2={100*l2(rows,('phi',)):.4f}% "
              f"gradient_L2={100*l2(rows,fields):.4f}% "
              f"max_abs_phi_difference={mx:.8g}")
print("Comparison only. No charge/current convergence or numerical acceptance claimed.")
PY
echo "COMPLETED V6 comparison: $OUT"

#!/usr/bin/env bash
# Safe candidate NA60+ P3 FEM weighting-field generation, NO particle transport.
# Reuses original C2 Gmsh/Elmer explicit cylinder geometry with corrected
# U/V stereo signs and 1.7mm strip width/pitch measured along ANODE wires.
# Refuses to run over an existing directory; never cleans source data.
# Usage: bash .../run_p3_2100_fem_field_only.sh baseline /absolute/new/path
#        bash .../run_p3_2100_fem_field_only.sh wire_refine /absolute/new/path
#        bash .../run_p3_2100_fem_field_only.sh curved /absolute/new/path
#        bash .../run_p3_2100_fem_field_only.sh curved_refine /absolute/new/path
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
MODE="${1:-}"
OUT="${2:-}"
if [[ "$MODE" != baseline && "$MODE" != wire_refine && "$MODE" != curved && "$MODE" != curved_refine ]] || [[ "$OUT" != /* ]] || [[ -e "$OUT" ]]; then
  echo "Usage: $0 {baseline|wire_refine|curved|curved_refine} /absolute/path/TO/NEW/OUTPUT_DIR" >&2
  echo "Refusing to overwrite any existing output path or choose output implicitly." >&2
  exit 2
fi
for cmd in python3 gmsh ElmerGrid ElmerSolver sha256sum; do
  command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 2; }
done
CONFIG="$HERE/p3_2100_stereo_width_along_w_1p7mm.json"
GENERATOR="$HERE/generate_weighting_model.py"
[[ -s "$CONFIG" && -s "$GENERATOR" ]] || {
  echo "Missing dedicated config or original FEM generator" >&2; exit 2;
}
mkdir -p -- "$OUT"
cp -- "$CONFIG" "$OUT/input_config.json"
if [[ "$MODE" == wire_refine || "$MODE" == curved || "$MODE" == curved_refine ]]; then
  # wire_refine: halve near-cylinder target mesh size with straight-sided P2.
  # curved: change only high-order geometry at the same baseline mesh size.
  # curved_refine: curved geometry and halved wire mesh size (second-stage test).
  python3 - "$OUT/input_config.json" "$MODE" <<'PY'
import json,sys,pathlib
p=pathlib.Path(sys.argv[1])
c=json.loads(p.read_text(encoding="utf-8"))
assert c["mesh"]["wire_size_mm"] == .012
mode=sys.argv[2]
if mode in ("wire_refine","curved_refine"):
    c["mesh"]["wire_size_mm"] = .006
if mode in ("curved","curved_refine"):
    c["mesh"]["curved_wire_surfaces"] = True
p.write_text(json.dumps(c,indent=2)+"\n",encoding="utf-8")
PY
fi
echo "OUTPUT=$OUT"
echo "MODE=$MODE"
echo "ONLY weighting fields (+ reference drift map required for ComponentElmer)"
echo "NO microscopic avalanche, no ion motion, no currents, NO electronics"
echo "No historical folders changed."
echo "input_config:"
cat "$OUT/input_config.json"
python3 "$GENERATOR" --config "$OUT/input_config.json" --output "$OUT/map" \
  > "$OUT/generate.log" 2>&1
(
  cd "$OUT/map"
  # New clean directory checked above. Explicit Gmsh/Elmer logs retained.
  gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh \
    > gmsh.log 2>&1
  ElmerGrid 14 2 chamber.msh -autoclean \
    > elmergrid.log 2>&1
  ElmerSolver drift.sif > drift.log 2>&1
  ElmerSolver weight_minus_0.sif > weight_minus.log 2>&1
  ElmerSolver weight_plus_0.sif > weight_plus.log 2>&1
)
for f in "$OUT/map/manifest.json" "$OUT/map/chamber/mesh.header" \
         "$OUT/map/chamber/mesh.nodes" "$OUT/map/chamber/mesh.elements" \
         "$OUT/map/chamber/drift.result" \
         "$OUT/map/chamber/weight_minus_0.result" \
         "$OUT/map/chamber/weight_plus_0.result" \
         "$OUT/map/dielectrics.dat"; do
  [[ -s "$f" ]] || { echo "Missing/empty output: $f" >&2; exit 3; }
done
(
  cd "$OUT"
  sha256sum input_config.json map/manifest.json \
       map/chamber/mesh.header map/chamber/mesh.nodes \
       map/chamber/mesh.elements map/chamber/drift.result \
       map/chamber/weight_minus_0.result \
       map/chamber/weight_plus_0.result > SHA256SUMS
)
cat > "$OUT/SCOPE.txt" <<'TXT'
3D Elmer/Gmsh FEM weighting fields for P3 reference: 9 grounded cylindrical
30um anodes at 3mm pitch, cathodes at v=+/-3mm.
Physical drift solve uses 2100V; weighting potentials are independent of HV.
Selected strip unit potential, all other electrode conductors grounded.
U_plus: w-0.1u, V_minus: w+0.1u, as in modern P3 proposed geometry.
Strip pitch and width 1.7mm measured along anode wires (projected for MATC).
No physical copper gaps, PCB dielectric, finite guard geometry or electronics.
Finite 3D FEM patch has artificial natural Neumann side walls.
A completed FEM solve is NOT validated strip-gradient convergence.
No ion drift, waveform, avalanche or induced-current calculation here.
TXT
echo "COMPLETED $MODE 3D FEM MAP GENERATION (not accuracy certification)"
du -sh "$OUT" "$OUT/map" 2>/dev/null || true
echo "Map: $OUT/map"
echo "Logs: $OUT/map/gmsh.log, elmergrid.log, drift.log, weight_minus.log, weight_plus.log"
echo "Integrity: (cd \"$OUT\" && sha256sum -c SHA256SUMS)"

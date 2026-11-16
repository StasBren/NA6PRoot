#!/usr/bin/env bash
# Isolated NA60+ P3 single-ion cathode-strip validation (Garfield++ + 3D Elmer).
# Geometry generation uses the repository's existing Gmsh/Elmer script.
# Physics, integration and checks: C++/Garfield++. Plots: ROOT exclusively.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
STAGE="$REPO/test/mwpc-stage-c"
CONFIG="${SINGLE_ION_CONFIG:-$HERE/p3_single_ion_fem.json}"
OUT="${SINGLE_ION_OUT:-$REPO/test_runs/mwpc_cathode_single_ion/p3_reference_test}"
BUILD="$STAGE/build-single-ion"
if [[ -e "$OUT" ]]; then
  echo "Refusing to overwrite existing run: $OUT" >&2
  echo "To repeat choose a fresh SINGLE_ION_OUT=/absolute/path" >&2
  exit 3
fi
if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
SETUP="$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
if [[ -f "$SETUP" ]]; then source "$SETUP"; fi
for command_name in python3 gmsh ElmerGrid ElmerSolver cmake root sha256sum; do
  command -v "$command_name" >/dev/null || {
    echo "Missing dependency: $command_name" >&2; exit 2;
  }
done
[[ -s "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 2; }
MOBILITY="${MWPC_ION_MOBILITY:-}"
if [[ -z "$MOBILITY" ]]; then
  MOBILITY="$(find "$HOME/na6p/install/garfieldpp" \
    -type f -name 'IonMobility_Ar+_Ar.txt' -print -quit 2>/dev/null || true)"
fi
[[ -s "$MOBILITY" ]] || {
  echo "Positive-ion mobility file not found." >&2
  echo "Set MWPC_ION_MOBILITY=/full/path/to/ion_mobility.txt" >&2
  exit 2
}
read_config() {
  python3 -c 'import json,sys
with open(sys.argv[1],encoding="utf8") as f: v=json.load(f)
for k in sys.argv[2].split("."): v=v[k]
print(v)' "$CONFIG" "$1"
}
HV="$(read_config electrostatics.anode_voltage_V)"
GMINUS="$(read_config geometry.gap_minus_mm)"
GPLUS="$(read_config geometry.gap_plus_mm)"
PITCH="$(read_config geometry.wire_pitch_mm)"
DIAM="$(read_config geometry.wire_diameter_um)"
HALF="$(read_config geometry.half_wires)"
echo "P3 test-only field: HV=$HV V, gaps=$GMINUS+$GPLUS mm, pitch=$PITCH mm"
echo "Source branch: $(git -C "$REPO" branch --show-current)"
echo "Mobility (PROVISIONAL): $MOBILITY"
mkdir -p "$OUT/map" "$OUT/results"
cp "$CONFIG" "$OUT/input_config.json"

python3 "$STAGE/numerical-weighting/generate_weighting_model.py" \
  --config "$OUT/input_config.json" --output "$OUT/map"
(
  cd "$OUT/map"
  gmsh chamber.geo -3 -order 2 -format msh2 -optimize -o chamber.msh \
    > gmsh.log 2>&1
  ElmerGrid 14 2 chamber.msh -autoclean > elmergrid.log 2>&1
  ElmerSolver drift.sif > drift.log 2>&1
  ElmerSolver weight_minus_0.sif > weight_minus.log 2>&1
  ElmerSolver weight_plus_0.sif > weight_plus.log 2>&1
)
for f in "$OUT/map/chamber/mesh.header" "$OUT/map/chamber/drift.result" \
         "$OUT/map/chamber/weight_minus_0.result" \
         "$OUT/map/chamber/weight_plus_0.result"; do
  [[ -s "$f" ]] || { echo "Missing/empty FEM file $f" >&2; exit 4; }
done

cmake -S "$STAGE" -B "$BUILD" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$BUILD" --target mwpc_stage_c_single_ion_stereo -j2
"$BUILD/mwpc_stage_c_single_ion_stereo" \
  --map-dir "$OUT/map" --out-dir "$OUT/results" --mobility "$MOBILITY" \
  --hv "$HV" --gap-minus-mm "$GMINUS" --gap-plus-mm "$GPLUS" \
  --wire-pitch-mm "$PITCH" --wire-diameter-um "$DIAM" \
  --half-wires "$HALF" --u-mm 0 --v-mm 0.04 --w-mm 0.2 \
  --tmax-us 250 --dt-ns 10,5,2 --steps-mm 0.04,0.02 \
  | tee "$OUT/results/physics.log"
root -l -b -q "$HERE/plot_single_ion.C(\"$OUT/results\")" \
  > "$OUT/results/root_plot.log" 2>&1
{
  echo "GIT_REV=$(git -C "$REPO" rev-parse HEAD)"
  echo "DATE_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "ROOT=$(root-config --version 2>/dev/null || echo unknown)"
  echo "GARFIELD_SETUP=$SETUP"
  echo "ION_MOBILITY=$MOBILITY"
  echo "NOTE=1800V is only a historical reference; future higher HV pending scan"
  echo "MESH=finite patch with natural side-wall Neumann conditions"
  echo "STATUS=local Garfield Ramo gates and ROOT plot completed"
} > "$OUT/manifest.txt"
(
  cd "$OUT"
  sha256sum input_config.json map/manifest.json \
    map/chamber/mesh.header map/chamber/drift.result \
    map/chamber/weight_minus_0.result map/chamber/weight_plus_0.result \
    results/single_ion_summary.csv results/physics.log > SHA256SUMS
)
echo "PASS: fresh isolated run; summary, waveforms and ROOT plots: $OUT/results"
echo "Verification: (cd \"$OUT\" && sha256sum -c SHA256SUMS)"

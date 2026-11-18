#!/usr/bin/env bash
# WSL: reuse existing real 3D stereo weighting/DriftLineRKF Stage-C benchmark.
# 2100V P3 9-anode geometry. Does not touch archived microscopic outputs.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CONFIG="$HERE/p3_single_ion_fem_2100.json"
OUT="${SINGLE_ION_OUT:-$REPO/test_runs/mwpc_stage_c/p3_2100_stereo_9wire_v1}"
[[ -s "$CONFIG" ]] || { echo "Missing pinned 2100V config: $CONFIG" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "Output exists: $OUT; refusing to overwrite" >&2;exit 3; }
[[ -f "$HOME/na6p/env.sh" ]] && source "$HOME/na6p/env.sh" || true
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
export SINGLE_ION_CONFIG="$CONFIG"
export SINGLE_ION_OUT="$OUT"
export MWPC_ION_MOBILITY="${MWPC_ION_MOBILITY:-$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt}"
[[ -s "$MWPC_ION_MOBILITY" ]] || {
  echo "Missing provisional Ar+ mobility table: $MWPC_ION_MOBILITY" >&2;
  echo "Set MWPC_ION_MOBILITY=/absolute/path/to/IonMobility_Ar+_Ar.txt" >&2;
  exit 2;
}
for cmd in python3 gmsh ElmerGrid ElmerSolver cmake root root-config sha256sum;do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Missing $cmd" >&2;exit 2; }
done
[[ -f "$HOME/na6p/install/garfieldpp/lib/cmake/Garfield/GarfieldConfig.cmake" ]] || {
  echo "Missing Garfield CMake config in $HOME/na6p/install/garfieldpp" >&2; exit 2;
}
echo "2100V P3 no-field-wire / nine anodes / both central stereo strips"
echo "Output: $OUT"
echo "Provisional mobility: $MWPC_ION_MOBILITY"
df -h "$HOME"
free -h || true
bash "$HERE/run_single_ion.sh"
echo "Real stereo-strip one-ion test has returned; inspect numerical gates and maps."
echo "Summary: $OUT/results/single_ion_summary.csv"
echo "Run SHA verification: (cd \"$OUT\" && sha256sum -c SHA256SUMS)"

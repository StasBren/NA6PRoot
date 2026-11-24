#!/usr/bin/env bash
# SYNTHETIC: C++ BirthMoments core vs numpy moments_from_births under a readout geometry (not NA60+ data).
#   bash tests/test_core_geometry.sh WORKDIR GEOMETRY.json
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; SA="$HERE/../spatial-ablation"
W="${1:?WORKDIR}"; G="${2:?GEOMETRY.json}"; [[ ! -e "$W" ]] || { echo "$W exists"; exit 2; }; mkdir -p "$W"
read -r A P < <(python3 -c "import json,sys; sys.path.insert(0,'$HERE'); import weighting_geometry as wg; g=wg.Geometry.load('$G'); print(repr(g.stereo_tan), repr(g.strip_pitch_mm))")
g++ -O2 -std=c++17 -Wall -Wextra -o "$W/core" "$SA/tests/birth_moments_core_test.cc"
python3 "$SA/tests/make_synthetic_births.py" "$W/syn" 7 > /dev/null
"$W/core" "$W/syn/births.csv" "$W/syn/frame.csv" "$W/default.csv"
"$W/core" "$W/syn/births.csv" "$W/syn/frame.csv" "$W/explicit_baseline.csv" 0.1 1.7
cmp "$W/default.csv" "$W/explicit_baseline.csv" && echo "PASS core: explicit baseline geometry (0.1, 1.7) byte-identical to the default"
"$W/core" "$W/syn/births.csv" "$W/syn/frame.csv" "$W/geom.csv" "$A" "$P"
python3 - "$W" "$A" "$P" "$SA/engine" <<'PY'
import csv, sys
import numpy as np
W, A, P, E = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
sys.path.insert(0, E); import cathode_engine as ce
b = list(csv.DictReader(open(f"{W}/syn/births.csv"))); fr = next(csv.DictReader(open(f"{W}/syn/frame.csv")))
col = lambda k: np.array([float(r[k]) for r in b])
m = ce.moments_from_births(col("primary").astype(int), col("u_cm"), col("v_cm"), col("w_cm"), col("t_ns"), col("weight"),
                           float(fr["patch_shift_cm"]), float(fr["u_in_mm"]), float(fr["w_in_mm"]), stereo_tan=A, strip_pitch_mm=P)
c = list(csv.DictReader(open(f"{W}/geom.csv"))); d = list(csv.DictReader(open(f"{W}/default.csv")))
assert (int(c[0]["base_U"]), int(c[0]["base_V"])) == (m.base_U, m.base_V), "base indices"
worst = 0.
for k in ("t_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A"):
    x = np.array([float(r[k]) for r in c]); y = m.cols[k]
    worst = max(worst, float(np.max(abs(x - y) / np.maximum(1, abs(y)))))
assert worst < 1e-12, worst
same = all(r1[k] == r2[k] for r1, r2 in zip(c, d) for k in ("t_ns", "N_fC", "r_eff_um", "C1", "S1", "C2", "S2", "A"))
moved = max(abs(float(r1["D_U_mm"]) - float(r2["D_U_mm"])) + abs(float(r1["D_V_mm"]) - float(r2["D_V_mm"])) for r1, r2 in zip(c, d))
print(f"PASS core vs numpy under stereo {A} pitch {P}: worst rel diff {worst:.2g}; base U/V {c[0]['base_U']}/{c[0]['base_V']} "
      f"(default {d[0]['base_U']}/{d[0]['base_V']}); transport-side moments (t, N, r, harmonics, A) unchanged: {same}; "
      f"D columns moved (max |dD_U|+|dD_V| {moved:.3g} mm)")
assert same
PY

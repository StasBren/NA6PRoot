#!/usr/bin/env python3
"""SYNTHETIC birth records + the export the FROZEN exporter would write from them (first vertex = birth, formulas of
p3_export_full_spectral_paths.cc), for test_birth_moments.sh. Not physics. Paths are 2-vertex stubs to the cathode
(only the first vertex enters the moments)."""
import sys
from pathlib import Path

import numpy as np

out = Path(sys.argv[1]); seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = np.random.default_rng(seed)
out.mkdir(parents=True, exist_ok=False)
u_in = rng.uniform(-20, 20)
# the patch is centred on the entrance up to the wire phase (as in the archives), so strip columns stay small
frame = dict(patch_shift_cm=0.1 * u_in + rng.uniform(-0.15, 0.15), u_in_mm=u_in, w_in_mm=rng.uniform(-50, 50))
uIn, wIn = .1 * frame["u_in_mm"], .1 * frame["w_in_mm"]
base_u = int(np.round((wIn - .1 * uIn) / .17)); base_v = int(np.round((wIn + .1 * uIn) / .17))
# python round is half-even; llround is half-away: avoid exact .5 in the synthetic frame
assert abs((wIn - .1 * uIn) / .17 % 1 - .5) > 1e-6 and abs((wIn + .1 * uIn) / .17 % 1 - .5) > 1e-6
B = []
for p in range(30):
    k = rng.choice([-1, 0, 1], p=[.2, .6, .2]); n_exact = int(rng.lognormal(7, 1.2)) + 1; n = min(256, n_exact)
    tp = rng.normal(26, 10); wtrack = rng.uniform(-0.05, 0.05)
    side = rng.choice([-1, 1])
    for i in range(n):
        th = rng.normal(side * np.pi / 2, 1.2)
        r = 15.0 + (1e-3 if i == 0 and p == 0 else rng.exponential(4.0))        # one ion 1 nm above the wire
        u = 0.3 * k + 1e-4 * r * np.cos(th); v = 1e-4 * r * np.sin(th)
        if v == 0: v = 1e-9
        if p in (1, 2) and i % 2:                      # a few ions near the cell boundary: clouds on two wires
            u = 0.3 * k + rng.choice([-1, 1]) * 0.1499 + rng.normal(0, 2e-4); v = rng.normal(0, 2e-3)
        B.append((p, u, v, wtrack + rng.normal(0, 2e-3), tp + rng.normal(0, 0.1), n_exact / n))
with open(out / "births.csv", "w") as f:
    f.write("primary,u_cm,v_cm,w_cm,t_ns,weight\n")
    for b in B: f.write(",".join([str(b[0])] + [repr(float(x)) for x in b[1:]]) + "\n")
(out / "frame.csv").write_text("patch_shift_cm,u_in_mm,w_in_mm\n" + ",".join(repr(float(frame[k])) for k in frame) + "\n")
ex = out / "frozen_export"; ex.mkdir()
V = []
with open(ex / "ions.csv", "w") as f:
    f.write("row,primary,weight,birth_ns,end_ns,first_vertex,vertices,upper_index,lower_index\n")
    for j, (p, u, v, w, t, wt) in enumerate(B):
        wu = 10. * (w + wIn - .1 * frame["patch_shift_cm"] - .17 * base_u)
        wv = 10. * (w + wIn + .1 * frame["patch_shift_cm"] - .17 * base_v)
        V += [(t, 10 * u, 10 * v, wu, wv), (t + 5e4, 10 * u, 3.0 * np.sign(v), wu, wv)]
        f.write(f"{j},{p},{float(wt)!r},{float(t)!r},{float(t + 5e4)!r},{2 * j},2,{base_u},{base_v}\n")
np.array(V, "<f8").tofile(ex / "vertices.f64")
tot = sum(b[5] for b in B)
(ex / "complete.txt").write_text(f"rows={len(B)}\nexact_ions={int(round(tot))}\nweight_sum={float(tot)!r}\nvertices={len(V)}\nstep_mm=0.005\n")
print(f"SYNTHETIC births={len(B)} primaries=30 base_U={base_u} base_V={base_v}")

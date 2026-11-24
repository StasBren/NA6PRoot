#!/usr/bin/env python3
"""SYNTHETIC test of geom_restrip_export.py (not NA60+ data).

Builds an export in the exporter's strip frame from per-ion longitudinal positions w_cm and a frame (exporter
formulas, frozen geometry), then
  1. identity: restrip with the baseline geometry -> vertices.f64 and ions.csv byte-identical to the source
  2. new geometry: restriped columns vs the exporter formula evaluated directly with the new geometry (<= 1e-12 mm),
     new base indices identical, transport columns (t, u, v) byte-identical
  3. birth path: cathode_engine.moments_from_births with the new geometry gives the D_U/D_V of the restriped export
     (charge-weighted mean of the first vertex) to <= 1e-12 mm
  test_restrip.py WORKDIR GEOMETRY.json
"""
import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GS = HERE.parent
sys.path.insert(0, str(GS))
sys.path.insert(0, str(GS.parent / "spatial-ablation" / "engine"))
import weighting_geometry as wg  # noqa: E402
from geom_restrip_export import llround  # noqa: E402


def exporter_cols(w_cm, frame, A, P_mm):
    Pc = .1 * P_mm
    uIn, wIn = .1 * frame["u_in_mm"], .1 * frame["w_in_mm"]
    bu, bv = llround((wIn - A * uIn) / Pc), llround((wIn + A * uIn) / Pc)
    wu = 10. * (w_cm + wIn - A * frame["patch_shift_cm"] - Pc * bu)
    wv = 10. * (w_cm + wIn + A * frame["patch_shift_cm"] - Pc * bv)
    return bu, bv, wu, wv


def main():
    work, geomp = Path(sys.argv[1]), Path(sys.argv[2])
    work.mkdir(parents=True, exist_ok=False)
    geom = wg.Geometry.load(geomp)
    subprocess.run([sys.executable, str(GS.parent / "timing-ablation" / "make_synthetic_export.py"), str(work / "raw"),
                    "--seed", "11"], check=True, stdout=subprocess.DEVNULL)
    rng = np.random.default_rng(5)
    frame = {"patch_shift_cm": 0.37, "u_in_mm": 2.3, "w_in_mm": -31.7}
    rows = list(csv.DictReader(open(work / "raw" / "ions.csv")))
    nv = int(dict(l.split("=") for l in (work / "raw" / "complete.txt").read_text().split())["vertices"])
    V = np.fromfile(work / "raw" / "vertices.f64", dtype="<f8").reshape(nv, 5)
    w_cm = rng.uniform(-.08, .08, len(rows))
    bu, bv, _, _ = exporter_cols(0., frame, .1, 1.7)
    frame.update(base_U=bu, base_V=bv, stereo_tan=.1, strip_pitch_mm=1.7)
    src = work / "src"; src.mkdir()
    for r, wc in zip(rows, w_cm):
        s, n = int(r["first_vertex"]), int(r["vertices"])
        _, _, wu, wv = exporter_cols(wc, frame, .1, 1.7)
        V[s:s + n, 3], V[s:s + n, 4] = wu, wv
        r["upper_index"], r["lower_index"] = str(bu), str(bv)
    V.astype("<f8").tofile(src / "vertices.f64")
    with open(src / "ions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (src / "complete.txt").write_text((work / "raw" / "complete.txt").read_text())
    (work / "frame.json").write_text(json.dumps(frame))
    R = GS / "geom_restrip_export.py"
    subprocess.run([sys.executable, str(R), str(src), str(work / "frame.json"), str(GS / "geometries" / "baseline.json"),
                    str(work / "id")], check=True)
    for f in ("vertices.f64", "ions.csv", "complete.txt"):
        assert (src / f).read_bytes() == (work / "id" / f).read_bytes(), f"identity restrip changed {f}"
    print("PASS identity restrip: byte-identical")
    subprocess.run([sys.executable, str(R), str(src), str(work / "frame.json"), str(geomp), str(work / "new")], check=True)
    Vn = np.fromfile(work / "new" / "vertices.f64", dtype="<f8").reshape(nv, 5)
    assert Vn[:, :3].tobytes() == V[:, :3].tobytes(), "transport columns changed"
    nb = exporter_cols(0., frame, geom.stereo_tan, geom.strip_pitch_mm)[:2]
    worst = 0.
    for r, wc in zip(csv.DictReader(open(work / "new" / "ions.csv")), w_cm):
        assert (int(r["upper_index"]), int(r["lower_index"])) == nb, "new base indices"
        s, n = int(r["first_vertex"]), int(r["vertices"])
        _, _, wu, wv = exporter_cols(wc, frame, geom.stereo_tan, geom.strip_pitch_mm)
        worst = max(worst, np.max(abs(Vn[s:s + n, 3] - wu)), np.max(abs(Vn[s:s + n, 4] - wv)))
    assert worst <= 1e-12, worst
    print(f"PASS restrip {geom.name}: base {(bu, bv)} -> {nb}, max |restriped - direct exporter formula| = {worst:.3g} mm,"
          " transport columns byte-identical")
    # birth path with the new geometry
    import cathode_engine as ce
    first = np.array([int(r["first_vertex"]) for r in rows])
    p = np.array([int(r["primary"]) for r in rows]); wt = np.array([float(r["weight"]) for r in rows])
    u_cm, v_cm, t = .1 * V[first, 1], .1 * V[first, 2], V[first, 0]
    m = ce.moments_from_births(p, u_cm, v_cm, w_cm, t, wt, frame["patch_shift_cm"], frame["u_in_mm"], frame["w_in_mm"],
                               stereo_tan=geom.stereo_tan, strip_pitch_mm=geom.strip_pitch_mm)
    assert (m.base_U, m.base_V) == nb
    k = np.rint(10 * u_cm / 3).astype(int)
    keys, cl = np.unique(np.stack([p, k], 1), axis=0, return_inverse=True); cl = cl.ravel()
    q = wt * ce.EFC
    du = np.bincount(cl, q * Vn[first, 3]) / np.bincount(cl, q); dv = np.bincount(cl, q * Vn[first, 4]) / np.bincount(cl, q)
    e = max(np.max(abs(du - m.cols["D_U_mm"])), np.max(abs(dv - m.cols["D_V_mm"])))
    assert e <= 1e-12, e
    print(f"PASS birth-path moments with {geom.name}: D_U/D_V agree with the restriped export to {e:.3g} mm")


if __name__ == "__main__":
    main()

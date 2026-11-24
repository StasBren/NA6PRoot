#!/usr/bin/env python3
"""Validate the Green-grid executable against the frozen exporter: same start point -> same trajectory.

  s_validate_transport.py points  EXPORT N POINTS.csv        pick N ion rows of a frozen/T1 export (spread over the
                                                             file), write their start points (u_mm, v_mm, t_ns)
  p3_green_grid_paths MOBILITY OUT --points POINTS.csv       drift them with the grid tool
  s_validate_transport.py compare EXPORT POINTS.csv OUT      vertex-by-vertex comparison

PASS requires identical vertex counts and |du|,|dv| < 1e-6 mm, |d(t - t_birth)| < 1e-6 ns on every vertex. The only
expected differences are rounding from the mm -> cm conversion of the start point (the frozen exporter starts from
the archive's cm value). Any larger difference means the grid tool does not reproduce the frozen transport.
"""
import argparse
import csv
from pathlib import Path

import numpy as np


def load(d):
    seal = dict(l.split("=", 1) for l in (d / "complete.txt").read_text().split() if "=" in l)
    rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
    V = np.memmap(d / "vertices.f64", dtype="<f8", mode="r", shape=(int(seal["vertices"]), 5))
    return rows, V


def path(rows, V, i):
    f, n = int(rows[i]["first_vertex"]), int(rows[i]["vertices"])
    return np.asarray(V[f:f + n])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("points"); p.add_argument("export", type=Path); p.add_argument("n", type=int); p.add_argument("out", type=Path)
    c = sub.add_parser("compare"); c.add_argument("export", type=Path); c.add_argument("points", type=Path); c.add_argument("grid", type=Path)
    a = ap.parse_args()
    if a.cmd == "points":
        rows, V = load(a.export)
        sel = np.unique(np.linspace(0, len(rows) - 1, min(a.n, len(rows))).astype(int))
        with open(a.out, "w") as f:
            f.write("u_mm,v_mm,t_ns\n")
            for i in sel:
                p0 = path(rows, V, i)[0]
                f.write(f"{p0[1]!r},{p0[2]!r},{p0[0]!r}\n")
        (a.out.with_suffix(".rows")).write_text("\n".join(map(str, sel)) + "\n")
        print(f"wrote {len(sel)} start points")
        return
    rows, V = load(a.export)
    grows, GV = load(a.grid)
    sel = [int(x) for x in a.points.with_suffix(".rows").read_text().split()]
    if len(sel) != len(grows):
        raise SystemExit("row count mismatch")
    worst = dict(n_mismatch=0, du=0.0, dv=0.0, dt=0.0)
    for j, i in enumerate(sel):
        pe, pg = path(rows, V, i), path(grows, GV, j)
        if len(pe) != len(pg):
            worst["n_mismatch"] += 1
            continue
        worst["du"] = max(worst["du"], float(np.max(np.abs(pe[:, 1] - pg[:, 1]))))
        worst["dv"] = max(worst["dv"], float(np.max(np.abs(pe[:, 2] - pg[:, 2]))))
        worst["dt"] = max(worst["dt"], float(np.max(np.abs((pe[:, 0] - pe[0, 0]) - (pg[:, 0] - pg[0, 0])))))
    ok = worst["n_mismatch"] == 0 and worst["du"] < 1e-6 and worst["dv"] < 1e-6 and worst["dt"] < 1e-6
    print(f"{'PASS' if ok else 'FAIL'} grid tool reproduces frozen export trajectories: {len(sel)} paths, "
          f"vertex-count mismatches {worst['n_mismatch']}, max |du| {worst['du']:.3g} mm, max |dv| {worst['dv']:.3g} mm, "
          f"max |d(t-t0)| {worst['dt']:.3g} ns")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()

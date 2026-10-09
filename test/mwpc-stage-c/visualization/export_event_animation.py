#!/usr/bin/env python3
"""
Build a ParaView time-series animation from the Stage-C2 event dump.

The animation is intentionally cumulative:
  * Heed clusters appear when their cluster time is reached.
  * primary electrons appear at their creation time.
  * microscopic electron trajectories grow point-by-point with their Garfield time.
  * avalanche ion-birth points accumulate as multiplication develops.

Detector geometry, FEM mesh and field slice remain static.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fv(row, key):
    return float(row[key])


def iv(row, key):
    return int(float(row[key]))


def write_points(path: Path, rows, scalar_defs=None):
    scalar_defs = scalar_defs or []
    pts = [(fv(r, "u_mm"), fv(r, "v_mm"), fv(r, "w_mm")) for r in rows]
    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write(path.stem + "\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(pts)} float\n")
        for x, y, z in pts:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        out.write(f"VERTICES {len(pts)} {2 * len(pts)}\n")
        for idx in range(len(pts)):
            out.write(f"1 {idx}\n")
        if pts and scalar_defs:
            out.write(f"POINT_DATA {len(pts)}\n")
            for name, key, typ in scalar_defs:
                out.write(f"SCALARS {name} {typ} 1\nLOOKUP_TABLE default\n")
                for row in rows:
                    out.write(
                        f"{iv(row, key)}\n" if typ == "int"
                        else f"{fv(row, key):.9g}\n"
                    )


def write_truncated_paths(path: Path, rows, tmax: float):
    groups = defaultdict(list)
    for row in rows:
        if fv(row, "t_ns") <= tmax:
            groups[(iv(row, "seed"), iv(row, "path"))].append(row)

    points = []
    lines = []
    line_meta = []
    for key, grows in groups.items():
        grows.sort(key=lambda r: iv(r, "point"))
        if len(grows) < 2:
            continue
        start = len(points)
        for row in grows:
            points.append((fv(row, "u_mm"), fv(row, "v_mm"), fv(row, "w_mm")))
        lines.append(list(range(start, start + len(grows))))
        line_meta.append(key)

    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("Growing microscopic electron paths\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(points)} float\n")
        for x, y, z in points:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        total = sum(len(line) + 1 for line in lines)
        out.write(f"LINES {len(lines)} {total}\n")
        for line in lines:
            out.write(str(len(line)) + " " + " ".join(map(str, line)) + "\n")
        if lines:
            out.write(f"CELL_DATA {len(lines)}\n")
            out.write("SCALARS seed int 1\nLOOKUP_TABLE default\n")
            for seed, _ in line_meta:
                out.write(f"{seed}\n")


def write_pvd(path: Path, entries):
    with path.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        out.write("  <Collection>\n")
        for time, filename in entries:
            out.write(
                f'    <DataSet timestep="{time:.9g}" group="" part="0" '
                f'file="{filename}"/>\n'
            )
        out.write("  </Collection>\n")
        out.write("</VTKFile>\n")


def quantile(values, q):
    if not values:
        return 0.0
    vals = sorted(values)
    pos = q * (len(vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    f = pos - lo
    return vals[lo] * (1 - f) + vals[hi] * f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event-dir", required=True, type=Path)
    ap.add_argument("--prefix", default="c2_event3d")
    ap.add_argument("--vtk-dir", required=True, type=Path)
    ap.add_argument("--frames", type=int, default=32)
    ap.add_argument("--tail-quantile", type=float, default=0.995)
    args = ap.parse_args()

    if args.frames < 4:
        raise SystemExit("--frames must be >= 4")

    p = args.prefix
    clusters = read_rows(args.event_dir / f"{p}_clusters.csv")
    primaries = read_rows(args.event_dir / f"{p}_primary_electrons.csv")
    paths = read_rows(args.event_dir / f"{p}_electron_paths.csv")
    ions = read_rows(args.event_dir / f"{p}_ion_births.csv")

    times = []
    times += [fv(r, "t_ns") for r in clusters]
    times += [fv(r, "t_ns") for r in primaries]
    times += [fv(r, "t_ns") for r in paths]
    times += [fv(r, "t_ns") for r in ions]
    if not times:
        raise SystemExit("No event times found.")

    t0 = min(times)
    # A tiny number of late Garfield path points can otherwise make the useful
    # avalanche evolution occupy only the first few frames.
    t1 = quantile(times, max(0.9, min(args.tail_quantile, 1.0)))
    if t1 <= t0:
        t1 = max(times)
    if t1 <= t0:
        t1 = t0 + 1.0

    anim = args.vtk_dir / "animation"
    anim.mkdir(parents=True, exist_ok=True)

    cluster_entries = []
    primary_entries = []
    path_entries = []
    ion_entries = []

    for iframe in range(args.frames):
        frac = iframe / (args.frames - 1)
        t = t0 + frac * (t1 - t0)

        crows = [r for r in clusters if fv(r, "t_ns") <= t]
        prows = [r for r in primaries if fv(r, "t_ns") <= t]
        irows = [r for r in ions if fv(r, "t_ns") <= t]

        cfile = f"clusters_{iframe:04d}.vtk"
        pfile = f"primaries_{iframe:04d}.vtk"
        efile = f"electron_paths_{iframe:04d}.vtk"
        ifile = f"ion_births_{iframe:04d}.vtk"

        write_points(
            anim / cfile, crows,
            [("n_primary", "n_primary", "int"), ("energy_eV", "energy_eV", "float")]
        )
        write_points(
            anim / pfile, prows,
            [("seed", "seed", "int"), ("energy_eV", "energy_eV", "float")]
        )
        write_truncated_paths(anim / efile, paths, t)
        write_points(
            anim / ifile, irows,
            [("seed", "seed", "int")]
        )

        rel = f"animation/{cfile}"
        cluster_entries.append((t, rel))
        primary_entries.append((t, f"animation/{pfile}"))
        path_entries.append((t, f"animation/{efile}"))
        ion_entries.append((t, f"animation/{ifile}"))

    write_pvd(args.vtk_dir / "clusters_time.pvd", cluster_entries)
    write_pvd(args.vtk_dir / "primary_electrons_time.pvd", primary_entries)
    write_pvd(args.vtk_dir / "electron_paths_time.pvd", path_entries)
    write_pvd(args.vtk_dir / "ion_births_time.pvd", ion_entries)

    script = args.vtk_dir / "open_animation.py"
    script.write_text(
f'''from paraview.simple import *
from pathlib import Path

base = Path(r"{args.vtk_dir}")

def legacy(name):
    return LegacyVTKReader(registrationName=name, FileNames=[str(base / name)])

# Static geometry and field.
cath = legacy("cathodes.vtk")
dc = Show(cath)
dc.Representation = "Surface"
dc.Opacity = 0.12
dc.DiffuseColor = [0.75, 0.75, 0.78]

wires = legacy("wire_axes.vtk")
tube = Tube(registrationName="Anode wires", Input=wires)
tube.Radius = 0.015
tube.NumberofSides = 16
dw = Show(tube)
dw.DiffuseColor = [0.15, 0.15, 0.15]

mu = legacy("muon_track.vtk")
mu_t = Tube(registrationName="Muon track", Input=mu)
mu_t.Radius = 0.035
mu_t.NumberofSides = 12
dm = Show(mu_t)
dm.DiffuseColor = [0.9, 0.15, 0.1]

field = legacy("field_slice_uv.vtk")
tri = Delaunay2D(registrationName="E-field slice", Input=field)
df = Show(tri)
df.Representation = "Surface"
df.Opacity = 0.30
ColorBy(df, ("POINTS", "E_mag_Vcm"))
df.RescaleTransferFunctionToDataRange(True, False)

# Temporal layers.
clusters = PVDReader(registrationName="Heed clusters (time)", FileName=str(base / "clusters_time.pvd"))
dcl = Show(clusters)
dcl.Representation = "Point Gaussian"
dcl.GaussianRadius = 0.075
dcl.DiffuseColor = [0.15, 0.55, 0.95]

prim = PVDReader(registrationName="Primary electrons (time)", FileName=str(base / "primary_electrons_time.pvd"))
dp = Show(prim)
dp.Representation = "Point Gaussian"
dp.GaussianRadius = 0.035
dp.DiffuseColor = [0.25, 0.85, 0.95]

paths = PVDReader(registrationName="Electron paths (time)", FileName=str(base / "electron_paths_time.pvd"))
dpa = Show(paths)
dpa.Representation = "Wireframe"
dpa.LineWidth = 1.0
dpa.DiffuseColor = [0.2, 0.45, 1.0]
dpa.Opacity = 0.55

ions = PVDReader(registrationName="Avalanche ion births (time)", FileName=str(base / "ion_births_time.pvd"))
di = Show(ions)
di.Representation = "Point Gaussian"
di.GaussianRadius = 0.018
di.DiffuseColor = [1.0, 0.55, 0.1]
di.Opacity = 0.50

mesh = XMLUnstructuredGridReader(
    registrationName="FEM mesh (hidden)", FileName=[str(base / "chamber.vtu")])
dmesh = Show(mesh)
dmesh.Representation = "Wireframe"
dmesh.Opacity = 0.08
Hide(mesh)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.PlayMode = "Snap To TimeSteps"

view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 1
view.AxesGrid.Visibility = 1
view.Background = [0.96, 0.96, 0.96]
ResetCamera()
Render()

print("Animated C2 scene loaded.")
print("Use the Play button in the Animation toolbar or press Space.")
print("Time range: {t0:.6g} to {t1:.6g} ns, frames: {args.frames}.")
''',
        encoding="utf-8",
    )

    print(f"Animation time range: {t0:.6g} .. {t1:.6g} ns")
    print(f"Frames: {args.frames}")
    print("Wrote:", script)
    print("Open with:")
    print(f'  paraview --script="{script}"')


if __name__ == "__main__":
    main()

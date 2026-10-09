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


def _write_data_array(out, name, values, vtk_type="Float64", ncomp=1):
    attrs = [f'type="{vtk_type}"', 'format="ascii"']
    if name:
        attrs.insert(1, f'Name="{name}"')
    if ncomp != 1:
        attrs.insert(2, f'NumberOfComponents="{ncomp}"')
    out.write("        <DataArray " + " ".join(attrs) + ">\n")
    if ncomp == 1:
        out.write("          " + " ".join(str(v) for v in values) + "\n")
    else:
        out.write(
            "          " +
            " ".join(" ".join(str(x) for x in value) for value in values) +
            "\n"
        )
    out.write("        </DataArray>\n")


def write_points(path: Path, rows, scalar_defs=None):
    """Write XML VTK PolyData (.vtp) point cloud.

    PVD collections are XML collection readers in ParaView and must reference
    XML VTK datasets (VTP/VTU), not legacy .vtk files.
    """
    scalar_defs = scalar_defs or []
    pts = [(fv(r, "u_mm"), fv(r, "v_mm"), fv(r, "w_mm")) for r in rows]
    n = len(pts)

    with path.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write(
            '<VTKFile type="PolyData" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        out.write("  <PolyData>\n")
        out.write(
            f'    <Piece NumberOfPoints="{n}" NumberOfVerts="{n}" '
            'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n'
        )

        out.write("      <PointData>\n")
        for name, key, typ in scalar_defs:
            vtk_type = "Int32" if typ == "int" else "Float64"
            vals = [iv(r, key) if typ == "int" else fv(r, key) for r in rows]
            _write_data_array(out, name, vals, vtk_type)
        out.write("      </PointData>\n")
        out.write("      <CellData/>\n")

        out.write("      <Points>\n")
        _write_data_array(out, None, pts, "Float64", 3)
        out.write("      </Points>\n")

        out.write("      <Verts>\n")
        _write_data_array(out, "connectivity", list(range(n)), "Int64")
        _write_data_array(out, "offsets", list(range(1, n + 1)), "Int64")
        out.write("      </Verts>\n")
        out.write("      <Lines/>\n")
        out.write("      <Strips/>\n")
        out.write("      <Polys/>\n")
        out.write("    </Piece>\n")
        out.write("  </PolyData>\n")
        out.write("</VTKFile>\n")


def write_truncated_paths(path: Path, rows, tmax: float):
    """Write cumulative microscopic electron trajectories as XML VTP."""
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

    connectivity = [idx for line in lines for idx in line]
    offsets = []
    running = 0
    for line in lines:
        running += len(line)
        offsets.append(running)

    with path.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write(
            '<VTKFile type="PolyData" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        out.write("  <PolyData>\n")
        out.write(
            f'    <Piece NumberOfPoints="{len(points)}" NumberOfVerts="0" '
            f'NumberOfLines="{len(lines)}" NumberOfStrips="0" '
            'NumberOfPolys="0">\n'
        )
        out.write("      <PointData/>\n")
        out.write("      <CellData>\n")
        _write_data_array(
            out, "seed", [seed for seed, _ in line_meta], "Int32"
        )
        _write_data_array(
            out, "path", [path_id for _, path_id in line_meta], "Int32"
        )
        out.write("      </CellData>\n")

        out.write("      <Points>\n")
        _write_data_array(out, None, points, "Float64", 3)
        out.write("      </Points>\n")

        out.write("      <Verts/>\n")
        out.write("      <Lines>\n")
        _write_data_array(out, "connectivity", connectivity, "Int64")
        _write_data_array(out, "offsets", offsets, "Int64")
        out.write("      </Lines>\n")
        out.write("      <Strips/>\n")
        out.write("      <Polys/>\n")
        out.write("    </Piece>\n")
        out.write("  </PolyData>\n")
        out.write("</VTKFile>\n")


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

        cfile = f"clusters_{iframe:04d}.vtp"
        pfile = f"primaries_{iframe:04d}.vtp"
        efile = f"electron_paths_{iframe:04d}.vtp"
        ifile = f"ion_births_{iframe:04d}.vtp"

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

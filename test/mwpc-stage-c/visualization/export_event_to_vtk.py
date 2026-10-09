#!/usr/bin/env python3
"""
Convert the C2 event-visualization CSV dump plus Gmsh chamber mesh into a small
set of VTK files that ParaView can load together.

Only meshio is required for the chamber .msh -> .vtu conversion.  All event
polydata are written directly in legacy ASCII VTK for maximum portability.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import meshio


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row, key):
    return float(row[key])


def i(row, key):
    return int(float(row[key]))


def write_points(path: Path, rows, point_columns, scalars=None):
    scalars = scalars or []
    pts = [(f(r, point_columns[0]), f(r, point_columns[1]), f(r, point_columns[2]))
           for r in rows]

    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write(path.stem + "\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(pts)} float\n")
        for x, y, z in pts:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")

        out.write(f"VERTICES {len(pts)} {2 * len(pts)}\n")
        for idx in range(len(pts)):
            out.write(f"1 {idx}\n")

        if pts and scalars:
            out.write(f"POINT_DATA {len(pts)}\n")
            for name, source, typ in scalars:
                out.write(f"SCALARS {name} {typ} 1\nLOOKUP_TABLE default\n")
                for row in rows:
                    if typ == "int":
                        out.write(f"{i(row, source)}\n")
                    else:
                        out.write(f"{f(row, source):.9g}\n")


def write_grouped_lines(path: Path, rows, group_keys, point_columns, cell_scalars=None):
    groups = defaultdict(list)
    for row in rows:
        key = tuple(i(row, k) for k in group_keys)
        groups[key].append(row)

    for key in groups:
        groups[key].sort(key=lambda r: i(r, "point"))

    points = []
    lines = []
    keys = []
    for key, grows in groups.items():
        start = len(points)
        for row in grows:
            points.append(
                (f(row, point_columns[0]), f(row, point_columns[1]), f(row, point_columns[2]))
            )
        if len(grows) >= 2:
            lines.append(list(range(start, start + len(grows))))
            keys.append(key)

    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write(path.stem + "\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(points)} float\n")
        for x, y, z in points:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")

        total = sum(len(line) + 1 for line in lines)
        out.write(f"LINES {len(lines)} {total}\n")
        for line in lines:
            out.write(str(len(line)) + " " + " ".join(map(str, line)) + "\n")

        if lines and cell_scalars:
            out.write(f"CELL_DATA {len(lines)}\n")
            for name, key_index in cell_scalars:
                out.write(f"SCALARS {name} int 1\nLOOKUP_TABLE default\n")
                for key in keys:
                    out.write(f"{key[key_index]}\n")


def write_muon_line(path: Path, rows):
    points = [(f(r, "u_mm"), f(r, "v_mm"), f(r, "w_mm")) for r in rows]
    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("Muon track\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(points)} float\n")
        for x, y, z in points:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        out.write(f"LINES 1 {len(points) + 1}\n")
        out.write(str(len(points)) + " " + " ".join(str(j) for j in range(len(points))) + "\n")


def write_field_points(path: Path, rows):
    valid = [r for r in rows if i(r, "valid") == 1]
    pts = [(f(r, "u_mm"), f(r, "v_mm"), f(r, "w_mm")) for r in valid]
    with path.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("C2 numerical field slice\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(pts)} float\n")
        for x, y, z in pts:
            out.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        out.write(f"VERTICES {len(pts)} {2 * len(pts)}\n")
        for idx in range(len(pts)):
            out.write(f"1 {idx}\n")

        out.write(f"POINT_DATA {len(pts)}\n")
        for name, source in [
            ("E_mag_Vcm", "E_mag_Vcm"),
            ("potential_V", "potential_V"),
            ("phi_minus_0", "phi_minus_0"),
            ("phi_plus_0", "phi_plus_0"),
        ]:
            out.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n")
            for row in valid:
                out.write(f"{f(row, source):.9g}\n")

        out.write("VECTORS E_Vcm float\n")
        for row in valid:
            out.write(
                f"{f(row, 'eu_Vcm'):.9g} {f(row, 'ev_Vcm'):.9g} {f(row, 'ew_Vcm'):.9g}\n"
            )


def write_detector_geometry(out_dir: Path, config):
    g = config["geometry"]
    u_half = float(g["u_half_span_mm"])
    w_half = float(g["w_half_span_mm"])
    gap_minus = float(g["gap_minus_mm"])
    gap_plus = float(g["gap_plus_mm"])
    pitch = float(g["wire_pitch_mm"])
    half_wires = int(g["half_wires"])

    cath = out_dir / "cathodes.vtk"
    with cath.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("Cathode planes\nASCII\nDATASET POLYDATA\n")
        pts = [
            (-u_half, -gap_minus, -w_half), (u_half, -gap_minus, -w_half),
            (u_half, -gap_minus, w_half), (-u_half, -gap_minus, w_half),
            (-u_half, gap_plus, -w_half), (u_half, gap_plus, -w_half),
            (u_half, gap_plus, w_half), (-u_half, gap_plus, w_half),
        ]
        out.write(f"POINTS {len(pts)} float\n")
        for p in pts:
            out.write(f"{p[0]} {p[1]} {p[2]}\n")
        out.write("POLYGONS 2 10\n")
        out.write("4 0 1 2 3\n")
        out.write("4 4 5 6 7\n")
        out.write("CELL_DATA 2\n")
        out.write("SCALARS cathode int 1\nLOOKUP_TABLE default\n0\n1\n")

    wires = out_dir / "wire_axes.vtk"
    points = []
    lines = []
    for k in range(-half_wires, half_wires + 1):
        idx = len(points)
        u = k * pitch
        points.append((u, 0.0, -w_half))
        points.append((u, 0.0, w_half))
        lines.append((idx, idx + 1, k))

    with wires.open("w", encoding="utf-8") as out:
        out.write("# vtk DataFile Version 3.0\n")
        out.write("Anode wire axes\nASCII\nDATASET POLYDATA\n")
        out.write(f"POINTS {len(points)} float\n")
        for p in points:
            out.write(f"{p[0]} {p[1]} {p[2]}\n")
        out.write(f"LINES {len(lines)} {3 * len(lines)}\n")
        for a, b, _ in lines:
            out.write(f"2 {a} {b}\n")
        out.write(f"CELL_DATA {len(lines)}\n")
        out.write("SCALARS wire_index int 1\nLOOKUP_TABLE default\n")
        for _, _, k in lines:
            out.write(f"{k}\n")


def write_paraview_script(out_dir: Path, wire_radius_mm: float):
    script = out_dir / "open_scene.py"
    content = f'''from paraview.simple import *
from pathlib import Path

DisableFirstRenderCameraReset()
base = Path(r"{out_dir}")

def reader(name):
    return LegacyVTKReader(registrationName=name, FileNames=[str(base / name)])

# Transparent detector envelope.
cath = reader("cathodes.vtk")
dc = Show(cath)
dc.Representation = "Surface"
dc.Opacity = 0.12
dc.DiffuseColor = [0.75, 0.75, 0.78]

wires = reader("wire_axes.vtk")
tube = Tube(registrationName="Anode wires", Input=wires)
tube.Radius = {wire_radius_mm}
tube.NumberofSides = 16
dw = Show(tube)
dw.DiffuseColor = [0.15, 0.15, 0.15]

mu = reader("muon_track.vtk")
mu_t = Tube(registrationName="Muon", Input=mu)
mu_t.Radius = 0.035
mu_t.NumberofSides = 12
dm = Show(mu_t)
dm.DiffuseColor = [0.9, 0.15, 0.1]

clusters = reader("clusters.vtk")
dcl = Show(clusters)
dcl.Representation = "Point Gaussian"
dcl.GaussianRadius = 0.075
dcl.DiffuseColor = [0.15, 0.55, 0.95]

prim = reader("primary_electrons.vtk")
dp = Show(prim)
dp.Representation = "Point Gaussian"
dp.GaussianRadius = 0.035
dp.DiffuseColor = [0.25, 0.85, 0.95]

paths = reader("electron_drift_lines.vtk")
dpa = Show(paths)
dpa.Representation = "Wireframe"
dpa.LineWidth = 1.0
dpa.DiffuseColor = [0.2, 0.45, 1.0]
dpa.Opacity = 0.45

ions = reader("avalanche_ion_births.vtk")
di = Show(ions)
di.Representation = "Point Gaussian"
di.GaussianRadius = 0.018
di.DiffuseColor = [1.0, 0.55, 0.1]
di.Opacity = 0.45

field = reader("field_slice_uv.vtk")
tri = Delaunay2D(registrationName="E-field slice", Input=field)
df = Show(tri)
df.Representation = "Surface"
df.Opacity = 0.42
ColorBy(df, ("POINTS", "E_mag_Vcm"))
df.RescaleTransferFunctionToDataRange(True, False)

# FEM mesh is available but starts hidden because it is visually dense.
mesh = XMLUnstructuredGridReader(
    registrationName="FEM mesh", FileName=[str(base / "chamber.vtu")])
dmesh = Show(mesh)
dmesh.Representation = "Wireframe"
dmesh.Opacity = 0.08
Hide(mesh)

view = GetActiveViewOrCreate("RenderView")
view.OrientationAxesVisibility = 1
view.AxesGrid.Visibility = 1
view.Background = [0.96, 0.96, 0.96]
ResetCamera()
Render()

print("C2 3D scene loaded.")
print("Use the Pipeline Browser to toggle FEM mesh, field slice, paths and points.")
'''
    script.write_text(content, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event-dir", required=True, type=Path)
    ap.add_argument("--map-dir", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--prefix", default="c2_event3d")
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(args.config.read_text(encoding="utf-8"))

    mesh = meshio.read(args.map_dir / "chamber.msh")
    meshio.write(args.output / "chamber.vtu", mesh)

    p = args.prefix
    muon = read_rows(args.event_dir / f"{p}_muon_track.csv")
    clusters = read_rows(args.event_dir / f"{p}_clusters.csv")
    primaries = read_rows(args.event_dir / f"{p}_primary_electrons.csv")
    paths = read_rows(args.event_dir / f"{p}_electron_paths.csv")
    ions = read_rows(args.event_dir / f"{p}_ion_births.csv")
    field = read_rows(args.event_dir / f"{p}_field_slice_uv.csv")

    write_muon_line(args.output / "muon_track.vtk", muon)
    write_points(
        args.output / "clusters.vtk", clusters, ("u_mm", "v_mm", "w_mm"),
        [("n_primary", "n_primary", "int"), ("energy_eV", "energy_eV", "float")],
    )
    write_points(
        args.output / "primary_electrons.vtk", primaries,
        ("u_mm", "v_mm", "w_mm"),
        [("seed", "seed", "int"), ("energy_eV", "energy_eV", "float")],
    )
    write_grouped_lines(
        args.output / "electron_drift_lines.vtk", paths,
        ("seed", "path"), ("u_mm", "v_mm", "w_mm"),
        [("seed", 0), ("path", 1)],
    )
    write_points(
        args.output / "avalanche_ion_births.vtk", ions,
        ("u_mm", "v_mm", "w_mm"),
        [("seed", "seed", "int")],
    )
    write_field_points(args.output / "field_slice_uv.vtk", field)
    write_detector_geometry(args.output, cfg)

    wire_radius_mm = 0.0005 * float(cfg["geometry"]["wire_diameter_um"])
    write_paraview_script(args.output, wire_radius_mm)

    print("Wrote ParaView scene files to", args.output)
    print("Open with:")
    print(f'  paraview --script="{args.output / "open_scene.py"}"')


if __name__ == "__main__":
    main()

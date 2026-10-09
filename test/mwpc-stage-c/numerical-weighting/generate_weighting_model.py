#!/usr/bin/env python3
"""
Generate a local 3-D Gmsh + Elmer weighting-field model for the MWPC.

The geometry is deliberately small and parameterised.  It contains:
  * gas between two cathodes,
  * explicit cylindrical anode wires,
  * one or more cathode-strip weighting solves.

The cathode segmentation is imposed as a position-dependent Dirichlet
boundary condition on the physical cathode plane.  This avoids geometrically
fragmenting the cathode surface while keeping the conductor boundary
conditions on the real anode wires.

Coordinates follow Stage C:
  Gmsh/Elmer x = local u  (across wires)
  Gmsh/Elmer y = local v  (chamber normal)
  Gmsh/Elmer z = local w  (along wires)
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable


def _need(mapping: Dict[str, Any], key: str) -> Any:
    if key not in mapping:
        raise KeyError(f"Missing required configuration key: {key}")
    return mapping[key]


def _fmt(x: float) -> str:
    return f"{x:.12g}"


def validate(cfg: Dict[str, Any]) -> None:
    g = _need(cfg, "geometry")
    e = _need(cfg, "electrostatics")
    m = _need(cfg, "mesh")
    w = _need(cfg, "weighting")

    positive = [
        ("gap_minus_mm", g["gap_minus_mm"]),
        ("gap_plus_mm", g["gap_plus_mm"]),
        ("wire_pitch_mm", g["wire_pitch_mm"]),
        ("wire_diameter_um", g["wire_diameter_um"]),
        ("u_half_span_mm", g["u_half_span_mm"]),
        ("w_half_span_mm", g["w_half_span_mm"]),
        ("strip_pitch_mm", g["strip_pitch_mm"]),
        ("strip_width_mm", g["strip_width_mm"]),
        ("bulk_size_mm", m["bulk_size_mm"]),
        ("cathode_size_mm", m["cathode_size_mm"]),
        ("cathode_refine_depth_mm", m["cathode_refine_depth_mm"]),
        ("wire_size_mm", m["wire_size_mm"]),
        ("wire_refine_radius_mm", m["wire_refine_radius_mm"]),
    ]
    for name, value in positive:
        if float(value) <= 0:
            raise ValueError(f"{name} must be > 0, got {value}")

    if int(g["half_wires"]) < 1:
        raise ValueError("half_wires must be >= 1")
    if float(g["strip_width_mm"]) > float(g["strip_pitch_mm"]) + 1e-12:
        raise ValueError("strip_width_mm must not exceed strip_pitch_mm")
    if float(g["tan_alpha"]) <= 0:
        raise ValueError("tan_alpha must be > 0")
    if float(e["gas_relative_permittivity"]) <= 0:
        raise ValueError("gas_relative_permittivity must be > 0")
    if int(m["wire_curvature_points"]) < 4:
        raise ValueError("wire_curvature_points must be >= 4")

    cathodes = list(w["cathodes"])
    allowed = {"minus", "plus"}
    if not cathodes or any(c not in allowed for c in cathodes):
        raise ValueError("weighting.cathodes must contain only 'minus'/'plus'")
    if not list(w["strip_indices"]):
        raise ValueError("weighting.strip_indices must not be empty")

    wire_pitch = float(g["wire_pitch_mm"])
    half_wires = int(g["half_wires"])
    min_u_half = (half_wires + 0.30) * wire_pitch
    if float(g["u_half_span_mm"]) <= min_u_half:
        raise ValueError(
            "u_half_span_mm is too small for requested wire bank: "
            f"need > {min_u_half:.3f} mm"
        )


def generate_geo(cfg: Dict[str, Any]) -> str:
    g = cfg["geometry"]
    m = cfg["mesh"]

    gm = float(g["gap_minus_mm"])
    gp = float(g["gap_plus_mm"])
    pitch = float(g["wire_pitch_mm"])
    radius = 0.5e-3 * float(g["wire_diameter_um"])  # um -> mm
    half_wires = int(g["half_wires"])
    u_half = float(g["u_half_span_mm"])
    w_half = float(g["w_half_span_mm"])

    bulk = float(m["bulk_size_mm"])
    cath = float(m["cathode_size_mm"])
    cath_depth = float(m["cathode_refine_depth_mm"])
    wire_size = float(m["wire_size_mm"])
    wire_depth = float(m["wire_refine_radius_mm"])
    curvature = int(m["wire_curvature_points"])

    vmin = -gm
    vmax = gp
    umin, umax = -u_half, u_half
    wmin, wmax = -w_half, w_half

    # Extend cutter wires beyond the local map in w so the wire holes are open
    # through the artificial end boundaries; this avoids unphysical end caps.
    overhang = max(1.0, 2.0 * bulk)
    eps = max(1e-5, radius * 0.1)

    lines = [
        'SetFactory("OpenCASCADE");',
        "",
        "// Units: mm. Coordinate mapping: (x,y,z)=(u,v,w).",
        f"uMin = {_fmt(umin)};",
        f"uMax = {_fmt(umax)};",
        f"vMin = {_fmt(vmin)};",
        f"vMax = {_fmt(vmax)};",
        f"wMin = {_fmt(wmin)};",
        f"wMax = {_fmt(wmax)};",
        f"wireR = {_fmt(radius)};",
        f"wirePitch = {_fmt(pitch)};",
        f"eps = {_fmt(eps)};",
        "",
        "// Gas box.",
        "gasBox = newv;",
        "Box(gasBox) = {uMin, vMin, wMin, uMax-uMin, vMax-vMin, wMax-wMin};",
        "",
        "// Explicit anode wires.  They are removed from the gas volume and",
        "// become true zero-potential conductor boundaries in weighting solves.",
        "wireVols[] = {};",
    ]

    for i in range(-half_wires, half_wires + 1):
        x = i * pitch
        lines += [
            "wv = newv;",
            (
                "Cylinder(wv) = "
                f"{{{_fmt(x)}, 0, {_fmt(wmin-overhang)}, "
                f"0, 0, {_fmt((wmax-wmin)+2*overhang)}, wireR}};"
            ),
            "wireVols[] += {wv};",
        ]

    lines += [
        "",
        "gas[] = BooleanDifference{ Volume{gasBox}; Delete; }"
        "{ Volume{wireVols[]}; Delete; };",
        "",
        "// Recover boundary surfaces by geometric location.",
        (
            "cathMinus[] = Surface In BoundingBox "
            "{uMin-eps, vMin-eps, wMin-eps, uMax+eps, vMin+eps, wMax+eps};"
        ),
        (
            "cathPlus[] = Surface In BoundingBox "
            "{uMin-eps, vMax-eps, wMin-eps, uMax+eps, vMax+eps, wMax+eps};"
        ),
        (
            "sideUMinus[] = Surface In BoundingBox "
            "{uMin-eps, vMin-eps, wMin-eps, uMin+eps, vMax+eps, wMax+eps};"
        ),
        (
            "sideUPlus[] = Surface In BoundingBox "
            "{uMax-eps, vMin-eps, wMin-eps, uMax+eps, vMax+eps, wMax+eps};"
        ),
        (
            "sideWMinus[] = Surface In BoundingBox "
            "{uMin-eps, vMin-eps, wMin-eps, uMax+eps, vMax+eps, wMin+eps};"
        ),
        (
            "sideWPlus[] = Surface In BoundingBox "
            "{uMin-eps, vMin-eps, wMax-eps, uMax+eps, vMax+eps, wMax+eps};"
        ),
        "",
        "wireSurfs[] = {};",
    ]

    for i in range(-half_wires, half_wires + 1):
        x = i * pitch
        lines += [
            (
                "ws[] = Surface In BoundingBox "
                f"{{{_fmt(x-radius-eps)}, {_fmt(-radius-eps)}, wMin-eps, "
                f"{_fmt(x+radius+eps)}, {_fmt(radius+eps)}, wMax+eps}};"
            ),
            "wireSurfs[] += {ws[]};",
        ]

    lines += [
        "",
        "// Stable physical IDs.  ElmerGrid -autoclean renumbers surfaces",
        "// independently in ascending order, yielding boundaries 1..4 below.",
        'Physical Volume("gas", 1) = {gas[]};',
        'Physical Surface("cathode_minus", 101) = {cathMinus[]};',
        'Physical Surface("cathode_plus", 102) = {cathPlus[]};',
        (
            'Physical Surface("side_walls", 103) = '
            "{sideUMinus[], sideUPlus[], sideWMinus[], sideWPlus[]};"
        ),
        'Physical Surface("anode_wires", 104) = {wireSurfs[]};',
        "",
        "// Mesh controls.",
        f"Mesh.CharacteristicLengthMin = {_fmt(min(wire_size, cath))};",
        f"Mesh.CharacteristicLengthMax = {_fmt(bulk)};",
        f"Mesh.MeshSizeFromCurvature = {curvature};",
        "Mesh.MeshSizeExtendFromBoundary = 0;",
        "",
        "// Fine mesh near the explicit wire surfaces.",
        "Field[1] = Distance;",
        "Field[1].SurfacesList = {wireSurfs[]};",
        "Field[1].Sampling = 100;",
        "Field[2] = Threshold;",
        "Field[2].InField = 1;",
        f"Field[2].SizeMin = {_fmt(wire_size)};",
        f"Field[2].SizeMax = {_fmt(bulk)};",
        "Field[2].DistMin = 0;",
        f"Field[2].DistMax = {_fmt(wire_depth)};",
        "",
        "// Resolve the strip-mask discontinuity on both cathode planes.",
        "Field[3] = Distance;",
        "Field[3].SurfacesList = {cathMinus[], cathPlus[]};",
        "Field[3].Sampling = 100;",
        "Field[4] = Threshold;",
        "Field[4].InField = 3;",
        f"Field[4].SizeMin = {_fmt(cath)};",
        f"Field[4].SizeMax = {_fmt(bulk)};",
        "Field[4].DistMin = 0;",
        f"Field[4].DistMax = {_fmt(cath_depth)};",
        "",
        "Field[5] = Min;",
        "Field[5].FieldsList = {2, 4};",
        "Background Field = 5;",
        "",
        "Mesh.Algorithm3D = 1;",
        "Mesh.Optimize = 1;",
        "Mesh.SecondOrderLinear = 0;",
        "Mesh.MshFileVersion = 2.2;",
        "",
    ]
    return "\n".join(lines)


def solver_common(mesh_dir: str, output_file: str, eps_r: float) -> str:
    return f"""Check Keywords Warn

Header
  Mesh DB "." "{mesh_dir}"
End

Simulation
  Coordinate System = Cartesian 3D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output File = "{output_file}"
  Post File = "{Path(output_file).stem}.ep"
End

Constants
  Permittivity Of Vacuum = 8.8541878128e-12
End

Body 1
  Name = "Gas"
  Equation = 1
  Material = 1
End

Material 1
  Name = "Gas"
  Relative Permittivity = {_fmt(eps_r)}
End

Equation 1
  Active Solvers(1) = 1
  Calculate Electric Energy = True
End

Solver 1
  Equation = "Stat Elec Solver"
  Variable = Potential
  Variable DOFs = 1
  Procedure = "StatElecSolve" "StatElecSolver"
  Calculate Electric Field = False
  Calculate Electric Flux = False
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 3000
  Linear System Abort Not Converged = True
  Linear System Convergence Tolerance = 1.0e-11
  Linear System Preconditioning = ILU1
  Steady State Convergence Tolerance = 1.0e-8
End
"""


def generate_drift_sif(cfg: Dict[str, Any], mesh_dir: str) -> str:
    e = cfg["electrostatics"]
    common = solver_common(
        mesh_dir, "drift.result", float(e["gas_relative_permittivity"])
    )
    hv = float(e["anode_voltage_V"])
    return common + f"""
! Boundary numbering after ElmerGrid -autoclean:
!   1 cathode_minus
!   2 cathode_plus
!   3 side_walls (natural Neumann; intentionally no BC)
!   4 anode_wires
Boundary Condition 1
  Name = "CathodeMinus"
  Target Boundaries(1) = 1
  Potential = 0.0
End

Boundary Condition 2
  Name = "CathodePlus"
  Target Boundaries(1) = 2
  Potential = 0.0
End

Boundary Condition 3
  Name = "AnodeWires"
  Target Boundaries(1) = 4
  Potential = {_fmt(hv)}
End
"""


def strip_expression(cfg: Dict[str, Any], cathode: str, index: int) -> str:
    g = cfg["geometry"]
    pitch = float(g["strip_pitch_mm"])
    width = float(g["strip_width_mm"])
    alpha = math.atan(float(g["tan_alpha"]))
    ca = math.cos(alpha)
    sa = math.sin(alpha)
    center = index * pitch

    # Variable list is Coordinate 1 (=u), Coordinate 3 (=w), hence
    # tx(0)=u and tx(1)=w in MATC.
    if cathode == "minus":
        projection = f"(tx(1)*{_fmt(ca)} - tx(0)*{_fmt(sa)})"
    elif cathode == "plus":
        projection = f"(tx(1)*{_fmt(ca)} + tx(0)*{_fmt(sa)})"
    else:
        raise ValueError(cathode)

    half_width = 0.5 * width
    return (
        f"if (abs({projection} - ({_fmt(center)})) <= {_fmt(half_width)}) "
        "1.0; else 0.0;"
    )


def generate_weighting_sif(
    cfg: Dict[str, Any], mesh_dir: str, cathode: str, index: int
) -> str:
    e = cfg["electrostatics"]
    name = result_stem(cathode, index)
    common = solver_common(
        mesh_dir, f"{name}.result", float(e["gas_relative_permittivity"])
    )
    expr = strip_expression(cfg, cathode, index)

    variable_bc = (
        "  Potential = Variable Coordinate 1, Coordinate 3\n"
        f'    Real MATC "{expr}"'
    )
    if cathode == "minus":
        minus = variable_bc
        plus = "  Potential = 0.0"
    else:
        minus = "  Potential = 0.0"
        plus = variable_bc

    return common + f"""
! Exact Shockley-Ramo conductor boundary conditions:
! selected strip = 1 V, every other conductor = 0 V.
! The selected strip is imposed as a coordinate mask on the cathode plane.
! Side walls retain the natural zero-normal-flux condition of the finite patch.
Boundary Condition 1
  Name = "CathodeMinus"
  Target Boundaries(1) = 1
{minus}
End

Boundary Condition 2
  Name = "CathodePlus"
  Target Boundaries(1) = 2
{plus}
End

Boundary Condition 3
  Name = "AnodeWires"
  Target Boundaries(1) = 4
  Potential = 0.0
End
"""


def result_stem(cathode: str, index: int) -> str:
    if index < 0:
        suffix = f"m{-index}"
    elif index > 0:
        suffix = f"p{index}"
    else:
        suffix = "0"
    return f"weight_{cathode}_{suffix}"


def write_all(cfg: Dict[str, Any], out: Path) -> None:
    validate(cfg)
    out.mkdir(parents=True, exist_ok=True)

    mesh_dir = "chamber"

    (out / "chamber.geo").write_text(generate_geo(cfg), encoding="utf-8")
    (out / "drift.sif").write_text(
        generate_drift_sif(cfg, mesh_dir), encoding="utf-8"
    )

    weighting_files = []
    for cathode in cfg["weighting"]["cathodes"]:
        for index in cfg["weighting"]["strip_indices"]:
            index = int(index)
            stem = result_stem(cathode, index)
            filename = f"{stem}.sif"
            (out / filename).write_text(
                generate_weighting_sif(cfg, mesh_dir, cathode, index),
                encoding="utf-8",
            )
            weighting_files.append(
                {
                    "cathode": cathode,
                    "strip_index": index,
                    "label": f"{cathode}_strip_{index}",
                    "sif": filename,
                    "result": f"{mesh_dir}/{stem}.result",
                }
            )

    # ComponentElmer material map: one gas material/body after -autoclean.
    eps_r = float(cfg["electrostatics"]["gas_relative_permittivity"])
    (out / "dielectrics.dat").write_text(
        f"1\n1 {_fmt(eps_r)}\n", encoding="utf-8"
    )

    g = cfg["geometry"]
    alpha = math.atan(float(g["tan_alpha"]))
    manifest = {
        "coordinate_mapping": {
            "elmer_x": "local_u_mm",
            "elmer_y": "local_v_mm",
            "elmer_z": "local_w_mm",
            "garfield_mapping": "(x,y,z)=(u,v,w)",
        },
        "boundary_ids_after_elmergrid_autoclean": {
            "1": "cathode_minus",
            "2": "cathode_plus",
            "3": "side_walls",
            "4": "anode_wires",
        },
        "geometry": cfg["geometry"],
        "electrostatics": cfg["electrostatics"],
        "mesh": cfg["mesh"],
        "weighting": cfg["weighting"],
        "derived": {
            "alpha_deg": alpha * 180.0 / math.pi,
            "wire_radius_mm": 0.5e-3 * float(g["wire_diameter_um"]),
            "full_gap_mm": float(g["gap_minus_mm"]) + float(g["gap_plus_mm"]),
        },
        "drift_result": f"{mesh_dir}/drift.result",
        "weighting_maps": weighting_files,
        "notes": [
            "Side walls use natural zero-normal-flux BC and are numerical truncation boundaries.",
            "Anode wires are explicit grounded conductor boundaries in every weighting solve.",
            "For strip_width < strip_pitch, the unselected cathode area is currently treated as grounded; dielectric PCB gaps are not yet included.",
        ],
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    cfg = json.loads(args.config.read_text(encoding="utf-8"))
    write_all(cfg, args.output)

    print(f"Generated numerical MWPC weighting model in {args.output}")
    print("  chamber.geo")
    print("  drift.sif")
    for cathode in cfg["weighting"]["cathodes"]:
        for index in cfg["weighting"]["strip_indices"]:
            print(f"  {result_stem(cathode, int(index))}.sif")
    print("  dielectrics.dat")
    print("  manifest.json")


if __name__ == "__main__":
    main()

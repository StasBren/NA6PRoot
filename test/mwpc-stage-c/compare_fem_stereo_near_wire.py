#!/usr/bin/env python3
"""Compare independently sampled U/V Elmer 3D near-anode weighting fields.

Input 1: completed 9-wire candidate baseline map near-wire sample CSV.
Input 2: independently regenerated 9-wire map with ONLY wire-surface mesh
         size changed (0.012 -> 0.006 mm), same geometry and stereo signs.

This is a finite-element MESH DISCRETIZATION test, not an independently
validated physical strip field, and not a current waveform simulation.
The script never deletes or overwrites archived FEM or muon results.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

KEY = ("side", "wire_index", "r_over_R", "theta_index", "w_mm")
POS = ("theta_rad", "u_mm", "v_mm", "h_mm")
GRAD = ("Eu_per_mm", "Ev_per_mm", "Ew_per_mm")
HALF = ("Eu_halfstep_per_mm", "Ev_halfstep_per_mm", "Ew_halfstep_per_mm")


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read(path: Path, expected_tag: str) -> dict:
    need(path.is_file(), f"missing CSV {path}")
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    need(bool(rows), f"empty FEM map samples: {path}")
    need(all(x.get("mesh_tag") == expected_tag for x in rows),
         f"unexpected mesh tag in {path}; expected {expected_tag}")
    found = {}
    for row in rows:
        key = tuple(row[col] for col in KEY)
        need(key not in found, f"duplicate sample key {key}: {path}")
        for col in POS + GRAD + HALF + ("phi", "gradient_step_difference_per_mm"):
            v = float(row[col])
            need(math.isfinite(v), f"non-finite {col}, key {key}, file {path}")
        found[key] = row
    return found


def check_config(base_path: Path, fine_path: Path, study: str) -> None:
    with base_path.open(encoding="utf-8") as f:
        a = json.load(f)
    with fine_path.open(encoding="utf-8") as f:
        b = json.load(f)
    for k in ("geometry", "weighting", "electrostatics"):
        need(a[k] == b[k], f"FEM map {k} differs; not a controlled comparison")
    ca, cb = dict(a["mesh"]), dict(b["mesh"])
    wa, wb = ca.pop("wire_size_mm"), cb.pop("wire_size_mm")
    curved_a = ca.pop("curved_wire_surfaces", False)
    curved_b = cb.pop("curved_wire_surfaces", False)
    need(ca == cb, "uncontrolled FEM mesh changes")
    need(math.isclose(wa, .012, abs_tol=1e-12),
         f"expected original 0.012mm baseline, got {wa}")
    need(curved_a is False, "baseline must use the original straight-sided mesh")
    if study == "curvature":
        need(math.isclose(wb, .012, abs_tol=1e-12) and curved_b is True,
             "curvature study requires unchanged 0.012mm wire mesh and curved-wire P2 option")
    elif study == "wire_refine":
        need(math.isclose(wb, .006, abs_tol=1e-12) and curved_b is False,
             "original wire refinement requires 0.012->0.006mm with straight-sided P2")
    else:
        raise ValueError(f"unknown controlled FEM comparison: {study}")
    need(a["weighting"].get("stereo_sign_convention") == "plus_negative_u",
         "unexpected legacy stereo channel convention")
    need(a["geometry"]["half_wires"] == 4 and
         a["geometry"]["wire_pitch_mm"] == 3 and
         a["geometry"]["gap_plus_mm"] == 3 and
         a["geometry"]["gap_minus_mm"] == 3,
         "not the nine-wire candidate geometry")


def norm2(v: tuple[float, ...]) -> float:
    return sum(x*x for x in v)


def norm(v: tuple[float, ...]) -> float:
    return math.sqrt(norm2(v))


def report(rows: list[tuple[dict, dict]], group: str, side: str, shell: str) -> dict:
    sq_delta, sq_fine, sq_base = 0., 0., 0.
    largest_grad, largest_rel, largest_pot = 0., 0., 0.
    worst_step_b, worst_step_f = 0., 0.
    min_phi = float("inf")
    for a, b in rows:
        va = tuple(float(a[k]) for k in GRAD)
        vb = tuple(float(b[k]) for k in GRAD)
        delta = tuple(x-y for x,y in zip(va,vb))
        d = norm(delta)
        sq_delta += d*d
        sq_fine += norm2(vb)
        sq_base += norm2(va)
        largest_grad = max(largest_grad, d)
        largest_rel = max(largest_rel, d/max(norm(vb), 1.e-30))
        largest_pot = max(largest_pot, abs(float(a["phi"])-float(b["phi"])))
        min_phi = min(min_phi, abs(float(b["phi"])))
        worst_step_b = max(worst_step_b, float(a["gradient_step_difference_per_mm"]))
        worst_step_f = max(worst_step_f, float(b["gradient_step_difference_per_mm"]))
    need(sq_fine > 0 and sq_base > 0,
         "zero norms in entire grouping; map may be invalid")
    return {
        "group": group,
        "side": side,
        "r_over_R": shell,
        "n_points": len(rows),
        "vector_gradient_relative_L2": math.sqrt(sq_delta/sq_fine),
        "max_local_relative_gradient_difference": largest_rel,
        "max_abs_gradient_difference_per_mm": largest_grad,
        "max_abs_potential_difference": largest_pot,
        "minimum_abs_refined_potential": min_phi,
        "max_baseline_fd_h_vs_h2_per_mm": worst_step_b,
        "max_refined_fd_h_vs_h2_per_mm": worst_step_f,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("wire_refine", "curvature"), default="wire_refine",
                        help="Only one controlled parameter may differ from original baseline")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--refined", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--refined-config", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    a = parser.parse_args()
    need(not a.out_summary.exists(),
         f"refusing to overwrite existing summary: {a.out_summary}")
    check_config(a.baseline_config, a.refined_config, a.study)
    coarse = read(a.baseline, "baseline")
    fine = read(a.refined, "curved" if a.study == "curvature" else "wire_refine")
    need(coarse.keys() == fine.keys(),
         "different sampled coordinates; cannot compare FEM meshes")
    need(len(coarse) == 480,
         f"expected 480 U/V points, got {len(coarse)}")
    by_side = defaultdict(list)
    by_side_shell = defaultdict(list)
    combined = []
    for key in sorted(coarse):
        x, y = coarse[key], fine[key]
        for col in POS:
            need(math.isclose(float(x[col]), float(y[col]), rel_tol=0.,
                              abs_tol=1.e-12),
                 f"sample coordinate mismatch {col}, key={key}")
        combined.append((x,y))
        by_side[x["side"]].append((x,y))
        by_side_shell[(x["side"],x["r_over_R"])].append((x,y))
    output = [report(combined, "all", "both", "all")]
    for side in sorted(by_side):
        output.append(report(by_side[side], "cathode",side,"all"))
    for (side,shell),items in sorted(
          by_side_shell.items(),key=lambda kv:(kv[0][0],float(kv[0][1]))):
        output.append(report(items,"radial_shell",side,shell))
    a.out_summary.parent.mkdir(parents=True,exist_ok=True)
    with a.out_summary.open("w",newline="",encoding="utf-8") as f:
        writer = csv.DictWriter(f,fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print(f"Study={a.study}: Compared {len(combined)} identical 3D near-wire U/V FEM samples")
    for x in output[:3]:
        print(x["group"],x["side"],
            "relative_grad_L2=",f'{x["vector_gradient_relative_L2"]:.6%}',
            "max_local=",f'{x["max_local_relative_gradient_difference"]:.6%}',
            "max_abs_potential=",f'{x["max_abs_potential_difference"]:.8g}')
    if a.study == "curvature":
        print("GEOMETRY SENSITIVITY: original straight P2 versus curved P2 at same nominal h.")
        print("This is not an h-convergence estimate; a refined CURVED mesh is required separately.")
    if output[0]["vector_gradient_relative_L2"] > .02:
        print("NOT_CONVERGED under 2% internal numerical screening criterion.")
    else:
        print("BELOW 2% internal screen for this TWO-MESH change only.")
    print("NOT independent field accuracy; cylinder BC, sidewalls and real PCB pending.")
    print("Summary:",a.out_summary)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Identities of a geometry-study run vs its baseline: what changed (weighting) and what must not have (transport).

  geom_manifest.py BASE_RUN GEOM_RUN TAG --out IDENTITIES.json

Records and checks:
  geometry_id         sha256 of the canonical geometry json (baseline and alternative)
  weighting_id        LUT sha256 (Green axes) - must DIFFER;  recomputed on the files (baseline LUT intact)
  transport_id        Green-grid vertices/ions sha256 - must be IDENTICAL (the grid paths are reused)
  per muon            restrip transport hash (t, u, v columns) of the source == of the restriped export; the baseline
                      T1 export's current vertices sha256 == the one recorded when it was restriped (baseline intact)
  kernel_id           sha256 over K.npy / KPHI0.npy / kernel.json sha256 (as the response engine), both geometries
Exit status 1 if any must-hold identity fails.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import weighting_geometry as wg  # noqa: E402
from geom_restrip_export import transport_sha  # noqa: E402


def kernel_id(d):
    h = hashlib.sha256()
    for f in ("K.npy", "KPHI0.npy", "kernel.json"):
        h.update(wg.sha256(d / f).encode())
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", type=Path); ap.add_argument("grun", type=Path); ap.add_argument("tag")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    import numpy as np
    g1 = wg.Geometry.load(a.grun / "geometry.json"); g0 = wg.Geometry.load(Path(__file__).parent / "geometries" / "baseline.json")
    ax0 = json.loads((a.base / f"green_{a.tag}" / "axes.json").read_text())
    ax1 = json.loads((a.grun / f"green_{a.tag}" / "axes.json").read_text())
    checks = {}
    checks["transport_id_identical"] = (ax0["grid_vertices_sha256"] == ax1["grid_vertices_sha256"] and
                                        ax0["grid_ions_sha256"] == ax1["grid_ions_sha256"] and ax0["grid"] == ax1["grid"])
    checks["transport_grid_files_intact"] = wg.sha256(Path(ax0["grid"]) / "vertices.f64") == ax0["grid_vertices_sha256"]
    checks["weighting_id_differs"] = ax0["lut_sha256"] != ax1["lut_sha256"]
    checks["baseline_lut_intact"] = wg.sha256(ax0["lut"]) == ax0["lut_sha256"]
    checks["alt_lut_intact"] = wg.sha256(ax1["lut"]) == ax1["lut_sha256"]
    checks["green_axes_geometry"] = (ax1["stereo_tan"] == g1.stereo_tan and ax1["pitch_mm"] == g1.strip_pitch_mm and
                                     ax0["stereo_tan"] == g0.stereo_tan and ax0["pitch_mm"] == g0.strip_pitch_mm)
    muons = {}
    for R in sorted(a.grun.glob("muon_*")):
        rs = json.loads((R / "T1_export" / "restrip.json").read_text())
        src = a.base / R.name / "T1_export"
        nv = int(dict(l.split("=", 1) for l in (src / "complete.txt").read_text().split() if "=" in l)["vertices"])
        V0 = np.memmap(src / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
        V1 = np.memmap(R / "T1_export" / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
        t0, t1 = transport_sha(V0), transport_sha(V1)
        muons[R.name] = {"transport_sha256_tuv_base": t0, "transport_sha256_tuv_alt": t1,
                         "base_T1_vertices_sha256_now": wg.sha256(src / "vertices.f64"),
                         "base_T1_vertices_sha256_at_restrip": rs["source_vertices_sha256"],
                         "new_frame": rs["new_frame"], "frozen_frame": rs["frozen_frame"]}
        checks[f"{R.name}_transport_identical"] = t0 == t1 == rs["transport_sha256_tuv"]
        checks[f"{R.name}_baseline_export_intact"] = muons[R.name]["base_T1_vertices_sha256_now"] == rs["source_vertices_sha256"]
    kid = {}
    for name, d in (("base", a.base / f"kernel_{a.tag}"), ("alt", a.grun / f"kernel_{a.tag}")):
        if (d / "K.npy").exists():
            kid[name] = kernel_id(d)
    out = {"geometry": {"base": {"id": g0.geometry_id, **g0.physics()}, "alt": {"id": g1.geometry_id, **g1.physics()}},
           "weighting_id": {"base": ax0["lut_sha256"], "alt": ax1["lut_sha256"]},
           "transport_id": {"base": ax0["grid_vertices_sha256"], "alt": ax1["grid_vertices_sha256"]},
           "kernel_id": kid, "muons": muons, "checks": checks, "pass": all(checks.values())}
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    for k, v in checks.items():
        print(f"  {'ok  ' if v else 'FAIL'} {k}")
    print(("PASS" if out["pass"] else "FAIL") + f" identities: geometry {g0.geometry_id[:12]} -> {g1.geometry_id[:12]}, "
          f"weighting {ax0['lut_sha256'][:12]} -> {ax1['lut_sha256'][:12]}, transport {ax1['grid_vertices_sha256'][:12]} unchanged")
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

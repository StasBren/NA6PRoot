#!/usr/bin/env python3
"""Re-express an archived ion-path export in the strip frame of another readout geometry. Transport is NOT touched.

  geom_restrip_export.py SRC_EXPORT FRAME.json GEOMETRY.json OUT_EXPORT

SRC_EXPORT: an export of the frozen exporter (or the T1 time-collapsed copy) in the frozen strip frame
(stereo 0.1, pitch 1.7 mm). FRAME.json: <moments>.frame.json written by p3_birth_moments for the SAME muon with the
DEFAULT geometry (patch_shift_cm, u_in_mm, w_in_mm, base_U, base_V read from the archive by the frozen reader).

What changes: vertex columns 3, 4 (strip-relative longitudinal coordinate of the U, V plane, mm; constant along each
path because the field is 2D) and ions.csv upper_index / lower_index. With the exporter's definitions
  idx_U = llround((w_in - A u_in) / P),  wU = 10 (w + w_in - A patch - P idx_U)       (cm inside, mm out)
  idx_V = llround((w_in + A u_in) / P),  wV = 10 (w + w_in + A patch - P idx_V)
the new columns are wU' = wU + 10 [(A0 - A1) patch - (P1 idx_U' - P0 idx_U)], wV' = wV + 10 [(A1 - A0) patch -
(P1 idx_V' - P0 idx_V)]: a constant per plane, exact up to rounding (~1e-15 mm), no ion-level information needed.
What does not change (checked, recorded as transport hashes): vertex columns 0-2 (t, u, v), row order, weights,
birth/end times, vertex counts, the seal. A path that is physically the same is reused; nothing is re-drifted.
Refuses to overwrite; read-only on the source.
"""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import weighting_geometry as wg  # noqa: E402

A0, P0C = wg.FROZEN["stereo_tan"], .1 * wg.FROZEN["strip_pitch_mm"]


def llround(x):
    return int(math.copysign(math.floor(abs(x) + .5), x))


def bases(frame, A, Pc):
    uIn, wIn = .1 * frame["u_in_mm"], .1 * frame["w_in_mm"]
    return llround((wIn - A * uIn) / Pc), llround((wIn + A * uIn) / Pc)


def shifts(frame, geom):
    """(idx_U', idx_V', dU_mm, dV_mm) for the new geometry relative to the frozen frame."""
    A1, P1C = geom.stereo_tan, .1 * geom.strip_pitch_mm
    u0, v0 = bases(frame, A0, P0C)
    u1, v1 = bases(frame, A1, P1C)
    ps = frame["patch_shift_cm"]
    du = 10. * ((A0 - A1) * ps - (P1C * u1 - P0C * u0))
    dv = 10. * ((A1 - A0) * ps - (P1C * v1 - P0C * v0))
    return u1, v1, du, dv


def transport_sha(V, chunk=1 << 20):
    h = hashlib.sha256()
    for s in range(0, len(V), chunk):
        h.update(np.ascontiguousarray(V[s:s + chunk, :3]).tobytes())
    return h.hexdigest()


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    src, framep, geomp, out = (Path(x) for x in sys.argv[1:])
    geom = wg.Geometry.load(geomp)
    frame = json.loads(framep.read_text())
    if frame.get("stereo_tan") != A0 or frame.get("strip_pitch_mm") != wg.FROZEN["strip_pitch_mm"]:
        raise SystemExit("FRAME.json must come from p3_birth_moments with the default (frozen) geometry")
    if out.exists():
        raise SystemExit(f"{out} exists")
    seal = dict(l.split("=", 1) for l in (src / "complete.txt").read_text().splitlines() if "=" in l)
    nv = int(seal["vertices"])
    if (src / "vertices.f64").stat().st_size != 40 * nv:
        raise SystemExit("vertices.f64 size does not match the seal")
    rows = list(csv.DictReader(open(src / "ions.csv", newline="")))
    fields = list(rows[0].keys())
    u0, v0 = bases(frame, A0, P0C)
    if (u0, v0) != (frame["base_U"], frame["base_V"]):
        raise SystemExit("frame.json base indices are not reproduced by the frozen formula")
    if any((int(r["upper_index"]), int(r["lower_index"])) != (u0, v0) for r in rows):
        raise SystemExit("export strip indices differ from the archive frame: wrong FRAME.json for this export?")
    u1, v1, du, dv = shifts(frame, geom)
    out.mkdir(parents=True)
    V = np.memmap(src / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
    with open(out / "vertices.f64", "wb") as f:
        for s in range(0, nv, 1 << 20):
            blk = np.array(V[s:s + (1 << 20)], dtype="<f8")
            if du != 0.0:
                blk[:, 3] += du
            if dv != 0.0:
                blk[:, 4] += dv
            f.write(blk.tobytes())
    with open(out / "ions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            r = dict(r); r["upper_index"] = str(u1); r["lower_index"] = str(v1); w.writerow(r)
    (out / "complete.txt").write_text((src / "complete.txt").read_text())
    Vo = np.memmap(out / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
    ts_in, ts_out = transport_sha(V), transport_sha(Vo)
    if ts_in != ts_out:
        raise SystemExit("transport columns changed: refusing")
    meta = {"tool": "geom_restrip_export.py", "source": str(src), "frame": frame,
            "geometry": {"name": geom.name, **geom.physics()}, "geometry_id": geom.geometry_id,
            "frozen_frame": {"stereo_tan": A0, "strip_pitch_mm": wg.FROZEN["strip_pitch_mm"], "base_U": u0, "base_V": v0},
            "new_frame": {"base_U": u1, "base_V": v1, "delta_U_mm": du, "delta_V_mm": dv},
            "transport_sha256_tuv": ts_out, "source_vertices_sha256": wg.sha256(src / "vertices.f64"),
            "source_ions_sha256": wg.sha256(src / "ions.csv"), "vertices_sha256": wg.sha256(out / "vertices.f64"),
            "identity": bool(du == 0.0 and dv == 0.0 and (u1, v1) == (u0, v0))}
    (out / "restrip.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"RESTRIPED {out}: base U {u0}->{u1} V {v0}->{v1}, dU {du:.17g} mm, dV {dv:.17g} mm, "
          f"transport(t,u,v) sha256 {ts_out[:16]} unchanged")


if __name__ == "__main__":
    main()

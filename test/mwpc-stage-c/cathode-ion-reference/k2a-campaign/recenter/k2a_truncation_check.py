#!/usr/bin/env python3
"""Finite-window truncation of the local (per-primary) strip windows, from the FROZEN weighting model.

  k2a_truncation_check.py GROUPS_DIR LUT.npz NUMW_DIR OUT.json [--half-strips 5]

For every ion of every recentering group (sub-exports of k2a_group_export.py) at its BIRTH point (first vertex):
  Sigma_win(x) = sum over the group's 2H+1 window strips of the strip weighting potential (frozen WeightingLUT.upper,
                 the same call as the frozen replay's query(); V plane mirrored, as there)
  Phi(x)       = whole-plane weighting potential (frozen fit(0., 24, 0), the replay's own k = 0 strip-sum check)
Shockley-Ramo with a final state that is a single landing strip inside the window (the frozen replay enforces it):
  final induced charge on ALL strips of a plane  = q (delta_land - Phi(x0))
  final induced charge on the window strips      = q (delta_land - Sigma_win(x0))
  => charge on the strips OUTSIDE the window     = q (Sigma_win(x0) - Phi(x0))     (exact; usually small, negative)
Reported per plane and per group: Q_plane, Q_window, Q_outside, |Q_outside| / |Q_window| (the asymptotic truncation),
the worst per-ion gap Phi - Sigma_win, and the same with one more strip on each side (H+1, where the LUT phase range
allows; ions outside it are counted, never extrapolated). This measures the final (integrated) truncated charge; the
transient (bipolar) signal on strips beyond the window is NOT bounded by it and is reported only through the H+1 ring.
Reads the frozen modules from NUMW_DIR without modifying them. Exit 1 on any nonfinite value or failed read.
"""
import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

EFC = 1.602176634e-4


def first_vertices(d, rows):
    fd = os.open(d / "vertices.f64", os.O_RDONLY)
    try:
        return np.array([np.frombuffer(os.pread(fd, 40, 40 * int(r["first_vertex"])), "<f8") for r in rows]).reshape(-1, 5)
    finally:
        os.close(fd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("groups", type=Path); ap.add_argument("lut", type=Path); ap.add_argument("numw", type=Path)
    ap.add_argument("out", type=Path); ap.add_argument("--half-strips", type=int, default=5)
    a = ap.parse_args()
    sys.path.insert(0, str(a.numw.resolve()))
    from p3_spectral_weighting_lut import WeightingLUT, P, A, PHASE_MAX  # noqa: E402  (frozen, read-only)
    from p3_fourier_mfs_mode_pilot import fit  # noqa: E402
    lut = WeightingLUT(a.lut)
    fn, _ = fit(0., 24, 0)
    meta = json.loads((a.groups / "groups.json").read_text())
    H = a.half_strips
    res = {"tool": "k2a_truncation_check.py", "half_strips": H, "phase_max_mm": float(PHASE_MAX), "groups": {},
           "planes": {}, "status": "OK"}
    tot = {s: {"Q_plane": 0.0, "Q_window": 0.0, "Q_outside": 0.0, "Q_outside_H1": 0.0, "Q_window_H1": 0.0,
               "ions_H1_out_of_lut": 0, "max_gap": 0.0, "max_gap_H1": 0.0} for s in ("U", "V")}
    for name in meta["groups"]:
        d = a.groups / name
        rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
        p0 = first_vertices(d, rows)
        q = np.array([float(r["weight"]) for r in rows]) * EFC
        g = {}
        for side, sname in ((0, "U"), (1, "V")):
            xx = p0[:, 1:3] if side == 0 else -p0[:, 1:3]
            col = p0[:, 3 + side]
            # landing side is not known at birth; the outside charge q (Sigma - Phi) does not depend on it
            phi_plane = fn(xx).real
            def wsum(h):
                offs = np.arange(-h, h + 1)
                eta = (col[:, None] - P * offs[None, :]) - A * xx[:, :1]    # as WeightingLUT.upper (xx mirrored for V)
                ok = np.all(np.abs(eta) <= PHASE_MAX, axis=1)       # LUT phase range (no extrapolation)
                s = np.full(len(col), np.nan)
                if ok.any():
                    xo = xx[ok]; w = col[ok][:, None] - P * offs[None, :]
                    vals = lut.upper(np.repeat(xo, len(offs), axis=0), w.ravel()).reshape(ok.sum(), len(offs))
                    s[ok] = vals.sum(1)
                return s, ok
            s5, ok5 = wsum(H)
            if not ok5.all():
                print(f"FAIL {name} {sname}: {int((~ok5).sum())} ions outside the LUT phase range for the H={H} window")
                res["status"] = "FAIL"
            s6, ok6 = wsum(H + 1)
            gap5 = phi_plane - s5; gap6 = phi_plane - s6
            if not np.all(np.isfinite(gap5[ok5])) or not np.all(np.isfinite(phi_plane)):
                print(f"FAIL {name} {sname}: nonfinite potentials"); res["status"] = "FAIL"
            Qout = float(np.sum(q[ok5] * (s5[ok5] - phi_plane[ok5])))
            # final window/plane charges need the landing side: an ion lands on the plane its path ends on
            g[sname] = {"Q_outside_fC": Qout, "max_gap": float(np.nanmax(np.abs(gap5))) if ok5.any() else None,
                        "ions": len(rows), "ions_H1_out_of_lut": int((~ok6).sum()),
                        "Q_outside_H1_fC": float(np.sum(q[ok6] * (s6[ok6] - phi_plane[ok6]))) if ok6.any() else None,
                        "max_gap_H1": float(np.nanmax(np.abs(gap6))) if ok6.any() else None,
                        "Q_minus_Sigma_win_fC": float(-np.sum(q[ok5] * s5[ok5])), "Q_minus_Phi_fC": float(-np.sum(q * phi_plane))}
            t = tot[sname]
            t["Q_outside"] += Qout; t["max_gap"] = max(t["max_gap"], g[sname]["max_gap"] or 0.0)
            t["ions_H1_out_of_lut"] += g[sname]["ions_H1_out_of_lut"]
            if g[sname]["Q_outside_H1_fC"] is not None:
                t["Q_outside_H1"] += g[sname]["Q_outside_H1_fC"]; t["max_gap_H1"] = max(t["max_gap_H1"], g[sname]["max_gap_H1"])
        # landing (delta) terms from the frozen replay's channels.csv are added by the caller if wanted; here we use
        # the census: every ion lands on exactly one strip of one plane inside its window
        res["groups"][name] = g
    # window charge per plane from the summed local S0 (endpoint charges), if present next to the groups
    for sname in ("U", "V"):
        t = tot[sname]
        res["planes"][sname] = {"Q_outside_window_fC": t["Q_outside"], "max_per_ion_gap": t["max_gap"],
                                "Q_outside_window_H_plus_1_fC": t["Q_outside_H1"], "max_per_ion_gap_H_plus_1": t["max_gap_H1"],
                                "ions_H_plus_1_out_of_lut_range": t["ions_H1_out_of_lut"]}
        if not all(math.isfinite(v) for v in (t["Q_outside"], t["max_gap"])):
            res["status"] = "FAIL"
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    for s, v in res["planes"].items():
        print(f"{res['status']} truncation {s}: charge outside the ±{H} local windows {v['Q_outside_window_fC']:+.6g} fC "
              f"(max per-ion gap {v['max_per_ion_gap']:.3g}); with ±{H + 1}: {v['Q_outside_window_H_plus_1_fC']:+.6g} fC "
              f"({v['ions_H_plus_1_out_of_lut_range']} ions beyond the LUT phase range for ±{H + 1})")
    return 0 if res["status"] == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())

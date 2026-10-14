#!/usr/bin/env python3
"""
Reconstruct 6-mm MWPC gas entry/exit segments from the historical MNP33-bore
station-plane crossings.

The LocalTrajectoryAudit study stored track states exactly at the nominal
station reference planes z = 300, 340, 550, 590, 810, 850 cm.

For the intrinsic response study below we place a nominal 6-mm gas layer around
that reference plane and extrapolate the saved GEANT state across it.  This is
NOT yet the exact per-chamber installed geometry: the current station layout
staggers individual chamber gas volumes in z.  Exact event-by-event chamber
entry/exit and gas Edep require replay through the current chamber geometry.

Stage-C Garfield coordinates:
  u = detector/global Y (across wires)
  v = chamber normal; wire plane at v = 0
  w = -detector/global X (along wires, installed +90 deg convention)

For the current 2+4 mm model:
  gas: v in [-2,+4] mm
  gas centre: v = +1 mm
Therefore the physical gas faces are +/-3 mm around the station reference
plane.  A straight segment is extrapolated from the saved GEANT state to those
two faces.

No GEANT gas energy deposition was stored in the historical LocalTrajectoryAudit
files.  edep_eV is therefore written as NaN.  This file is intended for the
geometric Poisson/binomial wire-sharing study; event-by-event gas Edep requires
a replay through the current MWPC gas geometry.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


STATION_Z_CM = {
    0: 300.0,
    1: 340.0,
    2: 550.0,
    3: 590.0,
    4: 810.0,
    5: 850.0,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracks", required=True, type=Path,
                    help="bore_muon_crossings.csv or c_pi_track_sample.csv")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--gap-minus-mm", type=float, default=2.0)
    ap.add_argument("--gap-plus-mm", type=float, default=4.0)
    ap.add_argument("--max-tracks", type=int, default=0)
    args = ap.parse_args()

    if args.gap_minus_mm <= 0 or args.gap_plus_mm <= 0:
        raise SystemExit("Gas gaps must be positive.")

    with args.tracks.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("No track rows found.")

    required = {
        "event", "track", "station",
        "x_cm", "y_cm", "z_cm",
        "px_GeV", "py_GeV", "pz_GeV",
    }
    missing = required - set(rows[0])
    if missing:
        raise SystemExit("Missing columns: " + ", ".join(sorted(missing)))

    if args.max_tracks > 0:
        rows = rows[:args.max_tracks]

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "event", "track_id", "det_id", "station",
        "channel", "job", "pdg",
        "x_in_cm", "y_in_cm", "z_in_cm",
        "x_out_cm", "y_out_cm", "z_out_cm",
        "u_in_mm", "v_in_mm", "w_in_mm",
        "u_out_mm", "v_out_mm", "w_out_mm",
        "px_in_GeV", "py_in_GeV", "pz_in_GeV",
        "pu_in_GeV", "pv_in_GeV", "pw_in_GeV",
        "path_length_mm", "edep_GeV", "edep_eV", "time_ns",
    ]

    max_plane_dz = 0.0
    n_backward = 0

    with args.output.open("w", newline="", encoding="utf-8") as f:
        wri = csv.DictWriter(f, fieldnames=fields)
        wri.writeheader()

        for r in rows:
            station = int(float(r["station"]))
            if station not in STATION_Z_CM:
                raise SystemExit(f"Unexpected station {station}")

            z0 = float(r["z_cm"])
            max_plane_dz = max(max_plane_dz, abs(z0 - STATION_Z_CM[station]))

            x0 = float(r["x_cm"])
            y0 = float(r["y_cm"])
            px = float(r["px_GeV"])
            py = float(r["py_GeV"])
            pz = float(r["pz_GeV"])
            if pz == 0:
                continue

            # For this intrinsic nominal-segment study we centre a 6-mm gas
            # layer on the historical station reference. In Garfield v the
            # wire plane is offset so the gas centre is +1 mm for -2/+4 mm.
            v_center = 0.5 * (args.gap_plus_mm - args.gap_minus_mm)
            v_a = -args.gap_minus_mm
            v_b = +args.gap_plus_mm

            # z displacement from gas centre in cm. Chamber normal is global Z.
            dz_a_cm = (v_a - v_center) / 10.0
            dz_b_cm = (v_b - v_center) / 10.0

            def point(dz_cm):
                scale = dz_cm / pz
                x = x0 + scale * px
                y = y0 + scale * py
                z = z0 + dz_cm
                # Installed chamber convention.
                u = 10.0 * y
                v = v_center + 10.0 * dz_cm
                ww = -10.0 * x
                return x, y, z, u, v, ww

            a = point(dz_a_cm)
            b = point(dz_b_cm)

            # Put entry first in the actual direction of motion.
            if pz > 0:
                pin, pout = a, b
                pu, pv, pw = py, pz, -px
            else:
                n_backward += 1
                pin, pout = b, a
                pu, pv, pw = py, pz, -px

            du = pout[3] - pin[3]
            dv = pout[4] - pin[4]
            dw = pout[5] - pin[5]
            length = math.sqrt(du * du + dv * dv + dw * dw)

            wri.writerow({
                "event": int(float(r["event"])),
                "track_id": int(float(r["track"])),
                "det_id": -1,
                "station": station,
                "channel": r.get("channel", ""),
                "job": r.get("job", ""),
                "pdg": r.get("pdg", ""),
                "x_in_cm": pin[0],
                "y_in_cm": pin[1],
                "z_in_cm": pin[2],
                "x_out_cm": pout[0],
                "y_out_cm": pout[1],
                "z_out_cm": pout[2],
                "u_in_mm": pin[3],
                "v_in_mm": pin[4],
                "w_in_mm": pin[5],
                "u_out_mm": pout[3],
                "v_out_mm": pout[4],
                "w_out_mm": pout[5],
                "px_in_GeV": px,
                "py_in_GeV": py,
                "pz_in_GeV": pz,
                "pu_in_GeV": pu,
                "pv_in_GeV": pv,
                "pw_in_GeV": pw,
                "path_length_mm": length,
                "edep_GeV": "nan",
                "edep_eV": "nan",
                "time_ns": "nan",
            })

    print("Historical GEANT station-plane -> nominal MWPC gas segments COMPLETE")
    print(f"  input tracks              : {len(rows)}")
    print(f"  backward pz tracks        : {n_backward}")
    print(f"  max |z-plane reference|   : {max_plane_dz:.6g} cm")
    print(f"  gas v range               : {-args.gap_minus_mm:g} .. {args.gap_plus_mm:g} mm")
    print(f"  output                    : {args.output}")
    print()
    print("NOTE: edep is NaN because the historical LocalTrajectoryAudit did not store gas Edep.")
    print("      These are nominal intrinsic 6-mm segments, not exact staggered chamber crossings.")
    print("      Wire-lattice phase uses global Y directly; a constant phase offset is still arbitrary.")


if __name__ == "__main__":
    main()

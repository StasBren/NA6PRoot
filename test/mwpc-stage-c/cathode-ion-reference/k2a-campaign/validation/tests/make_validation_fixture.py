#!/usr/bin/env python3
"""SYNTHETIC preflight tree for the validation-selector test (NOT NA60+ data).

  make_validation_fixture.py OUTDIR [--n 60]

OUTDIR/pf/inventory/slots.csv, OUTDIR/pf/domain/per_event.csv, OUTDIR/rec/per_event_recenter.csv,
OUTDIR/pf/moments_nonsealed/chunk_000/muon_NNNN.moments.csv (+ .event.json), OUTDIR/run/provenance_prepare_*.txt.
Slots 0000-0004 are 'prepared' (kernel muons), 0050-0059 sealed (no moments written: the selector must never need
them), a few previously_inspected / excluded / out-of-domain events, six stations, three channels, varied gain,
asymmetry, wire multiplicity and displacement.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

COLS = ["primary", "wire", "t_ns", "t_rms_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A",
        "min_clearance_um", "base_U", "base_V", "t_raw_rms_ns", "rows"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out", type=Path); ap.add_argument("--n", type=int, default=60)
    a = ap.parse_args()
    rng = np.random.default_rng(3)
    pf, rec, run = a.out / "pf", a.out / "rec", a.out / "run"
    md = pf / "moments_nonsealed" / "chunk_000"
    for d in (pf / "inventory", pf / "domain", rec, run, md):
        d.mkdir(parents=True)
    slots, dom, rr = [], [], []
    for i in range(a.n):
        s = f"{i:04d}"
        role = "sealed_validation" if 50 <= i < 60 else ("development" if i < 20 else ("previously_inspected" if i % 17 == 0 else "pool"))
        status = "EXCLUDED" if i == 33 else "USABLE"
        ch = ["jpsi", "phi", "omega"][i % 3]
        slots.append({"slot": s, "role": role, "status": status, "channel": ch})
        if i < 5:
            (run / f"provenance_prepare_{s}.txt").write_text("lut_sha256=x\n")
        if role == "sealed_validation":
            continue
        npri = int(rng.integers(3, 30)); rows = []
        disp = 0 if i % 4 else int(rng.integers(2, 6)) * (1 if i % 8 else -1)
        for p in range(npri):
            nw = 1 + int(i % 5 == 0) + int(i % 7 == 0)
            for w in range(nw):
                N = float(rng.gamma(2, 4 + 10 * (i % 9 == 0)))
                rows.append({"primary": p, "wire": w - nw // 2, "t_ns": float(rng.uniform(5, 300)), "t_rms_ns": 0.0, "N_fC": N,
                             "r_eff_um": float(rng.uniform(30, 800)), "D_U_mm": float(1.7 * disp * p / npri + rng.normal(0, .3)),
                             "D_V_mm": float(rng.normal(0, .3)), "C1": 0.1, "S1": 0.2, "C2": 0.0, "S2": 0.0,
                             "A": float(np.clip(rng.normal(0.1 * (i % 6), 0.3), -1, 1)), "min_clearance_um": 2.0,
                             "base_U": 10, "base_V": -3, "t_raw_rms_ns": 1.0, "rows": 100})
        with open(md / f"muon_{s}.moments.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)
        q = sum(r["N_fC"] for r in rows)
        wires = sorted({r["wire"] for r in rows})
        dom.append({"event": f"muon_{s}", "status": "IN_DOMAIN", "channel": ch, "station": i % 6,
                    "theta_u_deg": float(rng.normal(0, 8)), "theta_w_deg": float(rng.normal(0, 8)),
                    "n_primaries": npri + 2, "ion_births_exact": int(q / 1.602e-4), "rows": 256 * npri,
                    "clouds": len(rows), "n_wires_5pct": len(wires), "dominant_wire_share": 1.0 / len(wires),
                    "patch_edge_wire_charge_frac": 0.02 if i == 41 else 0.0})
        ood = i in (23, 47)
        rr.append({"event": f"muon_{s}", "status_after": "OUT_OF_DOMAIN" if ood else "IN_DOMAIN",
                   "causes_after": "R_EFF_HIGH: 1 clouds, 2% charge" if ood else "",
                   "charge_local_strip_outside_base_window": 0.3 if abs(disp) >= 4 else 0.0})
    for name, rows in ((pf / "inventory" / "slots.csv", slots), (pf / "domain" / "per_event.csv", dom), (rec / "per_event_recenter.csv", rr)):
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"SYNTHETIC preflight fixture in {a.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Read-only wire-domain probe of selected events (issue 1): where their charge sits in the nine-wire patch frame.

  k2a_wire_probe.py NEW_DIR --pf PF --kernel KERNEL_DIR --slots 0664,0599 [--scan DOMAIN_SCAN_DIR]

For each slot (never sealed; preflight moments PF/moments_nonsealed/**/muon_<slot>.moments.csv and its .event.json /
.frame.json written by p3_birth_moments):
  frame    patch_shift_cm (chamber u of the local patch origin), u_in_mm / w_in_mm (muon entry), theta_u / theta_w,
           path length, primaries, ion births; the entry offset from the patch origin u_in - 10 patch_shift and the
           expected u extent of the track across the 6 mm gap (6 tan theta_u)
  wires    per patch wire k = rint(u / 3 mm): clouds, charge and charge fraction, r_eff range (charge-weighted
           quantiles), D_U / D_V range, time range, primaries; in kernel? (kernel wires from kernel.json); |k| = 4 is
           the outermost wire of the nine-wire analytic patch (field-model edge)
  primaries spanning several wires; the charge on wires outside the kernel, and its split by |k|
  scan     (optional) the domain-scan record of the slot (status, blocked fraction) for cross-checking; a missing
           DIR/scan.csv is reported and skipped, never required
Writes NEW_DIR/wire_probe.json and wire_probe.md. Coordinate conventions (MomentAccumulator.hh / BirthMoments.hh):
u_mm = 10 u_patch_cm, k = rint(u_mm / 3) (round half to even), r and theta about wire k, D_U / D_V strip coordinates
WITHOUT the -/+ tan(alpha) u stereo term (the kernel adds it per wire), base_U / base_V from the muon entry.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def wq(x, w, qs):
    o = np.argsort(x); c = np.cumsum(w[o]) / w.sum()
    return [float(x[o][min(np.searchsorted(c, q), len(x) - 1)]) for q in qs]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--pf", type=Path, required=True); ap.add_argument("--kernel", type=Path, required=True)
    ap.add_argument("--slots", required=True); ap.add_argument("--scan", type=Path)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    slots = {r["slot"]: r for r in csv.DictReader(open(a.pf / "inventory" / "slots.csv"))}
    kw = json.loads((a.kernel / "kernel.json").read_text())["green_axes"]["wires"]
    scan, scan_note = {}, "no --scan given"
    if a.scan:
        if (a.scan / "scan.csv").is_file():
            scan = {r["slot"]: r for r in csv.DictReader(open(a.scan / "scan.csv"))}; scan_note = str(a.scan / "scan.csv")
        else:
            scan_note = f"{a.scan / 'scan.csv'} absent: scan cross-check skipped"
            print(f"note: {scan_note}")
    res = {"tool": "k2a_wire_probe.py", "kernel_wires": kw, "scan": scan_note, "events": {}}
    for s in [f"{int(x):04d}" for x in a.slots.split(",") if x.strip()]:
        if "sealed" in slots.get(s, {}).get("role", "sealed"):
            raise SystemExit(f"refusing: slot {s} is sealed or unknown")
        f = sorted((a.pf / "moments_nonsealed").rglob(f"muon_{s}.moments.csv"))
        if not f:
            res["events"][s] = {"error": "no preflight moments"}; continue
        f = f[0]
        ev = json.loads(Path(str(f) + ".event.json").read_text()) if Path(str(f) + ".event.json").is_file() else {}
        fr = json.loads(Path(str(f) + ".frame.json").read_text()) if Path(str(f) + ".frame.json").is_file() else {}
        rows = list(csv.DictReader(open(f)))
        X = {c: np.array([float(r[c]) for r in rows]) for c in ("primary", "wire", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "t_ns")}
        q = X["N_fC"]; Q = q.sum(); w = X["wire"].astype(int)
        per = {}
        for k in sorted(set(w)):
            m = w == k
            per[int(k)] = {"clouds": int(m.sum()), "charge_fC": float(q[m].sum()), "charge_frac": float(q[m].sum() / Q),
                           "in_kernel": int(k) in kw, "patch_edge": abs(int(k)) == 4,
                           "r_eff_um_q05_q50_q95": wq(X["r_eff_um"][m], q[m], (0.05, 0.5, 0.95)),
                           "D_U_mm": [float(X["D_U_mm"][m].min()), float(X["D_U_mm"][m].max())],
                           "D_V_mm": [float(X["D_V_mm"][m].min()), float(X["D_V_mm"][m].max())],
                           "t_ns": [float(X["t_ns"][m].min()), float(X["t_ns"][m].max())],
                           "primaries": int(len(set(X["primary"][m].astype(int))))}
        pw = defaultdict(set)
        for p_, k in zip(X["primary"].astype(int), w):
            pw[p_].add(int(k))
        multi = sum(1 for v in pw.values() if len(v) > 1)
        out_k = [k for k in per if not per[k]["in_kernel"]]
        tu = ev.get("theta_u_deg")
        frame = {k: ev.get(k, fr.get(k)) for k in ("patch_shift_cm", "u_in_mm", "w_in_mm")}
        rec = {"moments": str(f), "role": slots[s]["role"], "station": ev.get("station"), "channel": ev.get("source_channel"),
               "theta_u_deg": tu, "theta_w_deg": ev.get("theta_w_deg"), "path_length_mm": ev.get("path_length_mm"),
               "n_primaries": ev.get("n_primaries"), "ion_births_exact": ev.get("ion_births_exact"), "frame": frame,
               "entry_offset_from_patch_origin_mm": (frame["u_in_mm"] - 10 * frame["patch_shift_cm"])
               if frame.get("u_in_mm") is not None and frame.get("patch_shift_cm") is not None else None,
               "track_u_extent_across_6mm_gap_mm": 6 * math.tan(math.radians(abs(tu))) if tu is not None else None,
               "clouds": len(rows), "charge_fC": float(Q), "wires": per, "primaries_on_several_wires": multi,
               "charge_frac_outside_kernel_wires": float(sum(per[k]["charge_frac"] for k in out_k)),
               "outside_kernel_by_wire": {k: per[k]["charge_frac"] for k in out_k},
               "scan": {k: scan[s].get(k) for k in ("execution", "status_names", "blocked_frac", "blocked_causes")} if s in scan else None}
        res["events"][s] = rec
    a.out.mkdir(parents=True)
    (a.out / "wire_probe.json").write_text(json.dumps(res, indent=1) + "\n")
    L = ["# Wire-domain probe (read-only)", "", f"kernel wires {kw}; patch wires -4..4 (nine-wire analytic patch; |k| = 4 is its edge); scan: {scan_note}", ""]
    for s, r in res["events"].items():
        if "error" in r:
            L.append(f"## {s}: {r['error']}"); continue
        L += [f"## muon {s} ({r['role']}, station {r['station']}, {r['channel']})", "",
              f"theta_u {r['theta_u_deg']} deg, theta_w {r['theta_w_deg']} deg, path {r['path_length_mm']} mm; frame {r['frame']}; entry offset from patch "
              f"origin {r['entry_offset_from_patch_origin_mm']} mm; track u extent across the gap ~{r['track_u_extent_across_6mm_gap_mm']} mm; "
              f"{r['clouds']} clouds, {r['charge_fC']:.4g} fC; primaries on several wires {r['primaries_on_several_wires']}; "
              f"charge outside the kernel wires {100 * r['charge_frac_outside_kernel_wires']:.3g} % {r['outside_kernel_by_wire']}; scan {r['scan']}", "",
              "| wire | in kernel | patch edge | clouds | primaries | charge frac | r_eff q05/q50/q95 um | D_U mm | D_V mm | t ns |", "|---|---|---|---|---|---|---|---|---|---|"]
        for k, x in r["wires"].items():
            L.append(f"| {k} | {x['in_kernel']} | {x['patch_edge']} | {x['clouds']} | {x['primaries']} | {100 * x['charge_frac']:.3g} % | "
                     + "/".join(f"{v:.4g}" for v in x["r_eff_um_q05_q50_q95"]) + f" | {x['D_U_mm'][0]:.3g}..{x['D_U_mm'][1]:.3g} | "
                     f"{x['D_V_mm'][0]:.3g}..{x['D_V_mm'][1]:.3g} | {x['t_ns'][0]:.4g}..{x['t_ns'][1]:.4g} |")
        L.append("")
    (a.out / "wire_probe.md").write_text("\n".join(L) + "\n")
    print(f"wire probe -> {a.out / 'wire_probe.md'}")


if __name__ == "__main__":
    main()

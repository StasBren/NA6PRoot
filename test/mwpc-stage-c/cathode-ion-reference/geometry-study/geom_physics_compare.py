#!/usr/bin/env python3
"""Readout observables of the SAME muons (same archived trajectories, same electronics) under two strip geometries.

  geom_physics_compare.py --base BASE_RUN --alt GEOM_RUN --base-geom G0.json --alt-geom G1.json --out DIR
                          [--variants S0,K2Akern_TAG] [--tp 100,200] [--fracs 0.02,0.05,0.1,0.2,0.5] [--no-prompt]

For each muon present in both runs and each variant (S0 = exact frozen-replay reference; K2Akern_TAG = fast model):
shaped U/V waveforms with the t_shaped_compare shaper (PROVISIONAL CR-RC^n stand-in, unchanged), and per plane:
collected charge, U/V charge ratio, lead peak and time, threshold crossing time of the lead strip, multiplicity,
peak-weighted centroid. Thresholds are ABSOLUTE (fC-equivalent) and the same for both geometries: f x the lead peak
of the BASELINE S0 plane (same electronics, so a geometry change cannot move the threshold).
Coordinates: absolute strip indices are in the chamber frame, so with each geometry's own P and A
  w = P (cU + cV) / 2,   u_chamber = P (cV - cU) / (2 A)                         (mm)
and the charge-weighted ion-birth position from the export (w_phys = wU + 10 A patch + P idxU, u_ch = u + 10 patch)
is the per-muon truth. Writes per_muon.jsonl, summary.json, summary.md. Geometry effect = alt - base, reported for S0
and for the fast model; their difference is the fast model's error ON the geometry effect.
"""
import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import weighting_geometry as wg  # noqa: E402

spec = importlib.util.spec_from_file_location("tsc", HERE.parent / "timing-ablation" / "t_shaped_compare.py")
tsc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tsc)


def truth(export, geom, patch_shift_cm):
    seal = dict(l.split("=", 1) for l in (export / "complete.txt").read_text().splitlines() if "=" in l)
    V = np.memmap(export / "vertices.f64", dtype="<f8", mode="r", shape=(int(seal["vertices"]), 5))
    rows = list(csv.DictReader(open(export / "ions.csv", newline="")))
    first = np.array([int(r["first_vertex"]) for r in rows]); q = np.array([float(r["weight"]) for r in rows])
    iu = int(rows[0]["upper_index"])
    b = np.asarray(V[first])
    wphys = b[:, 3] + 10 * geom.stereo_tan * patch_shift_cm + geom.strip_pitch_mm * iu
    uch = b[:, 1] + 10 * patch_shift_cm
    return float(np.sum(q * wphys) / q.sum()), float(np.sum(q * uch) / q.sum())


def waves(d, prompt):
    lab, tc, cur = tsc.load_wave(d / "replay" / "waveforms.csv")
    if prompt:
        l2, tc2, pc = tsc.load_wave(d / "prompt" / "prompt_waveforms.csv")
        if l2 != lab or not np.array_equal(tc, tc2):
            raise SystemExit(f"{d}: prompt and replay grids differ")
        cur = cur + pc
    return lab, tc[1] - tc[0], cur


def observe(lab, dt_in, cur, h, dt, nb, geom, sgn, thr):
    Y = tsc.shape(cur[:nb], dt_in, h, dt)
    q = cur[:nb].sum(0) * dt_in
    o = {}
    for side in ("U", "V"):
        cols, idx = tsc.plane_cols(lab, side)
        p = tsc.plane_obs(Y[:, cols], idx, dt, sgn[side], thr[side], q[cols].sum())
        o[side] = {"charge_fC": p["charge"], "lead_index": p["lead_index"], "lead_peak": p["lead_peak"],
                   "lead_peak_time_ns": p["lead_peak_time"], "thr_time_ns": p["thr_time"], "mult": p["mult"],
                   "centroid_strip": p["centroid_strip"], "second_over_lead": p["second_over_lead"],
                   "amplitude": p["amplitude"]}
    P, A = geom.strip_pitch_mm, geom.stereo_tan
    o["w_mm"] = .5 * P * (o["U"]["centroid_strip"] + o["V"]["centroid_strip"])
    o["u_mm"] = P * (o["V"]["centroid_strip"] - o["U"]["centroid_strip"]) / (2 * A)
    o["UV_charge_ratio"] = o["U"]["charge_fC"] / o["V"]["charge_fC"]
    o["UV_amp_ratio"] = o["U"]["amplitude"] / o["V"]["amplitude"]
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True); ap.add_argument("--alt", type=Path, required=True)
    ap.add_argument("--base-geom", required=True); ap.add_argument("--alt-geom", required=True)
    ap.add_argument("--variants", default="S0"); ap.add_argument("--tp", default="100,200")
    ap.add_argument("--order", type=int, default=3); ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--window-ns", type=float, default=6000.0)
    ap.add_argument("--fracs", default="0.02,0.05,0.1,0.2,0.5"); ap.add_argument("--no-prompt", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    g0, g1 = wg.Geometry.load(a.base_geom), wg.Geometry.load(a.alt_geom)
    variants = a.variants.split(","); fracs = [float(x) for x in a.fracs.split(",")]
    if a.out.exists():
        raise SystemExit(f"{a.out} exists")
    muons = sorted(m.name for m in a.alt.glob("muon_*") if all((a.alt / m.name / v / "replay").is_dir() and
                   (a.base / m.name / v / "replay").is_dir() for v in variants))
    if not muons:
        raise SystemExit("no muon with all variants in both runs")
    a.out.mkdir(parents=True)
    rows, summ = [], {}
    for tp in [float(x) for x in a.tp.split(",")]:
        h, hlabel = tsc.impulse(tp, a.order, a.dt, max(12 * tp, 2000.0), None)
        for m in muons:
            rs = json.loads((a.alt / m / "T1_export" / "restrip.json").read_text())
            ps = rs["frame"]["patch_shift_cm"]
            t0 = truth(a.base / m / "T1_export", g0, ps); t1 = truth(a.alt / m / "T1_export", g1, ps)
            lab0, dt0, c0 = waves(a.base / m / "S0", not a.no_prompt)
            nb = int(round(a.window_ns / dt0))
            Y0 = tsc.shape(c0[:nb], dt0, h, a.dt)
            sgn, thr = {}, {}
            for side in ("U", "V"):
                cols, _ = tsc.plane_cols(lab0, side)
                ext = np.where(np.abs(Y0[:, cols].max(0)) >= np.abs(Y0[:, cols].min(0)), Y0[:, cols].max(0), Y0[:, cols].min(0))
                sgn[side] = 1.0 if ext[int(np.argmax(np.abs(ext)))] > 0 else -1.0
                lead = float(np.max(sgn[side] * Y0[:, cols]))
                thr[side] = {f"f{f:g}": f * lead for f in fracs}
            for v in variants:
                for gname, run, geom, tr in (("base", a.base, g0, t0), ("alt", a.alt, g1, t1)):
                    lab, dti, cur = waves(run / m / v, not a.no_prompt)
                    o = observe(lab, dti, cur, h, a.dt, nb, geom, sgn, thr)
                    o.update(muon=m, variant=v, geometry=gname, geometry_id=geom.geometry_id, Tp_ns=tp,
                             truth_w_mm=tr[0], truth_u_mm=tr[1], dw_mm=o["w_mm"] - tr[0], du_mm=o["u_mm"] - tr[1],
                             thresholds_fCeq=thr, shaper=hlabel)
                    rows.append(o)
            if abs(t0[0] - t1[0]) > 1e-9 or abs(t0[1] - t1[1]) > 1e-9:
                raise SystemExit(f"{m}: truth position differs between geometries ({t0} vs {t1}): frames inconsistent")
    with open(a.out / "per_muon.jsonl", "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    # geometry effect per muon: alt - base, for each variant; fast-model error on the effect
    keys = {"U charge [fC]": lambda o: o["U"]["charge_fC"], "V charge [fC]": lambda o: o["V"]["charge_fC"],
            "U/V charge ratio": lambda o: o["UV_charge_ratio"], "U lead peak": lambda o: o["U"]["lead_peak"],
            "V lead peak": lambda o: o["V"]["lead_peak"], "U second/lead": lambda o: o["U"]["second_over_lead"],
            "V second/lead": lambda o: o["V"]["second_over_lead"],
            "U t_thr f0.1 [ns]": lambda o: o["U"]["thr_time_ns"].get("f0.1", np.nan),
            "V t_thr f0.1 [ns]": lambda o: o["V"]["thr_time_ns"].get("f0.1", np.nan),
            "U mult f0.05": lambda o: o["U"]["mult"].get("f0.05", np.nan),
            "V mult f0.05": lambda o: o["V"]["mult"].get("f0.05", np.nan),
            "w - truth [mm]": lambda o: o["dw_mm"], "u - truth [mm]": lambda o: o["du_mm"]}
    L = ["# Strip-geometry effect on readout observables", "",
         f"base: {g0.name} (id {g0.geometry_id[:12]}), alt: {g1.name} (id {g1.geometry_id[:12]}); muons {', '.join(muons)}; "
         f"signal: {'ion term only' if a.no_prompt else 'ion term + instantaneous electron-step approximation'}; "
         "shaping PROVISIONAL; thresholds absolute, from the baseline S0 lead peak.", ""]
    for tp in sorted({r["Tp_ns"] for r in rows}):
        L += [f"## Tp = {tp:g} ns", "", "| muon | observable | " + " | ".join(f"{v} base | {v} alt" for v in variants)
              + (" | effect S0 | effect " + " | effect ".join(variants[1:]) + " | fast-model error on effect" if len(variants) > 1 else " | effect S0")
              + " |", "|---" * (2 + 2 * len(variants) + len(variants) + (1 if len(variants) > 1 else 0)) + "|"]
        for m in muons:
            sel = {(r["variant"], r["geometry"]): r for r in rows if r["muon"] == m and r["Tp_ns"] == tp}
            for k, fn in keys.items():
                vals = [fn(sel[(v, g)]) for v in variants for g in ("base", "alt")]
                eff = [fn(sel[(v, "alt")]) - fn(sel[(v, "base")]) for v in variants]
                cells = [f"{x:.5g}" for x in vals] + [f"{x:+.4g}" for x in eff]
                if len(variants) > 1:
                    cells.append(" / ".join(f"{e - eff[0]:+.3g}" for e in eff[1:]))
                L.append(f"| {m} | {k} | " + " | ".join(cells) + " |")
        L.append("")
    (a.out / "summary.md").write_text("\n".join(L) + "\n")
    (a.out / "summary.json").write_text(json.dumps({"base": g0.physics(), "alt": g1.physics(), "muons": muons,
                                                    "variants": variants}, indent=1) + "\n")
    print("\n".join(L[:40]))


if __name__ == "__main__":
    main()

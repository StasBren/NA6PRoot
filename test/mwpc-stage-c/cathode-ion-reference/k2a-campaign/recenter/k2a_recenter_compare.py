#!/usr/bin/env python3
"""Compare two replay-format variants on their COMMON absolute strips and time window; report what lies outside.

  k2a_recenter_compare.py REF_DIR CAND_DIR --out OUT.json [--window-ns 7000] [--tp 100,200] [--tol-l2 X] [--label TEXT]

REF_DIR / CAND_DIR: replay/waveforms.csv and prompt/prompt_waveforms.csv (labels U_<idx>, V_<idx>, absolute).
Per plane, ion-only and ion + electron step, over common labels and the common time prefix up to --window-ns:
  unshaped rel L2 of the current and of the cumulative charge, window charge rel diff, U/V charge ratio diff,
  lead strip (by window charge) identical?, max |per-strip window charge diff| / lead charge, shaped rel L2 at each Tp
  (t_shaped_compare CR-RC^3 stand-in, PROVISIONAL), and the window charge each side has on strips the other lacks
  (as a fraction of its plane charge) - e.g. signal outside the base-centred 11-strip output.
PASS/FAIL only if --tol-l2 is given (max over planes/terms of unshaped and shaped rel L2).
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location("tsc", HERE.parents[2] / "timing-ablation" / "t_shaped_compare.py")
tsc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tsc)


def load(d, prompt):
    lab, t, cur = tsc.load_wave(d / "replay" / "waveforms.csv")
    if prompt:
        l2, t2, pc = tsc.load_wave(d / "prompt" / "prompt_waveforms.csv")
        if l2 != lab or len(t2) != len(t) or not np.allclose(t2, t):
            raise SystemExit(f"{d}: prompt/replay grids differ")
        cur = cur + pc
    return lab, t, cur


def rel(a, b):
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref", type=Path); ap.add_argument("cand", type=Path); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--window-ns", type=float, default=7000.0); ap.add_argument("--tp", default="100,200")
    ap.add_argument("--tol-l2", type=float); ap.add_argument("--label", default="")
    a = ap.parse_args()
    res = {"ref": str(a.ref), "cand": str(a.cand), "label": a.label, "window_ns": a.window_ns, "terms": {}}
    worst = 0.0
    for term, prompt in (("ion_only", False), ("ion_plus_electron_step", True)):
        lr, tr, cr = load(a.ref, prompt); lc, tcand, cc = load(a.cand, prompt)
        dt = tr[1] - tr[0]
        if abs((tcand[1] - tcand[0]) - dt) > 1e-9 or abs(tcand[0] - tr[0]) > 1e-9:
            raise SystemExit("time grids differ")
        nb = min(len(tr), len(tcand), int(round(a.window_ns / dt)))
        T = {}
        for side in ("U", "V"):
            Lr = [l for l in lr if l[0] == side]; Lc = [l for l in lc if l[0] == side]
            common = [l for l in Lr if l in Lc]
            R = cr[:nb, [lr.index(l) for l in common]]; C = cc[:nb, [lc.index(l) for l in common]]
            qr = R.sum(0) * dt; qc = C.sum(0) * dt
            qr_all = cr[:nb, [lr.index(l) for l in Lr]].sum(0) * dt; qc_all = cc[:nb, [lc.index(l) for l in Lc]].sum(0) * dt
            only_r = [l for l in Lr if l not in Lc]; only_c = [l for l in Lc if l not in Lr]
            o = {"common_strips": len(common), "ref_only_strips": only_r, "cand_only_strips": only_c,
                 "unshaped_L2_rel": rel(C, R), "cumulative_L2_rel": rel(np.cumsum(C, 0), np.cumsum(R, 0)),
                 "window_charge_ref_fC": float(qr.sum()), "window_charge_cand_fC": float(qc.sum()),
                 "window_charge_rel_diff": float((qc.sum() - qr.sum()) / qr.sum()) if qr.sum() else float("nan"),
                 "lead_ref": common[int(np.argmax(np.abs(qr)))], "lead_cand": common[int(np.argmax(np.abs(qc)))],
                 "max_strip_charge_diff_over_lead": float(np.max(np.abs(qc - qr)) / np.max(np.abs(qr))),
                 "ref_charge_on_ref_only_strips_frac": float(qr_all[[Lr.index(l) for l in only_r]].sum() / qr_all.sum()) if only_r else 0.0,
                 "cand_charge_on_cand_only_strips_frac": float(qc_all[[Lc.index(l) for l in only_c]].sum() / qc_all.sum()) if only_c else 0.0}
            o["lead_identical"] = o["lead_ref"] == o["lead_cand"]
            for tp in (float(x) for x in a.tp.split(",")):
                h, _ = tsc.impulse(tp, 3, 1.0, max(12 * tp, 2000.0), None)
                o[f"shaped_L2_rel_Tp{tp:g}"] = rel(tsc.shape(C, dt, h, 1.0), tsc.shape(R, dt, h, 1.0))
            worst = max(worst, o["unshaped_L2_rel"], *(o[k] for k in o if k.startswith("shaped_L2_rel")))
            T[side] = o
        T["UV_charge_ratio_ref"] = T["U"]["window_charge_ref_fC"] / T["V"]["window_charge_ref_fC"]
        T["UV_charge_ratio_cand"] = T["U"]["window_charge_cand_fC"] / T["V"]["window_charge_cand_fC"]
        T["UV_charge_ratio_rel_diff"] = T["UV_charge_ratio_cand"] / T["UV_charge_ratio_ref"] - 1
        res["terms"][term] = T
    res["worst_L2_rel"] = worst
    if a.tol_l2 is not None:
        res["tol_l2"] = a.tol_l2; res["pass"] = worst <= a.tol_l2
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    for term, T in res["terms"].items():
        for s in ("U", "V"):
            o = T[s]
            print(f"{a.label} {term:22s} {s}: common {o['common_strips']:2d} L2 {o['unshaped_L2_rel']:.3e} "
                  f"shaped100 {o['shaped_L2_rel_Tp100']:.3e} Q rel {o['window_charge_rel_diff']:+.3e} lead "
                  f"{'same' if o['lead_identical'] else o['lead_ref'] + '->' + o['lead_cand']} "
                  f"ref-only {o['ref_charge_on_ref_only_strips_frac']:.3g} cand-only {o['cand_charge_on_cand_only_strips_frac']:.3g}")
        print(f"{a.label} {term:22s} U/V ratio rel diff {T['UV_charge_ratio_rel_diff']:+.3e}")
    if "pass" in res:
        print(("PASS" if res["pass"] else "FAIL") + f" {a.label}: worst rel L2 {worst:.3e} (tol {a.tol_l2:g})")
    return 0 if res.get("pass", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())

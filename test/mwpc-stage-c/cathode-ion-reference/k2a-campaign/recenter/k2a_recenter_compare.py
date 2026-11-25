#!/usr/bin/env python3
"""Compare two replay-format variants on common absolute strips; explicit status, never a silent pass.

  k2a_recenter_compare.py REF_DIR CAND_DIR --out OUT.json --kind KIND [--shift M] [--thresholds thresholds.json]
                          [--window-ns 7000] [--tp 100,200] [--label TEXT]

KIND
  e1          identity under an integer strip relabelling by --shift M: PASS / FAIL / INCONCLUSIVE
              (every metric <= the e1 limits, exact expected overlap 2H+1-|m| strips per plane, same lead strip)
  identity    e2 exact identity (e.g. K2A_local vs K2A base when every primary has m = 0): PASS / FAIL
  k2a_vs_s0   e2 diagnostic of the fast model against an exact reference: OK / WARN / FAIL / INCONCLUSIVE
  window      e2 diagnostic of the per-primary window choice (S0_local vs frozen base-centred S0): OK / WARN / ...
FAIL is reserved for numerical breakage (nonfinite values, identity violated, grids inconsistent) - or, for e1, any
deviation; WARN marks a documented threshold exceeded (thresholds.json); INCONCLUSIVE means there is nothing valid to
compare (missing files, zero/too few common strips, zero reference charge). Each check is listed with value, limit and
status. Metrics per plane (U, V) and term (ion only, ion + electron step), on common labels and the common time prefix:
unshaped rel L2 of the current and of the cumulative charge, shaped rel L2 per Tp (t_shaped_compare CR-RC^3 stand-in,
PROVISIONAL), window charge rel diff, U/V charge-ratio rel diff, lead strip, max per-strip charge diff / lead charge,
and the charge each side has on strips the other lacks (missing strip tails).
Exit: 0 PASS/OK, 1 FAIL, 2 WARN, 3 INCONCLUSIVE.
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location("tsc", HERE.parents[2] / "timing-ablation" / "t_shaped_compare.py")
tsc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tsc)
CODES = {"PASS": 0, "OK": 0, "FAIL": 1, "WARN": 2, "INCONCLUSIVE": 3}


class Inconclusive(Exception):
    pass


def load(d, prompt):
    f = d / "replay" / "waveforms.csv"
    if not f.is_file():
        raise Inconclusive(f"missing {f}")
    lab, t, cur = tsc.load_wave(f)
    if prompt:
        g = d / "prompt" / "prompt_waveforms.csv"
        if not g.is_file():
            raise Inconclusive(f"missing {g}")
        l2, t2, pc = tsc.load_wave(g)
        if l2 != lab or len(t2) != len(t) or not np.allclose(t2, t):
            raise ValueError(f"{d}: prompt/replay grids differ")
        cur = cur + pc
    return lab, t, cur


def rel(a, b):
    nb = float(np.linalg.norm(b))
    return float(np.linalg.norm(a - b)) / nb if nb > 0 else float("nan")


def metrics(lr, tr, cr, lc, tcand, cc, window_ns, tps):
    dt = tr[1] - tr[0]
    if abs((tcand[1] - tcand[0]) - dt) > 1e-9 or abs(tcand[0] - tr[0]) > 1e-9:
        raise ValueError("time grids differ")
    nb = min(len(tr), len(tcand), int(round(window_ns / dt)))
    out = {}
    for side in ("U", "V"):
        Lr = [l for l in lr if l[0] == side]; Lc = [l for l in lc if l[0] == side]
        common = [l for l in Lr if l in Lc]
        o = {"common_strips": len(common), "ref_strips": len(Lr), "cand_strips": len(Lc),
             "ref_only_strips": [l for l in Lr if l not in Lc], "cand_only_strips": [l for l in Lc if l not in Lr]}
        qr_all = cr[:nb, [lr.index(l) for l in Lr]].sum(0) * dt
        qc_all = cc[:nb, [lc.index(l) for l in Lc]].sum(0) * dt
        o["ref_plane_charge_fC"] = float(qr_all.sum()); o["cand_plane_charge_fC"] = float(qc_all.sum())
        o["ref_charge_on_ref_only_strips_frac"] = (float(qr_all[[Lr.index(l) for l in o["ref_only_strips"]]].sum() / qr_all.sum())
                                                   if o["ref_only_strips"] and qr_all.sum() else 0.0)
        o["cand_charge_on_cand_only_strips_frac"] = (float(qc_all[[Lc.index(l) for l in o["cand_only_strips"]]].sum() / qc_all.sum())
                                                     if o["cand_only_strips"] and qc_all.sum() else 0.0)
        if common:
            R = cr[:nb, [lr.index(l) for l in common]]; C = cc[:nb, [lc.index(l) for l in common]]
            qr, qc = R.sum(0) * dt, C.sum(0) * dt
            o["unshaped_L2_rel"] = rel(C, R)
            o["cumulative_L2_rel"] = rel(np.cumsum(C, 0), np.cumsum(R, 0))
            o["window_charge_ref_fC"] = float(qr.sum()); o["window_charge_cand_fC"] = float(qc.sum())
            o["window_charge_rel_diff"] = float((qc.sum() - qr.sum()) / qr.sum()) if qr.sum() else float("nan")
            o["lead_ref"] = common[int(np.argmax(np.abs(qr)))]; o["lead_cand"] = common[int(np.argmax(np.abs(qc)))]
            o["lead_identical"] = o["lead_ref"] == o["lead_cand"]
            mq = float(np.max(np.abs(qr)))
            o["max_strip_charge_diff_over_lead"] = float(np.max(np.abs(qc - qr)) / mq) if mq > 0 else float("nan")
            for tp in tps:
                h, _ = tsc.impulse(tp, 3, 1.0, max(12 * tp, 2000.0), None)
                o[f"shaped_L2_rel_Tp{tp:g}"] = rel(tsc.shape(C, dt, h, 1.0), tsc.shape(R, dt, h, 1.0))
        out[side] = o
    if out["U"].get("window_charge_ref_fC") and out["V"].get("window_charge_ref_fC"):
        a = out["U"]["window_charge_ref_fC"] / out["V"]["window_charge_ref_fC"]
        b = out["U"]["window_charge_cand_fC"] / out["V"]["window_charge_cand_fC"]
        out["UV_charge_ratio_ref"], out["UV_charge_ratio_cand"], out["UV_charge_ratio_rel_diff"] = a, b, b / a - 1
    return out


def check(checks, name, value, limit, bad_status, cmp="le"):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        checks.append({"check": name, "value": None if value is None else str(value), "limit": limit, "status": "FAIL",
                       "why": "nonfinite or missing value"})
        return
    ok = (abs(value) <= limit) if cmp == "le" else (value >= limit) if cmp == "ge" else bool(value) == bool(limit)
    checks.append({"check": name, "value": value, "limit": limit, "status": "OK" if ok else bad_status})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref", type=Path); ap.add_argument("cand", type=Path); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--kind", choices=("e1", "identity", "k2a_vs_s0", "window"), required=True)
    ap.add_argument("--shift", type=int, default=0); ap.add_argument("--half-strips", type=int, default=5)
    ap.add_argument("--thresholds", type=Path, default=HERE.parent / "thresholds.json")
    ap.add_argument("--window-ns", type=float, default=7000.0); ap.add_argument("--tp", default="100,200")
    ap.add_argument("--label", default="")
    a = ap.parse_args()
    th = json.loads(a.thresholds.read_text())
    tps = [float(x) for x in a.tp.split(",")]
    res = {"ref": str(a.ref), "cand": str(a.cand), "kind": a.kind, "label": a.label, "shift": a.shift,
           "window_ns": a.window_ns, "thresholds_file": str(a.thresholds), "thresholds": th, "terms": {}, "checks": []}
    checks = res["checks"]
    try:
        for term, prompt in (("ion_only", False), ("ion_plus_electron_step", True)):
            lr, tr, cr = load(a.ref, prompt); lc, tcd, cc = load(a.cand, prompt)
            if not (np.all(np.isfinite(cr)) and np.all(np.isfinite(cc))):
                raise ValueError(f"nonfinite waveform values ({term})")
            M = metrics(lr, tr, cr, lc, tcd, cc, a.window_ns, tps); res["terms"][term] = M
            for s in ("U", "V"):
                o = M[s]; p = f"{term}.{s}."
                if a.kind == "e1":
                    need = 2 * a.half_strips + 1 - abs(a.shift)
                    check(checks, p + "common_strips", o["common_strips"], need, "INCONCLUSIVE", "ge")
                    if o["common_strips"] != need and o["common_strips"] >= need:
                        checks.append({"check": p + "common_strips_exact", "value": o["common_strips"], "limit": need, "status": "FAIL"})
                else:
                    check(checks, p + "common_strips", o["common_strips"], th["e2"]["min_common_strips"], "INCONCLUSIVE", "ge")
                if not o["common_strips"]:
                    continue
                if not o.get("window_charge_ref_fC"):
                    checks.append({"check": p + "reference_charge", "value": o.get("window_charge_ref_fC"), "limit": ">0",
                                   "status": "INCONCLUSIVE"})
                    continue
                L = [o["unshaped_L2_rel"], o["cumulative_L2_rel"]] + [o[f"shaped_L2_rel_Tp{tp:g}"] for tp in tps]
                if a.kind == "e1":
                    E = th["e1"]
                    for nm, v in zip(["unshaped_L2_rel", "cumulative_L2_rel"] + [f"shaped_L2_rel_Tp{tp:g}" for tp in tps], L):
                        check(checks, p + nm, v, E["rel_L2_max"], "FAIL")
                    check(checks, p + "window_charge_rel_diff", o["window_charge_rel_diff"], E["window_charge_rel_max"], "FAIL")
                    check(checks, p + "max_strip_charge_diff_over_lead", o["max_strip_charge_diff_over_lead"],
                          E["max_strip_charge_diff_over_lead_max"], "FAIL")
                    check(checks, p + "lead_identical", o["lead_identical"], True, "FAIL", "eq")
                elif a.kind == "identity":
                    for nm, v in zip(["unshaped_L2_rel", "cumulative_L2_rel"] + [f"shaped_L2_rel_Tp{tp:g}" for tp in tps], L):
                        check(checks, p + nm, v, th["e2"]["identity_rel_L2_max"], "FAIL")
                    check(checks, p + "strip_sets_identical", not (o["ref_only_strips"] or o["cand_only_strips"]), True, "FAIL", "eq")
                else:
                    W = th["e2"]["K2A_vs_S0" if a.kind == "k2a_vs_s0" else "S0_window_choice"]
                    check(checks, p + "unshaped_L2_rel", o["unshaped_L2_rel"], W["unshaped_L2_rel_warn"], "WARN")
                    for tp in tps:
                        check(checks, p + f"shaped_L2_rel_Tp{tp:g}", o[f"shaped_L2_rel_Tp{tp:g}"], W["shaped_L2_rel_warn"], "WARN")
                    check(checks, p + "window_charge_rel_diff", o["window_charge_rel_diff"], W["window_charge_rel_warn"], "WARN")
                    check(checks, p + "max_strip_charge_diff_over_lead", o["max_strip_charge_diff_over_lead"],
                          W["max_strip_charge_diff_over_lead_warn"], "WARN")
                    check(checks, p + "lead_identical", o["lead_identical"], True, "WARN", "eq")
                    check(checks, p + "ref_charge_on_ref_only_strips_frac", o["ref_charge_on_ref_only_strips_frac"],
                          W["missing_strip_charge_frac_warn"], "WARN")
                    check(checks, p + "cand_charge_on_cand_only_strips_frac", o["cand_charge_on_cand_only_strips_frac"],
                          W["missing_strip_charge_frac_warn"], "WARN")
            if a.kind in ("k2a_vs_s0", "window") and "UV_charge_ratio_rel_diff" in M:
                W = th["e2"]["K2A_vs_S0" if a.kind == "k2a_vs_s0" else "S0_window_choice"]
                check(checks, f"{term}.UV_charge_ratio_rel_diff", M["UV_charge_ratio_rel_diff"], W["UV_ratio_rel_warn"], "WARN")
            elif a.kind == "e1" and "UV_charge_ratio_rel_diff" in M:
                check(checks, f"{term}.UV_charge_ratio_rel_diff", M["UV_charge_ratio_rel_diff"], th["e1"]["window_charge_rel_max"], "FAIL")
        st = [c["status"] for c in checks]
        if not st:
            status = "INCONCLUSIVE"
        elif "FAIL" in st:
            status = "FAIL"
        elif "INCONCLUSIVE" in st:
            status = "INCONCLUSIVE"
        elif "WARN" in st:
            status = "WARN"
        else:
            status = "PASS" if a.kind in ("e1", "identity") else "OK"
    except Inconclusive as e:
        status = "INCONCLUSIVE"; res["error"] = str(e)
    except (ValueError, KeyError, IndexError) as e:
        status = "FAIL"; res["error"] = f"{type(e).__name__}: {e}"
    res["status"] = status
    finite = [v for T in res["terms"].values() for s in ("U", "V") if s in T for k, v in T[s].items()
              if k.endswith("L2_rel") or k.startswith("shaped_L2_rel")]
    res["worst_L2_rel"] = max(finite) if finite and all(math.isfinite(x) for x in finite) else None
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, default=str) + "\n")
    bad = [c for c in checks if c["status"] != "OK"]
    print(f"{status:12s} {a.label} [{a.kind}] worst rel L2 {res['worst_L2_rel']}"
          + (f"; {len(bad)} non-OK checks, first: {bad[0]['check']} = {bad[0]['value']} (limit {bad[0]['limit']})" if bad else "")
          + (f"; error: {res['error']}" if "error" in res else ""))
    return CODES[status]


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""E2b analysis of one far-D muon, keeping four effects SEPARATE.

  k2a_far_analysis.py MUON_EXPERIMENT_DIR --kernel KERNEL_DIR --out NAME [--half-strips 5]

MUON_EXPERIMENT_DIR is recenter_experiment_*/muon_NNNN (far mode): prepared/T1_export, prepared/cloud_moments.csv,
groups/, S0groups/, S0_local/, K2A_local/, K2A_groups/, truncation.json, cmp_E2_K2Alocal_vs_S0local.json.

1. K2A approximation error  - K2A_local vs S0_local (same windows, same frozen weighting): total, and per recentering
                              group (K2A_groups/<g> vs S0groups/<g>), so a far group's error is not hidden by a near one
2. recovered signal         - charge and ion landing OUTSIDE the original entry-centred base +-H strips, in S0_local and
                              K2A_local: signed plane charge inside/outside, sum |q_strip| inside/outside, collected-ion
                              (landing) charge outside, lead strip and its distance from the entry strip
3. finite-window truncation - charge beyond each primary's own +-H window (truncation.json, exact final-state sum rule
                              from the frozen LUT), the +-(H+1) ring, and the frozen replay's own subset strip-sum gap
4. assumptions              - frozen-model limits shared by S0 and K2A (LUT cross-section |u| <= 13.5 mm, LUT phase
                              range margin of the local windows, outermost patch wires |k| = 4, longitudinal extent over
                              which w-invariance / infinite strips are assumed), and K2A-specific moment approximations
                              (per-cloud spread of D and of ln r around the cloud means used by the kernel)
Status: FAIL on missing/nonfinite inputs or a failed truncation check; WARN when a documented flag (thresholds.json
"e2b") is raised; otherwise OK. Writes NAME.json and NAME.md in the muon directory. Exit 0 OK, 1 FAIL, 2 WARN.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("rcmp", HERE / "k2a_recenter_compare.py")
rcmp = importlib.util.module_from_spec(spec); spec.loader.exec_module(rcmp)
EFC = 1.602176634e-4
A, P, UMAX, PHASE_MAX, RW = 0.1, 1.7, 13.5, 12.0, 0.015


def plane_split(d, base, H, window_ns):
    out = {}
    for term, prompt in (("ion_only", False), ("ion_plus_electron_step", True)):
        lab, t, cur = rcmp.load(d, prompt)
        dt = t[1] - t[0]; nb = min(len(t), int(round(window_ns / dt)))
        q = cur[:nb].sum(0) * dt
        T = {}
        for side, b in (("U", base[0]), ("V", base[1])):
            idx = [i for i, l in enumerate(lab) if l[0] == side]
            strips = np.array([int(lab[i].split("_")[1]) for i in idx]); qq = q[idx]
            inside = np.abs(strips - b) <= H
            lead = int(strips[np.argmax(np.abs(qq))])
            T[side] = {"strips": len(strips), "strip_range": [int(strips.min()), int(strips.max())], "base": int(b),
                       "Q_total_fC": float(qq.sum()), "Q_inside_base_window_fC": float(qq[inside].sum()),
                       "Q_outside_base_window_fC": float(qq[~inside].sum()),
                       "absQ_inside_fC": float(np.abs(qq[inside]).sum()), "absQ_outside_fC": float(np.abs(qq[~inside]).sum()),
                       "absQ_outside_frac": float(np.abs(qq[~inside]).sum() / np.abs(qq).sum()) if np.abs(qq).sum() else float("nan"),
                       "lead_strip": lead, "lead_distance_from_entry_strips": lead - int(b)}
        out[term] = T
    return out


def landing_outside(S0g, groups, base, H):
    res = {}
    for side, b in (("U", base[0]), ("V", base[1])):
        tot = out = 0.0
        for name in groups:
            for r in csv.DictReader(open(S0g / name / "replay" / "channels.csv", newline="")):
                if r["side"] != side:
                    continue
                c = float(r["collected_ion_charge_fC"]); tot += c
                if abs(int(r["strip_index"]) - b) > H:
                    out += c
        res[side] = {"collected_fC": tot, "collected_outside_base_window_fC": out,
                     "collected_outside_frac": out / tot if tot else 0.0}
    return res


def first_vertices(d):
    rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
    fd = os.open(d / "vertices.f64", os.O_RDONLY)
    try:
        V = np.array([np.frombuffer(os.pread(fd, 40, 40 * int(r["first_vertex"])), "<f8") for r in rows]).reshape(-1, 5)
    finally:
        os.close(fd)
    return rows, V


def wstat(x, w):
    m = np.sum(w * x) / w.sum()
    return m, math.sqrt(max(0.0, np.sum(w * (x - m) ** 2) / w.sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("E", type=Path); ap.add_argument("--kernel", type=Path, required=True)
    ap.add_argument("--out", default="far_analysis"); ap.add_argument("--half-strips", type=int, default=5)
    ap.add_argument("--window-ns", type=float, default=7000.0)
    ap.add_argument("--thresholds", type=Path, default=HERE / "thresholds.json")
    a = ap.parse_args()
    E, H = a.E, a.half_strips
    th = json.loads(a.thresholds.read_text()).get("e2b", {})
    res = {"tool": "k2a_far_analysis.py", "muon": E.name, "half_strips": H, "flags": [], "status": "OK"}
    try:
        T1 = E / "prepared" / "T1_export"
        rows, V = first_vertices(T1)
        base = (int(rows[0]["upper_index"]), int(rows[0]["lower_index"]))
        meta = json.loads((E / "groups" / "groups.json").read_text())
        groups = list(meta["groups"])
        res["groups"] = {g: {k: meta["groups"][g][k] for k in ("rows", "exact_ions", "m_U", "m_V")} for g in groups}
        res["entry_strips"] = {"U": base[0], "V": base[1]}
        # 1. K2A approximation error
        c = json.loads((E / "cmp_E2_K2Alocal_vs_S0local.json").read_text())
        one = {"status": c["status"], "worst_L2_rel": c.get("worst_L2_rel"), "per_plane": {}}
        for term, T in c["terms"].items():
            for s in ("U", "V"):
                o = T[s]
                one["per_plane"][f"{term}.{s}"] = {k: o.get(k) for k in ("unshaped_L2_rel", "cumulative_L2_rel", "shaped_L2_rel_Tp100",
                                                                         "shaped_L2_rel_Tp200", "window_charge_rel_diff",
                                                                         "max_strip_charge_diff_over_lead", "lead_identical")}
            one["per_plane"][f"{term}.UV_ratio_rel_diff"] = T.get("UV_charge_ratio_rel_diff")
        per_group = {}
        for g in groups:
            mm = {}
            for term, prompt in (("ion_only", False), ("ion_plus_electron_step", True)):
                lr, tr, cr = rcmp.load(E / "S0groups" / g, prompt); lc, tcd, cc = rcmp.load(E / "K2A_groups" / g, prompt)
                M = rcmp.metrics(lr, tr, cr, lc, tcd, cc, a.window_ns, [100.0, 200.0])
                mm[term] = {s: {k: M[s].get(k) for k in ("common_strips", "unshaped_L2_rel", "shaped_L2_rel_Tp100",
                                                         "shaped_L2_rel_Tp200", "window_charge_rel_diff")} for s in ("U", "V")}
            per_group[g] = mm
        one["per_group"] = per_group
        res["1_K2A_approximation_vs_local_S0"] = one
        # 2. recovered signal outside the entry-centred windows
        s0 = plane_split(E / "S0_local", base, H, a.window_ns); k2 = plane_split(E / "K2A_local", base, H, a.window_ns)
        res["2_signal_outside_entry_window"] = {"S0_local": s0, "K2A_local": k2,
                                               "S0_landing": landing_outside(E / "S0groups", groups, base, H),
                                               "K2A_minus_S0_absQ_outside_rel": {
                                                   f"{t}.{s}": (k2[t][s]["absQ_outside_fC"] - s0[t][s]["absQ_outside_fC"]) / s0[t][s]["absQ_outside_fC"]
                                                   if s0[t][s]["absQ_outside_fC"] else None for t in s0 for s in ("U", "V")}}
        # 3. truncation
        tr = json.loads((E / "truncation.json").read_text())
        gaps = {}
        for g in groups:
            sj = E / "S0groups" / g / "replay" / "summary.json"
            gaps[g] = json.loads(sj.read_text()).get("max_subset_strip_sum_gap") if sj.exists() else None
        land = res["2_signal_outside_entry_window"]["S0_landing"]
        res["3_window_truncation"] = {"status": tr["status"], "planes": tr["planes"], "replay_subset_strip_sum_gap": gaps,
                                      "Q_outside_over_collected": {s: (tr["planes"][s]["Q_outside_window_fC"] / land[s]["collected_fC"]
                                                                       if land[s]["collected_fC"] else None) for s in ("U", "V")}}
        if tr["status"] != "OK":
            res["status"] = "FAIL"; res["flags"].append(("FAIL", "truncation check failed"))
        # 4. assumptions
        q = np.array([float(r["weight"]) for r in rows]) * EFC
        u, v = V[:, 1], V[:, 2]
        k = np.rint(u / 3).astype(int); r = np.hypot(u - 3 * k, v)
        prim = np.array([int(r_["primary"]) for r_ in rows])
        shift = {int(p_): tuple(m) for p_, m in meta["primary_shift"].items()}
        mU = np.array([shift[p_][0] for p_ in prim]); mV = np.array([shift[p_][1] for p_ in prim])
        DU, DV = V[:, 3] - P * mU, V[:, 4] - P * mV
        offs = np.arange(-H, H + 1)
        etaU = np.abs((DU - A * u)[:, None] - P * offs[None, :]).max(1)
        etaV = np.abs((DV + A * u)[:, None] - P * offs[None, :]).max(1)
        cl = {}
        for i, key in enumerate(zip(prim, k)):
            cl.setdefault(key, []).append(i)
        sdD, sdlr = [], []
        for key, ii in cl.items():
            ii = np.array(ii); w = q[ii]
            sdD.append((max(wstat(V[ii, 3], w)[1], wstat(V[ii, 4], w)[1]), w.sum()))
            sdlr.append((wstat(np.log(r[ii]), w)[1], w.sum()))
        sd = np.array(sdD); sl = np.array(sdlr)
        kj = json.loads((a.kernel / "kernel.json").read_text())["green_axes"]
        res["4_assumptions"] = {
            "frozen_model_shared_by_S0_and_K2A": {
                "max_abs_u_patch_mm": float(np.abs(u).max()), "lut_cross_section_u_max_mm": UMAX,
                "outer_patch_wire_charge_frac": float(q[np.abs(k) >= 4].sum() / q.sum()),
                "max_local_window_phase_mm": float(max(etaU.max(), etaV.max())), "lut_phase_max_mm": PHASE_MAX,
                "longitudinal_extent_from_entry_mm": [float(min(V[:, 3].min(), V[:, 4].min())), float(max(V[:, 3].max(), V[:, 4].max()))],
                "note": "infinite contiguous strips and a w-invariant 2D field are assumed over this extent; strip ends, "
                        "PCB edges, frames and gaps are not modelled (neither in S0 nor K2A)"},
            "K2A_moment_approximation": {
                "cloud_D_rms_mm_charge_weighted_p50_p95_max": [float(np.quantile(sd[:, 0], .5)), float(np.quantile(sd[:, 0], .95)), float(sd[:, 0].max())],
                "cloud_lnr_rms_p50_p95_max": [float(np.quantile(sl[:, 0], .5)), float(np.quantile(sl[:, 0], .95)), float(sl[:, 0].max())],
                "kernel_r_um": [kj["r_um"][0], kj["r_um"][-1]], "clouds": len(cl)}}
        A4 = res["4_assumptions"]["frozen_model_shared_by_S0_and_K2A"]
        if PHASE_MAX - A4["max_local_window_phase_mm"] < th.get("phase_margin_warn_mm", 0.5):
            res["flags"].append(("WARN", f"local-window LUT phase margin {PHASE_MAX - A4['max_local_window_phase_mm']:.3f} mm"))
        if A4["max_abs_u_patch_mm"] > th.get("u_patch_warn_mm", 12.0):
            res["flags"].append(("WARN", f"ion births at |u| = {A4['max_abs_u_patch_mm']:.2f} mm, near the LUT cross-section edge"))
        if A4["outer_patch_wire_charge_frac"] > th.get("outer_wire_charge_frac_warn", 0.0):
            res["flags"].append(("WARN", f"{A4['outer_patch_wire_charge_frac']:.3%} of the charge on |k| = 4 (field-patch edge)"))
        tf = max(abs(v_) for v_ in res["3_window_truncation"]["Q_outside_over_collected"].values() if v_ is not None)
        if tf > th.get("truncation_frac_warn", 1e-3):
            res["flags"].append(("WARN", f"charge beyond the local ±{H} windows {tf:.3g} of the collected charge"))
        if res["status"] == "OK" and any(f[0] == "WARN" for f in res["flags"]):
            res["status"] = "WARN"
        finite = json.dumps(res, default=str)
        if "NaN" in finite or "Infinity" in finite:
            res["status"] = "FAIL"; res["flags"].append(("FAIL", "nonfinite value in the analysis"))
    except (FileNotFoundError, KeyError, ValueError, rcmp.Inconclusive) as e:
        res["status"] = "FAIL"; res["error"] = f"{type(e).__name__}: {e}"
    (E / f"{a.out}.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    L = [f"# E2b far-D analysis {E.name}: {res['status']}", ""]
    if "error" in res:
        L.append(f"error: {res['error']}")
    else:
        o = res["1_K2A_approximation_vs_local_S0"]
        L += [f"entry strips U {base[0]}, V {base[1]}; groups " + ", ".join(f"{g} ({res['groups'][g]['rows']} rows)" for g in groups), "",
              "## 1. K2A approximation error vs local S0 (same windows)", "",
              f"total: {o['status']}, worst rel L2 {o['worst_L2_rel']}", "",
              "| group | term | plane | strips | unshaped L2 | shaped L2 100 | shaped L2 200 | ΔQ/Q |", "|---|---|---|---|---|---|---|---|"]
        for g, mm in o["per_group"].items():
            for term, T in mm.items():
                for s, x in T.items():
                    L.append(f"| {g} | {term} | {s} | {x['common_strips']} | {x['unshaped_L2_rel']:.3g} | {x['shaped_L2_rel_Tp100']:.3g} | "
                             f"{x['shaped_L2_rel_Tp200']:.3g} | {x['window_charge_rel_diff']:+.3g} |")
        L += ["", "## 2. Signal outside the original entry-centred ±%d strips" % H, "",
              "| variant | term | plane | strips (range) | lead (dist) | Q in / out [fC] | Σ|q| outside frac | collected outside frac |",
              "|---|---|---|---|---|---|---|---|"]
        for var in ("S0_local", "K2A_local"):
            for term, T in res["2_signal_outside_entry_window"][var].items():
                for s, x in T.items():
                    lf = res["2_signal_outside_entry_window"]["S0_landing"][s]["collected_outside_frac"] if var == "S0_local" else None
                    L.append(f"| {var} | {term} | {s} | {x['strips']} ({x['strip_range'][0]}..{x['strip_range'][1]}) | {x['lead_strip']} "
                             f"({x['lead_distance_from_entry_strips']:+d}) | {x['Q_inside_base_window_fC']:.4g} / {x['Q_outside_base_window_fC']:.4g} | "
                             f"{x['absQ_outside_frac']:.3g} | {'' if lf is None else f'{lf:.3g}'} |")
        t3 = res["3_window_truncation"]
        L += ["", "## 3. Finite local-window truncation (frozen LUT sum rule at birth; final integrated charge)", ""]
        for s, x in t3["planes"].items():
            L.append(f"* {s}: charge beyond ±{H}: {x['Q_outside_window_fC']:+.4g} fC ({t3['Q_outside_over_collected'][s]:+.3g} of the collected), "
                     f"max per-ion gap {x['max_per_ion_gap']:.3g}; beyond ±{H + 1}: {x['Q_outside_window_H_plus_1_fC']:+.4g} fC "
                     f"({x['ions_H_plus_1_out_of_lut_range']} ions outside the LUT phase range at ±{H + 1})")
        L.append(f"* frozen replay subset strip-sum gap per group: {t3['replay_subset_strip_sum_gap']}")
        A4 = res["4_assumptions"]
        L += ["", "## 4. Model assumptions", "", "Frozen model (S0 and K2A alike):", ""] + [f"* {k}: {v}" for k, v in A4["frozen_model_shared_by_S0_and_K2A"].items()]
        L += ["", "K2A moment approximation:", ""] + [f"* {k}: {v}" for k, v in A4["K2A_moment_approximation"].items()]
    if res["flags"]:
        L += ["", "Flags:", ""] + [f"* {s}: {t}" for s, t in res["flags"]]
    (E / f"{a.out}.md").write_text("\n".join(L) + "\n")
    print(f"{res['status']} far analysis {E.name}" + (f": {res.get('error')}" if "error" in res else "")
          + ("; flags: " + "; ".join(t for _, t in res["flags"]) if res["flags"] else ""))
    return {"OK": 0, "FAIL": 1, "WARN": 2}[res["status"]]


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Read-only numerical probes of the FROZEN K2A kernel (and its Green grid paths). Nothing is rebuilt or modified.

  k2a_kernel_probe.py NEW_DIR --kernel KERNEL_DIR [--green-paths GREEN_PATHS_DIR] [--wire-pitch-mm 3.0]

1. TAU INTERPOLATION in the early bins (issue 3). The kernel stores cumulative induced charge Q(tau) on the tau grid
   [0, geomspace(0.01 ns, tau_max)]; the evaluator integrates a 10 ns bin as Q(edge - t) - Q(previous edge - t) with
   LINEAR interpolation. For every K0 row (side, wire, r, D) the first ten bin charges are formed for birth phases
   0.5, 1.5, ..., 9.5 ns inside a bin, once with linear and once with 4-point Lagrange interpolation in log tau; the
   difference relative to the row's largest early bin charge estimates the tau-interpolation error (an estimate: both
   are interpolants of the same nodes).
2. SPATIAL INTERPOLATION at fixed tau (10, 20, 50, 100, 1000 ns, 7 us): leave-one-out prediction of every interior
   radial node from its neighbours (linear in log r, as the evaluator) and of every interior D node (linear in D);
   LOO errors are reported raw and /4 (the mid-interval error of a smooth function interpolated on the actual spacing
   is about a quarter of the LOO error), relative to the row's max |Q|.
3. WIRE TRANSLATION (issue 1) between every pair of kernel wires (k, k0): is K[k](r, D, tau) equal to K[k0] at the
   stereo-shifted strip coordinate D -/+ tan(alpha) * pitch_wire * (k - k0) (U: minus, V: plus, from eta = w -/+ A u)?
   Reported per basis row (K0, Kc1, Ks1, KA2) as relative L2 over r, interior D and tau windows, for the predicted
   sign and the opposite sign (a convention check), next to the D-interpolation control of probe 2. A small error
   between interior wires supports a translated kernel ONLY for wires whose field is periodic; it says nothing about
   wires that are not in the table (edge wires of the nine-wire patch), which need their own Green entries to check.
4. TIME HORIZON (issue 2): tau_max of the table; with --green-paths, the drift (collection) time of every Green grid
   start point (grid.csv end_ns - birth): quantiles per wire and radius, i.e. the tau range an extended table must
   cover; whether the existing paths already reach the cathode (end_v = +-gap), in which case extending tau needs no
   new RKF paths, only re-tabulation.
Writes NEW_DIR/kernel_probe.json and kernel_probe.md.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def lagrange4(xn, yn, x):
    """4-point Lagrange interpolation of yn(xn) at x (xn increasing), vectorized over x (1D) and rows (yn: (..., n))."""
    n = len(xn)
    j = np.clip(np.searchsorted(xn, x) - 2, 0, n - 4)
    out = np.zeros(yn.shape[:-1] + (len(x),))
    for a in range(4):
        la = np.ones(len(x))
        for b in range(4):
            if b != a:
                la *= (x - xn[j + b]) / (xn[j + a] - xn[j + b])
        out += yn[..., j + a] * la
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--kernel", type=Path, required=True)
    ap.add_argument("--green-paths", type=Path); ap.add_argument("--wire-pitch-mm", type=float, default=3.0)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    kj = json.loads((a.kernel / "kernel.json").read_text()); ax = kj["green_axes"]
    K = np.load(a.kernel / "K.npy", mmap_mode="r")
    names = kj["basis"]; wires = ax["wires"]; rn = np.asarray(ax["r_um"]); D = np.asarray(ax["D_mm"]); tau = np.asarray(ax["tau_ns"])
    A = float(ax["stereo_tan"]); dD = float(ax["dD_mm"])
    res = {"tool": "k2a_kernel_probe.py", "kernel": str(a.kernel), "basis": names, "wires": wires, "r_um": [float(rn[0]), float(rn[-1])],
           "r_nodes": len(rn), "r_step_factor": float(rn[1] / rn[0]) if len(rn) > 1 else None, "D_mm": [float(D[0]), float(D[-1])],
           "dD_mm": dD, "tau_ns": [float(tau[0]), float(tau[1]), float(tau[-1])], "ntau": len(tau),
           "tau_rel_step": float(tau[2] / tau[1] - 1), "stereo_tan": A, "wire_pitch_mm_assumed": a.wire_pitch_mm}
    i0 = names.index("K0")
    lt = np.log(tau[1:])
    # 1. tau interpolation in the early bins
    tau_err = {}
    for si, s in enumerate(("U", "V")):
        errs = []
        for wk in range(len(wires)):
            Q = np.asarray(K[i0, si, wk], dtype=np.float64)                   # (r, D, tau)
            Q = Q.reshape(-1, len(tau))
            sel = np.argsort(-np.abs(Q[:, -1]))[: max(1, len(Q) // 20)]       # the 5 % rows with the largest total charge
            Qs = Q[sel]
            for ph in np.arange(0.5, 10.0, 1.0):
                x = np.r_[0.0, np.arange(1, 11) * 10.0 - ph]                  # birth at phase ph: first partial bin
                lin = np.stack([np.interp(x, tau, q, left=0.0) for q in Qs])
                xp = np.maximum(x, tau[1])
                cub = lagrange4(lt, Qs[:, 1:], np.log(xp)); cub[:, x < tau[1]] = lin[:, x < tau[1]]
                bl, bc = np.diff(lin, axis=1), np.diff(cub, axis=1)
                sc = np.max(np.abs(bc), axis=1)
                errs.append(np.max(np.abs(bl - bc), axis=1) / np.where(sc > 0, sc, np.nan))
        e = np.concatenate(errs); e = e[np.isfinite(e)]
        tau_err[s] = {"rows": int(len(e)), "median": float(np.median(e)), "p90": float(np.percentile(e, 90)), "max": float(e.max())}
    res["tau_interp_first10bins_rel_to_max_bin"] = tau_err
    # 2. spatial interpolation, leave-one-out
    tsel = [10.0, 20.0, 50.0, 100.0, 1000.0, 7000.0]
    sp = {}
    for si, s in enumerate(("U", "V")):
        rr, dd = {t: [] for t in tsel}, {t: [] for t in tsel}
        for wk in range(len(wires)):
            Q = np.asarray(K[i0, si, wk], dtype=np.float64)                   # (r, D, tau)
            Qt = np.stack([np.array([[np.interp(t, tau, Q[i, j]) for t in tsel] for j in range(Q.shape[1])]) for i in range(Q.shape[0])])  # (r, D, t)
            sc = np.max(np.abs(Qt), axis=(0, 1))                              # per t
            lr = np.log(rn)
            for i in range(1, len(rn) - 1):
                f = (lr[i] - lr[i - 1]) / (lr[i + 1] - lr[i - 1])
                pred = (1 - f) * Qt[i - 1] + f * Qt[i + 1]
                for k, t in enumerate(tsel):
                    rr[t].append(float(np.max(np.abs(pred[:, k] - Qt[i, :, k])) / sc[k]) if sc[k] > 0 else np.nan)
            pred = 0.5 * (Qt[:, :-2] + Qt[:, 2:])
            for k, t in enumerate(tsel):
                dd[t].append(float(np.max(np.abs(pred[:, :, k] - Qt[:, 1:-1, k])) / sc[k]) if sc[k] > 0 else np.nan)
        sp[s] = {f"{t:g}ns": {"r_LOO_max": float(np.nanmax(rr[t])), "r_mid_est": float(np.nanmax(rr[t]) / 4),
                              "D_LOO_max": float(np.nanmax(dd[t])), "D_mid_est": float(np.nanmax(dd[t]) / 4)} for t in tsel}
    res["spatial_interp_rel_to_max"] = sp
    # 3. wire translation
    tw = {"early_0_50ns": tau <= 50.0, "all": np.ones(len(tau), bool)}
    trans = []
    rows = [n for n in ("K0", "Kc1", "Ks1", "KA2") if n in names]
    for si, s in enumerate(("U", "V")):
        sgn = -1.0 if s == "U" else 1.0
        for ia, k in enumerate(wires):
            for ib, k0 in enumerate(wires):
                if k == k0:
                    continue
                for nm in rows:
                    b = names.index(nm)
                    Qk = np.asarray(K[b, si, ia], dtype=np.float64); Q0 = np.asarray(K[b, si, ib], dtype=np.float64)
                    rec = {"side": s, "wire": k, "from_wire": k0, "basis": nm}
                    for lab, sg in (("predicted_sign", sgn), ("opposite_sign", -sgn)):
                        shift = sg * A * a.wire_pitch_mm * (k - k0)
                        Ds = D + shift
                        ok = (Ds >= D[0]) & (Ds <= D[-1])
                        x = (Ds[ok] - D[0]) / dD; j = np.floor(x).astype(int); j = np.minimum(j, len(D) - 2); f = x - j
                        Q0s = (1 - f)[None, :, None] * Q0[:, j] + f[None, :, None] * Q0[:, j + 1]
                        for wl, m in tw.items():
                            num = np.linalg.norm(Qk[:, ok][..., m] - Q0s[..., m]); den = np.linalg.norm(Qk[:, ok][..., m])
                            rec[f"{lab}_{wl}_rel_L2"] = float(num / den) if den > 0 else None
                        rec[f"{lab}_shift_mm"] = shift
                    trans.append(rec)
    res["wire_translation"] = trans
    # 4. time horizon
    th = {"tau_max_ns": float(tau[-1])}
    if a.green_paths and (a.green_paths / "grid.csv").is_file():
        g = list(csv.DictReader(open(a.green_paths / "grid.csv")))
        num = lambda v: float(str(v).replace("np.float64(", "").rstrip(")"))
        ions = {int(r["row"]): r for r in csv.DictReader(open(a.green_paths / "ions.csv"))}
        dtv = np.array([num(r["end_ns"]) - float(ions[int(r["row"])]["birth_ns"]) for r in g])
        ev = np.array([num(r["end_v_mm"]) for r in g])
        th.update(start_points=len(g), drift_time_ns_quantiles={q: float(np.quantile(dtv, q)) for q in (0, 0.5, 0.9, 0.99, 1.0)},
                  frac_beyond_tau_max=float(np.mean(dtv > tau[-1])),
                  frac_landing_on_cathode=float(np.mean(np.abs(np.abs(ev) - np.max(np.abs(ev))) < 1e-6)),
                  cathode_v_mm=float(np.max(np.abs(ev))))
    res["time_horizon"] = th
    a.out.mkdir(parents=True)
    (a.out / "kernel_probe.json").write_text(json.dumps(res, indent=1) + "\n")
    p = lambda v: "n/a" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{100 * v:.3g} %"
    L = ["# Frozen K2A kernel probes (read-only)", "", f"kernel `{a.kernel}`; wires {wires}; r {rn[0]:.6g}-{rn[-1]:.6g} um ({len(rn)} log nodes, "
         f"step x{res['r_step_factor']:.4g}); D {D[0]:.4g}..{D[-1]:.4g} mm step {dD:.4g}; tau 0, {tau[1]:.3g} ... {tau[-1]:.6g} ns "
         f"({len(tau)} nodes, step {100 * res['tau_rel_step']:.3g} %); tan(stereo) {A}; wire pitch {a.wire_pitch_mm} mm (assumed)", "",
         "## 1. tau interpolation, first ten 10 ns bins (linear vs 4-point log-tau Lagrange; 5 % largest K0 rows; all birth phases)", ""]
    for s, x in tau_err.items():
        L.append(f"- {s}: median {p(x['median'])}, p90 {p(x['p90'])}, max {p(x['max'])} of the row's largest early bin charge ({x['rows']} row-phases)")
    L += ["", "## 2. spatial interpolation (leave-one-out max over rows; mid-interval estimate = LOO / 4), relative to max |Q|", "",
          "| side | tau | radial LOO | radial mid est. | D LOO | D mid est. |", "|---|---|---|---|---|---|"]
    for s, x in sp.items():
        for t, v in x.items():
            L.append(f"| {s} | {t} | {p(v['r_LOO_max'])} | {p(v['r_mid_est'])} | {p(v['D_LOO_max'])} | {p(v['D_mid_est'])} |")
    L += ["", "## 3. wire translation K[k](D) vs K[k0](D -/+ A p_w (k - k0)) (relative L2 over r, interior D)", "",
          "| side | wire <- from | basis | shift mm | predicted sign: early | all tau | opposite sign: early | all tau |", "|---|---|---|---|---|---|---|---|"]
    for r in trans:
        L.append(f"| {r['side']} | {r['wire']} <- {r['from_wire']} | {r['basis']} | {r['predicted_sign_shift_mm']:+.3g} | {p(r['predicted_sign_early_0_50ns_rel_L2'])} | "
                 f"{p(r['predicted_sign_all_rel_L2'])} | {p(r['opposite_sign_early_0_50ns_rel_L2'])} | {p(r['opposite_sign_all_rel_L2'])} |")
    L += ["", "## 4. time horizon", "", f"tau_max {th['tau_max_ns']:g} ns."]
    if "start_points" in th:
        L.append(f"Green grid paths: {th['start_points']} start points; drift time to the end of the path (ns) quantiles "
                 + ", ".join(f"q{k:g} {v:.4g}" for k, v in th["drift_time_ns_quantiles"].items())
                 + f"; {100 * th['frac_beyond_tau_max']:.3g} % of paths end after tau_max; {100 * th['frac_landing_on_cathode']:.3g} % end on the "
                 f"cathode (|v| = {th['cathode_v_mm']:g} mm): extending tau to the full drift needs " +
                 ("no new RKF paths (re-tabulation only)" if th["frac_landing_on_cathode"] > 0.999 else "checking: not every path reaches the cathode"))
    (a.out / "kernel_probe.md").write_text("\n".join(L) + "\n")
    print(f"kernel probe -> {a.out / 'kernel_probe.md'}")


if __name__ == "__main__":
    main()

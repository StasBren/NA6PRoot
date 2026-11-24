#!/usr/bin/env python3
"""S-ladder: U/V strip currents of one muon from the single-ion Green's table, for controlled spatial descriptions.

Every variant is a superposition of EXACT single-ion responses G (frozen RKF path from a grid start point + frozen
weighting LUT, s_green_table.py) with weights that sum to each ion's (or cloud's) charge. No charge is moved along a
path computed for another start point; a different spatial description means different grid start points, each with
its own trajectory. Timing is taken unchanged from the input export (use the T1 export: one time per primary).

Ion description from the export (first vertex of each path, patch frame, mm):
  wire k = round(u/3), r = |(u - 3k, v)|, theta = atan2(v, u - 3k), D_s = strip column s (U: col 3, V: col 4), t, q.
Cloud = (Heed primary, wire k).

Variants (comma list via --variants):
  G     closure: each ion at its own (r, theta, D); bilinear in log r and theta, linear in D   -> must reproduce S0
  S1    azimuthally symmetric: each ion spread uniformly in theta at its own r and D (radial + longitudinal kept)
  S2    each ion's angular delta truncated to harmonics n <= 1 (signed grid weights, total exactly 1); r, D per ion
        == the cloud's first angular moments resolved in r and D (radial weighting exact)
  S3    harmonics n <= 2, r and D per ion
  SR    compact radius only: every ion of a cloud at the cloud's r_eff = exp<log r>_q, theta and D per ion
  S2R   r_eff + n <= 1;   S3R  r_eff + n <= 2          (compact radial representation, Sec. S3 of the plan)
  SW    longitudinal: each cloud at its charge-weighted mean D_s, (r, theta) per ion          (S4, one bin per cloud)
  SWW   longitudinal: all ions on one wire at that wire's mean D_s (all primaries)             (S4, one bin per wire)
  K1    r_eff + n <= 1 + cloud mean D   = (N, <cos>, <sin>, r_eff, D) per cloud: the candidate kernel description
  K2    r_eff + n <= 2 + cloud mean D
Signed weights (S2, S3, ...) are a moment-matched surrogate distribution, not a positive charge density; each weight
still multiplies a physical single-ion response, and the total charge is exact.

Output per variant: OUT/<variant>/replay/waveforms.csv (ion term, fC/ns per 10 ns bin, same channel labels as the
reference) and OUT/<variant>/prompt/prompt_waveforms.csv (instantaneous electron-step APPROXIMATION), so
timing-ablation/t_shaped_compare.py compares them with S0 = the frozen replay of the same export.
"""
import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

EFC = 1.602176634e-4
DT = 10.0


def read_export(d):
    seal = dict(l.split("=", 1) for l in (d / "complete.txt").read_text().split() if "=" in l)
    rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
    nv = int(seal["vertices"])
    V = np.memmap(d / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
    first = np.array([int(r["first_vertex"]) for r in rows])
    p0 = np.asarray(V[first])
    w = np.array([float(r["weight"]) for r in rows])
    prim = np.array([int(r["primary"]) for r in rows])
    if not np.isclose(w.sum(), int(seal["exact_ions"]), rtol=1e-9):
        raise SystemExit("weight census mismatch")
    if np.any(np.abs(p0[:, 0] - np.array([float(r["birth_ns"]) for r in rows])) > 1e-7):
        raise SystemExit("first vertex time != birth_ns")
    return rows, p0, w, prim


def radial_weights(r, rn):
    lr, ln = np.log(r), np.log(rn)
    j = np.clip(np.searchsorted(ln, lr) - 1, 0, len(rn) - 2)
    f = (lr - ln[j]) / (ln[j + 1] - ln[j])
    return j, f


def angular_nodes(theta, M, nmax):
    """theta in rad. nmax=None: periodic linear interpolation (2 nodes); nmax=0: uniform; else Fourier truncation."""
    th_m = (np.arange(M) + 0.5) * 2 * np.pi / M
    if nmax is None:
        x = theta / (2 * np.pi / M) - 0.5
        j0 = np.floor(x).astype(int)
        f = x - j0
        return np.stack([j0 % M, (j0 + 1) % M], 1), np.stack([1 - f, f], 1)
    idx = np.broadcast_to(np.arange(M), (len(theta), M))
    wts = np.ones((len(theta), M))
    for n in range(1, nmax + 1):
        wts = wts + 2 * np.cos(n * (th_m[None, :] - theta[:, None]))
    return idx, wts / M


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("green", type=Path, help="s_green_table.py output directory")
    ap.add_argument("export", type=Path, help="path export of the muon (T1 timing), with vertices.f64")
    ap.add_argument("reference_waveforms", type=Path, help="S0 waveforms.csv (frozen replay of the same export)")
    ap.add_argument("out", type=Path)
    ap.add_argument("--variants", default="G,S1,S2,S3,SR,S2R,S3R,SW,SWW,K1,K2")
    ap.add_argument("--window-ns", type=float, default=7000.0)
    ap.add_argument("--suffix", default="", help="appended to every variant directory name (e.g. _r14t72)")
    a = ap.parse_args()
    ax = json.loads((a.green / "axes.json").read_text())
    G = np.load(a.green / "G.npy", mmap_mode="r")
    PHI0 = np.load(a.green / "PHI0.npy")
    wires, rn, D, tau = ax["wires"], np.array(ax["r_um"]), np.array(ax["D_mm"]), np.array(ax["tau_ns"])
    M, dD, dsub = len(ax["theta_deg"]), ax["dD_mm"], ax["d_sub"]
    with open(a.reference_waveforms) as f:
        labels = f.readline().strip().split(",")[1:]
    ref_cur = np.loadtxt(a.reference_waveforms, delimiter=",", skiprows=1, ndmin=2)[:, 1:]
    offs = {s: [int(l.split("_")[1]) for l in labels if l.startswith(s + "_")] for s in ("U", "V")}
    rows, p0, w, prim = read_export(a.export)
    base = {"U": int(rows[0]["upper_index"]), "V": int(rows[0]["lower_index"])}
    off = {s: np.array(offs[s]) - base[s] for s in ("U", "V")}
    nb = int(round(a.window_ns / DT))
    if tau[-1] < a.window_ns:
        raise SystemExit("Green table tau range shorter than the window")
    u, v, t = p0[:, 1], p0[:, 2], p0[:, 0]
    k = np.rint(u / 3.0).astype(int)
    du = u - 3.0 * k
    r = 1000.0 * np.hypot(du, v)
    th = np.mod(np.arctan2(v, du), 2 * np.pi)
    q = w * EFC
    if not set(np.unique(k)) <= set(wires):
        raise SystemExit(f"ions on wires {sorted(set(np.unique(k)) - set(wires))} not in the Green table")
    out_r = (r < rn[0]) | (r > rn[-1])
    clear = r - 15.0
    print(f"ion birth clearance from the 15 um wire surface: min {clear.min():.6g} um; grid r in [{float(rn[0])!r}, {float(rn[-1])!r}] um "
          f"(innermost node clearance {rn[0] - 15.0:.6g} um); ions outside the grid: {int(out_r.sum())}", flush=True)
    rep = {"rows": len(rows), "charge_fC": float(q.sum()), "frac_charge_r_outside_grid": float(q[out_r].sum() / q.sum()),
           "r_um_q01_q50_q99": np.percentile(r, [1, 50, 99]).tolist(), "min_clearance_um": float(clear.min()),
           "grid_r_um": [float(rn[0]), float(rn[-1])], "variants": {}}
    if out_r.any():
        raise SystemExit(f"{out_r.sum()} ions ({100 * rep['frac_charge_r_outside_grid']:.3f}% charge) outside the radial "
                         f"grid [{rn[0]}, {rn[-1]}] um: rebuild the grid (s_grid_spec.py)")
    kk = np.searchsorted(wires, k)
    cloud = np.unique(np.stack([prim, k], 1), axis=0, return_inverse=True)[1].ravel()
    wire_g = np.unique(k, return_inverse=True)[1].ravel()

    def gmean(x, grp, wt):
        return (np.bincount(grp, wt * x) / np.bincount(grp, wt))[grp]

    r_eff = np.exp(gmean(np.log(r), cloud, q))
    tgroups, tg = np.unique(np.round(t, 9), return_inverse=True)
    edges = np.arange(nb + 1) * DT
    a.out.mkdir(parents=True, exist_ok=True)
    for var in a.variants.split(","):
        if (a.out / (var + a.suffix)).exists():
            raise SystemExit(f"refusing to overwrite {a.out / (var + a.suffix)}")
    for var in a.variants.split(","):
        t0w = time.perf_counter()
        nmax = {"G": None, "S1": 0, "S2": 1, "S3": 2, "SR": None, "S2R": 1, "S3R": 2, "SW": None, "SWW": None,
                "K1": 1, "K2": 2}[var]
        rr = r_eff if var in ("SR", "S2R", "S3R", "K1", "K2") else r
        jr, fr = radial_weights(rr, rn)
        ja, fa = angular_nodes(th, M, nmax)
        total = np.zeros((nb, 2, len(off["U"])))
        prompt = np.zeros_like(total)
        for si, s in enumerate(("U", "V")):
            Ds = p0[:, 3 + si]
            if var in ("SW", "K1", "K2"):
                Ds = gmean(Ds, cloud, q)
            elif var == "SWW":
                Ds = gmean(Ds, wire_g, q)
            x = (Ds - D[0]) / dD
            jd = np.floor(x).astype(int)
            fd = x - jd
            lo = jd - dsub * off[s].max()
            hi = jd + 1 - dsub * off[s].min()
            if lo.min() < 0 or hi.max() >= len(D):
                raise SystemExit(f"{var}/{s}: strip-relative D outside the table; rebuild with larger --d-max-mm")
            # contributions: ion x (2 radial) x (angular nodes) x (2 D nodes)
            na = ja.shape[1]
            ion = np.repeat(np.arange(len(q)), 2 * na * 2)
            ir = (jr[:, None] + np.array([0, 1])[None, :])
            wr = np.stack([1 - fr, fr], 1)
            IR = np.broadcast_to(ir[:, :, None, None], (len(q), 2, na, 2)).ravel()
            IA = np.broadcast_to(ja[:, None, :, None], (len(q), 2, na, 2)).ravel()
            ID = np.broadcast_to((jd[:, None] + np.array([0, 1]))[:, None, None, :], (len(q), 2, na, 2)).ravel()
            WT = (q[:, None, None, None] * wr[:, :, None, None] * fa[:, None, :, None]
                  * np.stack([1 - fd, fd], 1)[:, None, None, :]).ravel()
            node = ((kk[ion] * len(rn) + IR) * M + IA) * len(D) + ID
            key = tg[ion].astype(np.int64) * (len(wires) * len(rn) * M * len(D)) + node
            uk, inv = np.unique(key, return_inverse=True)
            W = np.bincount(inv, WT)
            tgi = uk // (len(wires) * len(rn) * M * len(D))
            nd = uk % (len(wires) * len(rn) * M * len(D))
            Gs = G[si].reshape(-1, G.shape[-1])
            P0 = PHI0[si].reshape(-1)
            for g in np.unique(tgi):
                sel = tgi == g
                tp = tgroups[g]
                wg, ng = W[sel], nd[sel]
                b0 = int(np.floor(tp / DT))
                for oi, o in enumerate(off[s]):
                    rows_ = ng - dsub * o
                    Q = wg @ np.asarray(Gs[rows_], dtype=np.float64)
                    qe = np.interp(edges - tp, tau, Q, left=0.0, right=Q[-1])
                    total[:, si, oi] += np.diff(qe)
                    if b0 < nb:
                        prompt[b0, si, oi] += wg @ P0[rows_]
        tc = edges[:-1] + DT / 2
        lab = [f"U_{base['U'] + o}" for o in off["U"]] + [f"V_{base['V'] + o}" for o in off["V"]]
        if lab != labels:
            raise SystemExit("channel labels do not match the reference")
        # consistency checks against the reference (frozen replay of the same export): charge, polarity, mapping
        qref = ref_cur[:nb].reshape(nb, 2, -1).sum(0) * DT          # (side, strip) window charge, fC
        qv = total.sum(0)
        big = np.abs(qref) >= 0.01 * np.abs(qref).max()
        checks = {"plane_window_charge_rel_U": float(qv[0].sum() / qref[0].sum() - 1),
                  "plane_window_charge_rel_V": float(qv[1].sum() / qref[1].sum() - 1),
                  "strip_charge_sign_disagreements": int(np.sum(np.sign(qv[big]) != np.sign(qref[big]))),
                  "lead_strip_same_U": bool(np.argmax(np.abs(qv[0])) == np.argmax(np.abs(qref[0]))),
                  "lead_strip_same_V": bool(np.argmax(np.abs(qv[1])) == np.argmax(np.abs(qref[1])))}
        if var == "G" and (abs(checks["plane_window_charge_rel_U"]) > 0.01 or abs(checks["plane_window_charge_rel_V"]) > 0.01
                           or checks["strip_charge_sign_disagreements"] or not checks["lead_strip_same_U"]
                           or not checks["lead_strip_same_V"]):
            raise SystemExit(f"G closure consistency FAILED (charge, polarity or mapping): {checks}")
        for sub, arr, fn in (("replay", total, "waveforms.csv"), ("prompt", prompt, "prompt_waveforms.csv")):
            (a.out / (var + a.suffix) / sub).mkdir(parents=True)
            np.savetxt(a.out / (var + a.suffix) / sub / fn, np.column_stack((tc, arr.reshape(nb, -1) / DT)), delimiter=",",
                       header="time_center_ns," + ",".join(labels), comments="")
        rep["variants"][var] = {"wall_s": time.perf_counter() - t0w, "contributions": int(len(WT)),
                                "unique_nodes_last_side": int(len(uk)), "time_groups": int(len(tgroups)),
                                "window_charge_fC_U": float(total[:, 0].sum()), "window_charge_fC_V": float(total[:, 1].sum()),
                                "checks_vs_reference": checks}
        print(f"{var}{a.suffix}: {rep['variants'][var]['wall_s']:.1f} s, nodes {len(uk)}, time groups {len(tgroups)}, "
              f"window charge vs ref U {checks['plane_window_charge_rel_U']:+.2e} V {checks['plane_window_charge_rel_V']:+.2e}", flush=True)
    rep["green_axes_sha"] = {k_: ax[k_] for k_ in ("lut_sha256", "grid_ions_sha256", "grid_vertices_sha256")}
    rep["suffix"] = a.suffix
    (a.out / f"s_synthesize{a.suffix}_{time.strftime('%Y%m%d_%H%M%S')}.json").write_text(json.dumps(rep, indent=1) + "\n")


if __name__ == "__main__":
    main()

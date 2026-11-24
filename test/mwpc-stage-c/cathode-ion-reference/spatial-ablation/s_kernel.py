#!/usr/bin/env python3
"""Compact harmonic U/V kernels: evaluate K1/K2(/A) directly from per-cloud moments, no ion records at runtime.

  s_kernel.py build   GREEN OUT                          harmonic kernels from a Green table (s_green_table.py)
  s_kernel.py moments EXPORT OUT.csv                     per-cloud moments from an export (stand-in for the fast generator)
  s_kernel.py apply   KERNEL MOMENTS.csv REF_WAVEFORMS OUTDIR --model K2A [--suffix _TAG] [--repeat 20]

Identity used (exact, linear algebra): for a cloud of charge N at (wire k, r_eff, D_s) with charge-weighted angular
moments C_n = <cos n theta>, S_n = <sin n theta> (n <= nmax) and hemisphere asymmetry A = <sign(sin theta)>, the
K-variant of s_synthesize.py is N * sum_m a_m G[k, r_eff, theta_m, D_s, tau] with
  a_m = (1/M) [1 + 2 sum_n (C_n cos n theta_m + S_n sin n theta_m)] + beta s_perp_m,
  beta = (A - sum_m a0_m s_m) / (s_perp . s)          (hemisphere term only for the ...A models)
so it equals N [K0 + sum_n (C_n Kc_n + S_n Ks_n) + beta KA] with precomputed kernels
  K0 = (1/M) sum_m G_m,  Kc_n = (2/M) sum_m cos(n theta_m) G_m,  Ks_n = (2/M) sum_m sin(n theta_m) G_m,
  KA^(nmax) = sum_m s_perp^(nmax)_m G_m.
Runtime per cloud: 2 radial x 2 D interpolation corners x (1 + 2 nmax [+1]) kernel rows per strip, then a time shift.
`apply` must reproduce s_synthesize.py's K-variant waveforms to rounding (closure of the implementation), and reports
its wall time, which is compared with the frozen replay of the same muon (S0 summary.json wall_s).
The kernels are UNSHAPED (fC per unit charge vs tau); shaping stays a separate, swappable step.
"""
import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

EFC = 1.602176634e-4
DT = 10.0
R_WIRE_UM = 15.0
MODELS = {"K1": (1, False), "K2": (2, False), "K1A": (1, True), "K2A": (2, True)}


def s_perp(th_m, nmax):
    sm = np.sign(np.sin(th_m))
    sp = sm - sm.mean()
    for n in range(1, nmax + 1):
        for b in (np.cos(n * th_m), np.sin(n * th_m)):
            sp = sp - b * (sp @ b) / (b @ b)
    return sm, sp


def build(a):
    ax = json.loads((a.green / "axes.json").read_text())
    G = np.load(a.green / "G.npy", mmap_mode="r")          # (2, k, r, theta, D, tau)
    P0 = np.load(a.green / "PHI0.npy")                     # (2, k, r, theta, D)
    M = G.shape[3]
    th = (np.arange(M) + 0.5) * 2 * np.pi / M
    names, B = ["K0"], [np.full(M, 1.0 / M)]
    for n in (1, 2):
        names += [f"Kc{n}", f"Ks{n}"]; B += [2 / M * np.cos(n * th), 2 / M * np.sin(n * th)]
    for nm in (1, 2):
        names.append(f"KA{nm}"); B.append(s_perp(th, nm)[1])
    B = np.array(B)                                        # (nb, M)
    a.out.mkdir(parents=True, exist_ok=False)
    shape = (len(names),) + G.shape[:3] + G.shape[4:]
    K = np.lib.format.open_memmap(a.out / "K.npy", mode="w+", dtype=np.float32, shape=shape)
    t0 = time.perf_counter()
    for s in range(G.shape[0]):
        for k in range(G.shape[1]):
            for r in range(G.shape[2]):
                g = np.asarray(G[s, k, r], dtype=np.float64)          # (theta, D, tau)
                K[:, s, k, r] = np.tensordot(B, g, axes=(1, 0)).astype(np.float32)
    K.flush()
    np.save(a.out / "KPHI0.npy", np.tensordot(B, P0, axes=(1, 3)))   # (nb, 2, k, r, D)
    meta = {"basis": names, "green": str(a.green), "green_axes": ax, "build_s": time.perf_counter() - t0,
            "shape": list(shape), "bytes": int(K.nbytes)}
    (a.out / "kernel.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"PASS kernels {names} shape {shape} {K.nbytes / 2**20:.1f} MiB  {meta['build_s']:.1f} s")


def read_export(d):
    seal = dict(l.split("=", 1) for l in (d / "complete.txt").read_text().split() if "=" in l)
    rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
    V = np.memmap(d / "vertices.f64", dtype="<f8", mode="r", shape=(int(seal["vertices"]), 5))
    p0 = np.asarray(V[np.array([int(r["first_vertex"]) for r in rows])])
    return rows, p0, np.array([float(r["weight"]) for r in rows]), np.array([int(r["primary"]) for r in rows])


def moments(a):
    rows, p0, w, prim = read_export(a.export)
    u, v, t = p0[:, 1], p0[:, 2], p0[:, 0]
    k = np.rint(u / 3).astype(int)
    du = u - 3 * k
    r = 1000 * np.hypot(du, v)
    th = np.arctan2(v, du)
    q = w * EFC
    keys, cl = np.unique(np.stack([prim, k], 1), axis=0, return_inverse=True)
    cl = cl.ravel()
    Nq = np.bincount(cl, q)
    m = lambda x: np.bincount(cl, q * x) / Nq
    tmean = m(t)
    tspread = np.sqrt(np.maximum(m(t * t) - tmean ** 2, 0))
    cols = {"primary": keys[:, 0], "wire": keys[:, 1], "t_ns": tmean, "t_rms_ns": tspread, "N_fC": Nq,
            "r_eff_um": np.exp(m(np.log(r))), "D_U_mm": m(p0[:, 3]), "D_V_mm": m(p0[:, 4]),
            "C1": m(np.cos(th)), "S1": m(np.sin(th)), "C2": m(np.cos(2 * th)), "S2": m(np.sin(2 * th)),
            "A": m(np.sign(np.sin(th))), "min_clearance_um": np.array([np.min(r[cl == i]) - R_WIRE_UM for i in range(len(keys))]),
            "base_U": np.full(len(keys), int(rows[0]["upper_index"])), "base_V": np.full(len(keys), int(rows[0]["lower_index"]))}
    with open(a.out, "w") as f:
        f.write(",".join(cols) + "\n")
        for i in range(len(keys)):
            f.write(",".join(repr(float(cols[c][i])) if c not in ("primary", "wire", "base_U", "base_V") else str(int(cols[c][i])) for c in cols) + "\n")
    print(f"moments: {len(keys)} clouds from {len(rows)} ion rows; max within-cloud time rms {tspread.max():.3g} ns "
          f"(T1 input gives 0)")


def apply(a):
    meta = json.loads((a.kernel / "kernel.json").read_text())
    ax = meta["green_axes"]
    K = np.load(a.kernel / "K.npy", mmap_mode="r")
    KP = np.load(a.kernel / "KPHI0.npy")
    names = meta["basis"]
    nmax, hemi = MODELS[a.model]
    wires, rn, D, tau = ax["wires"], np.array(ax["r_um"]), np.array(ax["D_mm"]), np.array(ax["tau_ns"])
    M, dD, dsub = len(ax["theta_deg"]), ax["dD_mm"], ax["d_sub"]
    th = (np.arange(M) + 0.5) * 2 * np.pi / M
    sm, sp = s_perp(th, nmax)
    cm = {n: (np.cos(n * th) @ sm, np.sin(n * th) @ sm) for n in (1, 2)}
    with open(a.reference) as f:
        labels = f.readline().strip().split(",")[1:]
    mom = list(csv.DictReader(open(a.moments, newline="")))
    X = {c: np.array([float(r[c]) for r in mom]) for c in mom[0]}
    # strip-offset frames from the labels (same convention as s_synthesize: base read from the export)
    base = {"U": int(mom[0]["base_U"]), "V": int(mom[0]["base_V"])}
    off = {s: np.array([int(l.split("_")[1]) for l in labels if l.startswith(s + "_")]) - base[s] for s in ("U", "V")}
    nb = int(round(a.window_ns / DT))
    edges = np.arange(nb + 1) * DT
    wk = np.searchsorted(wires, X["wire"].astype(int))
    if not set(X["wire"].astype(int)) <= set(wires):
        raise SystemExit("cloud on a wire outside the kernel table")
    if X["r_eff_um"].min() < rn[0] or X["r_eff_um"].max() > rn[-1]:
        raise SystemExit("r_eff outside the radial grid (no extrapolation)")
    # basis coefficients per cloud
    coef = {"K0": np.ones(len(mom))}
    a0s = sm.sum() / M * np.ones(len(mom))
    for n in range(1, nmax + 1):
        coef[f"Kc{n}"], coef[f"Ks{n}"] = X[f"C{n}"], X[f"S{n}"]
        a0s = a0s + (2 / M) * (X[f"C{n}"] * cm[n][0] + X[f"S{n}"] * cm[n][1])
    if hemi:
        coef[f"KA{nmax}"] = (X["A"] - a0s) / (sp @ sm)
    bidx = [names.index(b) for b in coef]
    C = np.stack([coef[b] for b in coef], 1) * X["N_fC"][:, None]          # (clouds, nbasis)
    lr, ln = np.log(X["r_eff_um"]), np.log(rn)
    jr = np.clip(np.searchsorted(ln, lr) - 1, 0, len(rn) - 2)
    fr = (lr - ln[jr]) / (ln[jr + 1] - ln[jr])
    timings = []
    for rep in range(a.repeat):
        t0 = time.perf_counter()
        total = np.zeros((nb, 2, len(off["U"])))
        prompt = np.zeros_like(total)
        for si, s in enumerate(("U", "V")):
            x = (X[f"D_{s}_mm"] - D[0]) / dD
            jd = np.floor(x).astype(int); fd = x - jd
            if (jd - dsub * off[s].max()).min() < 0 or (jd + 1 - dsub * off[s].min()).max() >= len(D):
                raise SystemExit("D outside the kernel table")
            for c in range(len(mom)):
                Kc = K[bidx, si, wk[c]]                                       # (nbasis, r, D, tau) view
                b0 = int(np.floor(X["t_ns"][c] / DT))
                for oi, o in enumerate(off[s]):
                    d0 = jd[c] - dsub * o
                    Q = np.zeros(len(tau))
                    P = 0.0
                    for dr, wr_ in ((0, 1 - fr[c]), (1, fr[c])):
                        for dd, wd_ in ((0, 1 - fd[c]), (1, fd[c])):
                            wgt = wr_ * wd_
                            Q += wgt * (C[c] @ np.asarray(Kc[:, jr[c] + dr, d0 + dd], dtype=np.float64))
                            P += wgt * (C[c] @ KP[bidx, si, wk[c], jr[c] + dr, d0 + dd])
                    qe = np.interp(edges - X["t_ns"][c], tau, Q, left=0.0, right=Q[-1])
                    total[:, si, oi] += np.diff(qe)
                    if b0 < nb:
                        prompt[b0, si, oi] += P
        timings.append(time.perf_counter() - t0)
    name = f"{a.model}kern{a.suffix}"
    tc = edges[:-1] + DT / 2
    for sub, arr, fn in (("replay", total, "waveforms.csv"), ("prompt", prompt, "prompt_waveforms.csv")):
        (a.out / name / sub).mkdir(parents=True, exist_ok=False)
        np.savetxt(a.out / name / sub / fn, np.column_stack((tc, arr.reshape(nb, -1) / DT)), delimiter=",",
                   header="time_center_ns," + ",".join(labels), comments="")
    res = {"model": a.model, "clouds": len(mom), "basis_rows_per_cloud": len(bidx), "repeat": a.repeat,
           "wall_s_median": float(np.median(timings)), "wall_s_min": float(np.min(timings)),
           "implementation": "numpy, single thread, per-cloud Python loop (not optimized)"}
    (a.out / name / "kernel_apply.json").write_text(json.dumps(res, indent=1) + "\n")
    print(f"{name}: {len(mom)} clouds, {len(bidx)} kernel rows/cloud, median {res['wall_s_median'] * 1e3:.2f} ms "
          f"(min {res['wall_s_min'] * 1e3:.2f} ms over {a.repeat})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build"); b.add_argument("green", type=Path); b.add_argument("out", type=Path)
    m = sub.add_parser("moments"); m.add_argument("export", type=Path); m.add_argument("out", type=Path)
    p = sub.add_parser("apply"); p.add_argument("kernel", type=Path); p.add_argument("moments", type=Path)
    p.add_argument("reference", type=Path); p.add_argument("out", type=Path)
    p.add_argument("--model", choices=list(MODELS), default="K2A"); p.add_argument("--suffix", default="")
    p.add_argument("--window-ns", type=float, default=7000.0); p.add_argument("--repeat", type=int, default=5)
    a = ap.parse_args()
    {"build": build, "moments": moments, "apply": apply}[a.cmd](a)


if __name__ == "__main__":
    main()

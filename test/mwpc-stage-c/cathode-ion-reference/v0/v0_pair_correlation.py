#!/usr/bin/env python3
"""Between-primary residual correlation within muons, with explicit assumptions and admissibility checks.

Model (what rho means). For primaries i of one muon, with descriptor y_i (e.g. exact A), context c_i (signed start
cell x yield tercile) and weight w_i (ion yield):
    y_i = m(c_i) + s(c_i) eps_i + n_i,   Corr(eps_i, eps_j) = rho  (i != j, same muon),   eps independent across muons,
    n_i = reservoir (K256) sampling noise with known standard error se_i, independent of everything.
rho is the correlation left after conditioning on the start cell and yield: anything shared by a muon's primaries that
the start cell does not capture (track angle, cluster structure, physics not in the cell) shows up as rho > 0.

Estimator
    r_i   = y_i - m_LMO(c_i)                         m_LMO: cell mean computed WITHOUT the muon's own primaries
    s_c^2 = var_c(y) - mean_c(se^2)                  noise-corrected conditional variance (floored at 0)
    rho   = sum_m sum_{i!=j} w_i w_j r_i r_j / sum_m sum_{i!=j} w_i w_j s_i s_j
Why the numerator has no noise term: n_i, n_j are independent for i != j, so E[r_i r_j] = rho s_i s_j (+ a small
positive term var(m_LMO) for same-cell pairs, O(s^2 / n_cell)).

Admissibility. A common rho is a valid correlation only if every muon's weighted-sum variance is non-negative:
    A_m + rho B_m >= 0,  A_m = sum_i w_i^2 s_i^2,  B_m = sum_{i!=j} w_i w_j s_i s_j   ->   rho >= L_m = -A_m / B_m.
The binding bound is max_m L_m (the muon with the most effective primaries), not a median. Upper bound 1.
An estimate outside [max_m L_m, 1] is flagged INVALID; nothing is clipped.

Known limitations (reported, not corrected): the bootstrap resamples muons but keeps m_LMO and s_c fixed (CIs slightly
optimistic); m and s are cell averages (within-cell dependence of the mean on position leaks into rho as positive
correlation for nearby primaries of one track); yield weights are treated as fixed.

Run `python3 v0_pair_correlation.py --selftest` for the known-rho closure test.
"""
import argparse

import numpy as np
import pandas as pd


def pair_correlation(P, ycol, secol, cond_cols, wcol, nboot=2000, rng=None):
    rng = rng or np.random.default_rng(1)
    d = P.dropna(subset=[ycol, wcol] + cond_cols).copy()
    d = d[d[wcol] > 0]
    d["_ly"] = pd.qcut(np.log(d[wcol]), 3, labels=False, duplicates="drop")
    d["_c"] = d.groupby(cond_cols + ["_ly"]).ngroup()
    g = d.groupby("_c")
    sc, nc = g[ycol].transform("sum"), g[ycol].transform("size")
    gm = d.groupby(["_c", "muon_slot"])[ycol]
    scm, ncm = gm.transform("sum"), gm.transform("size")
    d["_r"] = d[ycol] - ((sc - scm) / (nc - ncm)).where(nc > ncm)
    var = g[ycol].transform("var")
    se2 = g[secol].transform(lambda s: np.nanmean(s ** 2)) if secol else 0.0
    d["_s"] = np.sqrt(np.maximum(var - se2, 0))
    d = d.dropna(subset=["_r", "_s"])
    per = []
    for _, x in d.groupby("muon_slot"):
        if len(x) < 2:
            continue
        w = x[wcol].to_numpy(float); r = x["_r"].to_numpy(); sd = x["_s"].to_numpy()
        wr, ws = w * r, w * sd
        A = (ws ** 2).sum()
        B = ws.sum() ** 2 - A
        per.append((wr.sum() ** 2 - (wr ** 2).sum(), B, A, w.sum() ** 2 / (w ** 2).sum()))
    per = np.array(per)
    est = per[:, 0].sum() / per[:, 1].sum()
    k = rng.integers(0, len(per), (nboot, len(per)))
    bs = per[k, 0].sum(1) / per[k, 1].sum(1)
    L = np.where(per[:, 1] > 0, -per[:, 2] / np.where(per[:, 1] > 0, per[:, 1], 1), -np.inf)
    bound = float(L.max())
    lo, hi = np.percentile(bs, [2.5, 97.5])
    ok = bound <= est <= 1
    return dict(y=ycol, est=float(est), ci95=(float(lo), float(hi)), bound=bound,
                bound_median=float(np.median(L)), frac_muons_violated=float(np.mean(est < L)),
                ci_inside_bounds=bool(lo >= bound and hi <= 1), neff_median=float(np.median(per[:, 3])),
                muons=len(per), primaries=len(d), cells=int(d["_c"].nunique()),
                flag="OK" if ok else "INVALID: outside admissible range (model or estimator wrong; do not use)")


def format_result(r):
    return (f"rho[{r['y']}] = {r['est']:+.4f}  CI95=[{r['ci95'][0]:+.4f},{r['ci95'][1]:+.4f}]  admissible >= {r['bound']:+.4f}"
            f" (binding muon; median muon {r['bound_median']:+.4f}), muons violating {r['frac_muons_violated']:.3f},"
            f" n_eff median {r['neff_median']:.1f}, muons={r['muons']} primaries={r['primaries']} cells={r['cells']}  {r['flag']}")


def synth(rho, n_mu, rng, noise=0.3):
    """Known-rho data: heteroscedastic cells, lognormal weights, reservoir noise with known se."""
    ncell = 48
    mc = rng.normal(0, 1, ncell); sc = rng.uniform(0.5, 1.5, ncell)
    rows = []
    for m in range(n_mu):
        n = int(rng.integers(10, 80))
        c = rng.integers(0, ncell, n)
        s = sc[c]
        C = rho * np.outer(s, s); np.fill_diagonal(C, s * s)
        e = np.linalg.cholesky(C) @ rng.normal(size=n)
        w = rng.lognormal(8, 1, n)
        se = noise * np.ones(n)
        y = mc[c] + e + se * rng.normal(size=n)
        for i in range(n):
            rows.append((m, c[i] // 8, c[i] % 8, w[i], y[i], se[i]))
    return pd.DataFrame(rows, columns=["muon_slot", "cu", "cv", "w", "y", "se"])


def selftest():
    rng = np.random.default_rng(7)
    print("known-rho closure (300 muons, 10-79 primaries each, 48 cells x yield terciles, heteroscedastic, se=0.3):")
    allok = True
    for rho in (-0.01, 0.0, 0.05, 0.2, 0.5):
        est = []
        for rep in range(24):
            D = synth(rho, 300, rng)
            r = pair_correlation(D, "y", "se", ["cu", "cv"], "w", nboot=300, rng=rng)
            est.append(r["est"])
        est = np.array(est)
        bias, sd = est.mean() - rho, est.std(ddof=1)
        ok = abs(bias) < max(3 * sd / np.sqrt(len(est)), 0.01)
        allok &= ok
        print(f"  true {rho:+.3f}: mean est {est.mean():+.4f}  bias {bias:+.4f}  sd {sd:.4f}  last bound {r['bound']:+.4f} {r['flag'][:7]}  {'PASS' if ok else 'FAIL'}")
    print(("PASS" if allok else "FAIL") + " pair-correlation closure (|bias| < max(3 sd/sqrt(24), 0.01))")
    return 0 if allok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    raise SystemExit(selftest() if a.selftest else 0)

#!/usr/bin/env python3
"""
C2 reconstruction study, v2.

Uses the completed 5-mm dense u scan and compares:
  1) nonlinear calibration of the current stereo CoG;
  2) full strip-pattern template matching with a pooled/shrunk covariance;
  3) ridge regression using CoG + the normalized strip pattern;
  4) ridge regression using CoG + low-order strip-shape moments.

All validation is grouped by random seed, so a microscopic seed is never used
both to train and test an estimator.

This is reconstruction-only: no Garfield or Elmer rerun is needed.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def ff(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return math.nan


def ii(row, key):
    return int(float(row[key]))


def discover_points(scan_dir: Path):
    out = []
    for d in sorted((scan_dir / "events").glob("u_*")):
        ev = d / "c2_u_dense_event_summary.csv"
        st = d / "c2_u_dense_strip_summary.csv"
        if ev.is_file() and st.is_file():
            out.append((d, ev, st))
    if not out:
        raise SystemExit("No completed dense-scan points found.")
    return out


def load_dataset(scan_dir: Path):
    records = []
    strip_ids = None

    for d, event_csv, strip_csv in discover_points(scan_dir):
        events = read_csv(event_csv)
        strip_rows = read_csv(strip_csv)

        by_event = defaultdict(lambda: {"minus": {}, "plus": {}})
        ids = set()
        for row in strip_rows:
            iev = ii(row, "event")
            side = row["side"].strip()
            k = ii(row, "strip")
            a = ff(row, "A_fC")
            if not math.isfinite(a):
                a = abs(ff(row, "Q_fC"))
            by_event[iev][side][k] = a
            ids.add(k)

        ids = tuple(sorted(ids))
        if strip_ids is None:
            strip_ids = ids
        elif ids != strip_ids:
            raise SystemExit(f"Strip ids differ in {d}.")

        for ev in events:
            if ii(ev, "readout_valid") != 1:
                continue
            iev = ii(ev, "event")
            if iev not in by_event:
                continue

            sides = {}
            ok = True
            for side in ("minus", "plus"):
                a = np.array(
                    [by_event[iev][side].get(k, math.nan) for k in strip_ids],
                    dtype=float,
                )
                if not np.all(np.isfinite(a)) or np.sum(a) <= 0:
                    ok = False
                    break
                sides[side] = a / np.sum(a)

            if not ok:
                continue

            feat = np.concatenate([sides["minus"], sides["plus"]])
            records.append({
                "event": iev,
                "seed": ii(ev, "random_seed"),
                "u_true": ff(ev, "u0_mm"),
                "u_raw": ff(ev, "u_stereo_mm"),
                "active_wires": ii(ev, "active_wires"),
                "wire_cog_u": ff(ev, "wire_cog_u_mm"),
                "features": feat,
                "minus": sides["minus"],
                "plus": sides["plus"],
            })

    if not records:
        raise SystemExit("No valid events loaded.")
    return records, strip_ids


def seed_folds(records, nfolds):
    seeds = sorted({r["seed"] for r in records})
    nfolds = max(2, min(nfolds, len(seeds)))
    folds = [set() for _ in range(nfolds)]
    for i, s in enumerate(seeds):
        folds[i % nfolds].add(s)
    return folds


def mean_curve(records, getter, step):
    by_u = defaultdict(list)
    for r in records:
        by_u[r["u_true"]].append(np.atleast_1d(getter(r)).astype(float))

    nodes = np.array(sorted(by_u), dtype=float)
    means = np.vstack([np.mean(np.vstack(by_u[u]), axis=0) for u in nodes])

    lo, hi = float(nodes[0]), float(nodes[-1])
    n = int(round((hi - lo) / step)) + 1
    grid = np.linspace(lo, hi, n)
    mean_grid = np.column_stack(
        [np.interp(grid, nodes, means[:, j]) for j in range(means.shape[1])]
    )
    return nodes, means, grid, mean_grid


def baseline_predictor(train, step):
    _, _, grid, mg = mean_curve(train, lambda r: [r["u_raw"]], step)

    # A single global residual scale is deliberately used; per-u variance made
    # the previous template unnecessarily unstable.
    def predict(r):
        d = mg[:, 0] - r["u_raw"]
        return float(grid[int(np.argmin(d * d))])

    return predict


def pooled_template_predictor(train, step, shrink=0.75):
    nodes, means, grid, mean_grid = mean_curve(
        train, lambda r: r["features"], step
    )

    node_mean = {float(u): means[i] for i, u in enumerate(nodes)}
    residuals = np.vstack([
        r["features"] - node_mean[float(r["u_true"])]
        for r in train
    ])
    cov = np.cov(residuals, rowvar=False, ddof=1)
    diag = np.diag(np.diag(cov))
    scale = float(np.trace(cov) / cov.shape[0])
    cov_reg = (
        (1.0 - shrink) * cov
        + shrink * diag
        + max(scale * 1.0e-4, 1.0e-10) * np.eye(cov.shape[0])
    )
    inv = np.linalg.pinv(cov_reg)

    def predict(r):
        d = mean_grid - r["features"][None, :]
        score = np.einsum("ij,jk,ik->i", d, inv, d)
        return float(grid[int(np.argmin(score))])

    return predict


def strip_moments(r, strip_ids):
    k = np.asarray(strip_ids, dtype=float)
    out = []
    for side in ("minus", "plus"):
        p = r[side]
        m = float(np.dot(p, k))
        c = k - m
        v = float(np.dot(p, c * c))
        s3 = float(np.dot(p, c * c * c))
        skew = s3 / max(v ** 1.5, 1.0e-12)
        center = float(p[np.argmin(np.abs(k))])
        near = float(np.sum(p[np.abs(k) <= 1]))
        outer = float(np.sum(p[np.abs(k) >= 2]))
        out += [m, v, skew, center, near, outer]
    return np.asarray(out, dtype=float)


def full_ridge_features(r):
    u = r["u_raw"]
    return np.concatenate([
        np.asarray([u, u * u, u * u * u], dtype=float),
        r["features"],
    ])


def moment_ridge_features(r, strip_ids):
    u = r["u_raw"]
    return np.concatenate([
        np.asarray([u, u * u, u * u * u], dtype=float),
        strip_moments(r, strip_ids),
    ])


def fit_ridge(train, getter, lam):
    X = np.vstack([getter(r) for r in train]).astype(float)
    y = np.asarray([r["u_true"] for r in train], dtype=float)

    mu = np.mean(X, axis=0)
    sd = np.std(X, axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd < 1.0e-12)] = 1.0
    Z = (X - mu) / sd

    y0 = float(np.mean(y))
    yc = y - y0
    a = Z.T @ Z + lam * np.eye(Z.shape[1])
    b = Z.T @ yc
    beta = np.linalg.solve(a, b)

    def predict(r):
        z = (getter(r) - mu) / sd
        return float(y0 + z @ beta)

    return predict


def choose_lambda(train, getter, candidates):
    # Inner seed-grouped CV. Three folds are enough for the tiny lambda search.
    folds = seed_folds(train, 3)
    best = None
    for lam in candidates:
        sq = []
        for test_seeds in folds:
            tr = [r for r in train if r["seed"] not in test_seeds]
            te = [r for r in train if r["seed"] in test_seeds]
            if not tr or not te:
                continue
            pred = fit_ridge(tr, getter, lam)
            for r in te:
                e = pred(r) - r["u_true"]
                sq.append(e * e)
        rmse = math.sqrt(float(np.mean(sq))) if sq else math.inf
        if best is None or rmse < best[0]:
            best = (rmse, lam)
    return best[1]


def metric(rows, key):
    e = np.asarray([r[key] - r["u_true"] for r in rows], dtype=float)
    ae = np.abs(e)
    return {
        "n": len(rows),
        "bias_mm": float(np.mean(e)),
        "rmse_mm": float(np.sqrt(np.mean(e * e))),
        "mae_mm": float(np.mean(ae)),
        "q68_abs_error_mm": float(np.quantile(ae, 0.68)),
        "q95_abs_error_mm": float(np.quantile(ae, 0.95)),
        "outlier_gt_0p5_fraction": float(np.mean(ae > 0.5)),
    }


def write_rows(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", type=Path, required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--grid-step-mm", type=float, default=0.005)
    args = ap.parse_args()

    records, strip_ids = load_dataset(args.scan_dir)
    outer = seed_folds(records, args.folds)
    lambdas = [1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0, 10.0, 100.0]

    preds = []
    fold_info = []

    for ifold, test_seeds in enumerate(outer):
        train = [r for r in records if r["seed"] not in test_seeds]
        test = [r for r in records if r["seed"] in test_seeds]

        p_cog = baseline_predictor(train, args.grid_step_mm)
        p_tpl = pooled_template_predictor(
            train, args.grid_step_mm, shrink=0.75
        )

        full_getter = full_ridge_features
        moment_getter = lambda r: moment_ridge_features(r, strip_ids)

        lam_full = choose_lambda(train, full_getter, lambdas)
        lam_mom = choose_lambda(train, moment_getter, lambdas)
        p_full = fit_ridge(train, full_getter, lam_full)
        p_mom = fit_ridge(train, moment_getter, lam_mom)

        fold_info.append({
            "fold": ifold,
            "lambda_full_strip_ridge": lam_full,
            "lambda_moment_ridge": lam_mom,
        })

        for r in test:
            preds.append({
                **r,
                "fold": ifold,
                "u_cog_cal_mm": p_cog(r),
                "u_pooled_template_mm": p_tpl(r),
                "u_full_strip_ridge_mm": p_full(r),
                "u_moment_ridge_mm": p_mom(r),
            })

    preds.sort(key=lambda r: (r["u_true"], r["seed"]))

    estimator_keys = [
        ("nonlinear_cog_calibration", "u_cog_cal_mm"),
        ("pooled_strip_template", "u_pooled_template_mm"),
        ("cog_plus_full_strip_ridge", "u_full_strip_ridge_mm"),
        ("cog_plus_strip_moment_ridge", "u_moment_ridge_mm"),
    ]

    pred_rows = []
    for r in preds:
        pred_rows.append({
            "u_true_mm": r["u_true"],
            "random_seed": r["seed"],
            "event": r["event"],
            "fold": r["fold"],
            "active_wires": r["active_wires"],
            "wire_cog_u_mm": r["wire_cog_u"],
            "u_raw_mm": r["u_raw"],
            **{key: r[key] for _, key in estimator_keys},
        })

    write_rows(
        args.scan_dir / "c2_strip_reco_v2_predictions.csv",
        pred_rows,
        list(pred_rows[0]),
    )

    regions = {
        "all": lambda r: True,
        "core_abs_u_le_1p5": lambda r: abs(r["u_true"]) <= 1.5 + 1e-12,
        "edge_abs_u_gt_1p5": lambda r: abs(r["u_true"]) > 1.5 + 1e-12,
        "single_wire": lambda r: r["active_wires"] == 1,
        "multiwire": lambda r: r["active_wires"] > 1,
    }

    mrows = []
    for region, selector in regions.items():
        sub = [r for r in preds if selector(r)]
        if not sub:
            continue
        for name, key in estimator_keys:
            mrows.append({"region": region, "estimator": name, **metric(sub, key)})

    metric_file = args.scan_dir / "c2_strip_reco_v2_metrics.csv"
    write_rows(
        metric_file,
        mrows,
        [
            "region", "estimator", "n", "bias_mm", "rmse_mm", "mae_mm",
            "q68_abs_error_mm", "q95_abs_error_mm",
            "outlier_gt_0p5_fraction",
        ],
    )

    byu = []
    for u in sorted({r["u_true"] for r in preds}):
        sub = [r for r in preds if r["u_true"] == u]
        for name, key in estimator_keys:
            byu.append({"u_true_mm": u, "estimator": name, **metric(sub, key)})

    byu_file = args.scan_dir / "c2_strip_reco_v2_by_u.csv"
    write_rows(
        byu_file,
        byu,
        [
            "u_true_mm", "estimator", "n", "bias_mm", "rmse_mm", "mae_mm",
            "q68_abs_error_mm", "q95_abs_error_mm",
            "outlier_gt_0p5_fraction",
        ],
    )

    fold_file = args.scan_dir / "c2_strip_reco_v2_folds.csv"
    write_rows(
        fold_file,
        fold_info,
        ["fold", "lambda_full_strip_ridge", "lambda_moment_ridge"],
    )

    print("C2 reconstruction v2: 5-mm baseline")
    print(f"  valid events : {len(records)}")
    print(f"  strip ids    : {strip_ids}")
    print()
    for row in mrows:
        if row["region"] == "all":
            print(
                f"  {row['estimator']:30s} "
                f"RMSE={row['rmse_mm']:.4f} mm  "
                f"MAE={row['mae_mm']:.4f} mm  "
                f"68%={row['q68_abs_error_mm']:.4f} mm  "
                f"95%={row['q95_abs_error_mm']:.4f} mm"
            )
    print()
    print("Wrote:")
    print(f"  {metric_file}")
    print(f"  {byu_file}")
    print(f"  {fold_file}")
    print(f"  {args.scan_dir / 'c2_strip_reco_v2_predictions.csv'}")


if __name__ == "__main__":
    main()

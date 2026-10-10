#!/usr/bin/env python3
"""
C2 reconstruction study, v3: symmetry-aware edge recovery.

Motivation from v2:
  * nonlinear stereo-CoG calibration is best in the central part of the cell;
  * full strip-pattern template matching becomes substantially better near the
    positive cell edge, but is asymmetric and less stable on the negative side;
  * unconstrained linear ridge models can extrapolate catastrophically on rare
    multi-wire/tail events.

This study therefore does three deliberately conservative things:
  1) keeps the nonlinear CoG calibration as the baseline/core estimator;
  2) enforces the expected u -> -u detector symmetry on the strip templates by
     augmenting training data with mirrored events (minus/plus cathodes swap);
  3) builds a hybrid estimator that switches to the strip template only when
     the template itself indicates an edge-like event.

The hybrid switching threshold is selected inside each outer training fold using
an inner seed-grouped cross-validation.  Therefore the final test event's seed
never influences either the calibration/template or the switch threshold.

No Garfield/Elmer rerun is required.
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
    except (KeyError, TypeError, ValueError):
        return math.nan


def ii(row, key):
    return int(float(row[key]))


def discover_points(scan_dir: Path):
    points = []
    for d in sorted((scan_dir / "events").glob("u_*")):
        ev = d / "c2_u_dense_event_summary.csv"
        st = d / "c2_u_dense_strip_summary.csv"
        if ev.is_file() and st.is_file():
            points.append((d, ev, st))
    if not points:
        raise SystemExit("No completed dense-scan points found.")
    return points


def load_dataset(scan_dir: Path):
    records = []
    strip_ids = None

    for d, event_csv, strip_csv in discover_points(scan_dir):
        events = read_csv(event_csv)
        strips = read_csv(strip_csv)

        by_event = defaultdict(lambda: {"minus": {}, "plus": {}})
        ids = set()
        for row in strips:
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
            raise SystemExit(f"Strip-id mismatch in {d}.")

        for ev in events:
            if ii(ev, "readout_valid") != 1:
                continue
            iev = ii(ev, "event")
            if iev not in by_event:
                continue

            side_vectors = {}
            ok = True
            for side in ("minus", "plus"):
                a = np.asarray(
                    [by_event[iev][side].get(k, math.nan) for k in strip_ids],
                    dtype=float,
                )
                if not np.all(np.isfinite(a)) or float(np.sum(a)) <= 0:
                    ok = False
                    break
                side_vectors[side] = a / np.sum(a)
            if not ok:
                continue

            records.append({
                "event": iev,
                "seed": ii(ev, "random_seed"),
                "u_true": ff(ev, "u0_mm"),
                "u_raw": ff(ev, "u_stereo_mm"),
                "active_wires": ii(ev, "active_wires"),
                "wire_cog_u": ff(ev, "wire_cog_u_mm"),
                "minus": side_vectors["minus"],
                "plus": side_vectors["plus"],
                "features": np.concatenate(
                    [side_vectors["minus"], side_vectors["plus"]]
                ),
            })

    if not records:
        raise SystemExit("No valid dense-scan events loaded.")
    return records, strip_ids


def mirror_record(r):
    # At w=0 in the symmetric 2.5+2.5 mm chamber, u -> -u exchanges the two
    # stereo coordinates.  The measured templates verify that same-index
    # minus/plus exchange is substantially closer to the symmetry than
    # additionally reversing the strip index.
    minus = np.array(r["plus"], copy=True)
    plus = np.array(r["minus"], copy=True)
    return {
        **r,
        "u_true": -r["u_true"],
        "minus": minus,
        "plus": plus,
        "features": np.concatenate([minus, plus]),
        "_mirrored": True,
    }


def augment_symmetry(records):
    out = []
    for r in records:
        rr = dict(r)
        rr["_mirrored"] = False
        out.append(rr)
        out.append(mirror_record(rr))
    return out


def seed_folds(records, nfolds):
    seeds = sorted({r["seed"] for r in records})
    nfolds = max(2, min(nfolds, len(seeds)))
    folds = [set() for _ in range(nfolds)]
    for i, seed in enumerate(seeds):
        folds[i % nfolds].add(seed)
    return folds


def mean_curve(records, getter, step):
    by_u = defaultdict(list)
    for r in records:
        by_u[r["u_true"]].append(np.atleast_1d(getter(r)).astype(float))

    nodes = np.asarray(sorted(by_u), dtype=float)
    means = np.vstack(
        [np.mean(np.vstack(by_u[u]), axis=0) for u in nodes]
    )

    lo, hi = float(nodes[0]), float(nodes[-1])
    n = int(round((hi - lo) / step)) + 1
    grid = np.linspace(lo, hi, n)
    mean_grid = np.column_stack(
        [np.interp(grid, nodes, means[:, j]) for j in range(means.shape[1])]
    )
    return nodes, means, grid, mean_grid


def fit_cog_calibration(train, step):
    _, _, grid, mean_grid = mean_curve(
        train, lambda r: np.asarray([r["u_raw"]]), step
    )

    def predict(r):
        d = mean_grid[:, 0] - r["u_raw"]
        return float(grid[int(np.argmin(d * d))])

    return predict


def fit_symmetry_template(train, step, shrink=0.85):
    aug = augment_symmetry(train)
    nodes, means, grid, mean_grid = mean_curve(
        aug, lambda r: r["features"], step
    )

    node_mean = {float(u): means[i] for i, u in enumerate(nodes)}
    residuals = np.vstack([
        r["features"] - node_mean[float(r["u_true"])]
        for r in aug
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


def choose_hybrid_threshold(
    train,
    grid_step,
    threshold_candidates,
    inner_folds=3,
):
    folds = seed_folds(train, inner_folds)
    scores = []

    for threshold in threshold_candidates:
        sq = []
        ae = []
        for test_seeds in folds:
            tr = [r for r in train if r["seed"] not in test_seeds]
            te = [r for r in train if r["seed"] in test_seeds]
            if not tr or not te:
                continue

            pc = fit_cog_calibration(tr, grid_step)
            pt = fit_symmetry_template(tr, grid_step)

            for r in te:
                uc = pc(r)
                ut = pt(r)
                pred = ut if abs(ut) >= threshold else uc
                e = pred - r["u_true"]
                sq.append(e * e)
                ae.append(abs(e))

        rmse = math.sqrt(float(np.mean(sq))) if sq else math.inf
        mae = float(np.mean(ae)) if ae else math.inf
        scores.append((rmse, mae, threshold))

    scores.sort()
    return scores[0][2], scores


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
    outer_folds = seed_folds(records, args.folds)
    threshold_candidates = np.arange(1.10, 1.86, 0.05)

    predictions = []
    fold_rows = []

    for ifold, test_seeds in enumerate(outer_folds):
        train = [r for r in records if r["seed"] not in test_seeds]
        test = [r for r in records if r["seed"] in test_seeds]

        threshold, inner_scores = choose_hybrid_threshold(
            train,
            args.grid_step_mm,
            threshold_candidates,
        )

        pc = fit_cog_calibration(train, args.grid_step_mm)
        pt = fit_symmetry_template(train, args.grid_step_mm)

        best_inner = min(inner_scores)
        fold_rows.append({
            "fold": ifold,
            "hybrid_threshold_mm": threshold,
            "inner_rmse_mm": best_inner[0],
            "inner_mae_mm": best_inner[1],
        })

        for r in test:
            uc = pc(r)
            ut = pt(r)
            uh = ut if abs(ut) >= threshold else uc
            predictions.append({
                **r,
                "fold": ifold,
                "hybrid_threshold_mm": threshold,
                "u_cog_cal_mm": uc,
                "u_sym_strip_template_mm": ut,
                "u_hybrid_mm": uh,
                "hybrid_used_strip": int(abs(ut) >= threshold),
            })

    predictions.sort(key=lambda r: (r["u_true"], r["seed"]))

    estimator_keys = [
        ("nonlinear_cog_calibration", "u_cog_cal_mm"),
        ("symmetry_strip_template", "u_sym_strip_template_mm"),
        ("hybrid_cog_plus_strip", "u_hybrid_mm"),
    ]

    pred_rows = []
    for r in predictions:
        pred_rows.append({
            "u_true_mm": r["u_true"],
            "random_seed": r["seed"],
            "event": r["event"],
            "fold": r["fold"],
            "active_wires": r["active_wires"],
            "wire_cog_u_mm": r["wire_cog_u"],
            "u_raw_mm": r["u_raw"],
            "hybrid_threshold_mm": r["hybrid_threshold_mm"],
            "hybrid_used_strip": r["hybrid_used_strip"],
            **{key: r[key] for _, key in estimator_keys},
        })

    pred_file = args.scan_dir / "c2_strip_reco_v3_predictions.csv"
    write_rows(pred_file, pred_rows, list(pred_rows[0]))

    regions = {
        "all": lambda r: True,
        "core_abs_u_le_1p5": lambda r: abs(r["u_true"]) <= 1.5 + 1e-12,
        "edge_abs_u_gt_1p5": lambda r: abs(r["u_true"]) > 1.5 + 1e-12,
        "single_wire": lambda r: r["active_wires"] == 1,
        "multiwire": lambda r: r["active_wires"] > 1,
    }

    metric_rows = []
    for region, selector in regions.items():
        sub = [r for r in predictions if selector(r)]
        if not sub:
            continue
        for name, key in estimator_keys:
            metric_rows.append({
                "region": region,
                "estimator": name,
                **metric(sub, key),
            })

    metric_file = args.scan_dir / "c2_strip_reco_v3_metrics.csv"
    write_rows(
        metric_file,
        metric_rows,
        [
            "region", "estimator", "n", "bias_mm", "rmse_mm", "mae_mm",
            "q68_abs_error_mm", "q95_abs_error_mm",
            "outlier_gt_0p5_fraction",
        ],
    )

    by_u_rows = []
    for u in sorted({r["u_true"] for r in predictions}):
        sub = [r for r in predictions if r["u_true"] == u]
        strip_fraction = float(np.mean([r["hybrid_used_strip"] for r in sub]))
        for name, key in estimator_keys:
            row = {
                "u_true_mm": u,
                "estimator": name,
                **metric(sub, key),
            }
            if name == "hybrid_cog_plus_strip":
                row["hybrid_strip_fraction"] = strip_fraction
            else:
                row["hybrid_strip_fraction"] = math.nan
            by_u_rows.append(row)

    by_u_file = args.scan_dir / "c2_strip_reco_v3_by_u.csv"
    write_rows(
        by_u_file,
        by_u_rows,
        [
            "u_true_mm", "estimator", "n", "bias_mm", "rmse_mm", "mae_mm",
            "q68_abs_error_mm", "q95_abs_error_mm",
            "outlier_gt_0p5_fraction", "hybrid_strip_fraction",
        ],
    )

    fold_file = args.scan_dir / "c2_strip_reco_v3_folds.csv"
    write_rows(
        fold_file,
        fold_rows,
        ["fold", "hybrid_threshold_mm", "inner_rmse_mm", "inner_mae_mm"],
    )

    print("C2 reconstruction v3: symmetry-aware hybrid")
    print(f"  valid events : {len(records)}")
    print(f"  strip ids    : {strip_ids}")
    print()
    for row in metric_rows:
        if row["region"] == "all":
            print(
                f"  {row['estimator']:28s} "
                f"RMSE={row['rmse_mm']:.4f} mm  "
                f"MAE={row['mae_mm']:.4f} mm  "
                f"68%={row['q68_abs_error_mm']:.4f} mm  "
                f"95%={row['q95_abs_error_mm']:.4f} mm"
            )
    print()
    print("Outer-fold hybrid thresholds:")
    for row in fold_rows:
        print(
            f"  fold {row['fold']}: {row['hybrid_threshold_mm']:.2f} mm "
            f"(inner RMSE {row['inner_rmse_mm']:.4f} mm)"
        )
    print()
    print("Wrote:")
    print(f"  {metric_file}")
    print(f"  {by_u_file}")
    print(f"  {fold_file}")
    print(f"  {pred_file}")


if __name__ == "__main__":
    main()

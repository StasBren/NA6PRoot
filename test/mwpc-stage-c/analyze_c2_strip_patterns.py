#!/usr/bin/env python3
"""
Compare the current stereo-CoG estimator with a calibration that uses the full
two-cathode strip-charge pattern from an already completed dense C2 scan.

No Garfield/Elmer rerun is required.  The script reads the per-u0
c2_u_dense_event_summary.csv and c2_u_dense_strip_summary.csv files written by
run_c2_u_cell_scan_long.sh.

The comparison is seed-group cross-validated: the same microscopic seed never
appears in both the template training sample and its test sample.  This matters
because the dense scan intentionally reuses the same seed block at every u0.

Observables used by the strip-template estimator:
  a_{side,k} = |Q_{side,k}| / sum_j |Q_{side,j}|

Thus total avalanche gain is removed and only the charge-sharing shape is used.
No truth, ion-birth, gain-weighted seed coordinate, or wire assignment enters
the estimator.

Outputs:
  c2_strip_estimator_predictions.csv
  c2_strip_estimator_metrics.csv
  c2_strip_estimator_by_u.csv
  c2_strip_templates.csv
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


def as_float(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def as_int(row, key):
    return int(float(row[key]))


def discover_points(scan_dir: Path):
    event_root = scan_dir / "events"
    found = []
    for point_dir in sorted(event_root.glob("u_*")):
        event_csv = point_dir / "c2_u_dense_event_summary.csv"
        strip_csv = point_dir / "c2_u_dense_strip_summary.csv"
        if event_csv.is_file() and strip_csv.is_file():
            found.append((point_dir, event_csv, strip_csv))
    if not found:
        raise SystemExit(
            f"No complete dense-scan point directories found under {event_root}"
        )
    return found


def load_dataset(scan_dir: Path):
    records = []
    strip_ids = None

    for point_dir, event_csv, strip_csv in discover_points(scan_dir):
        events = read_csv(event_csv)
        strips = read_csv(strip_csv)

        by_event = defaultdict(lambda: {"minus": {}, "plus": {}})
        ids = set()
        for row in strips:
            iev = as_int(row, "event")
            side = row["side"].strip()
            k = as_int(row, "strip")
            amp = as_float(row, "A_fC")
            if not math.isfinite(amp):
                amp = abs(as_float(row, "Q_fC"))
            by_event[iev][side][k] = amp
            ids.add(k)

        ids = tuple(sorted(ids))
        if strip_ids is None:
            strip_ids = ids
        elif ids != strip_ids:
            raise SystemExit(
                f"Strip-id mismatch in {point_dir}: {ids} vs {strip_ids}"
            )

        for event in events:
            if as_int(event, "readout_valid") != 1:
                continue
            iev = as_int(event, "event")
            if iev not in by_event:
                continue

            feat = []
            ok = True
            for side in ("minus", "plus"):
                amps = np.array(
                    [by_event[iev][side].get(k, math.nan) for k in strip_ids],
                    dtype=float,
                )
                if not np.all(np.isfinite(amps)):
                    ok = False
                    break
                denom = float(np.sum(amps))
                if denom <= 0:
                    ok = False
                    break
                feat.extend((amps / denom).tolist())

            if not ok:
                continue

            records.append(
                {
                    "point_dir": str(point_dir),
                    "event": iev,
                    "seed": as_int(event, "random_seed"),
                    "u_true": as_float(event, "u0_mm"),
                    "u_raw": as_float(event, "u_stereo_mm"),
                    "active_wires": as_int(event, "active_wires"),
                    "wire_cog_u": as_float(event, "wire_cog_u_mm"),
                    "features": np.asarray(feat, dtype=float),
                }
            )

    if not records:
        raise SystemExit("No valid events with strip patterns were loaded.")
    return records, strip_ids


def build_template_grid(
    records,
    feature_getter,
    grid_step,
    variance_floor_fraction=0.05,
):
    by_u = defaultdict(list)
    for rec in records:
        by_u[rec["u_true"]].append(np.atleast_1d(feature_getter(rec)).astype(float))

    u_nodes = np.array(sorted(by_u), dtype=float)
    means = np.vstack(
        [np.mean(np.vstack(by_u[u]), axis=0) for u in u_nodes]
    )
    variances = np.vstack(
        [
            np.var(np.vstack(by_u[u]), axis=0, ddof=1)
            if len(by_u[u]) > 1
            else np.zeros(means.shape[1])
            for u in u_nodes
        ]
    )

    lo = float(u_nodes[0])
    hi = float(u_nodes[-1])
    ngrid = int(round((hi - lo) / grid_step)) + 1
    u_grid = np.linspace(lo, hi, ngrid)

    mean_grid = np.column_stack(
        [np.interp(u_grid, u_nodes, means[:, j]) for j in range(means.shape[1])]
    )
    var_grid = np.column_stack(
        [np.interp(u_grid, u_nodes, variances[:, j]) for j in range(variances.shape[1])]
    )

    positive = variances[variances > 0]
    if positive.size:
        variance_floor = variance_floor_fraction * float(np.median(positive))
    else:
        variance_floor = 1.0e-8
    variance_floor = max(variance_floor, 1.0e-8)

    return {
        "u_nodes": u_nodes,
        "means": means,
        "variances": variances,
        "u_grid": u_grid,
        "mean_grid": mean_grid,
        "var_grid": var_grid,
        "variance_floor": variance_floor,
    }


def predict_from_template(x, template, use_variance=True):
    x = np.atleast_1d(x).astype(float)
    diff = template["mean_grid"] - x[None, :]
    if use_variance:
        scale2 = template["var_grid"] + template["variance_floor"]
        score = np.sum(diff * diff / scale2, axis=1)
    else:
        score = np.sum(diff * diff, axis=1)
    return float(template["u_grid"][int(np.argmin(score))])


def make_seed_folds(records, nfolds):
    seeds = sorted({r["seed"] for r in records})
    if len(seeds) < 2:
        raise SystemExit("Need at least two independent seeds for cross-validation.")
    nfolds = max(2, min(nfolds, len(seeds)))
    folds = [set() for _ in range(nfolds)]
    for i, seed in enumerate(seeds):
        folds[i % nfolds].add(seed)
    return folds


def metrics(rows, pred_key):
    err = np.array([r[pred_key] - r["u_true"] for r in rows], dtype=float)
    ae = np.abs(err)
    return {
        "n": len(rows),
        "bias_mm": float(np.mean(err)),
        "rmse_mm": float(np.sqrt(np.mean(err * err))),
        "mae_mm": float(np.mean(ae)),
        "q68_abs_error_mm": float(np.quantile(ae, 0.68)),
        "q95_abs_error_mm": float(np.quantile(ae, 0.95)),
        "outlier_gt_0p5_fraction": float(np.mean(ae > 0.5)),
    }


def write_dict_rows(path: Path, rows, fields):
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

    if args.grid_step_mm <= 0:
        raise SystemExit("--grid-step-mm must be > 0")

    records, strip_ids = load_dataset(args.scan_dir)
    folds = make_seed_folds(records, args.folds)

    predictions = []
    for ifold, test_seeds in enumerate(folds):
        train = [r for r in records if r["seed"] not in test_seeds]
        test = [r for r in records if r["seed"] in test_seeds]

        # Baseline: use only the current raw stereo-CoG coordinate, but allow
        # the same data-driven nonlinear calibration machinery.
        base_template = build_template_grid(
            train,
            lambda r: np.array([r["u_raw"]]),
            args.grid_step_mm,
        )

        # Candidate improvement: complete normalized strip-sharing vector.
        strip_template = build_template_grid(
            train,
            lambda r: r["features"],
            args.grid_step_mm,
        )

        for r in test:
            pred_base = predict_from_template(
                np.array([r["u_raw"]]), base_template, use_variance=True
            )
            pred_strip = predict_from_template(
                r["features"], strip_template, use_variance=True
            )
            predictions.append(
                {
                    **r,
                    "fold": ifold,
                    "u_cog_cal_mm": pred_base,
                    "u_strip_template_mm": pred_strip,
                }
            )

    predictions.sort(key=lambda r: (r["u_true"], r["seed"]))

    pred_rows = []
    for r in predictions:
        row = {
            "u_true_mm": r["u_true"],
            "random_seed": r["seed"],
            "event": r["event"],
            "fold": r["fold"],
            "active_wires": r["active_wires"],
            "wire_cog_u_mm": r["wire_cog_u"],
            "u_raw_mm": r["u_raw"],
            "u_cog_cal_mm": r["u_cog_cal_mm"],
            "u_strip_template_mm": r["u_strip_template_mm"],
            "err_cog_cal_mm": r["u_cog_cal_mm"] - r["u_true"],
            "err_strip_template_mm": r["u_strip_template_mm"] - r["u_true"],
        }
        pred_rows.append(row)

    pred_file = args.scan_dir / "c2_strip_estimator_predictions.csv"
    write_dict_rows(pred_file, pred_rows, list(pred_rows[0]))

    regions = {
        "all": lambda r: True,
        "core_abs_u_le_1p5": lambda r: abs(r["u_true"]) <= 1.5 + 1e-12,
        "edge_abs_u_gt_1p5": lambda r: abs(r["u_true"]) > 1.5 + 1e-12,
        "single_wire": lambda r: r["active_wires"] == 1,
        "multiwire": lambda r: r["active_wires"] > 1,
    }

    metric_rows = []
    for region, selector in regions.items():
        subset = [r for r in predictions if selector(r)]
        if not subset:
            continue
        for estimator, key in (
            ("nonlinear_cog_calibration", "u_cog_cal_mm"),
            ("full_strip_template", "u_strip_template_mm"),
        ):
            m = metrics(subset, key)
            metric_rows.append(
                {"region": region, "estimator": estimator, **m}
            )

    metric_file = args.scan_dir / "c2_strip_estimator_metrics.csv"
    write_dict_rows(
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
        subset = [r for r in predictions if r["u_true"] == u]
        for estimator, key in (
            ("nonlinear_cog_calibration", "u_cog_cal_mm"),
            ("full_strip_template", "u_strip_template_mm"),
        ):
            m = metrics(subset, key)
            by_u_rows.append(
                {"u_true_mm": u, "estimator": estimator, **m}
            )

    by_u_file = args.scan_dir / "c2_strip_estimator_by_u.csv"
    write_dict_rows(
        by_u_file,
        by_u_rows,
        [
            "u_true_mm", "estimator", "n", "bias_mm", "rmse_mm", "mae_mm",
            "q68_abs_error_mm", "q95_abs_error_mm",
            "outlier_gt_0p5_fraction",
        ],
    )

    # Finally build a descriptive full-data template table.  This is not used
    # for the cross-validation numbers above; it is a compact calibration
    # object for inspecting how the strip pattern evolves with u0.
    full_template = build_template_grid(
        records, lambda r: r["features"], args.grid_step_mm
    )
    template_rows = []
    nstrip = len(strip_ids)
    for iu, u in enumerate(full_template["u_nodes"]):
        row = {"u_true_mm": float(u)}
        for j, k in enumerate(strip_ids):
            row[f"minus_k{k}_mean_frac"] = float(full_template["means"][iu, j])
            row[f"minus_k{k}_std_frac"] = float(
                math.sqrt(max(full_template["variances"][iu, j], 0.0))
            )
        for j, k in enumerate(strip_ids):
            jj = nstrip + j
            row[f"plus_k{k}_mean_frac"] = float(full_template["means"][iu, jj])
            row[f"plus_k{k}_std_frac"] = float(
                math.sqrt(max(full_template["variances"][iu, jj], 0.0))
            )
        template_rows.append(row)

    template_file = args.scan_dir / "c2_strip_templates.csv"
    write_dict_rows(template_file, template_rows, list(template_rows[0]))

    print("C2 full-strip estimator comparison")
    print(f"  valid events     : {len(records)}")
    print(f"  strip ids        : {strip_ids}")
    print(f"  seed folds       : {len(folds)}")
    print(f"  interpolation    : {args.grid_step_mm:.6g} mm")
    print()
    for row in metric_rows:
        if row["region"] != "all":
            continue
        print(
            f"  {row['estimator']:28s} "
            f"RMSE={row['rmse_mm']:.4f} mm  "
            f"MAE={row['mae_mm']:.4f} mm  "
            f"68%={row['q68_abs_error_mm']:.4f} mm  "
            f"95%={row['q95_abs_error_mm']:.4f} mm"
        )
    print()
    print("Wrote:")
    print(f"  {pred_file}")
    print(f"  {metric_file}")
    print(f"  {by_u_file}")
    print(f"  {template_file}")


if __name__ == "__main__":
    main()

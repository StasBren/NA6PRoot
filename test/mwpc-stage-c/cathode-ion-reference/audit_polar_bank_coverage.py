#!/usr/bin/env python3
"""Read-only audit of *intrinsic 2D* training coverage for the 2100V P3 whole planes.

Compare the frozen 256 positive-v source bank to ONE already examined,
independent weighted muon (0301). Avoid inspecting held-out muons 0408/0446.
No Garfield, new RKF, electronics, FEM, or new avalanche simulation.

This is a geometric/linear-algebra diagnostic, NOT a response accuracy test.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from study_response_rank_and_moments import load_bank, features
from study_event_moment_replay import load_ions
from study_regularized_moment_replay import positive_features

PITCH_CM = .3
WIRE_RADIUS_CM = .0015
FEATURES = ("dipole", "dipole_radial", "quadrupole", "separatrix")


def polar(u, v):
    wire = np.rint(u / PITCH_CM).astype(int)
    local = u - wire * PITCH_CM
    r = np.hypot(local, v)
    if np.any(~np.isfinite(r)) or np.any(r <= WIRE_RADIUS_CM):
        raise ValueError("Invalid birth radius or metal-wire birth")
    return wire, r * 10000., np.arctan2(np.abs(v), local), np.log(r / WIRE_RADIUS_CM)


def quantiles(x, weights):
    idx = np.argsort(x)
    v, w = x[idx], weights[idx]
    c = np.cumsum(w) - .5 * w
    c /= np.sum(w)
    return {f"p{p:g}": float(np.interp(p / 100., c, v))
            for p in (0, 5, 25, 50, 75, 95, 100)}


def closeness(train_xy, event_xy, weights):
    scales = np.maximum(train_xy.std(axis=0), (.02, .08))
    a = train_xy / scales
    b = event_xy / scales
    distance2 = np.sum((a[:, None, :] - a[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(distance2, np.inf)
    tr_knn = np.sqrt(np.min(distance2, axis=1))
    ev_distance2 = np.sum((b[:, None, :] - a[None, :, :]) ** 2, axis=2)
    ev_knn = np.sqrt(np.min(ev_distance2, axis=1))
    q95 = float(np.quantile(tr_knn, .95))
    return {
        "scales_log_r_theta": scales.tolist(),
        "source_leave_one_out_neighbor_distance": quantiles(tr_knn, np.ones(len(tr_knn))),
        "target_nearest_source_neighbor_distance": quantiles(ev_knn, weights),
        "source_95pct_nearest_neighbor_threshold": q95,
        "target_weight_fraction_farther_than_source_95pct":
            float(weights[ev_knn > q95].sum() / weights.sum())
    }


def feature_diagnostics(train_u, train_v, target_u, target_v, weights):
    out = []
    for model in FEATURES:
        A = positive_features(train_u, train_v, model)
        B = positive_features(target_u, target_v, model)
        sigma = np.maximum(np.sqrt(np.mean(A * A, axis=0)), 1e-12)
        X = A / sigma
        Z = B / sigma
        lo, hi = A.min(axis=0), A.max(axis=0)
        margin = .01 * np.maximum(hi - lo, 1e-8)
        outside = np.any((B[:, 1:] < lo[None, 1:] - margin[None, 1:]) |
                         (B[:, 1:] > hi[None, 1:] + margin[None, 1:]), axis=1)
        gram = X.T @ X / len(X)
        cond0 = float(np.linalg.cond(gram))
        regularization = np.eye(len(sigma)) * .001
        regularization[0, 0] = 0
        cond1 = float(np.linalg.cond(gram + regularization))
        out.append(dict(
            model=model,
            n_features=A.shape[1],
            cond_normalized_gram=cond0,
            cond_normalized_gram_with_lambda_1e_minus3=cond1,
            target_weight_fraction_outside_marginal_feature_ranges=float(
                weights[outside].sum() / weights.sum()),
            max_standardized_training_feature=float(np.abs(X).max()),
            max_standardized_target_feature=float(np.abs(Z).max()),
            max_absolute_source_prediction_feature=float(np.abs(A).max()),
            max_absolute_target_prediction_feature=float(np.abs(B).max()),
        ))
    # The alleged "singular separatrix feature" is a smooth bounded gate.
    # Explicitly test max raw gate and radial multipliers separately.
    train_phi = polar(train_u, train_v)[2]
    target_phi = polar(target_u, target_v)[2]
    for phi, label in ((train_phi, "train"), (target_phi, "target")):
        g = np.exp(-.5 * (np.sin(phi) / .15) ** 2)
        if not np.all((g >= 0) & (g <= 1)):
            raise AssertionError("Gaussian separatrix gate is not bounded in [0,1]")
        out.append(dict(model="bounded_gate_" + label,
                        n_features=1,
                        max_absolute_source_prediction_feature=float(np.max(g)),
                        max_absolute_target_prediction_feature=float(np.min(g))))
    return out


def angle_hist(theta, w, bin_count=16):
    bins = np.linspace(0, 180, bin_count + 1)
    hist, _ = np.histogram(np.degrees(theta), bins=bins, weights=w)
    return hist / w.sum()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bank-root", type=Path, required=True)
    parser.add_argument("--ion-file", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    a = parser.parse_args()
    if a.out_dir.exists():
        parser.error("Refusing to overwrite existing analysis directory: " + str(a.out_dir))
    tu, tv, tw, original_primary, waves, dt, n_bins = load_bank(a.bank_root)
    meta, ions, exact = load_ions(a.ion_file)
    if len(tu) != 256 or not np.all(tv > 0) or abs(dt - 10) > 1e-8 or n_bins != 10000:
        raise ValueError("Unexpected 256-site training bank geometry/time grid")
    eu = np.asarray(ions["u_patch_cm"], dtype=float)
    ev = np.asarray(ions["v_absolute_cm"], dtype=float)
    ew = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    if not np.isclose(ew.sum(), exact.sum(), rtol=1e-8):
        raise ValueError("Source weights differ from exact original ion census")
    source_wire, source_radius, source_theta, source_logr = polar(tu, tv)
    event_wire, event_radius, event_theta, event_logr = polar(eu, ev)
    source_xy = np.column_stack((source_logr, source_theta))
    event_xy = np.column_stack((event_logr, event_theta))
    hist_train = angle_hist(source_theta, np.ones(len(tu)))
    hist_event = angle_hist(event_theta, ew)
    train_radial_range = (float(np.min(source_radius)), float(np.max(source_radius)))
    train_angle_range = (float(np.min(source_theta)), float(np.max(source_theta)))
    intrinsic_outside = ((event_radius < train_radial_range[0]) |
                         (event_radius > train_radial_range[1]) |
                         (event_theta < train_angle_range[0]) |
                         (event_theta > train_angle_range[1]))
    wire_rows = []
    for k in sorted(set(source_wire.tolist()) | set(event_wire.tolist())):
        train_n = int(np.sum(source_wire == k))
        source_w = float(ew[event_wire == k].sum())
        wire_rows.append(dict(local_wire_index=k, n_bank_sites=train_n,
                              n_source_sample_rows=int(np.sum(event_wire == k)),
                              source_weighted_ions=source_w,
                              source_weight_fraction=source_w / ew.sum()))
    result = dict(
        scope="Ideal symmetric nine-wire whole-cathode positive-ion P3 response",
        training_bank=str(a.bank_root.resolve()),
        test_ion_archive=str(a.ion_file.resolve()),
        crossing_id=meta["crossing_id"],
        n_bank_sites=len(tu),
        n_bank_original_primaries=int(len(set(original_primary.tolist()))),
        n_test_samples=len(eu),
        n_true_weighted_ions=float(ew.sum()),
        n_test_negative_v=int(np.sum(ev < 0)),
        bank_radius_um=quantiles(source_radius, np.ones(len(tu))),
        target_radius_um=quantiles(event_radius, ew),
        bank_theta_degrees=quantiles(np.degrees(source_theta), np.ones(len(tu))),
        target_theta_degrees=quantiles(np.degrees(event_theta), ew),
        target_fraction_outside_intrinsic_radius_or_theta_range=float(
            ew[intrinsic_outside].sum() / ew.sum()),
        source_vs_bank_2d_nearest_neighbors=closeness(source_xy, event_xy, ew),
        near_saddle_proxy_fraction_sin_theta_lt_015={
            "bank": float(np.mean(np.sin(source_theta) < .15)),
            "source_weighted": float(ew[np.sin(event_theta) < .15].sum() / ew.sum())
        },
        caution=[
            "Training sites come from six original primaries of a single source muon",
            "2D neighbor metrics are diagnostic, not a generalization theorem",
            "No proof of a saddle crossing from sin(theta) alone",
            "Wire-relative folding is NOT exact under finite nine-wire patch edge effects",
            "No unobserved full-ion populations inferred from K256 weighted sample",
            "Do not analyze locked independent events 0408 or 0446 yet"
        ])
    a.out_dir.mkdir(parents=True, exist_ok=False)
    (a.out_dir / "coverage_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with (a.out_dir / "feature_conditioning.csv").open("w", newline="") as f:
        diagnostics = feature_diagnostics(tu, tv, eu, np.abs(ev), ew)
        fields = sorted({key for row in diagnostics for key in row})
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(diagnostics)
    with (a.out_dir / "wire_index.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(wire_rows[0]))
        w.writeheader(); w.writerows(wire_rows)
    with (a.out_dir / "angle_coverage.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("theta_deg_min", "theta_deg_max",
                    "bank_fraction", "source_weight_fraction"))
        for j in range(len(hist_train)):
            w.writerow((180*j/len(hist_train), 180*(j+1)/len(hist_train),
                        hist_train[j], hist_event[j]))
    print("PASS read-only 2D birth-coverage audit")
    print("Source:", result["crossing_id"])
    print("Original bank 256 positive-v; source macroions:", len(eu),
          "negative-v source ions:", int(np.sum(ev < 0)))
    print("Out of intrinsic r/theta min-max range (weighted):",
          round(result["target_fraction_outside_intrinsic_radius_or_theta_range"], 5))
    print("Outside source-bank LOO 95th-percentile nearest-neighbor distance:",
          round(result["source_vs_bank_2d_nearest_neighbors"][
              "target_weight_fraction_farther_than_source_95pct"], 5))
    print("See files:", a.out_dir)


if __name__ == "__main__":
    main()

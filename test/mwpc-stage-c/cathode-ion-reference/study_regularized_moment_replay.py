#!/usr/bin/env python3
"""Gate 3: training-only grouped cross-validation and ridge-stabilized muon moment replay.

Scientific scope: positive-ion 10 ns waveform, *ideal whole-cathode planes*,
symmetric 9-wire 2100V fixed field, provisional Ar+ mobility. Neither actual
U/V strips nor electron prompt nor an avalanche generator are modeled.

Never use the held-out muon to tune ridge penalties or feature definitions.
Train and choose lambda ONLY from the frozen single-ion response bank.
CV groups are original Heed primary IDs (6 groups in known bank).
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from study_response_rank_and_moments import features, load_bank
from study_event_moment_replay import (
    BINS, DT_NS, MODELS, load_ions, check_direct, load_reference,
    load_baseline, deposit_moments, synthesize, metrics
)

LAMBDAS = (0.0, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)


def positive_features(u, v, model):
    if np.any(v <= 0):
        raise ValueError("This training-basis function requires positive v")
    array = features(u, v, model)
    return array[:, : array.shape[1] // 2]


def fit_ridge(F, Y, lam):
    """Relative lambda: normalized feature columns; intercept is unpenalized."""
    sigma = np.maximum(np.sqrt(np.mean(F * F, axis=0)), 1e-12)
    X = F / sigma
    if lam == 0:
        B = np.linalg.pinv(X, rcond=1e-10) @ Y
    else:
        n = len(F)
        gram = X.T @ X / n
        cross = X.T @ Y / n
        penalty = np.eye(F.shape[1]) * lam
        penalty[0, 0] = 0  # monopole/constant remains exact at all penalties
        try:
            B = np.linalg.solve(gram + penalty, cross)
        except np.linalg.LinAlgError as e:
            raise ValueError("Bad ridge linear system") from e
    return sigma, B


def rel_l2_rows(pred, y):
    delta = pred - y
    return np.sqrt(np.sum(delta * delta, axis=1) /
                   np.maximum(np.sum(y * y, axis=1), 1e-60))


def cross_validate(F, Y, primary, model, group_min=2):
    distinct = sorted(set(int(x) for x in primary))
    distinct = [g for g in distinct if np.count_nonzero(primary == g) >= group_min]
    if len(distinct) < 3:
        raise ValueError("Need >=3 populated original-Heed-primary groups for LOPO CV")
    summary = []
    for lam in (LAMBDAS if model != "monopole" else (0.0,)):
        group_scores = []
        row_scores = []
        for g in distinct:
            valid = primary == g
            fit = ~valid
            if np.count_nonzero(fit) <= F.shape[1] + 3:
                raise ValueError("Insufficient training sites in leave-one-primary-out CV")
            scale, B = fit_ridge(F[fit], Y[fit], lam)
            predict = (F[valid] / scale) @ B
            row_errors = rel_l2_rows(predict, Y[valid])
            group_scores.append(float(np.mean(row_errors)))
            row_scores.extend(row_errors.tolist())
        summary.append(dict(
            model=model, ridge_lambda_relative=lam, n_groups=len(distinct),
            cv_equal_primary_mean_L2=float(np.mean(group_scores)),
            cv_site_mean_L2=float(np.mean(row_scores)),
            cv_site_median_L2=float(np.median(row_scores)),
            cv_site_p90_L2=float(np.quantile(row_scores, .9)),
            cv_worst_primary_mean_L2=float(np.max(group_scores)),
            cv_site_max_L2=float(np.max(row_scores)),
        ))
    best = min(summary, key=lambda r: (r["cv_equal_primary_mean_L2"],
                                       r["ridge_lambda_relative"]))
    return summary, best


def diagnostics(F_train, F_event, weights):
    """Very conservative per-column min/max extrapolation proxy, not a proof."""
    lo = F_train.min(axis=0)
    hi = F_train.max(axis=0)
    # Ignore the constant and use a small 1% dynamic-range numerical margin.
    margin = .01 * np.maximum(hi - lo, 1e-8)
    outside = (F_event[:, 1:] < (lo[1:] - margin[1:])) | (
        F_event[:, 1:] > (hi[1:] + margin[1:]))
    mask = np.any(outside, axis=1)
    return dict(
        fraction_rows_outside_training_column_ranges=float(np.mean(mask)),
        fraction_weighted_ions_outside_training_ranges=float(
            weights[mask].sum() / weights.sum()),
        n_columns_with_outside_rows=int(np.sum(np.any(outside, axis=0))),
    )


def analysis(bank_path, source_path, ref_wave, ref_summary, out_dir, baseline=None):
    if out_dir.exists():
        raise ValueError("Output directory exists; refusing overwrite: " + str(out_dir))
    u, v, w, primaries, Y, dt, n_bins = load_bank(bank_path)
    if n_bins != BINS or abs(dt - DT_NS) > 1e-6 or len(u) != 256:
        raise ValueError("Require original frozen 256-site 10ns 100us training bank")
    if not np.all(v > 0):
        raise ValueError("Frozen bank unexpectedly contains non-positive-v starts")
    meta, ions, exact = load_ions(source_path)
    check_direct(meta, ref_summary)
    direct = load_reference(ref_wave)
    event_u = np.asarray(ions["u_patch_cm"], dtype=float)
    event_v = np.asarray(ions["v_absolute_cm"], dtype=float)
    weights = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    shifted = np.asarray(ions["t_birth_ns"], dtype=float) / DT_NS
    starts = np.floor(shifted).astype(int)
    fractions = shifted - starts
    if np.any(starts < 0):
        raise ValueError("Negative birth-time bin")
    rows = [metrics("REFERENCE_DIRECT_RKF", direct, direct)]
    cv_rows = []
    selection_rows = []
    output = {"direct_RKF": direct}
    if baseline is not None:
        old = load_baseline(baseline, direct)
        rows.append(metrics("EXISTING_IDW4_event", old, direct))
        output["existing_IDW4"] = old
    for model in MODELS:
        X = positive_features(u, v, model)
        cv, selected = cross_validate(X, Y, primaries, model)
        cv_rows.extend(cv)
        lam = float(selected["ridge_lambda_relative"])
        scale, B = fit_ridge(X, Y, lam)
        # Compare with the original unregularized model on the same untouched event.
        F = positive_features(event_u, np.abs(event_v), model)
        info = diagnostics(X, F, weights)
        selection_rows.append(dict(
            model=model, selected_lambda_relative=lam,
            selected_cv_equal_primary_mean_L2=selected["cv_equal_primary_mean_L2"],
            **info,
        ))
        pos = deposit_moments(starts, fractions, weights,
                              F / scale, event_v > 0)
        neg = deposit_moments(starts, fractions, weights,
                              F / scale, event_v < 0)
        pred = synthesize(pos, neg, B)
        if not np.all(np.isfinite(pred)):
            raise ValueError("Predicted waveform has nonfinite entries: " + model)
        label = "ridge_CV_" + model
        rows.append(metrics(label, pred, direct))
        output[label] = pred
    provenance = dict(
        scope="Grouped source-primary CV; independent original K256 muon",
        source_muon=meta, input_bank=str(bank_path.resolve()),
        direct_reference=str(ref_wave.resolve()),
        direct_summary=str(ref_summary.resolve()),
        existing_IDW4=None if baseline is None else str(baseline.resolve()),
        n_training_sites=len(u),
        training_primaries=sorted(set(primaries.tolist())),
        n_source_primaries=len(exact),
        source_true_ions=int(exact.sum()),
        source_saved_ions=int(len(event_u)),
        time_bin_ns=DT_NS, window_us=100,
        model_selection="Leave-one-original-primary-out on muon0 only, equal group weight; NEVER target muon",
        limits=[
            "Whole-plane 2100V symmetric noFW analytic model only",
            "No actual U/V strips or electronics",
            "No avalanche generator; true weighted K256 birth positions retained",
            "No prompt electrons or Heed primary ions",
            "Uncalibrated Ar+ ion mobility",
            "Original 256 templates cover v>0 only: v<0 exact ideal reflection",
            "Covariate diagnostic is per-feature range proxy, not mathematical extrapolation guarantee",
            "Training only six original primary groups; LOPO may not represent other muons",
            "Selection lambdas and feature families are provisional engineering tests",
            "Original K256 sampling uncertainty separate from regression error",
            "Windowed 100us charge not complete endpoint for late ions"
        ]
    )
    out_dir.mkdir(parents=True, exist_ok=False)
    for name, data in (("cv_grid.csv", cv_rows),
                       ("selected_models.csv", selection_rows),
                       ("event_metrics.csv", rows)):
        with (out_dir / name).open("w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=list(data[0]))
            wr.writeheader(); wr.writerows(data)
    (out_dir / "metadata.json").write_text(json.dumps(provenance, indent=2) + "\n")
    # Save first 100ns currents, then selected 1000/10000 time grid metrics CSV
    with (out_dir / "first_100ns_waveforms.csv").open("w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["time_center_ns"] + [
            name + "_" + side + "_fC_per_ns"
            for name in output for side in ("plus", "minus")
        ])
        for t in range(10):
            wr.writerow([(t+.5)*DT_NS] + [
                float(output[name][t, c]) for name in output for c in (0,1)
            ])
    print("PASS frozen-bank grouped cross-validation and independent-muon ridge replay")
    print("Source", meta["crossing_id"], "saved", meta["n_macroions"],
          "true ions", meta["n_true_ions"])
    for selection in selection_rows:
        print("SELECTED", selection["model"],
              "lambda", selection["selected_lambda_relative"],
              "CV", round(selection["selected_cv_equal_primary_mean_L2"], 6),
              "weighted_source_outside_train", round(
                selection["fraction_weighted_ions_outside_training_ranges"], 4))
    for row in rows:
        print(f'{row["model"]:26s} L2={row["combined_L2"]:.6f}'
              f' plus={row["plus_L2_floor"]:.6f} minus={row["minus_L2_floor"]:.6f}'
              f' early100_err=({row["plus_first100_error_fC"]:.4f},'
              f'{row["minus_first100_error_fC"]:.4f})fC')
    print("No target-muon metric is used in lambda selection; compare with prior unregularized run.")


def self_test():
    # Known linear data with a rank-deficient/unstable extra feature:
    # Ridge must remain finite; zero-lambda should match pseudoinverse.
    rng = np.random.default_rng(311)
    n = 60
    x = rng.uniform(-1, 1, n)
    primary = np.repeat(np.arange(6), n//6)
    F = np.column_stack((np.ones(n), x, x*x, np.exp(-200*(x-.4)**2)))
    Y = F[:, :3] @ rng.normal(size=(3, 40))
    all_cv, best = cross_validate(F, Y, primary, "test")
    if not np.isfinite(best["cv_equal_primary_mean_L2"]):
        raise AssertionError("LOPO cross-validation non-finite")
    scale, beta = fit_ridge(F, Y, 0.0)
    error = np.max(np.abs((F / scale) @ beta - Y))
    if error > 1e-10:
        raise AssertionError("Unregularized synthetic exact-rank fit incorrect: "+str(error))
    scale, beta = fit_ridge(F, Y, 1e-2)
    if not np.all(np.isfinite(beta)):
        raise AssertionError("Ridge produced nonfinite coefficients")
    print("PASS synthetic grouped CV, pseudoinverse and ridge stability")
    print("SYNTHETIC", json.dumps(best, sort_keys=True))
    print("NOTE: synthetic correctness is not evidence of physics accuracy.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--bank-root", type=Path)
    p.add_argument("--ion-file", type=Path)
    p.add_argument("--reference-waveforms", type=Path)
    p.add_argument("--reference-summary", type=Path)
    p.add_argument("--baseline-waveforms", type=Path)
    p.add_argument("--out-dir", type=Path)
    a = p.parse_args()
    if a.self_test:
        self_test(); return
    for name in ("bank_root", "ion_file", "reference_waveforms",
                 "reference_summary", "out_dir"):
        if getattr(a, name) is None:
            p.error("Required option: --" + name.replace("_", "-"))
    analysis(a.bank_root, a.ion_file, a.reference_waveforms,
             a.reference_summary, a.out_dir, a.baseline_waveforms)


if __name__ == "__main__":
    main()

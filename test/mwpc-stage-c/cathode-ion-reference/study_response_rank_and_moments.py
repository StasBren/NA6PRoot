#!/usr/bin/env python3
"""Read-only single-ion whole-plane response rank and physical-moment pilot.

The input is the *existing* 256-site, 10 ns training_bank ROOT tree. Each row
is a directly simulated single-ion RKF current waveform, NOT a full avalanche.
The operator SVD is an oracle compression test: held-out response coefficients
are obtained from held-out true waveforms. They are not fast predictions.
The physical-feature models instead predict held-out responses from ion
starting coordinates, and their features can be summed into avalanche moments
only when birth times (and any spatial phases) are dealt with consistently.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

R_CM = 0.0015
PITCH_CM = 0.3
E_FC = 1.602176634e-4


def features(u, v, model):
    """Linear physical-coordinate bases. Sum rows to form ion-cloud moments."""
    wire_u = np.rint(u / PITCH_CM) * PITCH_CM
    du = u - wire_u
    r = np.hypot(du, v)
    if np.any(~np.isfinite(r)) or np.any(r <= R_CM):
        raise ValueError("Ion starts must lie outside the 15 um anode radius")
    theta = np.arctan2(v, du)
    x = np.log(r / R_CM)
    c, s = np.cos(theta), np.sin(theta)
    ones = np.ones_like(x)
    cols = [ones]
    if model in ("dipole", "dipole_radial", "quadrupole", "separatrix"):
        cols += [c, s]
    if model in ("dipole_radial", "quadrupole", "separatrix"):
        cols += [x, x * x, x * c, x * s]
    if model in ("quadrupole", "separatrix"):
        cols += [np.cos(2 * theta), np.sin(2 * theta),
                 x * np.cos(2 * theta), x * np.sin(2 * theta)]
    if model == "separatrix":
        # The width is a hypothesis, not a physical separatrix calibration.
        gate = np.exp(-0.5 * (np.sin(theta) / 0.15) ** 2)
        cols += [gate, gate * x, gate * c]
    # Distinct half-plane drift and induction families: allow independent
    # coefficients without inventing symmetry for a one-sided bank.
    base = np.column_stack(cols)
    return np.concatenate((base * (v[:, None] > 0),
                           base * (v[:, None] < 0)), axis=1)


def load_bank(path):
    try:
        import uproot
        import awkward as ak
    except ImportError as exc:
        raise SystemExit("Install python3 -m pip install --user numpy uproot awkward") from exc
    with uproot.open(path) as root:
        for key in ("training_bank", "average_plus", "average_minus"):
            if key not in root:
                raise ValueError("Missing ROOT object: " + key)
        tree = root["training_bank"]
        required = ["u_cm", "v_cm", "w_cm", "original_primary",
                    "waveform_plus_then_minus"]
        for key in required:
            if key not in tree:
                raise ValueError("Missing training_bank branch: " + key)
        u = np.asarray(tree["u_cm"].array(library="np"), dtype=np.float64)
        v = np.asarray(tree["v_cm"].array(library="np"), dtype=np.float64)
        w = np.asarray(tree["w_cm"].array(library="np"), dtype=np.float64)
        primary = np.asarray(tree["original_primary"].array(library="np"), dtype=int)
        raw = tree["waveform_plus_then_minus"].array(library="ak")
        waves = np.asarray(ak.to_numpy(ak.to_regular(raw)), dtype=np.float64)
        edges = root["average_plus"].axis().edges()
        minus_edges = root["average_minus"].axis().edges()
    n_bins = len(edges) - 1
    if waves.ndim != 2 or waves.shape != (len(u), 2 * n_bins):
        raise ValueError("Not a two-plane equal-length waveform bank")
    if len(u) < 32 or not np.all(np.diff(edges) > 0):
        raise ValueError("Need >=32 templates and positive time-bin widths")
    if not np.allclose(np.diff(edges), edges[1] - edges[0], rtol=0, atol=1e-7):
        raise ValueError("Nonuniform time grid")
    if not np.array_equal(edges, minus_edges):
        raise ValueError("Plane time axes differ")
    if np.any(~np.isfinite(waves)) or np.any(~np.isfinite(u + v + w)):
        raise ValueError("Non-finite archive values")
    if np.any(v == 0):
        raise ValueError("Ion born on v=0 separatrix; explicit handling required")
    # In the validated whole-plane analytic field, w invariance is exact.
    return u, v, w, primary, waves, float(edges[1] - edges[0]), n_bins


def split_indices(v, fraction, seed):
    rng = np.random.default_rng(seed)
    train, test = [], []
    for mask in (np.flatnonzero(v > 0), np.flatnonzero(v < 0)):
        if len(mask) == 0:
            continue
        if len(mask) < 5:
            raise ValueError("Insufficient training/held-out samples for a populated hemisphere")
        mask = rng.permutation(mask)
        n_test = max(2, min(len(mask) - 3, int(round(fraction * len(mask)))))
        test.extend(mask[:n_test])
        train.extend(mask[n_test:])
    train = np.asarray(sorted(train), dtype=int)
    test = np.asarray(sorted(test), dtype=int)
    if len(test) < 8 or len(train) < 24:
        raise ValueError("Insufficient train/test data")
    return train, test


def low_rank(Y, rank_cap):
    # Uncentered response-operator SVD: zero ion population -> zero signal.
    # Avoid an enormous (2*N_time) square SVD.
    gram = Y @ Y.T
    eig, U = np.linalg.eigh(gram)
    idx = np.argsort(eig)[::-1]
    eig = np.maximum(eig[idx], 0.0)
    scale = float(eig[0])
    valid = np.flatnonzero(eig > scale * 1e-12)
    if scale <= 0 or len(valid) == 0:
        raise ValueError("Zero or numerically degenerate response operator")
    k = min(rank_cap, len(valid))
    V = (Y.T @ U[:, idx[:k]]) / np.sqrt(eig[:k])[None, :]
    return eig, V


def idw_predict(u, v, Y, train, test, neighbors=4):
    r = np.hypot(u - np.rint(u / PITCH_CM) * PITCH_CM, v)
    # Exact old-pilot 30 um u/v and log(r) metric, updated to local-wire radius.
    pred = np.zeros((len(test), Y.shape[1]), dtype=np.float64)
    for j, source in enumerate(test):
        ok = train[(v[train] > 0) == (v[source] > 0)]
        if len(ok) < neighbors:
            raise ValueError("Too few training sites in the same hemisphere")
        d = ((u[ok] - u[source]) / 0.003) ** 2
        d += ((v[ok] - v[source]) / 0.003) ** 2
        d += (np.log(r[ok] / r[source]) / 0.4) ** 2
        nearest = np.argsort(d)[:neighbors]
        if d[nearest[0]] < 1e-14:
            weights = np.zeros(len(nearest)); weights[0] = 1
        else:
            weights = 1 / np.maximum(d[nearest], 1e-14)
            weights /= weights.sum()
        pred[j] = weights @ Y[ok[nearest]]
    return pred


def evaluate(truth, pred, n_bins, dt_ns):
    delta = pred - truth
    denom = np.sum(truth ** 2, axis=1)
    err = np.sqrt(np.sum(delta ** 2, axis=1) / np.maximum(denom, 1e-60))
    # Weak-plane floors prevent ratios of near-zero bipolar net signals
    # masquerading as a catastrophic physical relative error.
    plane_err = []
    for k in (0, 1):
        d = delta[:, k * n_bins:(k + 1) * n_bins]
        y = truth[:, k * n_bins:(k + 1) * n_bins]
        floor = np.maximum(np.sum(truth ** 2, axis=1) * 1e-4, 1e-60)
        plane_err.append(np.sqrt(np.sum(d * d, axis=1) /
                                 np.maximum(np.sum(y * y, axis=1), floor)))
    first = min(n_bins, max(1, int(round(100.0 / dt_ns))))
    early = np.stack([np.abs(np.sum(delta[:, k * n_bins:k * n_bins + first],
                                    axis=1) * dt_ns) for k in (0, 1)], axis=1)
    recorded = np.stack([np.abs(np.sum(delta[:, k * n_bins:(k + 1) * n_bins],
                                       axis=1) * dt_ns) for k in (0, 1)], axis=1)
    return dict(mean_l2=float(err.mean()), median_l2=float(np.median(err)),
                p90_l2=float(np.quantile(err, .90)),
                max_l2=float(err.max()),
                max_plus_l2_floor=float(plane_err[0].max()),
                max_minus_l2_floor=float(plane_err[1].max()),
                median_abs_early100_fC=float(np.median(early)),
                max_abs_recorded_charge_fC=float(recorded.max()))


def run(u, v, waves, dt_ns, n_bins, fraction, seed, max_rank):
    train, test = split_indices(v, fraction, seed)
    T, V = waves[train], waves[test]
    eig, basis = low_rank(T, max_rank)
    results = []
    ranks = [r for r in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64)
             if r <= basis.shape[1]]
    if ranks[-1] != basis.shape[1]:
        ranks.append(basis.shape[1])
    def add(label, pred, dimension, oracle=False):
        r = evaluate(V, pred, n_bins, dt_ns)
        r.update(model=label, n_coeff=dimension, oracle=oracle)
        results.append(r)

    add("IDW4_coordinate_prediction", idw_predict(u, v, waves, train, test), 4)
    for rank in ranks:
        B = basis[:, :rank]
        add("SVD_ORACLE_heldout_waveform", (V @ B) @ B.T, rank, True)
    for name in ("monopole", "dipole", "dipole_radial", "quadrupole", "separatrix"):
        F = features(u, v, name)
        # Predict from birth position. For an event with common birth time,
        # sums of the input physical features are sufficient for this model.
        Ft, Fv = F[train], F[test]
        # Scale nonzero features to avoid arbitrary radial units conditioning.
        sigma = np.maximum(np.sqrt(np.mean(Ft * Ft, axis=0)), 1e-12)
        fit = np.linalg.pinv(Ft / sigma, rcond=1e-10) @ T
        pred = (Fv / sigma) @ fit
        add("physical_"+name+"_coordinate_prediction", pred,
            int(np.linalg.matrix_rank(Ft)))
    spec = []
    for i, val in enumerate(eig):
        spec.append(dict(mode=i + 1, eigenvalue=float(val),
                         cumulative_energy=float(np.sum(eig[:i + 1]) /
                                                 np.sum(eig))))
    return results, spec, len(train), len(test)


def synthetic_selftest():
    rng = np.random.default_rng(42)
    n, bins = 96, 60
    theta = np.concatenate((rng.uniform(.4, 2.74, n // 2),
                            rng.uniform(3.54, 5.88, n // 2)))
    r = rng.uniform(.002, .004, n)
    u, v = r * np.cos(theta), r * np.sin(theta)
    t = np.arange(bins) * 10.0
    pulse0 = np.exp(-t / 70)
    pulse1 = np.exp(-t / 200) * (1 - t / 250)
    pulse2 = np.exp(-t / 30) * (1 + t / 50)
    base = np.stack([np.r_[pulse0, .7 * pulse0],
                     np.r_[pulse1, -.5 * pulse1],
                     np.r_[pulse2, .1 * pulse2]])
    W = (np.column_stack([np.ones(n), np.cos(theta), np.sin(theta)]) @ base)
    results, _, _, _ = run(u, v, W, 10.0, bins, .25, 123, 12)
    oracle = next(x for x in results if x["model"] == "SVD_ORACLE_heldout_waveform"
                  and x["n_coeff"] == 3)
    dipole = next(x for x in results if x["model"] == "physical_dipole_coordinate_prediction")
    if oracle["max_l2"] > 1e-7 or dipole["max_l2"] > 1e-7:
        raise AssertionError(f"Low-rank/moment synthetic recovery failed: {oracle} {dipole}")
    print("PASS synthetic rank-3 recovery and dipole moment prediction")
    print(json.dumps(dict(oracle=oracle, dipole=dipole), indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bank-root", type=Path, help="existing 10ns ..._training.root")
    p.add_argument("--out-dir", type=Path, help="new output directory (must not exist)")
    p.add_argument("--test-fraction", type=float, default=.25)
    p.add_argument("--seed", type=int, default=210021)
    p.add_argument("--max-rank", type=int, default=32)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        synthetic_selftest()
        return
    if not a.bank_root or not a.out_dir:
        p.error("--bank-root and --out-dir are required")
    if not 0.1 <= a.test_fraction <= .45 or a.max_rank < 3:
        p.error("invalid test fraction or maximum rank")
    if a.out_dir.exists():
        p.error("Output directory already exists; refusing overwrite: " + str(a.out_dir))
    u, v, w, primary, wave, dt, bins = load_bank(a.bank_root)
    if abs(dt - 10.0) > 1e-6:
        p.error("This stage is 10ns-only; refuse old 100ns response bank")
    if bins != 10000:
        p.error("Expected complete 100us = 10000 bins; do not mix windows")
    rows, spectrum, n_train, n_test = run(
        u, v, wave, dt, bins, a.test_fraction, a.seed, a.max_rank)
    # Protect data: the only writes are into a completely new directory.
    a.out_dir.mkdir(parents=True, exist_ok=False)
    for name, data in (("results.csv", rows), ("spectrum.csv", spectrum)):
        with (a.out_dir / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    meta = dict(scope="10ns whole-plane positive-ion single-start response only",
                source=str(a.bank_root.resolve()),
                entries=len(u), n_train=n_train, n_test=n_test,
                n_primaries_in_bank=int(len(set(primary.tolist()))),
                time_bins=bins, dt_ns=dt, n_plus=int(np.sum(v > 0)),
                n_minus=int(np.sum(v < 0)), seed=a.seed,
                test_fraction=a.test_fraction,
                study_limits=["No physical U/V strips", "No electron contribution",
                              "No avalanche-level sufficiency verified",
                              "SVD_ORACLE uses held-out true waveforms",
                              "No new ion transport or microscopic avalanche",
                              "Ion mobility remains provisional",
                              "All rows are from existing original training bank",
                              "Independent-muon 4.81 percent result is a separate test"])
    (a.out_dir / "metadata.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("PASS read-only rank/moment pilot; outputs:", a.out_dir)
    print("Source:", a.bank_root, "train:", n_train, "held-out sites:", n_test)
    for row in rows:
        if row["model"].startswith("physical_") or (
            row["model"] == "SVD_ORACLE_heldout_waveform"
            and row["n_coeff"] in (2, 4, 8, 16, 32)):
            print(f'{row["model"]:43s} k={row["n_coeff"]:2d}'
                  f' median_L2={row["median_l2"]:.5f}'
                  f' p90_L2={row["p90_l2"]:.5f}'
                  f' worst_L2={row["max_l2"]:.5f}')
    print("NOTE: SVD oracle is a compression bound, NOT a deployable fast model.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Read-only held-out full-muon prediction using weighted avalanche-ion moments.

Fit a physical-coordinate single-ion response operator on the frozen 256-site
10 ns training bank from muon0. Sum the *weighted real ion-birth coordinates*
of independent muons into per-time-bin moment histograms, convolve against
the learned response, and compare with existing direct-RKF whole-plane
waveforms of the identical original sampled ion populations.

No new Garfield/RKF/FEM, no microscopic avalanche synthesis, no electronics,
no physical U/V strips. This tests *response compression given true ion births*
and is NOT a generative avalanche model.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np

from study_response_rank_and_moments import features, load_bank

MODELS = ("monopole", "dipole", "dipole_radial", "quadrupole", "separatrix")
BINS = 10000
DT_NS = 10.0


def only_row(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle))
    if len(records) != 1:
        raise ValueError(f"Expected exactly one CSV row in {path}")
    return records[0]


def load_reference(path):
    with Path(path).open(encoding="utf-8") as f:
        header = f.readline().strip()
    if header != ("time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
                  "cumulative_plus_fC,cumulative_minus_fC"):
        raise ValueError("Unexpected direct-RKF waveform CSV header")
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    if a.shape != (BINS, 5):
        raise ValueError("Expected precisely 10000 x 5 direct waveform values")
    if not np.allclose(a[:, 0], (np.arange(BINS) + .5) * DT_NS, rtol=0, atol=1e-6):
        raise ValueError("Direct waveform time axis is not 10 ns, 100 us")
    if not np.all(np.isfinite(a)):
        raise ValueError("Nonfinite direct waveform")
    return a[:, 1:3]


def load_baseline(path, direct):
    with Path(path).open(encoding="utf-8") as f:
        header = f.readline().strip()
    if not header.startswith("time_center_ns,plus_kernel_fC_per_ns,minus_kernel_fC_per_ns,"):
        raise ValueError("Unexpected prior IDW baseline CSV header")
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    if a.shape != (BINS, 9) or not np.allclose(
            a[:, 0], (np.arange(BINS) + .5) * DT_NS, rtol=0, atol=1e-6):
        raise ValueError("Bad prior IDW baseline dimensions/time")
    if not np.allclose(a[:, 3:5], direct, rtol=3e-6, atol=2e-10):
        raise ValueError("Prior baseline and direct RKF are not the same event")
    return a[:, 1:3]


def root_scalar(tree, key):
    return tree[key].array(library="np")[0].item() if isinstance(
        tree[key].array(library="np")[0], np.generic
    ) else tree[key].array(library="np")[0]


def load_ions(path):
    try:
        import uproot
    except ImportError as exc:
        raise SystemExit("Install 'uproot awkward' into the study venv") from exc
    with uproot.open(path) as f:
        mu, pri, ion = (f[k] for k in ("muons", "primaries", "ions"))
        meta = {
            "crossing_id": str(root_scalar(mu, "crossing_id")),
            "source_channel": str(root_scalar(mu, "source_channel")),
            "station": int(root_scalar(mu, "station")),
            "det_id": int(root_scalar(mu, "det_id")),
            "storage_mode": str(root_scalar(mu, "storage_mode")),
            "reservoir_per_primary": int(root_scalar(mu, "reservoir_per_primary")),
        }
        if meta["storage_mode"] != "weighted":
            raise ValueError("Only original per-primary weighted ion archives are allowed")
        names = ("primary_index", "u_patch_cm", "v_absolute_cm",
                 "w_relative_cm", "t_birth_ns", "inverse_inclusion_weight")
        ions = {k: np.asarray(ion[k].array(library="np")) for k in names}
        exact = np.asarray(pri["ion_births_exact"].array(library="np"), dtype=np.int64)
        saved = np.asarray(pri["ion_births_saved"].array(library="np"), dtype=np.int64)
        good = np.asarray(pri["successful"].array(library="np"), dtype=int)
    n = len(ions["primary_index"])
    if n == 0 or len(exact) == 0 or len(good) != len(exact) or np.any(good != 1):
        raise ValueError("Missing/unsuccessful microscopic primaries")
    p = np.asarray(ions["primary_index"], dtype=np.int64)
    w = np.asarray(ions["inverse_inclusion_weight"], dtype=np.float64)
    if np.any((p < 0) | (p >= len(exact))) or np.any(w <= 0) or not np.all(np.isfinite(w)):
        raise ValueError("Invalid ion primary assignment or sampling weight")
    counts = np.bincount(p, minlength=len(exact))
    totals = np.bincount(p, weights=w, minlength=len(exact))
    if not np.array_equal(counts, saved) or not np.array_equal(saved, np.minimum(
            exact, meta["reservoir_per_primary"])):
        raise ValueError("Stored sample census is inconsistent with primaries")
    if not np.allclose(totals, exact, rtol=1e-8, atol=1e-5):
        raise ValueError("Inverse-inclusion ion weights do not close exact birth counts")
    for name in names:
        if not np.all(np.isfinite(ions[name])):
            raise ValueError("Nonfinite source ions: " + name)
    if np.any(ions["t_birth_ns"] < 0) or np.any(ions["v_absolute_cm"] == 0):
        raise ValueError("Negative birth time or unsupported v=0")
    if np.any(np.abs(ions["u_patch_cm"]) >= 1.35) or np.any(np.abs(ions["v_absolute_cm"]) >= .3):
        raise ValueError("Ion starts outside original 9-wire analytic field")
    meta.update(n_primaries=len(exact), n_true_ions=int(exact.sum()),
                n_macroions=n, n_plus=int(np.sum(ions["v_absolute_cm"] > 0)),
                n_minus=int(np.sum(ions["v_absolute_cm"] < 0)))
    return meta, ions, exact


def check_direct(meta, path):
    record = only_row(path)
    expectations = dict(crossing_id=meta["crossing_id"],
                        archive_mode="weighted", response_mode="direct",
                        true_ion_births=str(meta["n_true_ions"]),
                        retained_ion_rows=str(meta["n_macroions"]),
                        station=str(meta["station"]),
                        det_id=str(meta["det_id"]))
    for key, expected in expectations.items():
        if record.get(key) != expected:
            raise ValueError(f"Reference {key} differs: got {record.get(key)!r}, expected {expected!r}")
    if not np.isclose(float(record["dt_ns"]), DT_NS) or not np.isclose(
            float(record["window_us"]), 100.0):
        raise ValueError("Reference summary does not use 10ns/100us")
    return record


def train_physical_models(bank_path):
    u, v, w, primary, wave, dt, bins = load_bank(bank_path)
    if bins != BINS or abs(dt - DT_NS) > 1e-7:
        raise ValueError("Frozen operator bank must be exactly 10ns / 100us")
    # Training archive is positive-v only. This is enforced deliberately:
    # negative-v extension is by the *validated whole-plane* mirror identity.
    if not np.all(v > 0):
        raise ValueError("Expected frozen single-positive-hemisphere 256-site training bank")
    if len(u) != 256:
        raise ValueError("Expected exact frozen 256 template sites")
    coeffs = {}
    for name in MODELS:
        F = features(u, v, name)
        n_feat = F.shape[1] // 2
        F = F[:, :n_feat]  # Only the positive-v columns are populated.
        scale = np.maximum(np.sqrt(np.mean(F * F, axis=0)), 1e-12)
        beta = np.linalg.pinv(F / scale, rcond=1e-10) @ wave
        coeffs[name] = (scale, beta)
    return coeffs


def deposit_moments(birth_bins, frac, weight, feats, mask, bins=BINS):
    """Fractional-bin birth-time shift, identical to original C++ shift/add."""
    moments = np.zeros((bins, feats.shape[1]), dtype=np.float64)
    start = np.flatnonzero(mask & (birth_bins < bins))
    if len(start):
        np.add.at(moments, birth_bins[start],
                  (weight[start] * (1.0 - frac[start]))[:, None] * feats[start])
    later = np.flatnonzero(mask & (birth_bins + 1 < bins))
    if len(later):
        np.add.at(moments, birth_bins[later] + 1,
                  (weight[later] * frac[later])[:, None] * feats[later])
    return moments


def synthesize(moments_pos, moments_neg, beta, n_bins=BINS):
    """Linear moment accumulation and signal convolution on two planes."""
    nfft = 1 << (2 * n_bins - 1).bit_length()
    up = np.zeros(n_bins)
    down = np.zeros(n_bins)
    for side, moments in ((1, moments_pos), (-1, moments_neg)):
        freq = np.fft.rfft(moments, n=nfft, axis=0)
        for c in (0, 1):
            original_c = c if side == 1 else 1 - c
            response = beta[:, original_c * n_bins:(original_c + 1) * n_bins]
            filt = np.fft.rfft(response, n=nfft, axis=1)
            predicted = np.fft.irfft(np.sum(freq * filt.T, axis=1), n=nfft)[:n_bins]
            if c == 0:
                up += predicted
            else:
                down += predicted
    return np.column_stack((up, down))


def metrics(name, predicted, reference):
    d = predicted - reference
    total = float(np.sum(reference ** 2))
    per = []
    for c in (0, 1):
        denom = max(float(np.sum(reference[:, c] ** 2)), 1e-4 * total, 1e-50)
        per.append(float(np.sqrt(np.sum(d[:, c] ** 2) / denom)))
    hundred = int(100 / DT_NS)
    q0 = DT_NS * np.sum(reference[:hundred], axis=0)
    dq0 = DT_NS * np.sum(d[:hundred], axis=0)
    qwin = DT_NS * np.sum(reference, axis=0)
    dqwin = DT_NS * np.sum(d, axis=0)
    return dict(model=name,
                combined_L2=float(np.sqrt(np.sum(d * d) / max(total, 1e-50))),
                plus_L2_floor=per[0], minus_L2_floor=per[1],
                plus_first100_direct_fC=float(q0[0]),
                minus_first100_direct_fC=float(q0[1]),
                plus_first100_error_fC=float(dq0[0]),
                minus_first100_error_fC=float(dq0[1]),
                plus_window100us_direct_fC=float(qwin[0]),
                minus_window100us_direct_fC=float(qwin[1]),
                plus_window_error_fC=float(dqwin[0]),
                minus_window_error_fC=float(dqwin[1]))


def source_descriptors(ions, exact):
    p = np.asarray(ions["primary_index"], dtype=int)
    weight = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    u = np.asarray(ions["u_patch_cm"], dtype=float)
    v = np.asarray(ions["v_absolute_cm"], dtype=float)
    r = np.hypot(u - np.rint(u / .3) * .3, v)
    theta = np.arctan2(v, u - np.rint(u / .3) * .3)
    t = np.asarray(ions["t_birth_ns"], dtype=float)
    w = np.asarray(ions["w_relative_cm"], dtype=float)
    out = []
    for j, n in enumerate(exact):
        mask = p == j
        if n == 0:
            out.append(dict(primary=j, exact_ion_births=0, n_samples=0))
            continue
        q = weight[mask]; z = float(n)
        mean_w = float(np.sum(q * w[mask]) / z)
        mean_t = float(np.sum(q * t[mask]) / z)
        mean_r = float(np.sum(q * r[mask]) / z)
        sep = np.abs(np.sin(theta[mask])) < .15
        out.append(dict(primary=j, exact_ion_births=int(n),
                        n_samples=int(np.sum(mask)), mean_birth_t_ns=mean_t,
                        mean_w_cm=mean_w,
                        sd_w_cm=float(np.sqrt(np.sum(q * (w[mask] - mean_w) ** 2) / z)),
                        mean_r_cm=mean_r,
                        cos1=float(np.sum(q * np.cos(theta[mask]))),
                        sin1=float(np.sum(q * np.sin(theta[mask]))),
                        cos2=float(np.sum(q * np.cos(2 * theta[mask]))),
                        sin2=float(np.sum(q * np.sin(2 * theta[mask]))),
                        separatrix_sector_fraction=float(np.sum(q * sep) / z)))
    return out


def run_event(bank_path, ion_path, direct_path, direct_summary, out_dir,
              baseline=None, do_write=True):
    if out_dir.exists() and do_write:
        raise ValueError(f"Output directory already exists: {out_dir}")
    meta, ions, exact = load_ions(ion_path)
    check_direct(meta, direct_summary)
    direct = load_reference(direct_path)
    models = train_physical_models(bank_path)
    p = np.asarray(ions["primary_index"], dtype=int)
    u = np.asarray(ions["u_patch_cm"], dtype=float)
    v = np.asarray(ions["v_absolute_cm"], dtype=float)
    wt = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    starts = np.asarray(ions["t_birth_ns"], dtype=float) / DT_NS
    bins = np.floor(starts).astype(int)
    frac = starts - bins
    if np.max(bins) >= BINS:
        meta["late_ion_birth_rows_outside_window"] = int(np.sum(bins >= BINS))
    else:
        meta["late_ion_birth_rows_outside_window"] = 0
    results = [metrics("REFERENCE_DIRECT_RKF", direct, direct)]
    output = {"direct_RKF": direct}
    if baseline is not None:
        b = load_baseline(baseline, direct)
        results.append(metrics("EXISTING_IDW4_event", b, direct))
        output["existing_IDW4"] = b
    for name, (scale, beta) in models.items():
        # u is unchanged under the exact v reflection of this *ideal* model;
        # swapping the +/- weighting electrodes handles negative v.
        raw = features(u, np.abs(v), name)
        feat = raw[:, :len(scale)] / scale
        pos = deposit_moments(bins, frac, wt, feat, v > 0)
        neg = deposit_moments(bins, frac, wt, feat, v < 0)
        y = synthesize(pos, neg, beta)
        if not np.all(np.isfinite(y)):
            raise ValueError("Nonfinite reconstructed waveform: " + name)
        results.append(metrics("moments_" + name, y, direct))
        output[name] = y
    meta.update(
        training_bank=str(bank_path.resolve()), test_ion_archive=str(ion_path.resolve()),
        reference_csv=str(direct_path.resolve()), reference_summary=str(direct_summary.resolve()),
        baseline_csv=None if baseline is None else str(baseline.resolve()),
        dt_ns=DT_NS, window_us=100.0,
        prediction="physical moment operators fit to muon0 independent ion sites",
        phase="fixed-field positive-ion response; NO new avalanche generation",
        symmetry="v mirror and whole-plane swapping ONLY for symmetric 3+3mm noFW",
        limits=["Sampled K256/K512 ions, not full ion cloud",
                "No electron prompt; no Heed primary-ion current",
                "No real stereo U/V fields or VMM frontend",
                "Late arrivals beyond 100us are truncated consistently with direct reference",
                "Different source muons; same provisional mobility and analytic field",
                "Physical feature model trained on 256 starts from one reference muon",
                "No direct performance or compression of avalanche generation established"])
    if not do_write:
        return results, meta
    out_dir.mkdir(parents=True, exist_ok=False)
    with (out_dir / "event_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader(); writer.writerows(results)
    descriptors = source_descriptors(ions, exact)
    with (out_dir / "per_primary_sampled_descriptors.csv").open("w", newline="") as f:
        keys = list(next((d for d in descriptors if d["exact_ion_births"] > 0),
                         descriptors[0]))
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for record in descriptors:
            writer.writerow(record)
    with (out_dir / "waveforms.csv").open("w", newline="") as f:
        keys = ["time_center_ns"] + [
            label + "_" + plane + "_fC_per_ns"
            for label in output for plane in ("plus", "minus")]
        writer = csv.writer(f); writer.writerow(keys)
        for t in range(BINS):
            writer.writerow([(t + .5) * DT_NS] + [
                float(output[label][t, c]) for label in output for c in (0, 1)])
    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"PASS independent-muon moments replay {meta['crossing_id']}")
    print("Source K rows:", meta["n_macroions"], "total weighted ions:",
          meta["n_true_ions"], "plus/minus:", meta["n_plus"], meta["n_minus"])
    for r in results:
        print(f'{r["model"]:25s} L2={r["combined_L2"]:.5f}'
              f' +/-=({r["plus_L2_floor"]:.5f}, {r["minus_L2_floor"]:.5f})'
              f' first100 err fC=({r["plus_first100_error_fC"]:.4g},'
              f' {r["minus_first100_error_fC"]:.4g})')
    print("NOTE: separate ion-waveform compression from new avalanche synthesis.")
    return results, meta


def self_test():
    rng = np.random.default_rng(7)
    n, n_feat, n_ions = 64, 4, 20
    response = rng.normal(size=(n_feat, 2 * n))
    feat = rng.normal(size=(n_ions, n_feat))
    time = rng.uniform(0, 13, n_ions)
    b = np.floor(time).astype(int)
    frac = time - b
    wt = rng.uniform(.3, 2, n_ions)
    v = np.where(np.arange(n_ions) % 2, 1, -1)
    pos = deposit_moments(b, frac, wt, feat, v > 0, n)
    neg = deposit_moments(b, frac, wt, feat, v < 0, n)
    fast = synthesize(pos, neg, response, n)
    slow = np.zeros((n, 2))
    for j in range(n_ions):
        ion = wt[j] * (feat[j] @ response).reshape(2, n)
        for channel in (0, 1):
            src = channel if v[j] > 0 else 1 - channel
            for t in range(n):
                a = b[j] + t
                if a < n:
                    slow[a, channel] += (1 - frac[j]) * ion[src, t]
                if a + 1 < n:
                    slow[a + 1, channel] += frac[j] * ion[src, t]
    err = np.max(np.abs(fast - slow))
    if err > 1e-9:
        raise AssertionError(f"Fractional time shift/hemisphere FFT failed: {err}")
    print("PASS synthetic exact moment summation, v mirror and fractional 10 ns shifts:",
          err)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--bank-root", type=Path)
    parser.add_argument("--ion-file", type=Path)
    parser.add_argument("--reference-waveforms", type=Path)
    parser.add_argument("--reference-summary", type=Path)
    parser.add_argument("--baseline-waveforms", type=Path)
    parser.add_argument("--out-dir", type=Path)
    a = parser.parse_args()
    if a.self_test:
        self_test()
        return
    required = ("bank_root", "ion_file", "reference_waveforms",
                "reference_summary", "out_dir")
    for k in required:
        if getattr(a, k) is None:
            parser.error("--" + k.replace("_", "-") + " is required")
    run_event(a.bank_root, a.ion_file, a.reference_waveforms,
              a.reference_summary, a.out_dir, a.baseline_waveforms)


if __name__ == "__main__":
    main()

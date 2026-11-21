#!/usr/bin/env python3
"""Gate 4: direct bilinear interpolation of deterministic ion RKF waveforms.

Read a new deterministic polar bank for the SAME symmetric ideal whole-plane
P3 field and an EXISTING independent muon K256 reservoir + direct reference.
No new microscopic avalanche, ion drift, FEM, electrons, or stereo strips
are executed here. No SVD, ridge or other response compression is applied.
The goal is spatial-coverage error alone (including the 10ns time-shift rule).
"""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
import numpy as np

from study_event_moment_replay import (
    BINS, DT_NS, load_ions, check_direct, load_reference, load_baseline,
    metrics
)

NFFT = 1 << (2 * BINS - 1).bit_length()
R_WIRE_UM = 15.


def load_grid(filename):
    try:
        import uproot
        import awkward as ak
    except ImportError as e:
        raise SystemExit("This analysis requires uproot and awkward") from e
    with uproot.open(filename) as f:
        if "polar_bank" not in f or "scope" not in f:
            raise ValueError("Not a validated deterministic polar RKF bank")
        scope = str(f["scope"].title)
        if ("Two ideal whole cathode planes" not in scope or
                "positive-v" not in scope):
            raise ValueError("Unexpected grid scope; refuse non-whole-plane bank")
        t = f["polar_bank"]
        keys = ("i_radius", "i_angle", "r_um", "theta_deg", "u_cm", "v_cm",
                "late_after_100us", "end_time_ns",
                "waveform_plus_then_minus")
        for key in keys:
            if key not in t:
                raise ValueError("Missing grid branch: " + key)
        ia = np.asarray(t["i_angle"].array(library="np"), dtype=int)
        ir = np.asarray(t["i_radius"].array(library="np"), dtype=int)
        radii = np.asarray(t["r_um"].array(library="np"), dtype=float)
        angles = np.asarray(t["theta_deg"].array(library="np"), dtype=float)
        u = np.asarray(t["u_cm"].array(library="np"), dtype=float)
        v = np.asarray(t["v_cm"].array(library="np"), dtype=float)
        late = np.asarray(t["late_after_100us"].array(library="np"), dtype=bool)
        end_time = np.asarray(t["end_time_ns"].array(library="np"), dtype=float)
        raw = t["waveform_plus_then_minus"].array(library="ak")
        waves = np.asarray(ak.to_numpy(ak.to_regular(raw)), dtype=float)
    nr, nt = int(np.max(ir)) + 1, int(np.max(ia)) + 1
    if nr < 3 or nt < 6 or nr * nt > 2048 or nr * nt != len(ia):
        raise ValueError("Grid incomplete or exceeds allowed dimensions")
    if not np.array_equal(ir, np.repeat(np.arange(nr), nt)) or not np.array_equal(
            ia, np.tile(np.arange(nt), nr)):
        raise ValueError("Grid rows out of original radius-major order")
    rs = radii.reshape(nr, nt)[:, 0]
    ts = angles.reshape(nr, nt)[0, :]
    if not np.allclose(radii.reshape(nr, nt), rs[:, None], atol=1e-8) or not np.allclose(
            angles.reshape(nr, nt), ts[None, :], atol=1e-8):
        raise ValueError("Grid is not structured in (radius, theta)")
    if np.any(np.diff(rs) <= 0) or np.any(np.diff(ts) <= 0):
        raise ValueError("Grid axes are not strictly increasing")
    if not (rs[0] > R_WIRE_UM and ts[0] > 0 and ts[-1] < 180):
        raise ValueError("Grid includes wire metal or exact separatrix")
    reconstructed_r = np.hypot(u, v) * 10000
    reconstructed_angle = np.degrees(np.arctan2(v, u))
    if not np.allclose(reconstructed_r, radii, atol=1e-8) or not np.allclose(
            reconstructed_angle, angles, atol=1e-8):
        raise ValueError("Cartesian and polar start coordinates inconsistent")
    if waves.shape != (nr * nt, 2 * BINS) or not np.all(np.isfinite(waves)):
        raise ValueError("Expected 2 x 10000 finite 10ns currents for every node")
    if np.any(~np.isfinite(end_time)):
        raise ValueError("Nonfinite ion arrival time")
    return dict(nradius=nr, nangle=nt, radius_um=rs, angle_deg=ts,
                waves=waves, late_rows=int(np.sum(late)),
                max_arrival_ns=float(np.max(end_time)),
                scope=scope)


def interpolants(grid, u, v):
    # Central wire ONLY, by explicit first-stage design.
    if np.any(np.rint(u/.3).astype(int) != 0):
        raise ValueError("Target muon contains a noncentral wire; refuse folding")
    radius = np.hypot(u, v) * 10000
    angle = np.degrees(np.arctan2(np.abs(v), u))
    if np.any(radius <= R_WIRE_UM) or np.any(angle <= 0) or np.any(angle >= 180):
        raise ValueError("A target ion lies in wire metal or exact v=0")
    radii = grid["radius_um"]
    angles = grid["angle_deg"]
    margin = 1e-10
    outside = (radius < radii[0] - margin) | (radius > radii[-1] + margin) | (
        angle < angles[0] - margin) | (angle > angles[-1] + margin)
    if np.any(outside):
        raise ValueError(
            f"Polar bank does not cover {int(outside.sum())}/{len(outside)} macroions; "
            "no extrapolation, clipping, or data-based grid fitting permitted. "
            f"Observed radius [{radius.min():.8g}, {radius.max():.8g}] um; "
            f"theta [{angle.min():.8g}, {angle.max():.8g}] deg; grid r "
            f"[{radii[0]:.8g}, {radii[-1]:.8g}], theta "
            f"[{angles[0]:.8g}, {angles[-1]:.8g}].")
    # Important: interpolate in log(radius) but do NOT pretend that
    # ln(r) is a field-line coordinate. Near separatrix may remain nonsmooth.
    logr = np.log(radius)
    grid_logr = np.log(radii)
    i = np.clip(np.searchsorted(grid_logr, logr, side="right") - 1,
                0, len(radii) - 2)
    j = np.clip(np.searchsorted(angles, angle, side="right") - 1,
                0, len(angles) - 2)
    a = (logr - grid_logr[i]) / (grid_logr[i + 1] - grid_logr[i])
    b = (angle - angles[j]) / (angles[j + 1] - angles[j])
    if np.any(a < -1e-10) or np.any(a > 1+1e-10) or np.any(b < -1e-10) or np.any(b > 1+1e-10):
        raise AssertionError("Invalid bilinear fraction; outside grid")
    # Four corners in radius-major flattened grid.
    nt = len(angles)
    indices = np.stack([i * nt + j, (i+1) * nt + j,
                        i * nt + j + 1, (i+1) * nt + j + 1], axis=1)
    weights = np.stack([(1-a)*(1-b), a*(1-b),
                        (1-a)*b, a*b], axis=1)
    if not np.allclose(weights.sum(axis=1), 1.0, rtol=0, atol=1e-12):
        raise AssertionError("Interpolated node coefficients do not sum to 1")
    return indices, weights, radius, angle


def group_replay(grid, ions, n_bins=BINS):
    n_sites = grid["waves"].shape[0]
    u = np.asarray(ions["u_patch_cm"], dtype=float)
    v = np.asarray(ions["v_absolute_cm"], dtype=float)
    w = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    t = np.asarray(ions["t_birth_ns"], dtype=float) / DT_NS
    offset = np.floor(t).astype(int)
    frac = t - offset
    if np.any(t < 0) or np.any(offset >= n_bins):
        raise ValueError("Birth outside existing reference 100us window")
    indices, interp, radius, theta = interpolants(grid, u, v)
    streams = []
    for sign in (1, -1):
        moment = np.zeros((n_bins, n_sites), dtype=float)
        ids = np.flatnonzero((v > 0) if sign == 1 else (v < 0))
        for k in range(4):
            ix = indices[ids, k]
            coefficients = interp[ids, k] * w[ids]
            np.add.at(moment, (offset[ids], ix), coefficients * (1-frac[ids]))
            valid = offset[ids] + 1 < n_bins
            np.add.at(moment, (offset[ids][valid] + 1, ix[valid]),
                      coefficients[valid] * frac[ids][valid])
        streams.append(moment)
    output = np.zeros((n_bins, 2), dtype=float)
    # FFT convolution is EXACTLY equivalent to the previous template
    # fractional-shift deposition, up to numerical FFT roundoff.
    for sign, stream in zip((1, -1), streams):
        F = np.fft.rfft(stream, n=(1 << (2*n_bins - 1).bit_length()), axis=0)
        for channel in (0, 1):
            original_channel = channel if sign > 0 else 1-channel
            template = grid["waves"][:, original_channel*n_bins:(original_channel+1)*n_bins]
            T = np.fft.rfft(template, n=(1 << (2*n_bins - 1).bit_length()), axis=1)
            prediction = np.fft.irfft(
                np.einsum("fn,nf->f", F, T, optimize=True),
                n=(1 << (2*n_bins - 1).bit_length()))[:n_bins]
            output[:, channel] += prediction
    if not np.all(np.isfinite(output)):
        raise ValueError("Nonfinite replay signal")
    return output, dict(
        min_r_um=float(radius.min()), max_r_um=float(radius.max()),
        min_theta_deg=float(theta.min()), max_theta_deg=float(theta.max()),
        exact_four_corner_sum=True, out_of_grid_weight_fraction=0.0)


def synthetic_selftest():
    rng = np.random.default_rng(20261121)
    nr, nt, n_bins = 3, 6, 64
    r = np.array([15.00001, 20, 80.])
    a = np.linspace(.001, 179.999, nt)
    W = np.empty((nr*nt, 2*n_bins))
    q = np.exp(-np.arange(n_bins) / 8.)
    for i in range(nr):
        for j in range(nt):
            x = np.log(r[i]); y = a[j]/180.
            W[i*nt+j] = np.r_[q*(.3 + .2*x + .1*y + .03*x*y),
                               q*(-.4 + .1*x - .3*y + .05*x*y)]
    grid = dict(radius_um=r, angle_deg=a, waves=W)
    theta = np.array([2., 75., 110., 176.]) * np.pi / 180
    rr = np.array([15.1, 25., 50., 79.])
    u = rr * np.cos(theta) / 10000
    v = rr * np.sin(theta) / 10000
    v[1] = -v[1]
    v[3] = -v[3]
    wt = rng.uniform(.2, 3., 4)
    time_ns = np.array([2.1, 22.7, 49.2, 105.3])
    ions = dict(u_patch_cm=u, v_absolute_cm=v,
                inverse_inclusion_weight=wt, t_birth_ns=time_ns)
    fast, _ = group_replay(grid, ions, n_bins)
    idx, alpha, _, _ = interpolants(grid, u, v)
    slow = np.zeros_like(fast)
    for ion in range(len(u)):
        signal = alpha[ion] @ W[idx[ion]]
        offset = int(np.floor(time_ns[ion]/10.))
        frac = (time_ns[ion]/10.) - offset
        for c in range(2):
            src = c if v[ion] > 0 else 1-c
            for j in range(n_bins):
                current = wt[ion] * signal[src*n_bins+j]
                if offset+j < n_bins:
                    slow[offset+j,c] += (1-frac)*current
                if offset+j+1 < n_bins:
                    slow[offset+j+1,c] += frac*current
    difference = float(np.max(np.abs(fast-slow)))
    if difference > 1e-10:
        raise AssertionError("Bilinear+fractional timing FFT identity failed: "+str(difference))
    print("PASS synthetic bilinear grid, positive/negative-v plane swapping, "
          "weighted 10ns fractional convolution; max_abs_difference=", difference)


def write_result(out, grid_path, input_path, ref_wave, ref_summary, baseline, grid, meta, rows, prediction):
    out.mkdir(parents=True, exist_ok=False)
    with (out/"event_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with (out/"first_100ns_waveforms.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("time_center_ns",
                         "polar_plus_fC_per_ns", "polar_minus_fC_per_ns",
                         "direct_plus_fC_per_ns", "direct_minus_fC_per_ns"))
        direct = load_reference(ref_wave)
        for j in range(10):
            writer.writerow(((j+.5)*DT_NS, prediction[j,0], prediction[j,1],
                             direct[j,0], direct[j,1]))
    note = dict(
        scope="Independent original K256 muon, deterministic direct RKF grid",
        grid_path=str(grid_path.resolve()), ion_archive=str(input_path.resolve()),
        direct_waveform=str(ref_wave.resolve()), direct_summary=str(ref_summary.resolve()),
        IDW_baseline=None if baseline is None else str(baseline.resolve()),
        grid_nodes_r=int(grid["nradius"]), grid_nodes_theta=int(grid["nangle"]),
        radius_grid_um=grid["radius_um"].tolist(),
        theta_grid_deg=grid["angle_deg"].tolist(),
        grid_late_rows=grid["late_rows"],
        grid_max_arrival_ns=grid["max_arrival_ns"],
        grid_scope=grid["scope"],
        source=meta,
        models=["direct RKF source/reference", "optional unchanged prior IDW4",
                "direct uncompressed bilinear(logr,theta) interpolation"],
        caveats=[
            "Whole cathode planes only; cannot infer physical U/V readout",
            "No SVD, polynomial ridge moments or predicted generative avalanche",
            "One independent muon already inspected: development event only",
            "Exact ion-reservoir K256 weighted sampling shared with direct reference",
            "v reflection only; other wire-index folding explicitly forbidden",
            "Grid interpolation errors can be large near v=0 separatrix",
            "No claim that all field lines have time shifts representable at 10ns",
            "Late arrival after 100us gives finite-window current only",
            "No prompt electron, electronics, calibrated ion mixture mobility"
        ])
    (out/"metadata.json").write_text(json.dumps(note, indent=2)+"\n")
    print("PASS uncompressed deterministic polar grid replay, source", meta["crossing_id"])
    for x in rows:
        print(f"{x['model']:25s} relative_L2={x['combined_L2']:.7f}"
              f" first100 errors=({x['plus_first100_error_fC']:.6g},"
              f"{x['minus_first100_error_fC']:.6g}) fC")
    print("NOTE: no SVD and no new microscopic simulation; original source ROOT is read-only.")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--grid-root", type=Path)
    parser.add_argument("--ion-file", type=Path)
    parser.add_argument("--reference-waveforms", type=Path)
    parser.add_argument("--reference-summary", type=Path)
    parser.add_argument("--baseline-waveforms", type=Path)
    parser.add_argument("--out-dir", type=Path)
    a=parser.parse_args()
    if a.self_test:
        synthetic_selftest(); return
    for name in ("grid_root", "ion_file", "reference_waveforms", "reference_summary",
                 "out_dir"):
        if getattr(a, name) is None:
            parser.error("Required option: --"+name.replace("_","-"))
    if a.out_dir.exists():
        parser.error("Refusing to overwrite existing output directory: "+str(a.out_dir))
    start=time.perf_counter()
    grid=load_grid(a.grid_root)
    meta, ions, exact=load_ions(a.ion_file)
    check_direct(meta, a.reference_summary)
    direct=load_reference(a.reference_waveforms)
    rows=[metrics("REFERENCE_DIRECT_RKF", direct, direct)]
    if a.baseline_waveforms:
        baseline=load_baseline(a.baseline_waveforms, direct)
        rows.append(metrics("EXISTING_IDW4", baseline, direct))
    calc_start=time.perf_counter()
    predicted, coverage=group_replay(grid, ions)
    replay_sec=time.perf_counter()-calc_start
    meta.update(coverage)
    meta["pure_replay_elapsed_seconds"]=replay_sec
    meta["total_elapsed_seconds_before_files"]=time.perf_counter()-start
    rows.append(metrics("polar_grid_bilinear_uncompressed", predicted, direct))
    write_result(a.out_dir, a.grid_root, a.ion_file, a.reference_waveforms,
                 a.reference_summary, a.baseline_waveforms,
                 grid, meta, rows, predicted)


if __name__=="__main__":
    main()

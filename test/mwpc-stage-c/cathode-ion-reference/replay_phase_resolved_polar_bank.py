#!/usr/bin/env python3
"""Deterministic 3-coordinate (ln r, theta, 10ns birth-phase) ion operator.

All data are existing muon0301 2100V K256 ion birth samples, EXACT same archived
direct-RKF reference, whole planes only. Direct single-ion precomputation uses
original Garfield++ RKF field and provisional Ar+ mobility. This script uses
no Garfield, no FEM, no electron, no TrackHeed, no microscopic avalanches,
no SVD and no regression.

Relative to original polar replay, birth-phase is interpolated in genuinely
recomputed Garfield Sensor waveforms at 0,2.5,5,7.5 ns. The phase 10 ns
endpoint is exactly the phase-0 response shifted by one integer 10 ns bin.
Spatial coordinates use the already-tested bilinear (ln r,theta) scheme.
"""
import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from study_event_moment_replay import (
    BINS, DT_NS, load_ions, check_direct, load_reference,
    load_baseline, metrics
)
from replay_deterministic_polar_bank import interpolants

NFFT = 1 << (2 * BINS - 1).bit_length()


def load_phase_grid(filename, n_bins=BINS):
    try:
        import uproot
        import awkward as ak
    except ImportError as e:
        raise SystemExit("uproot and awkward are required") from e
    with uproot.open(filename) as root:
        if "polar_bank" not in root or "scope" not in root or "phase_nodes" not in root:
            raise ValueError("Not a phase-resolved deterministic polar grid")
        scope = str(root["scope"].member("fTitle"))
        if "positive-v" not in scope or "Two ideal whole cathode planes" not in scope:
            raise ValueError("Input is not the expected two ideal cathode planes")
        tree = root["polar_bank"]
        names = ("i_radius", "i_angle", "i_phase", "phase_ns", "r_um",
                 "theta_deg", "u_cm", "v_cm", "late_after_100us",
                 "end_time_ns", "waveform_plus_then_minus")
        if any(n not in tree for n in names):
            raise ValueError("Phase grid missing one or more required ROOT branches")
        ir = np.asarray(tree["i_radius"].array(library="np"), dtype=int)
        ia = np.asarray(tree["i_angle"].array(library="np"), dtype=int)
        ip = np.asarray(tree["i_phase"].array(library="np"), dtype=int)
        phase_ns = np.asarray(tree["phase_ns"].array(library="np"), dtype=float)
        radius = np.asarray(tree["r_um"].array(library="np"), dtype=float)
        angle = np.asarray(tree["theta_deg"].array(library="np"), dtype=float)
        uu = np.asarray(tree["u_cm"].array(library="np"), dtype=float)
        vv = np.asarray(tree["v_cm"].array(library="np"), dtype=float)
        end_time = np.asarray(tree["end_time_ns"].array(library="np"), dtype=float)
        late = np.asarray(tree["late_after_100us"].array(library="np"), dtype=bool)
        waveform = np.asarray(ak.to_numpy(ak.to_regular(
            tree["waveform_plus_then_minus"].array(library="ak"))), dtype=np.float64)

    nr, nt, np_ = int(np.max(ir))+1, int(np.max(ia))+1, int(np.max(ip))+1
    if (np_ != 4 or nr < 3 or nt < 6 or
            len(ir) != nr*nt*np_ or nr*nt*np_ > 2048):
        raise ValueError("Requires complete four-phase polar grid, <=2048 RKF starts")
    if not np.array_equal(ir, np.repeat(np.arange(nr), nt*np_)):
        raise ValueError("Radius-major ordering is invalid")
    if not np.array_equal(ia, np.tile(np.repeat(np.arange(nt), np_), nr)):
        raise ValueError("Angle-major ordering is invalid")
    if not np.array_equal(ip, np.tile(np.arange(np_), nr*nt)):
        raise ValueError("Phase-major ordering is invalid")
    rr = radius.reshape(nr,nt,np_)[:,0,0]
    th = angle.reshape(nr,nt,np_)[0,:,0]
    expected_phase = np.arange(np_)*DT_NS/np_
    if not np.allclose(phase_ns, np.tile(expected_phase, nr*nt), rtol=0, atol=1e-10):
        raise ValueError("Required 0,2.5,5,7.5ns phases are missing")
    if (np.any(np.diff(rr) <= 0) or np.any(np.diff(th) <= 0) or
            rr[0] <= 15. or th[0] <= 0. or th[-1] >= 180.):
        raise ValueError("Invalid radial/angular support, or wire/separatrix node")
    if not np.allclose(radius.reshape(nr,nt,np_), rr[:,None,None], atol=1e-8):
        raise ValueError("Radial nodes differ across phases or angles")
    if not np.allclose(angle.reshape(nr,nt,np_), th[None,:,None], atol=1e-8):
        raise ValueError("Angular nodes differ across phases or radii")
    if (not np.allclose(np.hypot(uu,vv)*10000, radius, atol=1e-8) or
            not np.allclose(np.degrees(np.arctan2(vv,uu)), angle, atol=1e-8)):
        raise ValueError("ROOT Cartesian starts inconsistent with polar coordinates")
    if waveform.shape != (nr*nt*np_, 2*n_bins) or not np.all(np.isfinite(waveform)):
        raise ValueError("Wrong waveform size/nonfinite RKF impulse")
    if not np.all(np.isfinite(end_time)):
        raise ValueError("Bad RKF endpoint time")
    # Site-major, phase-next, channel-next, time-next.
    waves = waveform.reshape(nr*nt,np_,2,n_bins)
    return dict(radius_um=rr, angle_deg=th, waves=waves,
                nradius=nr, nangle=nt, nphase=np_,
                late_rows=int(late.sum()), max_arrival_ns=float(end_time.max()),
                scope=scope)


def phase_weights(t_ns, nphase, n_bins):
    t = np.asarray(t_ns, dtype=float) / DT_NS
    offset = np.floor(t).astype(int)
    if np.any(~np.isfinite(t)) or np.any(offset < 0) or np.any(offset >= n_bins):
        raise ValueError("Source ion birth outside original current window")
    f = t-offset
    x = f*nphase
    k = np.minimum(nphase-1, np.floor(x).astype(int))
    frac = x-k
    if np.any(frac < 0) or np.any(frac >= 1.0+1e-10):
        raise AssertionError("Bad sub-bin timing fraction")
    # Phase endpoint nPhase means phase0 + ONE integer bin, never fractional bin.
    low_ph = k
    high_ph = np.where(k+1 == nphase, 0, k+1)
    low_bin = offset
    high_bin = offset+(k+1 == nphase)
    return low_ph, high_ph, low_bin, high_bin, 1-frac, frac


def replay_fft(grid, ions, n_bins=BINS):
    u = np.asarray(ions["u_patch_cm"], dtype=float)
    v = np.asarray(ions["v_absolute_cm"], dtype=float)
    weights = np.asarray(ions["inverse_inclusion_weight"], dtype=float)
    indexes, spatial, radii, angles = interpolants(grid,u,v)
    lowph, highph, lowbin, highbin, loww, highw = phase_weights(
        ions["t_birth_ns"], grid["nphase"], n_bins)
    ns = grid["waves"].shape[0]
    nf = 1 << (2*n_bins-1).bit_length()
    reconstructed = np.zeros((n_bins,2), dtype=np.float64)
    total_in = float(np.sum(weights))
    total_coeff = 0.0

    for sign in (1,-1):
        ids_sign = np.flatnonzero((v > 0) if sign == 1 else (v < 0))
        # Stream one phase at a time: keeps FFT memory bounded for 12x24x4.
        for p in range(grid["nphase"]):
            moment = np.zeros((n_bins,ns), dtype=float)
            for pp, bb, ww in ((lowph, lowbin, loww),
                                (highph, highbin, highw)):
                valid = ids_sign[(pp[ids_sign] == p) &
                                 (bb[ids_sign] < n_bins) &
                                 (ww[ids_sign] > 0)]
                if valid.size == 0:
                    continue
                for corner in range(4):
                    coefficient = weights[valid]*spatial[valid,corner]*ww[valid]
                    np.add.at(moment, (bb[valid], indexes[valid,corner]), coefficient)
                    total_coeff += float(np.sum(coefficient))
            if not np.any(moment):
                continue
            fm = np.fft.rfft(moment, n=nf, axis=0)
            del moment
            for channel in (0,1):
                original_c = channel if sign > 0 else 1-channel
                fw = np.fft.rfft(grid["waves"][:,p,original_c,:], n=nf, axis=1)
                freq = np.einsum("fn,nf->f",fm,fw,optimize=True)
                reconstructed[:,channel] += np.fft.irfft(freq,n=nf)[:n_bins]
                del fw, freq
            del fm
    if not np.all(np.isfinite(reconstructed)):
        raise ValueError("Nonfinite predicted signal")
    if abs(total_coeff-total_in) > 1e-8 * max(1, total_in):
        # The check is of coefficients, not total integrated current.
        raise AssertionError("Weighted source mass was lost in deposition")
    return reconstructed, dict(
        n_sites=ns, nphase=grid["nphase"], no_extrapolation=True,
        radius_min_um=float(radii.min()), radius_max_um=float(radii.max()),
        angle_min_deg=float(angles.min()), angle_max_deg=float(angles.max()),
        n_target_macroions=len(u), weighted_ion_population=total_in
    )


def brute_deposit(grid, ions, n_bins):
    """Slow direct eight-corner scalar reference for small synthetic self-test."""
    u=np.asarray(ions["u_patch_cm"],float)
    v=np.asarray(ions["v_absolute_cm"],float)
    weight=np.asarray(ions["inverse_inclusion_weight"],float)
    idx,alpha,*_=interpolants(grid,u,v)
    pl,ph,bl,bh,wl,wh=phase_weights(ions["t_birth_ns"],grid["nphase"],n_bins)
    out=np.zeros((n_bins,2),float)
    for ion in range(len(u)):
        for corner in range(4):
            for p,b,ph_w in ((pl[ion],bl[ion],wl[ion]),(ph[ion],bh[ion],wh[ion])):
                if ph_w <= 0 or b>=n_bins: continue
                for channel in (0,1):
                    original_c=channel if v[ion] > 0 else 1-channel
                    out[b:,channel] += (
                        weight[ion]*alpha[ion,corner]*ph_w*
                        grid["waves"][idx[ion,corner],p,original_c,:n_bins-b])
    return out


def selftest():
    rng=np.random.default_rng(20261121)
    n_bins=64
    nr,nt,nphase=3,6,4
    radii=np.array([15.00001,22.,80.])
    angles=np.linspace(.001,179.999,nt)
    waves=np.zeros((nr*nt,nphase,2,n_bins))
    tt=np.arange(n_bins)
    for ir in range(nr):
        for it in range(nt):
            for ip in range(nphase):
                f=(.5+.15*np.log(radii[ir])+angles[it]/1000.)*(1+.08*ip)
                waves[ir*nt+it,ip,0,:] = f*np.exp(-tt/(7.+ip))*(tt<40)
                waves[ir*nt+it,ip,1,:] = -f*np.exp(-tt/(9.+ip))*(tt<50)
    grid=dict(radius_um=radii, angle_deg=angles,
              waves=waves, nphase=nphase)
    theta=np.array([1.,51.,110.,176.,160.,75.])*np.pi/180
    r=np.array([15.3,21.,38.,79.,57.,15.0002])
    u=r*np.cos(theta)/10000
    v=r*np.sin(theta)/10000
    v[[1,4]]=-v[[1,4]]
    ions=dict(u_patch_cm=u,v_absolute_cm=v,
              inverse_inclusion_weight=rng.uniform(.2,8.,len(u)),
              t_birth_ns=np.array([.001,2.4,4.999,19.91,37.5,99.9]))
    quick,info=replay_fft(grid,ions,n_bins)
    slow=brute_deposit(grid,ions,n_bins)
    error=float(np.max(np.abs(quick-slow)))
    if error > 3e-12:
        raise AssertionError("Phase/spatial/negative-v FFT did not match brute replay: "
                             + str(error))
    print("PASS 3D bilinear(logr,theta)+4 phase FFT replay vs scalar brute "
          "including upper/lower swap, boundary phase, weights and birth-time offset; "
          "max_abs=",error)


def compatible_trace(path,direct,kind):
    a=np.loadtxt(path,delimiter=",",skiprows=1)
    if kind=="phase4":
        if a.shape != (BINS,9):
            raise ValueError("Expected 10000x9 exact-spatial phase4 C++ comparison")
        old_direct=a[:,7:9]
        prediction=a[:,5:7]
    elif kind=="oldpolar":
        if a.shape != (BINS,5):
            raise ValueError("Expected 10000x5 old polar prediction")
        old_direct=a[:,3:5]
        prediction=a[:,1:3]
    else:
        raise AssertionError("Unknown trace kind")
    if not np.allclose(a[:,0],(np.arange(BINS)+.5)*DT_NS,rtol=0,atol=1e-6):
        raise ValueError("Different current bins")
    if not np.allclose(old_direct,direct,rtol=2e-10,atol=1e-11):
        raise ValueError("Previous response is not the same direct RKF event")
    return prediction


def window_results(comparisons,direct):
    bands=(("0to100ns",0,10),("100ns_to_1us",10,100),
           ("1us_to_10us",100,1000),("10us_to_100us",1000,BINS))
    results=[]
    for label,model in comparisons:
        residual=model-direct
        for band,lo,hi in bands:
            base=direct[lo:hi]
            err=residual[lo:hi]
            results.append(dict(
                model=label,window=band,
                combined_L2_window=float(np.sqrt(
                    np.sum(err*err)/max(np.sum(base*base),1e-50))),
                plus_charge_error_fC=float(DT_NS*np.sum(err[:,0])),
                minus_charge_error_fC=float(DT_NS*np.sum(err[:,1]))
            ))
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--grid-root",type=Path)
    p.add_argument("--ion-file",type=Path)
    p.add_argument("--reference-waveforms",type=Path)
    p.add_argument("--reference-summary",type=Path)
    p.add_argument("--baseline-waveforms",type=Path)
    p.add_argument("--old-polar-waveforms",type=Path)
    p.add_argument("--exact-spatial-phase-comparison",type=Path)
    p.add_argument("--out-dir",type=Path)
    a=p.parse_args()
    if a.self_test:
        selftest();return
    for key in ("grid_root","ion_file","reference_waveforms","reference_summary",
                "old_polar_waveforms","exact_spatial_phase_comparison","out_dir"):
        if getattr(a,key) is None:
            p.error("Required parameter --"+key.replace("_","-"))
    if a.out_dir.exists():
        p.error("Refusing to overwrite "+str(a.out_dir))
    begin=time.perf_counter()
    grid=load_phase_grid(a.grid_root)
    meta,ions,exact=load_ions(a.ion_file)
    check_direct(meta,a.reference_summary)
    direct=load_reference(a.reference_waveforms)
    exact_phase=compatible_trace(a.exact_spatial_phase_comparison,direct,"phase4")
    prior_polar=compatible_trace(a.old_polar_waveforms,direct,"oldpolar")
    controls=[("FROZEN_DIRECT_RKF",direct)]
    if a.baseline_waveforms:
        controls.append(("IDW4_original",load_baseline(a.baseline_waveforms,direct)))
    controls.extend([("polar_12x24_original_10ns_shift",prior_polar),
                     ("exact_spatial_phase4",exact_phase)])
    time0=time.perf_counter()
    prediction,geometry=replay_fft(grid,ions)
    elapsed=time.perf_counter()-time0
    models=controls+[("polar_12x24_phase4",prediction)]
    rows=[metrics(name,x,direct) for name,x in models]
    window=window_results(models,direct)
    spatial=prediction-exact_phase
    timing=exact_phase-direct
    total=prediction-direct
    denom=float(np.sum(direct*direct))
    spatial_e2=float(np.sum(spatial*spatial)/denom)
    timing_e2=float(np.sum(timing*timing)/denom)
    cross=float(2*np.sum(spatial*timing)/denom)
    total_e2=float(np.sum(total*total)/denom)
    if not np.isclose(total_e2,spatial_e2+timing_e2+cross,rtol=1e-10,atol=1e-14):
        raise AssertionError("Error decomposition failed")
    decomposition=dict(
        phase4_spatial_L2_relative_to_direct=float(np.sqrt(spatial_e2)),
        exact_spatial_phase4_timing_L2_relative_to_direct=float(np.sqrt(timing_e2)),
        cross_term_in_relative_squared_L2=cross,
        total_phase4_polar_L2=float(np.sqrt(total_e2)),
        replay_elapsed_seconds=elapsed,
        full_analysis_elapsed_before_files_seconds=time.perf_counter()-begin
    )
    a.out_dir.mkdir(parents=True,exist_ok=False)
    for filename,records in (("event_metrics.csv",rows),
                             ("window_metrics.csv",window)):
        with (a.out_dir/filename).open("x",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader();writer.writerows(records)
    with (a.out_dir/"complete_waveforms.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("time_center_ns","polar_phase4_plus_fC_per_ns",
                    "polar_phase4_minus_fC_per_ns",
                    "exact_spatial_phase4_plus_fC_per_ns",
                    "exact_spatial_phase4_minus_fC_per_ns",
                    "original_direct_plus_fC_per_ns",
                    "original_direct_minus_fC_per_ns"))
        for j in range(BINS):
            w.writerow(((j+.5)*DT_NS,*prediction[j],*exact_phase[j],*direct[j]))
    info=dict(
        scope="Phase-resolved polar ion response only: whole-plane model, no U/V",
        source=meta,
        source_archive=str(a.ion_file.resolve()),
        exact_ion_births=int(exact.sum()),
        grid_file=str(a.grid_root.resolve()),
        grid_radial_nodes=grid["radius_um"].tolist(),
        grid_angular_nodes=grid["angle_deg"].tolist(),
        phase_ns=[0.,2.5,5.,7.5],
        grid_late_rows=grid["late_rows"],
        grid_max_arrival_ns=grid["max_arrival_ns"],
        numerical_contribution_decomposition=decomposition,
        coverage=geometry,
        limitations=[
            "Only development muon0301 was evaluated; 0408/0446 remain sealed",
            "Original 4096 sampled ion coordinates and weights are GIVEN",
            "No demonstration that ion clouds can be generated from few descriptors",
            "No SVD: this operator is not yet compressed",
            "Two ideal whole cathode planes, no stereo U/V strips or electronics",
            "Uncalibrated Ar+ ion mobility surrogate",
            "Not a proof of event-population generalization",
            "Four-phase error includes its own phase interpolation approximation"
        ]
    )
    (a.out_dir/"metadata.json").write_text(json.dumps(info,indent=2)+"\n")
    print("PASS phase-aware deterministic polar grid 12x24x4, existing muon 0301")
    for r in rows:
        print("%-40s L2=%0.8f  Q100=(%+.6g,%+.6g) fC"
              %(r["model"],r["combined_L2"],
                r["plus_first100_error_fC"],r["minus_first100_error_fC"]))
    print("DECOMPOSITION",json.dumps(decomposition,indent=2))
    print("No new microscopic avalanches, no spatial extrapolation, "
          "no new independent holdout evaluated.")


if __name__=="__main__":
    main()

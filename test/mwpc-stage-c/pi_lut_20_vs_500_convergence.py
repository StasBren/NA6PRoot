#!/usr/bin/env python3
"""Read-only microscopic-LUT convergence: 20 vs 500 trials/cell at p=3mm.

Uses exactly the same native Geant4 hit segments and deterministic Heed
births as the original Stage-C fold, but stores only a 900 x (nwire+2)
birth-occupancy grid. Rather than re-sampling one LUT result for every
Heed electron, integrate over ALL independent microscopic trajectories in
each occupied LUT cell, preserving the finite physical anode aperture.

Per-phi-bin Garfield uncertainty uses independent drift trials PER SPATIAL
LUT CELL, not the millions of reused Heed draws or reflected symmetry
copies. Reuse of each microscopic response at many physical anode indices
is treated as correlated (within that LUT cell). The error bars are
conditional on sampled Heed births and do not include Heed/Geant4 event
statistics or systematic LUT grid error.

Outputs: JSON, CSV, fully self-contained HTML/SVG. No C++/Garfield/Geant4
or original frozen 6M folds are modified or started.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from pi_pitch_fieldwire_pilot import tag
from pi_six_geometry_native_timing import (
    MAPS, HEED_TRACKS, HEED_CLUSTERS, birth_sample, first_indices,
    load_map, n_anodes, stable_seed,
)
from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord
from run_pi_per_chamber_response import read_heed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT500 = REPO / "test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/pilot"
OUT = REPO / "test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_v1"
EDGES = np.arange(-180, 181, 10)
CENTERS = (EDGES[:-1]+EDGES[1:])/2
NCELL = 900
PITCH = 3.
CONFIGS = (0., 75.)


def occupancies(station, chunks, max_hits, tracks, clusters):
    """Heed birth occupancy keyed by (cell, physical nearest wire).

    Physical nearest wire can be -1 or n_wires in the finite gas margin.
    The Garfield response may migrate those electrons to a valid anode.
    """
    nwire = n_anodes(station, PITCH)
    out = np.zeros((NCELL, nwire+2), dtype=np.int64)
    seen = sampled = born = zero_paths = 0
    native_sha = []
    for chunk in range(1, chunks+1):
        path, provenance = crossing_provenance("Jpsi", chunk)
        native_sha.append(provenance["csv_sha256"])
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            for row in reader:
                if int(row["station"]) != station or abs(int(row["pdg"])) != 13:
                    continue
                seen += 1
                if max_hits and sampled >= max_hits:
                    continue
                length = float(row["path_length_mm"])
                if length == 0.:
                    # Match the immutable native-fold exception exactly; an
                    # inconsistent zero chord is an input error, not silently
                    # discarded detector response.
                    if not coincident_zero_chord(row):
                        raise ValueError("Invalid zero-chord native gas crossing")
                    zero_paths += 1
                    sampled += 1
                    continue
                if length < 0. or not math.isfinite(length):
                    raise ValueError(f"Invalid native gas chord: {length}")
                seed = stable_seed(
                    "native-six-pitch-shared-Heed-v1", "jpsi40-bore",
                    f"chunk{chunk:04d}", row["event"], row["hit_index"],
                    row["track_id"])
                pos, _ = birth_sample(row, tracks, clusters, seed)
                sampled += 1
                if not len(pos):
                    continue
                u, v = pos[:, 0], pos[:, 1]
                nearest, cell = first_indices(u, PITCH, nwire)
                if np.any((nearest < -1) | (nearest > nwire)):
                    raise ValueError("Primary electron outside finite anode margin")
                iu = np.clip(np.floor((cell+PITCH/2)/.1).astype(int), 0, 29)
                iv = np.clip(np.floor((v+3.)/.2).astype(int), 0, 29)
                np.add.at(out, (iu*30+iv, nearest+1), 1)
                born += len(pos)
        print(f"Native J/psi chunk {chunk}/{chunks}: "
              f"selected station gas crossings so far={sampled:,}, "
              f"Heed births={born:,}", flush=True)
    if born == 0:
        raise ValueError("No Heed births selected; increase input exposure")
    if int(out.sum()) != born:
        raise AssertionError("Heed occupancy conservation failed")
    return out, {
        "station": station, "chunks": chunks, "source": "jpsi40-bore",
        "selected_muon_gas_crossings": sampled,
        "all_muon_gas_crossings_seen": seen,
        "zero_length_crossings": zero_paths,
        "heed_primaries": born, "native_csv_sha256": native_sha,
        "physical_anodes": nwire, "geometrical_lut_cells": NCELL,
        "full_source_exposure": chunks == 40 and not max_hits,
    }


def microscopic_pdf(occupancy, pool):
    """Weighted all-trial mean and delta-method MC SE with cell-level reuse.

    W[b,nearest] counts actual Heed births. For each Garfield response trial
    in box b, w[b,k] is the total birth population for WHICH the returned
    relative wire index points to a finite physical anode. All repeats of a
    Garfield drift outcome at different physical anode positions are then
    correlated by construction and handled via their shared w[b,k].
    """
    cells, K, features = pool.shape
    if cells != NCELL or features != 5 or K < 2:
        raise ValueError("Wrong LUT schema or insufficient independent trials")
    nwire = occupancy.shape[1]-2
    collected = pool[:, :, 0] == 1.
    offset = pool[:, :, 1].astype(np.int32)
    phi = pool[:, :, 4]
    weight = np.zeros((NCELL, K), dtype=np.float64)
    for shift in np.unique(offset[collected]):
        shift = int(shift)
        # Actual nearest wires range -1..nwire; packed column nearest+1.
        lower, upper = max(-1, -shift), min(nwire, nwire-1-shift)
        if lower <= upper:
            possible = occupancy[:, lower+1:upper+2].sum(axis=1)
            mask = collected & (offset == shift)
            weight[mask] = np.broadcast_to(possible[:, None], mask.shape)[mask]
    finite = collected & np.isfinite(phi) & (phi >= -180.) & (phi <= 180.)
    if np.any((weight > 0) & ~finite):
        raise ValueError("Collected Garfield drift has undefined impact azimuth")
    weight[~finite] = 0.
    j = np.clip(np.floor((np.where(finite,phi,-180.)+180.)/10.).astype(int),0,35)
    H = np.zeros((NCELL, 36), dtype=np.float64)
    H2 = np.zeros((NCELL, 36), dtype=np.float64)
    for b in range(NCELL):
        ok = weight[b] > 0
        if not np.any(ok):
            continue
        H[b] = np.bincount(j[b,ok], weights=weight[b,ok], minlength=36)
        H2[b] = np.bincount(j[b,ok], weights=weight[b,ok]**2,
                             minlength=36)
    N = H.sum(axis=1)
    N2 = H2.sum(axis=1)
    total = N.sum()/K
    if total <= 0:
        raise ValueError("Zero finite anode collection from weighted LUT")
    p = H.sum(axis=0)/(K*total)
    if not np.isclose(p.sum(), 1., rtol=1.e-10, atol=1.e-10):
        raise AssertionError("Weighted angular probability does not normalize")
    # For each cell b and azimuth bin j, microscopic trial variables:
    # X_bkj = w_bk*(I(phi_bk==j) - p_j*I(collected_bk)).
    # Unbiased sample variance of X, divided by K, then summed across
    # independent cells; normalization derivative supplies 1/total².
    moment = H2*(1-2*p[None,:]) + N2[:,None]*p[None,:]**2
    center = (H-N[:,None]*p[None,:])**2/K
    variance = np.maximum(0., (moment-center).sum(axis=0) /
                          (K*(K-1)*total**2))
    se = np.sqrt(variance)
    q2 = float(np.dot(p, np.cos(2*np.deg2rad(CENTERS))))
    return {
        "trials_per_cell": K,
        "heeds_in_scope": int(occupancy.sum()),
        "expected_collected_heed_electrons": float(total),
        "collection_fraction": float(total/occupancy.sum()),
        "phi_probabilities": p.tolist(),
        "garfield_mc_standard_errors": se.tolist(),
        "Q2_binned": q2,
        "occupied_cells": int(np.count_nonzero(occupancy.sum(axis=1))),
        "independent_microscopic_trials_total": int(cells*K),
        "uncertainty_scope": "Finite Garfield trial SE at fixed Heed cell "
                             "occupancies, with common LUT response across "
                             "physical wires treated as correlated",
    }


def svg_plot(name, result20, result500):
    p20 = np.asarray(result20["phi_probabilities"])
    p500 = np.asarray(result500["phi_probabilities"])
    e20 = np.asarray(result20["garfield_mc_standard_errors"])
    e500 = np.asarray(result500["garfield_mc_standard_errors"])
    m = max(float(np.max(p20+2*e20)),float(np.max(p500+2*e500)),1.e-9)
    W,H = 1050,320
    left, right, top, bottom = 64, 24, 24, 49
    w,h=W-left-right,H-top-bottom
    x=lambda i: left+(i+.5)*w/36
    y=lambda v: top+h*(1.-v/m)
    parts=[f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
           'role="img" aria-label="LUT 20 vs 500 angular probabilities">']
    for frac in (0.,.25,.5,.75,1.):
        yy=y(frac*m)
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{W-right}" '
                     f'y2="{yy:.2f}" stroke="#aaa" stroke-opacity=".25"/>')
        parts.append(f'<text x="{left-7}" y="{yy+4:.2f}" '
                     f'text-anchor="end" font-size="11">{frac*m*100:.2f}%</text>')
    for k in range(0,36,3):
        parts.append(f'<text x="{x(k):.2f}" y="{H-18}" '
                     f'text-anchor="middle" font-size="11">{CENTERS[k]:g}°</text>')
    def poly(p,color):
        line=" ".join(f"{x(i):.2f},{y(float(z)):.2f}" for i,z in enumerate(p))
        parts.append(f'<polyline points="{line}" fill="none" '
                     f'stroke="{color}" stroke-width="2.1" stroke-linejoin="round"/>')
    poly(p20,"#dc7c43")
    poly(p500,"#267fbb")
    for i in range(36):
        xx=x(i)
        parts.append(f'<line x1="{xx:.2f}" x2="{xx:.2f}" '
                     f'y1="{y(max(0,p500[i]-1.96*e500[i])):.2f}" '
                     f'y2="{y(p500[i]+1.96*e500[i]):.2f}" '
                     'stroke="#267fbb" stroke-opacity=".65" stroke-width="1.2"/>')
    parts.append(f'<text x="{left+w/2}" y="{H-2}" text-anchor="middle" '
                 'font-size="12">φ_hit [degrees], 36 × 10° bins</text>')
    parts.append('</svg>')
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--station",type=int,choices=range(6),default=0)
    ap.add_argument("--chunks",type=int,default=2)
    ap.add_argument("--max-hits",type=int,default=0,
                    help="Total selected gas crossings; 0 means all in chunks")
    ap.add_argument("--root20",type=Path,default=MAPS)
    ap.add_argument("--root500",type=Path,default=ROOT500)
    ap.add_argument("--out",type=Path,default=OUT)
    args=ap.parse_args()
    if not 1 <= args.chunks <= 40 or args.max_hits < 0:
        ap.error("chunks must be in 1..40, max-hits nonnegative")
    old,new={},{}
    lut_identity={}
    for fw in CONFIGS:
        key=tag(PITCH,fw)
        old[key],oldmeta=load_map(PITCH,fw,args.root20)
        new[key],newmeta=load_map(PITCH,fw,args.root500)
        if (old[key].shape!=(NCELL,20,5) or
            new[key].shape!=(NCELL,500,5)):
            raise ValueError(f"Expected original 20 and new 500 trial LUTs: {key}")
        lut_identity[key]={"LUT20":oldmeta,"LUT500":newmeta}
        if (oldmeta["csv_sha256"]==newmeta["csv_sha256"]):
            raise ValueError("New LUT must be independently computed")
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    occ,scope=occupancies(args.station,args.chunks,args.max_hits,lam,marks)
    results={}
    args.out.mkdir(parents=True,exist_ok=True)
    html=['<!doctype html><html lang="en"><meta charset="utf-8">',
          '<title>Garfield 20 vs 500 convergence | NA60+</title>',
          '<style>body{max-width:1100px;margin:25px auto;padding:0 18px;'
          'font:15px system-ui,sans-serif;line-height:1.5}'
          'svg{max-width:100%;height:auto}section{border:1px solid #ddd;'
          'padding:15px;margin:20px 0;border-radius:12px}'
          'small{color:#555}</style>',
          '<h1>NA60+ MWPC: Garfield LUT 20 vs 500</h1>',
          f'<p>J/ψ, MS{args.station}, first {args.chunks} native 50k-parent chunks, '
          f'{scope["selected_muon_gas_crossings"]:,} real muon gas crossings, '
          f'{scope["heed_primaries"]:,} shared deterministic Heed births.</p>',
          '<p>Orange = 20 microscopic trials/cell; blue = 500. '
          'Blue error bars ±1.96 × LUT-only standard error. '
          'Expected histograms integrate all Garfield trials with the '
          'same Heed spatial-cell occupancy, not a 171M electron resampling.</p>']
    with (args.out/"phi_comparison.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(("geometry","phi_lower_deg","phi_upper_deg",
                         "p20","se20","p500","se500","delta_p",
                         "delta_p_over_combined_se"))
        for fw in CONFIGS:
            key=tag(PITCH,fw)
            a=microscopic_pdf(occ,old[key])
            b=microscopic_pdf(occ,new[key])
            pa,pb=np.asarray(a["phi_probabilities"]),np.asarray(b["phi_probabilities"])
            ea,eb=np.asarray(a["garfield_mc_standard_errors"]),np.asarray(b["garfield_mc_standard_errors"])
            ed=np.sqrt(ea**2+eb**2)
            valid=(pb>1.e-3)&(ed>0)
            delta=pb-pa
            metrics={
                "total_variation_20_vs_500":float(np.abs(delta).sum()/2),
                "rms_relative_change_in_populated_bins":float(
                    np.sqrt(np.mean((delta[valid]/pb[valid])**2))) if np.any(valid) else None,
                "rms_standardized_bin_shift_descriptive":float(
                    np.sqrt(np.mean((delta[valid]/ed[valid])**2))) if np.any(valid) else None,
                "populated_bins":int(np.count_nonzero(valid)),
                "independence_note":"Bins are correlated by normalization; "
                    "RMS standardized shift is descriptive, not a chi-squared test.",
            }
            for i in range(36):
                writer.writerow((key,EDGES[i],EDGES[i+1],
                    pa[i],ea[i],pb[i],eb[i],delta[i],
                    delta[i]/ed[i] if ed[i]>0 else ""))
            results[key]={"LUT20":a,"LUT500":b,"comparison":metrics}
            print(f"{key}: TV(20,500)={metrics['total_variation_20_vs_500']:.6g}; "
                  f"RMS relative difference in populated bins="
                  f"{metrics['rms_relative_change_in_populated_bins']:.3%}; "
                  f"Q2: {a['Q2_binned']:+.5f} -> {b['Q2_binned']:+.5f}",
                  flush=True)
            html.extend([
                f'<section><h2>{key}</h2>',
                f'<p>TV = {metrics["total_variation_20_vs_500"]:.5f}; '
                f'Q₂ = {a["Q2_binned"]:+.4f} (20) vs '
                f'{b["Q2_binned"]:+.4f} (500); '
                f'expected anode collection = '
                f'{a["collection_fraction"]:.2%} (20) vs '
                f'{b["collection_fraction"]:.2%} (500).</p>',
                svg_plot(key,a,b),'</section>'])
    html.extend([
        '<p><small>The displayed uncertainty is Monte Carlo error from '
        'finite independent Garfield trajectories in each spatial cell. '
        'It does NOT estimate Geant4 parent-event uncertainty, systematic '
        'field/gas errors, spatial LUT discretization, or charge avalanche. '
        'Mirrored trajectories are not counted as independent.</small></p>',
        '</html>'])
    (args.out/"compare.html").write_text("\n".join(html),encoding="utf-8")
    (args.out/"comparison.json").write_text(
        json.dumps({"schema":"NA60_GARFIELD_LUT_20_VS_500_V1",
                    "scope":scope,"lut_identity":lut_identity,"results":results},indent=2,
                   ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print(f"PASS: 20 vs 500 read-only LUT convergence data: {args.out}",flush=True)
    print("Open compare.html; detailed 36-bin results in phi_comparison.csv",flush=True)


if __name__=="__main__":
    main()

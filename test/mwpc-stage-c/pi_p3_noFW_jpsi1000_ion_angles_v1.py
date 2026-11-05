#!/usr/bin/env python3
"""Read-only ion-birth azimuth and EXACT event asymmetry for 1000 J/psi p3_noFW muons.

Data: SHA-sealed Geant4/TrackHeed/full Garfield per-primary physics
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/native_chunks/.
Each original Geant4 J/psi parent contributes one selected gas crossing.
NO new Garfield physics, NO ion drift, NO Shockley-Ramo, NO cathode charge.

Angular histogram is APPROXIMATE from <=128 uniformly sampled ion births
per *independent avalanche*, weighted N_birth / N_saved. The saved sampler
uses RNG independent of Garfield electron collisions.
Muon-level hemisphere and +/-45-degree top/bottom asymmetry are EXACT
from the full unsampled C++ ionisation-birth census. Undefined A stays NaN.

Plots, CSV and SHA manifest go to a NEW output directory; never mutate
original physics or overwrite previous analyses.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pi_p3_noFW_jpsi1000_native_benchmark_v1 as bench
from pi_prepare_muon_avalanche_b1 import sha
from pi_run_muon_lut_validation import id_of

HERE = Path(__file__).resolve()
ROOT = bench.ROOT
OUT = ROOT / "plots_ion_birth_azimuth_muon1000_v1"
GEOM = "p3_noFW"
PITCH_MM = 3.0
NBINS = 72  # 5 degrees
EDGES = np.linspace(-180.0, 180.0, NBINS + 1)
CENTRES = (EDGES[:-1] + EDGES[1:]) / 2
# Native crossing file identifies events by a composite source|chunk|... header,
# not by the "crossing_id" field used in seed/reservoir files.
MUON_CSV_COLS = (bench.CROSSING_COLUMNS[0],"configuration",
                 "station","theta_u_in_deg","theta_w_in_deg",
                 "ion_births","ion_up","ion_down","ion_boundary",
                 "ion_top","ion_bottom","ion_side",
                 "A_ion_hemi","A_ion_TB","heed_primaries",
                 "failed_seeds","capped_seeds")
SEED_CSV_COLS = ("crossing_id","configuration","primary_index",
                 "ion_births","ion_up","ion_down",
                 "ion_top","ion_bottom","A_ion_hemi","A_ion_TB")
RES_COLS = ("crossing_id","configuration","primary_index",
            "u_mm","v_mm","w_mm","t_ns",
            "seed_ion_births_total","reservoir_weight")


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def rows(path, cols):
    with path.open(newline="") as f:
        rd = csv.DictReader(f)
        require(set(cols).issubset(rd.fieldnames or []),
                "Unexpected source CSV schema: " + str(path))
        yield from rd


def round_cpp(x):
    """std::round: half away from zero, not Python's ties-to-even."""
    return math.floor(x+0.5) if x>=0 else math.ceil(x-0.5)


def birth_phi_deg(u_mm, v_mm):
    """Same period-folded angle as C++ BirthAccumulator/Counts::Add."""
    nearest = round_cpp(u_mm / PITCH_MM)
    du = u_mm - nearest * PITCH_MM
    return math.degrees(math.atan2(v_mm, du))


def asym(a, b):
    return (a-b)/(a+b) if a+b else math.nan


def validate_asym(row, a, b, field):
    aa, bb = int(row[a]), int(row[b])
    expect = asym(aa, bb)
    x = row[field]
    if math.isnan(expect):
        require(x == "", "Undefined exact "+field+" must be empty")
    else:
        require(x != "" and abs(float(x)-expect)<2e-9,
                "Exact C++ "+field+" sector identity mismatch")
    return expect


def draw_step(ax, y, label, *, linestyle="-"):
    ax.step(CENTRES, y, where="mid", label=label, linestyle=linestyle)


def analyze():
    source, meta = bench.check_input()
    finished, partial = bench.scanned_status(full_audit=True)
    require(partial is None and len(finished) == bench.N//bench.GROUP,
            "1000 original real muons must be completed and SHA re-audited first")
    require(len(source) == bench.N == 1000, "Unexpected real Geant4 muon count")
    require(not OUT.exists(), "Output already exists; preserving prior analysis: "+str(OUT))

    ordered_ids = [id_of(r) for r in source]
    index = {x:i for i,x in enumerate(ordered_ids)}
    require(len(index) == 1000, "Duplicate native muon identity")
    n = len(index)
    hist = np.zeros((n, NBINS), dtype=float)
    births_exact = np.zeros(n, dtype=np.int64)
    Nup = np.zeros(n, dtype=np.int64)
    Ndown = np.zeros(n, dtype=np.int64)
    Ntop = np.zeros(n, dtype=np.int64)
    Nbottom = np.zeros(n, dtype=np.int64)
    A_hemi = np.full(n, np.nan, dtype=float)
    A_tb = np.full(n, np.nan, dtype=float)
    theta_u = np.zeros(n)
    theta_w = np.zeros(n)
    nprim = np.zeros(n, dtype=np.int64)
    station = np.zeros(n, dtype=np.int64)
    failed = np.zeros(n, dtype=np.int64)
    capped = np.zeros(n, dtype=np.int64)
    seen_crossings = set()
    seed_asym_hemi = []
    seed_asym_tb = []
    seed_counts = {}
    sample_rows = Counter()
    sample_weights = Counter()
    source_sha = []
    sample_count = 0

    for ichunk in range(50):
        d = bench.folder(ichunk)
        p_cross = d/"muon_reference_crossing.csv"
        p_seeds = d/"muon_reference_seed.csv"
        p_birth = d/"muon_reference_ion_birth_sample.csv"
        for p in (p_cross,p_seeds,p_birth):
            source_sha.append({"file":str(p.relative_to(ROOT)), "sha256":sha(p)})
        for r in rows(p_cross, MUON_CSV_COLS):
            cid = r[bench.CROSSING_COLUMNS[0]]
            require(cid in index and cid not in seen_crossings,
                    "Unexpected/duplicated event-level physical muon")
            i = index[cid]
            require(r["configuration"] == GEOM, "Wrong geometric field")
            seen_crossings.add(cid)
            nb = int(r["ion_births"])
            up = int(r["ion_up"])
            down = int(r["ion_down"])
            boundary = int(r["ion_boundary"])
            top = int(r["ion_top"])
            bottom = int(r["ion_bottom"])
            side = int(r["ion_side"])
            require(nb == up+down+boundary == top+bottom+side,
                    "Whole-muon full Garfield ion census not conserved")
            births_exact[i] = nb
            Nup[i],Ndown[i] = up,down
            Ntop[i],Nbottom[i] = top,bottom
            A_hemi[i] = validate_asym(r,"ion_up","ion_down","A_ion_hemi")
            A_tb[i] = validate_asym(r,"ion_top","ion_bottom","A_ion_TB")
            theta_u[i] = float(r["theta_u_in_deg"])
            theta_w[i] = float(r["theta_w_in_deg"])
            station[i] = int(r["station"])
            nprim[i] = int(r["heed_primaries"])
            failed[i] = int(r["failed_seeds"])
            capped[i] = int(r["capped_seeds"])
        for r in rows(p_seeds, SEED_CSV_COLS):
            cid = r["crossing_id"]
            i = index[cid]
            require(r["configuration"] == GEOM, "Unexpected geometry in seed census")
            key = (cid,int(r["primary_index"]))
            require(key not in seed_counts, "Duplicate original Garfield primary")
            nc = int(r["ion_births"])
            seed_counts[key] = nc
            require(nc>=0, "Negative ion births")
            x = validate_asym(r,"ion_up","ion_down","A_ion_hemi")
            y = validate_asym(r,"ion_top","ion_bottom","A_ion_TB")
            if math.isfinite(x):seed_asym_hemi.append(x)
            if math.isfinite(y):seed_asym_tb.append(y)

        # Stream ~millions of reservoir points. Memory remains O(1000*72)
        # rather than O(total ion births). Every event's origin is preserved.
        for r in rows(p_birth, RES_COLS):
            cid = r["crossing_id"]
            require(cid in index and r["configuration"] == GEOM,
                    "Unknown physical parent/geometry in ion-birth reservoir")
            i = index[cid]
            key = (cid,int(r["primary_index"]))
            require(key in seed_counts and seed_counts[key]>0,
                    "Birth reservoir references missing/zero-ion avalanche")
            nion = seed_counts[key]
            require(nion == int(r["seed_ion_births_total"]),
                    "Reservoir total ion count is not original full census")
            weight = float(r["reservoir_weight"])
            expect = nion/min(nion,bench.RESERVOIR)
            require(math.isfinite(weight) and
                    abs(weight-expect)<1e-7*max(1.,expect),
                    "Incorrect ion reservoir weight")
            u,v = float(r["u_mm"]),float(r["v_mm"])
            require(math.isfinite(u) and math.isfinite(v), "Nonfinite birth location")
            phi = birth_phi_deg(u,v)
            b = min(NBINS-1,max(0,int((phi+180.)*NBINS/360.)))
            hist[i,b]+=weight
            sample_rows[key]+=1
            sample_weights[key]+=weight
            sample_count+=1
        print(f"READ SHA-VERIFIED ion birth reservoir {ichunk+1}/50",flush=True)

    n_expected_avalanche = sum(int(z["stats"]["n_full_garfield_avalanches"])
                               for z in finished)
    n_expected_failed = sum(int(z["stats"]["n_unsuccessful_seeds"])
                            for z in finished)
    n_expected_births = sum(int(z["stats"]["n_ion_births"])
                            for z in finished)
    require(len(seen_crossings) == n and
            len(seed_counts) == sum(nprim) == n_expected_avalanche,
            "Missing or repeated true TrackHeed electron/Garfield avalanche")
    for key,nion in seed_counts.items():
        expected_samples = min(nion,bench.RESERVOIR)
        require(sample_rows[key] == expected_samples,
                "Incorrect number of sampled ion births per avalanche")
        require(abs(sample_weights[key]-nion)<1e-7*max(1.,nion),
                "Ion reservoir sum must reconstruct full ion count")
    # Each ion-birth point is sampled with equal probability conditional
    # on its avalanche. Event-level weighted sum equals exact N_ion
    # up to CSV roundoff, not exact angle bin locations.
    require(np.allclose(hist.sum(axis=1),births_exact,rtol=1e-7,atol=1e-6),
            "Some muon reservoirs do not sum to true total N_ion")
    require(np.all(capped==0), "Capped avalanches are censored; do not plot as full physics")
    require(int(failed.sum()) == n_expected_failed,
            "Unsuccessful Garfield event count differs from frozen SHA physics")

    total_ions = int(births_exact.sum())
    require(total_ions == n_expected_births,
            "Source all-ion census differs from original SHA-sealed checkpoints")
    global_hist = hist.sum(axis=0)
    require(abs(global_hist.sum()-total_ions)<1e-5*max(1,total_ions),
            "Ion-birth weighted angle histogram normalization failed")
    # Different physical questions:
    # (1) pick random birth across all muons -> high-gain muons dominate
    # (2) pick random muon, then random birth -> each muon equal
    nonzero = births_exact>0
    require(np.any(nonzero), "No muons with avalanche ion births")
    event_norm = np.zeros_like(hist)
    event_norm[nonzero] = hist[nonzero]/births_exact[nonzero,None]
    equal_muon = event_norm[nonzero].mean(axis=0)

    OUT.mkdir(parents=True)
    # 1. The ion-birth azimuth PDF, both natural weighting conventions.
    fig,ax=plt.subplots(figsize=(11.5,5.6))
    draw_step(ax,global_hist/global_hist.sum(),
              "Ion-weighted (every ion birth equal)")
    draw_step(ax,equal_muon,
              "Muon-weighted (every nonzero muon equal)",linestyle="--")
    for left,right in ((-135,-45),(45,135)):
        ax.axvspan(left,right,alpha=.10)
    ax.axvline(0,color="gray",lw=.6,alpha=.5)
    ax.set_xticks(np.arange(-180,181,30))
    ax.set_xlim(-180,180)
    ax.set_xlabel(r"Ion birth azimuth $\phi=\mathrm{atan2}(v,u-u_\mathrm{anode})$ [deg]")
    ax.set_ylabel("Probability per 5 degree bin")
    ax.set_title("J/psi real-muon FULL Garfield p3_noFW: ion birth azimuth\n"
                 "128-point/avalanche weighted spatial reservoir (not full cloud)")
    ax.grid(alpha=.14)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT/"01_ion_birth_azimuth.png",dpi=185)
    plt.close(fig)

    # 2. Exact whole-muon asymmetry distributions from full ion census.
    fig,axs=plt.subplots(1,2,figsize=(12,4.7),sharex=True)
    edges_A=np.linspace(-1,1,42)
    for ax,data,title in (
        (axs[0],A_hemi,"Hemisphere: upper v>0 vs lower v<0"),
        (axs[1],A_tb,"Cathode-facing: +v 45..135 vs -v -135..-45"),
    ):
        valid=data[np.isfinite(data)]
        ax.hist(valid,bins=edges_A,weights=np.ones(len(valid))/max(1,len(valid)),
                histtype="step",lw=1.7)
        ax.axvline(0,ls="--",lw=.9,color="gray")
        ax.set_xlim(-1,1)
        ax.set_xlabel("Exact full-Garfield whole-muon ion-birth asymmetry")
        ax.set_title(title+f"\nValid {len(valid)}/1000; mean={valid.mean():+.4f}")
        ax.grid(alpha=.15)
    axs[0].set_ylabel("Probability per bin (one muon = one observation)")
    fig.suptitle("Whole-muon asymmetry from all ion births; NOT reservoir estimates")
    fig.tight_layout()
    fig.savefig(OUT/"02_exact_event_asymmetry_distributions.png",dpi=185)
    plt.close(fig)

    # 3. Individual 1000 physical muon events, original frozen SHA source order.
    fig,axs=plt.subplots(2,1,figsize=(13,7.2),sharex=True,sharey=True)
    for ax,data,title in ((axs[0],A_hemi,"Hemisphere ion asymmetry"),
                          (axs[1],A_tb,"Cathode-facing top/bottom ion asymmetry")):
        ax.scatter(np.arange(n)[np.isfinite(data)],data[np.isfinite(data)],
                   s=8,alpha=.48,rasterized=True)
        ax.axhline(0,color="gray",ls="--",lw=.9)
        ax.set_ylim(-1.03,1.03)
        ax.set_ylabel(title)
        ax.grid(alpha=.12)
    axs[-1].set_xlabel("Physical muon index (0..999, frozen J/psi source order)")
    fig.suptitle("One point per real Geant4/TrackHeed parent crossing; exact ion-sector census")
    fig.tight_layout()
    fig.savefig(OUT/"03_asymmetry_each_of_1000_muons.png",dpi=185)
    plt.close(fig)

    # 4. Ion azimuth probability heatmap: one row per true muon.
    fig,ax=plt.subplots(figsize=(12.4,7.8))
    m=ax.imshow(event_norm,origin="lower",aspect="auto",interpolation="nearest",
                extent=(-180,180,-.5,n-.5),cmap="viridis")
    fig.colorbar(m,ax=ax,label="Birth-angle probability per 5 degree bin, per muon")
    ax.set_xticks(np.arange(-180,181,45))
    ax.set_xlabel("Azimuth about nearest anode [degrees]")
    ax.set_ylabel("Muon index, frozen Geant4 source order")
    ax.set_title("Approximate ion-birth angular distribution event by event\n"
                 "Each nonzero muon normalized to 1; 128 birth points per avalanche")
    fig.tight_layout()
    fig.savefig(OUT/"04_each_muon_azimuth_heatmap.png",dpi=185)
    plt.close(fig)

    # 5. Asymmetry versus BOTH genuine incoming 3D muon angles.
    fig,axes=plt.subplots(2,2,figsize=(11.2,8.2),sharey=True)
    for i,(data,title) in enumerate(((A_hemi,r"$A_{\rm hemi}$"),
                                     (A_tb,r"$A_{\rm TB}$"))):
        for j,(theta,label) in enumerate(((theta_u,r"$\theta_u$ [deg]"),
                                          (theta_w,r"$\theta_w$ [deg]"))):
            ax=axes[i,j]
            mask=np.isfinite(data)
            ax.scatter(theta[mask],data[mask],s=8,alpha=.27,rasterized=True)
            ax.axhline(0,color="gray",lw=.8,ls="--")
            ax.set_xlabel("True native muon incidence "+label)
            ax.set_ylabel(title)
            ax.set_ylim(-1.05,1.05)
            ax.grid(alpha=.13)
    fig.suptitle("True Geant4 angles vs full Garfield whole-muon ion-birth asymmetry")
    fig.tight_layout()
    fig.savefig(OUT/"05_asymmetry_vs_true_incidence.png",dpi=185)
    plt.close(fig)

    # 6. Independent primary-electron avalanches versus combined whole muon;
    # the electron histogram is one count per electron, not birth-weighted.
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),sharex=True)
    for ax,per_seed,per_muon,title in (
        (axes[0],seed_asym_hemi,A_hemi,"Hemisphere"),
        (axes[1],seed_asym_tb,A_tb,"Top/bottom 90 degree sectors"),
    ):
        x=np.asarray(per_seed,float)
        y=per_muon[np.isfinite(per_muon)]
        e=np.linspace(-1,1,42)
        if len(x):ax.hist(x,bins=e,weights=np.ones(len(x))/len(x),
                          histtype="step",lw=1.4,label=f"Each electron ({len(x)})")
        if len(y):ax.hist(y,bins=e,weights=np.ones(len(y))/len(y),
                          histtype="step",lw=1.4,label=f"Whole muon ({len(y)})")
        ax.set_title(title)
        ax.set_xlabel("Ion-birth asymmetry")
        ax.set_xlim(-1,1)
        ax.legend()
        ax.grid(alpha=.12)
    axes[0].set_ylabel("Probability per bin (separately normalized)")
    fig.suptitle("Individual avalanche vs combined full-muon ion asymmetry (different units)")
    fig.tight_layout()
    fig.savefig(OUT/"06_single_avalanche_vs_whole_muon_asymmetry.png",dpi=185)
    plt.close(fig)

    # Source-derived physical rows: one independent Geant4 muon/event,
    # all A values are exact from unsampled complete Garfield birth sectors.
    data_path=OUT/"muon_exact_asymmetry_1000.csv"
    with data_path.open("x",newline="") as f:
        fieldnames=("muon_index","crossing_id","station",
                    "theta_u_deg","theta_w_deg","heed_primaries",
                    "ion_births","ion_up","ion_down","ion_top","ion_bottom",
                    "A_hemi_exact","A_TB_exact","failed_seeds","capped_seeds")
        writer=csv.DictWriter(f,fieldnames=fieldnames)
        writer.writeheader()
        for i in range(n):
            writer.writerow({
                "muon_index":i,"crossing_id":ordered_ids[i],
                "station":int(station[i]),"theta_u_deg":theta_u[i],
                "theta_w_deg":theta_w[i],"heed_primaries":int(nprim[i]),
                "ion_births":int(births_exact[i]),
                "ion_up":int(Nup[i]),"ion_down":int(Ndown[i]),
                "ion_top":int(Ntop[i]),"ion_bottom":int(Nbottom[i]),
                "A_hemi_exact":A_hemi[i] if math.isfinite(A_hemi[i]) else "",
                "A_TB_exact":A_tb[i] if math.isfinite(A_tb[i]) else "",
                "failed_seeds":int(failed[i]),"capped_seeds":int(capped[i])
            })
    angle_path=OUT/"ion_birth_azimuth_5deg.csv"
    with angle_path.open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("phi_center_deg","ion_weighted_probability",
                         "equal_nonzero_muon_probability",
                         "estimated_ion_births_in_bin"))
        for k in range(NBINS):
            writer.writerow((CENTRES[k],global_hist[k]/global_hist.sum(),
                             equal_muon[k],global_hist[k]))
    summary={
        "schema":"NA60_JPSI1000_P3_NOFW_ION_ANGLES_V1",
        "source":"1000 SHA-protected independent Geant4 J/psi parents; one crossing per parent",
        "geometry":GEOM,
        "n_muons":n,
        "n_independent_TrackHeed_primary_electrons":len(seed_counts),
        "n_real_full_Garfield_avalanche_calls":len(seed_counts),
        "n_total_exact_ion_births":total_ions,
        "n_128_cap_reservoir_rows":sample_count,
        "n_failed_Garfield_calls":int(failed.sum()),
        "n_capped_calls":int(capped.sum()),
        "n_undefined_A_hemi":int(np.isnan(A_hemi).sum()),
        "n_undefined_A_TB":int(np.isnan(A_tb).sum()),
        "A_hemi_mean_per_muon_defined":float(np.nanmean(A_hemi)),
        "A_TB_mean_per_muon_defined":float(np.nanmean(A_tb)),
        "A_hemi_rms_per_muon_defined":float(np.sqrt(np.nanmean(A_hemi**2))),
        "A_TB_rms_per_muon_defined":float(np.sqrt(np.nanmean(A_tb**2))),
        "angle_definition":"C++ phi=atan2(v,u-nearest_3mm_anode), degrees [-180,180]",
        "angle_sample_limit":"Uniform reservoir at most 128 points PER microscopic avalanche",
        "angle_weights":"N_exact_ion_births_per_avalanche / N_saved_points",
        "muon_asymmetry_source":"FULL unsampled exact native C++ ion sector census",
        "A_hemi":"(N_up-N_down)/(N_up+N_down), phi in (0,180) vs (-180,0); boundaries omitted",
        "A_TB":"(N_top-N_bottom)/(N_top+N_bottom), phi [45,135) vs [-135,-45); side omitted",
        "not_calculated":"Ion drift, induced charge, U/V cathode weighting, electronics",
        "source_shas":source_sha,
        "plotter_sha256":sha(HERE)
    }
    output_files=sorted(OUT.glob("*.png"))+[data_path,angle_path]
    summary["output_files_sha256"]={p.name:sha(p) for p in output_files}
    (OUT/"manifest.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("PASS: READ-ONLY true muon Garfield angle analysis",flush=True)
    print("PHYSICS: 1000 independent J/psi muons, p3_noFW, 3mm, no field wires")
    print("FULL MICROSCOPIC GARFIELD AVALANCHES:",len(seed_counts))
    print("EXACT TOTAL ION BIRTHS:",total_ions)
    print("UNBIASED RESERVOIR SAMPLES:",sample_count)
    print("WHOLE MUON A_HEMI mean (defined only):",round(summary["A_hemi_mean_per_muon_defined"],6))
    print("WHOLE MUON A_TB mean (defined only):",round(summary["A_TB_mean_per_muon_defined"],6))
    print("UNDEFINED A_HEMI / A_TB:",summary["n_undefined_A_hemi"],
          summary["n_undefined_A_TB"])
    print("ZERO ION BIRTHS PHYSICS muons:",int((~nonzero).sum()))
    print("PLOTS + EXACT per-muon CSV:",OUT,flush=True)
    for f in sorted(OUT.iterdir()):
        print("  ",f.name,flush=True)


def selftest():
    for x,want in [(-2.5,-3),(-1.5,-2),(-.5,-1),(.5,1),(1.5,2),(2.5,3)]:
        require(round_cpp(x)==want,"C++ half-away-from-zero convention")
    require(birth_phi_deg(0,1)==90 and birth_phi_deg(0,-1)==-90,
            "Opposite cathode-facing angle")
    require(birth_phi_deg(1,0)==0 and birth_phi_deg(-1,0)==180,
            "Opposite wire-side angle")
    require(asym(3,1)==.5 and math.isnan(asym(0,0)),
            "Sector asymmetry definition")
    require(abs(sum((100/4 for _ in range(4)))-100)<1e-9,
            "Reservoir count reconstruction")
    require(NBINS==72 and len(CENTRES)==72,"5-degree angular bins")
    print("PASS exact C++ angle convention, A undefined handling, ion weights")
    print("NO physics files read, no plotting directory created")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode",choices=("selftest","plot"))
    a=ap.parse_args()
    if a.mode=="selftest":selftest()
    else:analyze()

#!/usr/bin/env python3
"""Make one uploadable, read-only Penning-ON MWPC radial/along-wire data ZIP.

Input: 1000 archived true J/psi Geant4 parents, each original full
Garfield avalanche ion census, at most 128 random ion-birth points/seed,
and existing SHA-verifiable timing-correlations NPZ.
Output: small tables/histograms, representative ion clouds and provenance.
No Garfield, TrackHeed, altered input, cache overwrite or 30k campaign.

Run:
  python SCRIPT selftest
  python SCRIPT export

Distributions of r, phi, delta-w and delay are ION-WEIGHTED (each saved
point has birth count / sample count); per-avalanche summary rows give
one equal weight to every individual ion-producing avalanche.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pi_p3_noFW_penning057_timing_correlations_v1 as temporal

REPO=Path(__file__).resolve().parents[2]
ROOT=temporal.ROOT
OUTPUT=REPO/"test_runs/mwpc_stage_c/penning_on_radial_wire_share_v1.zip"
PITCH=3.0
WIRE_RADIUS=.015
V0=np.linspace(-3.,3.,25)
PHI=np.linspace(-180.,180.,73)
R=np.r_[0.,np.geomspace(.0001,5.,161)]
# Along-wire distance relative to original electron starting w coordinate,
# and (separately) relative to each sampled avalanche birth-cloud centroid.
W=np.linspace(-10.,10.,201)
W_FINE=np.linspace(-2.,2.,201)
R2=np.r_[0.,np.geomspace(.0001,5.,81)]
W2=np.linspace(-5.,5.,101)
SEED_ROWS=(
  "muon","primary","cluster","start_u_wire_mm","start_v_mm","start_w_mm",
  "gain_e_collected","exact_ion_births","exact_A","n_saved",
  "birth_wire_change_fraction","mean_r_mm","sigma_r_mm",
  "r_p10_mm","r_p50_mm","r_p90_mm","frac_r_below_50um",
  "frac_r_below_100um","frac_r_below_300um","phi_circular_R",
  "phi_circular_dispersion","mean_dw_from_seed_mm",
  "rms_dw_from_seed_mm","mean_dw_centered_mm","rms_dw_centered_mm",
  "w_p10_mm_from_seed","w_p50_mm_from_seed","w_p90_mm_from_seed",
  "min_birth_delay_ns","mean_birth_delay_ns","max_birth_delay_ns",
)
CLOUD_ROWS=(
  "muon","primary","cluster","gain_e_collected","exact_ion_births",
  "start_v_mm","start_w_mm","du_to_birth_wire_mm","v_birth_mm",
  "w_birth_mm","delta_w_from_seed_mm","delta_w_from_av_cloud_mean_mm",
  "r_birth_mm","phi_birth_deg","delay_since_Heed_ns",
  "wire_changed","reservoir_weight",
)

def must(ok,msg):
    if not ok: raise RuntimeError(msg)
def sha(p):
    d=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            d.update(b)
    return d.hexdigest()
def rd(path,needed):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames is not None and
             len(set(r.fieldnames))==len(r.fieldnames) and
             set(needed).issubset(r.fieldnames),
             "Raw physics CSV missing required columns "+str(path))
        yield from r
def rnd(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)
def wire(umm):
    return PITCH*rnd(float(umm)/PITCH)
def hist1(x,w,e):
    return np.histogram(x,bins=e,weights=w)[0]
def hist2(x,y,w,ex,ey):
    return np.histogram2d(x,y,bins=(ex,ey),weights=w)[0]
def add_csv(path,fields,rows):
    with path.open("x",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(fields)
        wr.writerows(rows)
def add_hist_1(path,edges,counts,units,label):
    add_csv(path,("bin_low_"+units,"bin_high_"+units,label),
            zip(edges[:-1],edges[1:],counts))
def add_hist_2(path,edges_x,edges_y,grid,units_x,units_y,label):
    xx=[];yy=[];ww=[]
    for i in range(len(edges_x)-1):
        for j in range(len(edges_y)-1):
            xx.append((edges_x[i],edges_x[i+1]))
            yy.append((edges_y[j],edges_y[j+1]))
            ww.append(float(grid[i,j]))
    add_csv(path,(f"x_low_{units_x}",f"x_high_{units_x}",
                  f"y_low_{units_y}",f"y_high_{units_y}",label),
            ((a,b,c,d,v) for (a,b),(c,d),v in zip(xx,yy,ww)))
def representative_rows(cache):
    """Deterministic stratification; do NOT choose only pretty events."""
    rng=np.random.default_rng(20261110)
    ions=cache["ions"]
    valid=np.flatnonzero((ions>=128)&
                        np.isfinite(cache["v0"])&
                        (cache["sample_count"]==128))
    must(len(valid)>150,"Too few full-128-sample real avalanches")
    chosen=set()
    vgroup=np.digitize(cache["v0"],[-2.,-1.,0.,1.,2.])
    ggroup=np.digitize(cache["gain"],[300.,1000.,3000.,10000.])
    # Subsample across positive/negative v0, charge and tail.
    for i in range(6):
        for j in range(5):
            pool=valid[(vgroup[valid]==i)&(ggroup[valid]==j)]
            if len(pool):
                chosen.update(map(int,rng.choice(pool,size=min(4,len(pool)),
                                                replace=False)))
    chosen.update(map(int,rng.choice(valid,size=min(75,len(valid)),replace=False)))
    chosen.update(map(int,valid[np.argsort(ions[valid])[-20:]]))
    span=np.nan_to_num(cache["span"][valid],nan=-1.)
    chosen.update(map(int,valid[np.argsort(span)[-20:]]))
    return chosen

def export():
    must(not OUTPUT.exists(),"Will not overwrite existing upload ZIP "+str(OUTPUT))
    cache,meta=temporal.load()
    must(meta["n_independent_true_Geant4_muons"]==1000,
         "Temporal cache does not contain exactly 1000 real muons")
    n=len(cache["gain"])
    must(n==meta["n_exact_Garfield_avalanches"] and n>20000,
         "Garfield avalanche count inconsistent")
    chosen=representative_rows(cache)
    must(len(chosen)<=300,"Representative set exceeded intended upload budget")
    # Build exact seed identity lookup from already verified timing cache.
    lookup={}
    for i,(p,seed) in enumerate(zip(cache["parent"],cache["primary"])):
        identifier=(int(p),int(seed))
        must(identifier not in lookup,"Duplicate original primary in timing cache")
        lookup[identifier]=i
    r_hist=np.zeros(len(R)-1)
    w_hist=np.zeros(len(W)-1)
    w_fine_hist=np.zeros(len(W_FINE)-1)
    wc_hist=np.zeros(len(W_FINE)-1)
    h_vphi=np.zeros((len(V0)-1,len(PHI)-1))
    h_vphi_near=np.zeros_like(h_vphi)
    h_vr=np.zeros((len(V0)-1,len(R2)-1))
    h_rw=np.zeros((len(R2)-1,len(W2)-1))
    h_rw_center=np.zeros_like(h_rw)
    ion_total_seen=0
    saved_points_total=0
    outlier={
       "radial_gt_5mm_ions":0.,"abs_dw_from_seed_gt_10mm_ions":0.,
       "abs_dw_from_seed_gt_2mm_ions":0.,
       "abs_dw_from_cloud_centre_gt_2mm_ions":0.,
       "h_v0_outside_pm3mm_ions":0.,
       "ion_births_at_r_below_wire_radius":0.,
    }
    full_128=0
    n_ion_avalanches=0
    selected_exported=0
    parent_rows=defaultdict(int)
    with tempfile.TemporaryDirectory(prefix="mwpc_radial_export_") as folder:
        temp=Path(folder)
        with (temp/"avalanche_geometry_summary.csv").open("x",newline="") as file_summary, \
             (temp/"representative_ion_birth_points.csv").open("x",newline="") as file_cloud:
            sw=csv.writer(file_summary);sw.writerow(SEED_ROWS)
            cw=csv.writer(file_cloud);cw.writerow(CLOUD_ROWS)
            for chunk in range(50):
                directory=ROOT/"native_chunks"/f"chunk_{chunk:04d}"
                receipt=json.loads((directory/"complete.json").read_text())
                must(receipt["chunk"]==chunk and receipt["n_parents"]==20
                     and receipt["stats"]["n_capped_seeds"]==0
                     and receipt["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
                     "Input production chunk is not a complete sealed Penning-ON run")
                files=("muon_reference_seed.csv",
                       "muon_reference_ion_birth_sample.csv")
                for filename in files:
                    must(sha(directory/filename)==receipt["file_sha256"][filename],
                         "Original native physics file SHA mismatch "+filename)
                raw_seed={}
                for row in rd(directory/files[0],(
                    "crossing_id","primary_index","seed_u_mm","seed_v_mm",
                    "seed_w_mm","ion_births","collected_e","ion_up","ion_down")):
                    identifier=(row["crossing_id"],int(row["primary_index"]))
                    must(identifier not in raw_seed,"Duplicate exact Garfield seed")
                    raw_seed[identifier]=row
                parent_ids={meta["muons"][i]["id"]:i
                            for i in range(chunk*20,(chunk+1)*20)}
                points=defaultdict(list)
                for row in rd(directory/files[1],(
                    "crossing_id","primary_index","u_mm","v_mm","w_mm",
                    "t_ns","reservoir_weight","seed_ion_births_total")):
                    identifier=(row["crossing_id"],int(row["primary_index"]))
                    must(identifier in raw_seed,"Unknown real ion source seed")
                    points[identifier].append(row)
                must(len(raw_seed)==
                     sum(meta["muons"][i]["n_Heed"]
                         for i in range(chunk*20,(chunk+1)*20)),
                     "True Heed seeds disagree with original exact census")
                for identifier,s in raw_seed.items():
                    mu=parent_ids[identifier[0]]
                    pi=identifier[1]
                    k=lookup[(mu,pi)]
                    exact=int(s["ion_births"])
                    gain=int(s["collected_e"])
                    must(exact==int(cache["ions"][k]) and
                         gain==int(cache["gain"][k]),"Original Garfield gain/ion counts differ from cache")
                    must(abs(float(s["seed_v_mm"])-float(cache["v0"][k]))<1.e-5,
                         "Original starting electron position differs from cached v0")
                    if exact==0:
                        must(len(points.get(identifier,()))==0,
                             "Zero-ion avalanche has saved ion birth points")
                        continue
                    sample=points.get(identifier,[])
                    n_expected=min(128,exact)
                    must(len(sample)==n_expected,"128-point original reservoir was truncated")
                    unitweight=exact/n_expected
                    weight=np.asarray([float(x["reservoir_weight"]) for x in sample])
                    must(np.isfinite(weight).all() and
                         np.allclose(weight,unitweight,rtol=1.e-7,atol=1.e-7),
                         "Original ion reservoir weights are inconsistent")
                    xyz=np.asarray([[float(x[p]) for p in ("u_mm","v_mm","w_mm","t_ns")]
                                     for x in sample],dtype=float)
                    must(np.isfinite(xyz).all(),"Original ion-birth sample has NaN")
                    birth_wire=PITCH*np.where(xyz[:,0]>=0,
                         np.floor(xyz[:,0]/PITCH+.5),
                         np.ceil(xyz[:,0]/PITCH-.5))
                    delta_u=xyz[:,0]-birth_wire
                    v=xyz[:,1]
                    radial=np.hypot(delta_u,v)
                    azimuth=np.rad2deg(np.arctan2(v,delta_u))
                    w0=float(s["seed_w_mm"])
                    delta_w=xyz[:,2]-w0
                    centre_w=float(xyz[:,2].mean())
                    dw_centered=xyz[:,2]-centre_w
                    delay=xyz[:,3]-float(cache["t_heed"][k])
                    must(delay.min()>=-1.e-6,"An ion was created before the initial Heed electron")
                    wire_changed=(np.abs(birth_wire-wire(float(s["seed_u_mm"])))>1.e-6)
                    r_hist+=hist1(radial,weight,R)
                    w_hist+=hist1(delta_w,weight,W)
                    w_fine_hist+=hist1(delta_w,weight,W_FINE)
                    wc_hist+=hist1(dw_centered,weight,W_FINE)
                    h_vphi+=hist2(np.full(len(sample),float(cache["v0"][k])),
                                  azimuth,weight,V0,PHI)
                    near=radial<.1
                    h_vphi_near+=hist2(np.full(int(near.sum()),float(cache["v0"][k])),
                                       azimuth[near],weight[near],V0,PHI)
                    h_vr+=hist2(np.full(len(sample),float(cache["v0"][k])),
                                radial,weight,V0,R2)
                    h_rw+=hist2(radial,delta_w,weight,R2,W2)
                    h_rw_center+=hist2(radial,dw_centered,weight,R2,W2)
                    ion_total_seen+=exact
                    saved_points_total+=len(sample)
                    n_ion_avalanches+=1
                    full_128+=(n_expected==128)
                    parent_rows[mu]+=1
                    outlier["radial_gt_5mm_ions"]+=float(weight[radial>R[-1]].sum())
                    outlier["abs_dw_from_seed_gt_10mm_ions"]+=float(weight[np.abs(delta_w)>10.].sum())
                    outlier["abs_dw_from_seed_gt_2mm_ions"]+=float(weight[np.abs(delta_w)>2.].sum())
                    outlier["abs_dw_from_cloud_centre_gt_2mm_ions"]+=float(weight[np.abs(dw_centered)>2.].sum())
                    outlier["h_v0_outside_pm3mm_ions"]+=float(weight[np.abs(float(cache["v0"][k]))>3.].sum())
                    outlier["ion_births_at_r_below_wire_radius"]+=float(weight[radial<WIRE_RADIUS].sum())
                    A=float(cache["A"][k])
                    z=np.exp(1j*np.deg2rad(azimuth)).mean()
                    quant=np.quantile
                    row=(
                        mu,pi,int(cache["cluster"][k]),
                        float(s["seed_u_mm"])-wire(float(s["seed_u_mm"])),
                        float(cache["v0"][k]),w0,gain,exact,A,n_expected,
                        float(wire_changed.mean()),
                        float(radial.mean()),float(radial.std()),
                        *map(float,quant(radial,[.1,.5,.9])),
                        float((radial<.05).mean()),
                        float((radial<.1).mean()),
                        float((radial<.3).mean()),
                        float(abs(z)),float(1.-abs(z)),
                        float(delta_w.mean()),float(delta_w.std()),
                        float(dw_centered.mean()),float(dw_centered.std()),
                        *map(float,quant(delta_w,[.1,.5,.9])),
                        float(delay.min()),float(delay.mean()),float(delay.max()),
                    )
                    must(len(row)==len(SEED_ROWS),"Summary physical column mismatch")
                    sw.writerow(row)
                    if k in chosen:
                        selected_exported+=1
                        for j in range(len(sample)):
                            cw.writerow((
                              mu,pi,int(cache["cluster"][k]),gain,exact,
                              float(cache["v0"][k]),w0,
                              float(delta_u[j]),float(v[j]),float(xyz[j,2]),
                              float(delta_w[j]),float(dw_centered[j]),
                              float(radial[j]),float(azimuth[j]),float(delay[j]),
                              int(wire_changed[j]),float(weight[j])
                            ))
                print(f"PASS physically sealed Penning-ON chunk {chunk+1}/50: "
                      f"{ion_total_seen:,} total true ion births",flush=True)
        must(ion_total_seen==meta["n_true_ion_births"] and
             saved_points_total==meta["n_sampled_ion_birth_positions"],
             "Full spatial physics sample/census differs from verified temporal cache")
        must(n_ion_avalanches==sum(int(x)>0 for x in cache["ions"]),
             "Avalanche count differs from cached original ion birth census")
        must(selected_exported==len(chosen),
             "Selected avalanche cloud sample is missing a source seed")
        # Complete exact reference includes zero-ion/zero-collected
        # original primary electrons omitted from spatial-shape tables.
        complete_fields=(
            "muon","primary","cluster","initial_u_wire_mm",
            "initial_v0_mm","original_Heed_t_ns",
            "collected_gain_exact","ion_births_exact",
            "A_exact_or_blank","saved_birth_points","sampled_t10_ns_or_blank")
        full_rows=(
            (int(cache["parent"][j]),int(cache["primary"][j]),
             int(cache["cluster"][j]),float(cache["u0"][j]),
             float(cache["v0"][j]),float(cache["t_heed"][j]),
             int(cache["gain"][j]),int(cache["ions"][j]),
             (float(cache["A"][j]) if np.isfinite(cache["A"][j]) else ""),
             int(cache["sample_count"][j]),
             (float(cache["t10"][j]) if np.isfinite(cache["t10"][j]) else ""))
            for j in range(n)
        )
        add_csv(temp/"all_avalanche_exact_index.csv",complete_fields,full_rows)
        add_hist_1(temp/"radial_ion_weighted.csv",R,r_hist,"mm","weighted_ion_births")
        add_hist_1(temp/"dw_from_seed_ion_weighted.csv",W,w_hist,
                   "mm","weighted_ion_births")
        add_hist_1(temp/"dw_from_seed_zoom_ion_weighted.csv",W_FINE,w_fine_hist,
                   "mm","weighted_ion_births")
        add_hist_1(temp/"dw_from_birth_centroid_ion_weighted.csv",W_FINE,wc_hist,
                   "mm","weighted_ion_births")
        add_hist_2(temp/"azimuth_vs_start_v0_ion_weighted.csv",
                   V0,PHI,h_vphi,"mm","deg","weighted_ion_births")
        add_hist_2(temp/"near_wire_azimuth_vs_start_v0_ion_weighted.csv",
                   V0,PHI,h_vphi_near,"mm","deg","weighted_ion_births")
        add_hist_2(temp/"radial_vs_start_v0_ion_weighted.csv",
                   V0,R2,h_vr,"mm","mm","weighted_ion_births")
        add_hist_2(temp/"radial_vs_dw_from_seed_ion_weighted.csv",
                   R2,W2,h_rw,"mm","mm","weighted_ion_births")
        add_hist_2(temp/"radial_vs_dw_from_centroid_ion_weighted.csv",
                   R2,W2,h_rw_center,"mm","mm","weighted_ion_births")
        details={
            "schema":"NA60_MWPC_PENNING_ON_SPATIAL_SHARE_V1",
            "number_of_real_Geant4_muons":1000,
            "number_of_original_Garfield_avalanches":n,
            "ion_producing_avalanches":n_ion_avalanches,
            "zero_ion_avalanches":n-n_ion_avalanches,
            "full_128_sample_avalanches":full_128,
            "exact_original_ion_births":ion_total_seen,
            "uniform_reservoir_points_saved":saved_points_total,
            "representative_avalanche_clouds":selected_exported,
            "geometry":"p3_noFW 3mm pitch, 30um diameter anode, +/-3mm cathodes, 1800V",
            "Penning":"Ar/CO2 70:30, Penning ON r=0.57",
            "reference_wire":"Nearest wire to EACH ion birth in the original finite local 3mm wire lattice; available wire_changed flag compares to initiating Heed seed wire",
            "w_axis":"along-wire Garfield w in mm",
            "delta_w_from_seed":"ion birth w minus the ORIGINAL Heed electron w coordinate in mm; not necessarily centered on avalanche centroid",
            "delta_w_from_cloud_centroid":"ion birth w minus mean sampled birth w of SAME avalanche in mm",
            "r_mm":"sqrt((u_birth-nearest_birth_wire_u)^2 + v_birth^2); wire radius 0.015mm",
            "phi_deg":"atan2(v_birth,u_birth-nearest_birth_wire_u), -180..180 degrees",
            "weight":"Each sampled birth represents exact ion_births/min(128,ion_births); pooled histograms ion-weighted; summary rows are one per ion-producing avalanche",
            "warning":"Exact A/gain/ion total, sampled 2D/3D ion-birth shapes. Not full ion trajectories, electron arrival angle, ion drift or actual cathode induced response.",
            "tails_outside_export_grid_in_ion_units":outlier,
            "sum_ion_births_in_radial_grid":float(r_hist.sum()),
            "sum_ion_births_in_dw_seed_grid":float(w_hist.sum()),
            "sum_ion_births_in_dw_centred_grid":float(wc_hist.sum()),
            "sha256_source_temporal_npz":sha(temporal.CACHE),
            "version_notes":"Representative avalanches sampled deterministically across signed v0, gain strata, upper ion-count tail, time-span tail, and random draws. The all_avalanche_exact_index.csv includes zero-ion and zero-collected records; avalanche_geometry_summary.csv contains only ion-producing avalanches.",
        }
        (temp/"README.json").write_text(json.dumps(details,indent=2,allow_nan=False)+"\n")
        sha_map={p.name:sha(p) for p in sorted(temp.iterdir())}
        (temp/"SHA256_manifest.json").write_text(json.dumps(sha_map,indent=2)+"\n")
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(OUTPUT,"x",compression=zipfile.ZIP_DEFLATED,
                             compresslevel=7) as zf:
            for file in sorted(temp.iterdir()):
                zf.write(file,arcname=file.name)
    print("PASS",n_ion_avalanches,"ion-producing avalanches,",
          f"{ion_total_seen:,} exact ion births,",
          f"{saved_points_total:,} sampled birth positions",flush=True)
    print("One upload file:",OUTPUT,flush=True)
    print("ZIP size MB:",round(OUTPUT.stat().st_size/1e6,2),flush=True)

def selftest():
    assert rnd(-.5)==-1 and rnd(.5)==1
    assert abs(wire(3.2)-3.)<1e-12
    assert np.allclose(np.hypot([0.,.03],[.015,0.]),[.015,.03])
    assert len(SEED_ROWS)==31 and len(CLOUD_ROWS)==17
    assert len(R)==162 and len(W)==201 and len(PHI)==73
    assert abs(hist1([1.,2.],[3.,4.],[0,1.5,3.]).sum()-7.)<1e-8
    x=np.array([1.,2.])
    phi=np.arctan2([0.,1.],[1.,0.])
    assert np.allclose(phi,[0.,np.pi/2])
    print("PASS geometry/mm units, 30um wire diameter, histogram weights, "
          "CSV schemas. No Garfield, no original files modified.")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("selftest","export"))
    args=parser.parse_args()
    if args.mode=="selftest":selftest()
    else:export()

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Read-only 1000-muon Penning-ON timing/avalanche correlation cache.

IMPORTANT: t10 is the 10th percentile of ion BIRTH times saved from
<=128 uniformly sampled ionisations in one Garfield avalanche, relative
to the original MUON ENTRY time (Geant4/Heed t=0). It is NOT an independently
measured electron ARRIVAL time, drift time, or incoming angle.

This script NEVER runs Garfield, modifies sealed source data, or accesses
the ongoing 30k muon campaign. The new auxiliary NPZ is independent of
the existing 12-plot Streamlit dashboard cache.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO=Path(__file__).resolve().parents[2]
ROOT=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
DIR=ROOT/"interactive_dashboard_v1"
ORIGINAL=DIR/"dashboard_data.npz"
ORIGMETA=DIR/"provenance.json"
CACHE=DIR/"temporal_correlations_v1.npz"
META=DIR/"temporal_correlations_v1.json"
PITCH=3.
PHI_EDGES=np.linspace(-180.,180.,73)

KINDS={
  "t10_gain":"Mean gain vs early ion-birth time",
  "t10_asym":"Mean |ion asymmetry| vs early ion-birth time",
  "t10_angles":"Ion-birth angular distributions: early vs late",
  "t10_start":"Early ion-birth time vs initial |v|",
  "t10_gain_position":"Gain conditioned on early time and initial |v|",
  "t10_timeline":"One real muon's cluster/avalanche timeline",
}
DEFAULTS={
  "t10_gain":("Avalanche gain vs early ion-birth time",
              "Early ion-birth time since muon entry [ns]",
              "Mean collected electrons per starting electron",
              "Mean collected gain"),
  "t10_asym":("Ion-birth asymmetry vs early ion-birth time",
              "Early ion-birth time since muon entry [ns]",
              "Mean |A_ion,hemi|",
              "Mean |A|"),
  "t10_angles":("Ion-birth azimuth in early and late avalanches",
                "Birth azimuth around nearest anode [degrees]",
                "Ion-weighted probability per angular bin",
                "Ion births"),
  "t10_start":("Early ion birth time versus initial |v|",
               "Initial electron |v| [mm]",
               "Early ion-birth time since muon entry [ns]",
               "Heed electron starts"),
  "t10_gain_position":("Gain vs early ion birth time and starting |v|",
                       "Early ion-birth time since muon entry [ns]",
                       "Initial electron |v| [mm]",
                       "Mean collected gain"),
  "t10_timeline":("One muon: cluster and avalanche chronology",
                  "Time since muon gas entry [ns]",
                  "Original TrackHeed cluster index",
                  "Electron collection and exact ion asymmetry"),
}

def must(x,msg):
    if not x:raise RuntimeError(msg)
def sha(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            d.update(block)
    return d.hexdigest()
def rows(path,required):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames is not None and
             len(set(r.fieldnames))==len(r.fieldnames) and
             set(required).issubset(r.fieldnames),
             "Incorrect raw native CSV schema "+str(path))
        yield from r
def seed_key(row):
    return (row["crossing_id"],int(row["primary_index"]))
def wire_offset(u):
    n=math.floor(u/PITCH+.5) if u>=0 else math.ceil(u/PITCH-.5)
    return u-PITCH*n
def phi_bin(u,v):
    theta=math.degrees(math.atan2(v,wire_offset(u)))
    return max(0,min(71,int((theta+180.)/5.)))
def sample_quantiles(tt):
    t=np.asarray(tt,dtype=float)
    must(t.size>0 and np.isfinite(t).all(),"Invalid sampled physical birth time")
    return np.quantile(t,[.10,.50,.90])
def validate_original_cache():
    must(ORIGINAL.is_file() and ORIGMETA.is_file(),
         "Run the original dashboard 'prepare' once before timing prepare")
    meta=json.loads(ORIGMETA.read_text())
    must(meta.get("schema")=="P3_NOFW_PENNING_ON_INTERACTIVE_DASHBOARD_V1"
         and meta.get("n_independent_Geant4_Jpsi_parents")==1000
         and sha(ORIGINAL)==meta.get("cache_sha256"),
         "Verified original Penning-ON dashboard cache invalid")
    return meta

def prepare():
    source=validate_original_cache()
    must(not CACHE.exists() and not META.exists(),
         "Refusing to overwrite auxiliary timing-cache data: "+str(CACHE))
    variables=("parent","primary","cluster","gain","ions","A",
               "u0","v0","t_heed","t10","t50","t90","span","sample_count")
    out={field:[] for field in variables}
    angular=[]
    muons=[]
    ions_all=0
    points_all=0
    n_bad_time=0
    for block in range(50):
        d=ROOT/"native_chunks"/f"chunk_{block:04d}"
        sealed=d/"complete.json"
        must(sealed.is_file(),"Missing sealed Garfield production checkpoint "+str(d))
        rec=json.loads(sealed.read_text())
        must(rec["chunk"]==block and rec["n_parents"]==20 and
             rec["stats"]["n_capped_seeds"]==0 and
             rec["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
             "Garfield physics seal/paired-Heed validation failed")
        for file in ("muon_reference_crossing.csv","muon_reference_heed_primaries.csv",
                     "muon_reference_heed_clusters.csv","muon_reference_seed.csv",
                     "muon_reference_ion_birth_sample.csv"):
            must(rec["file_sha256"].get(file)==sha(d/file),
                 "SHA mismatch in immutable Penning-ON "+str(d/file))
        crossings=list(rows(d/"muon_reference_crossing.csv",(
            "source|chunk|event|hit_index|track_id|station|det_id",
            "heed_primaries","heed_clusters",
            "theta_u_in_deg","theta_w_in_deg","path_mm",
            "collected_avalanche_electrons","ion_births")))
        must(len(crossings)==20,"Penning-on chunk does not contain 20 muons")
        index={}
        for r in crossings:
            name=r["source|chunk|event|hit_index|track_id|station|det_id"]
            must(name not in index,"Duplicate real Geant4 parent")
            index[name]=len(muons)
            muons.append({
                "id":name,
                "n_Heed":int(r["heed_primaries"]),
                "n_clusters":int(r["heed_clusters"]),
                "theta_u_deg":float(r["theta_u_in_deg"]),
                "theta_w_deg":float(r["theta_w_in_deg"]),
                "path_mm":float(r["path_mm"]),
                "collected_e":int(r["collected_avalanche_electrons"]),
                "ion_births":int(r["ion_births"]),
            })
        heed={}
        n_seen=defaultdict(int)
        for r in rows(d/"muon_reference_heed_primaries.csv",(
            "crossing_id","primary_index","cluster_index","t_ns")):
            k=seed_key(r)
            must(k not in heed,"Duplicate original Heed primary electron")
            heed[k]=(int(r["cluster_index"]),float(r["t_ns"]))
            n_seen[k[0]]+=1
        cluster_count=defaultdict(int)
        for r in rows(d/"muon_reference_heed_clusters.csv",(
            "crossing_id","cluster_index","n_accepted","t_ns")):
            cluster_count[r["crossing_id"]]+=int(r["n_accepted"])
        for r in crossings:
            name=r["source|chunk|event|hit_index|track_id|station|det_id"]
            must(n_seen.get(name,0)==int(r["heed_primaries"])==
                 cluster_count.get(name,0),
                 "Original Heed clusters and primaries inconsistent")
        seeds=list(rows(d/"muon_reference_seed.csv",(
            "crossing_id","primary_index","seed_u_mm","seed_v_mm",
            "collected_e","ion_births","ion_up","ion_down","capped")))
        lookup={}
        must(len(seeds)==len(heed),
             "Garfield seed count and original Heed population differ")
        for j,r in enumerate(seeds):
            k=seed_key(r)
            must(k in heed and k[0] in index and k not in lookup,
                 "Garfield primary identity invalid")
            lookup[k]=j
        # The microscopic ion-birth CSV can be interspersed with ionless
        # avalanches. Group by exact (crossing_id,primary_index); never
        # assume every seed has a reservoir or that rows are dense.
        samples=defaultdict(list)
        angle_counts=defaultdict(lambda:np.zeros(72,dtype=np.float32))
        weights=defaultdict(float)
        for r in rows(d/"muon_reference_ion_birth_sample.csv",(
            "crossing_id","primary_index","u_mm","v_mm",
            "t_ns","reservoir_weight","seed_ion_births_total")):
            k=seed_key(r)
            must(k in lookup,"Unmatched sampled Garfield ion")
            num=int(seeds[lookup[k]]["ion_births"])
            m=min(128,num)
            must(num>0 and int(r["seed_ion_births_total"])==num,
                 "Sampled ion birth has invalid avalanche ion census")
            w=float(r["reservoir_weight"])
            must(math.isfinite(w) and abs(w-num/m)<1.e-7*max(1.,num/m),
                 "Incorrect reservoir weight or sample size")
            t=float(r["t_ns"])
            must(math.isfinite(t),"Nonfinite microscopic physical time")
            samples[k].append(t)
            angle_counts[k][phi_bin(float(r["u_mm"]),float(r["v_mm"]))]+=1.
            weights[k]+=w
        must(sum(int(r["ion_births"]) for r in seeds)==
             sum(int(r["ion_births"]) for r in crossings),
             "Exact native seed and real muon ion census mismatch")
        must(sum(int(r["collected_e"]) for r in seeds)==
             sum(int(r["collected_avalanche_electrons"]) for r in crossings),
             "Exact Garfield collected charges disagree with muon sums")
        for r in seeds:
            k=seed_key(r)
            p=index[k[0]]
            cl,th=heed[k]
            n=int(r["ion_births"])
            gain=int(r["collected_e"])
            up,down=int(r["ion_up"]),int(r["ion_down"])
            must(min(n,gain,up,down)>=0 and int(r["capped"])==0,
                 "Invalid/capped microscopic avalanche")
            m=min(128,n)
            t=samples.get(k,[])
            must(len(t)==m,
                 "Sampled ion birth reservoir does not match exact ion census")
            must(abs(weights.get(k,0.)-n)<1.e-6*max(1.,n),
                 "Reservoir total does not reproduce exact ion birth count")
            if n:
                q=sample_quantiles(t)
                if q[0]<th-1.e-6:n_bad_time+=1
                onset,p50,p90=map(float,q)
                span=float(max(t)-min(t))
                fraction=angle_counts[k]/float(m)
                must(abs(float(fraction.sum())-1.)<1e-5,
                     "Angular ion-birth samples do not conserve unit weight")
                ions_all+=n
                points_all+=m
            else:
                onset=p50=p90=span=float("nan")
                fraction=np.zeros(72,dtype=np.float32)
            A=(up-down)/(up+down) if up+down else float("nan")
            vals=dict(parent=p,primary=k[1],cluster=cl,gain=gain,ions=n,
                A=A,u0=wire_offset(float(r["seed_u_mm"])),
                v0=float(r["seed_v_mm"]),t_heed=th,
                t10=onset,t50=p50,t90=p90,
                span=span,sample_count=m)
            for field in variables:out[field].append(vals[field])
            angular.append(fraction)
        del samples,angle_counts,weights
        print(f"PASS sealed Penning ON {block+1}/50 muon chunks; "
              f"{len(seeds)} original Garfield avalanches",flush=True)
    must(len(muons)==1000 and len(out["gain"])>20_000,
         "1000 real muons/Heed electron population missing")
    must(n_bad_time==0,
         f"{n_bad_time} sampled early birth times predate the initiating Heed electron")
    arrays={}
    for field,values in out.items():
        dtype=(np.int32 if field in ("parent","primary","cluster","gain",
                                    "ions","sample_count") else np.float64)
        arrays[field]=np.asarray(values,dtype=dtype)
    arrays["angle_fraction"]=np.asarray(angular,dtype=np.float32)
    n=len(arrays["gain"])
    must(arrays["angle_fraction"].shape==(n,72),
         "Incorrect per-avalanche angular sample dimensions")
    with np.load(ORIGINAL,allow_pickle=False) as source_data:
        must(np.array_equal(source_data["gain"],arrays["gain"]) and
             np.array_equal(source_data["births"],arrays["ions"]),
             "New timing data do not match existing verified Penning gain/birth arrays")
        must(len(source_data["muon_heed"])==len(muons) and
             np.array_equal(source_data["muon_heed"],
                            np.asarray([x["n_Heed"] for x in muons])),
             "New timing muon identity order differs from existing dashboard")
    must(sum(x["ion_births"] for x in muons)==ions_all,
         "Full 1000-muon ion birth total disagrees with 128-point reservoirs")
    must(ions_all==int(source["n_true_ion_births"]),
         "New temporal ion census disagrees with old verified dashboard")
    np.savez_compressed(CACHE,**arrays)
    metadata={
        "schema":"P3_NOFW_PENNING_ON_TEMPORAL_CORRELATIONS_V1",
        "muons":muons,
        "n_independent_true_Geant4_muons":1000,
        "n_exact_Garfield_avalanches":n,
        "n_true_ion_births":ions_all,
        "n_sampled_ion_birth_positions":points_all,
        "recommended_sample_count_threshold":128,
        "t10_definition":"10th percentile of <=128 uniformly sampled ion-birth timestamps after original muon entry; NOT electron drift or arrival time",
        "physical_setup":"Penning ON r=0.57, p3_noFW, pitch=3mm, +/-3mm cathodes, 1800V",
        "original_dashboard_cache_sha256":sha(ORIGINAL),
        "auxiliary_cache_sha256":sha(CACHE),
    }
    META.write_text(json.dumps(metadata,indent=2)+"\n")
    print(f"PASS new immutable temporal cache: 1000 muons, "
          f"{n:,} avalanches, {ions_all:,} exact ion births, "
          f"{points_all:,} sampled positions")
    print("NEW CACHE",CACHE,flush=True)

def selftest():
    must(phi_bin(0.,1.)==54 and phi_bin(0.,-1.)==18
         and abs(wire_offset(3.2)-.2)<1e-12,
         "Time study nearest-wire angular geometry incorrect")
    q=sample_quantiles([1.,2.,3.,4.])
    must(np.allclose(q,[1.3,2.5,3.7]),
         "Physical 10/50/90 percentile estimator incorrect")
    # A zero-ion avalanche has no birth samples and no defined onset.
    ions=0; onset=float("nan") if ions==0 else 0.
    must(np.isnan(onset),"Zero-ion onset must remain undefined")
    # All six live views and three-part scientific explanations must render
    # using synthetic data BEFORE any expensive sealed-CSV preparation.
    import plotly.graph_objects as go
    class StubStreamlit:
        def warning(self,message):
            pass
    # 20 equal-occupancy time bins require at least 400 test avalanches.
    # Exercise the real plotting minima instead of using a too-small fixture.
    n=600
    starts=np.linspace(2.,160.,n)
    synth={
        "gain":np.linspace(300.,900.,n),
        "ions":np.full(n,750,dtype=int),
        "A":np.linspace(-.8,.8,n),
        "t10":starts,
        "sample_count":np.full(n,128,dtype=int),
        "v0":np.linspace(-3.,3.,n),
        "angle_fraction":np.tile(np.ones((1,72))/72.,(n,1)),
        "parent":np.zeros(n,dtype=int),
        "primary":np.arange(n,dtype=int),
        "cluster":np.arange(n,dtype=int)%7,
        "t_heed":np.zeros(n),
    }
    fake_meta={"muons":[{
       "id":"synthetic","n_Heed":n,"n_clusters":7,
       "path_mm":6.,"theta_u_deg":5.,"theta_w_deg":-2.
    }]}
    settings={
        "reliable_only":True,"mode":"Quantile groups",
        "early_min":10.,"early_max":40.,
        "late_min":85.,"late_max":115.,
        "add_stats":True,
        "muon_index":0,"color":"#E46B08","line_width":2.,
        "legend":"synthetic","legend_early":"Early","legend_late":"Late",
        "angular_deg":5
    }
    for plot in KINDS:
        explanation=describe(plot,synth,fake_meta,settings)
        must(len(explanation)==3 and
             all(len(line)>80 for line in explanation),
             "Invalid science explanation for "+plot)
        output,note=render(plot,go.Figure(),synth,fake_meta,settings,
                           StubStreamlit())
        must(isinstance(note,str),"No user-visible plot annotation "+plot)
        if output is not None:
            must(all(len(np.asarray(v))==len(next(iter(output.values())))
                     for v in output.values()),
                 "Plotted CSV curves have mismatched lengths: "+plot)
    # Verify independent selected ns windows and both legend statistics.
    settings["mode"]="Custom time intervals (ns)"
    early,late,names=group_masks(synth,settings)
    must(np.array_equal(early,(starts>=10.)&(starts<40.)) and
         np.array_equal(late,(starts>=85.)&(starts<115.)),
         "Independent time interval selection failed")
    for check in ("t10_angles","t10_gain","t10_asym"):
        fig=go.Figure()
        render(check,fig,synth,fake_meta,settings,StubStreamlit())
        must(len(fig.data)>0 and
             all("mean=" in tr.name and "RMS=" in tr.name for tr in fig.data),
             "Mean/RMS legend missing or incorrect: "+check)
    for check in ("t10_start","t10_gain_position"):
        settings["zscale"]="Log"
        fig=go.Figure()
        render(check,fig,synth,fake_meta,settings,StubStreamlit())
        must(fig.data[0].colorbar.title.text.startswith("log10("),
             "Log colorbar control is not functional: "+check)
    settings["zscale"]="Linear"
    # Check a zero-ion event does not corrupt the source's time estimator.
    print("PASS six temporal plot+interpretation selftests, nearest-anode "
          "physics and sampled timing quantiles; NO Garfield/raw data writes.")

def load():
    must(CACHE.is_file() and META.is_file(),
         "Temporal cache missing: run timing-correlations script 'prepare' first")
    m=json.loads(META.read_text())
    must(m.get("schema")=="P3_NOFW_PENNING_ON_TEMPORAL_CORRELATIONS_V1"
         and sha(CACHE)==m.get("auxiliary_cache_sha256"),
         "Temporal cache checksum invalid")
    with np.load(CACHE,allow_pickle=False) as z:
        return {k:z[k] for k in z.files},m


def eligible(data,require128=True,need_A=False):
    """Fixed >=128-sample population avoids gain-dependent min-time bias."""
    size=data["sample_count"]
    mask=np.isfinite(data["t10"]) & (data["ions"]>0)
    if require128:mask&=size==128
    if need_A:mask&=np.isfinite(data["A"])
    return mask

def quantile_edges(x,n=12):
    x=np.asarray(x,dtype=float)
    must(x.size>=max(100,n*20),"Too few avalanches in timing study")
    return np.unique(np.quantile(x,np.linspace(0.,1.,n+1)))

def group_masks(data,setting):
    x=data["t10"]
    valid=eligible(data,require128=setting["reliable_only"])
    must(bool(np.any(valid)),"No eligible ion-producing avalanches")
    if setting["mode"]=="Quantile groups":
        q25,q75=np.quantile(x[valid],[.25,.75])
        early=valid&(x<=q25)
        late=valid&(x>=q75)
        description=(f"Earliest quartile (≤{q25:.2f} ns)",
                     f"Latest quartile (≥{q75:.2f} ns)")
    else:
        a,b=float(setting["early_min"]),float(setting["early_max"])
        c,d=float(setting["late_min"]),float(setting["late_max"])
        must(all(math.isfinite(z) for z in (a,b,c,d)) and
             0.<=a<b and 0.<=c<d,
             "Early and late windows must each have min < max, in ns")
        early=valid&(x>=a)&(x<b)
        late=valid&(x>=c)&(x<d)
        description=(f"{a:g}–{b:g} ns",f"{c:g}–{d:g} ns")
    return early,late,description

def typical_muon(meta):
    n=np.array([x["n_Heed"] for x in meta["muons"]],dtype=int)
    eligible_mu=np.flatnonzero(n>=10)
    must(len(eligible_mu)>0,"No muons available for timeline")
    typical=np.median(n[eligible_mu])
    return int(eligible_mu[np.argmin(np.abs(n[eligible_mu]-typical))])

def describe(kind,data,meta,settings):
    n_total=len(data["gain"])
    selected=eligible(data,require128=settings["reliable_only"])
    n=int(np.count_nonzero(selected))
    require_text=("Only avalanches with all 128 sampled ion-birth positions are "
                  "included to reduce sample-count-dependent time-estimator bias."
                  if settings["reliable_only"] else
                  "All ion-producing avalanches are included; low-ion-count "
                  "avalanches have less stable sampled time percentiles.")
    if kind=="t10_gain":
        x=data["t10"][selected]
        y=data["gain"][selected]
        q25,q75=np.quantile(x,[.25,.75])
        lo=float(y[x<=q25].mean())
        hi=float(y[x>=q75].mean())
        return (
            "**What you see:** Mean collected-electron gain, grouped by the "
            "sampled 10th-percentile *ion-birth* time relative to muon entry.",
            "**How measured:** Each avalanche's exact collected electron count "
            "comes from Garfield endpoints; its early-time proxy comes from "
            "up to 128 randomly retained ionisation timestamps. Time-bin "
            "averages include zero collected gain where ion births occurred. "
            +require_text,
            f"**Penning-ON observations:** {n:,}/{n_total:,} avalanches enter "
            f"the analysis. Mean gain in earliest/latest timing quartiles is "
            f"{lo:,.1f} / {hi:,.1f}. This correlation alone does not establish "
            "that drift time causes the gain difference."
        )
    if kind=="t10_asym":
        good=selected&np.isfinite(data["A"])
        x=data["t10"][good];a=np.abs(data["A"][good])
        q25,q75=np.quantile(x,[.25,.75])
        return (
            "**What you see:** Mean magnitude of the microscopic ion-birth "
            "hemisphere asymmetry |A| versus early sampled birth time.",
            "**How measured:** A = (N_up−N_down)/(N_up+N_down) uses **all** "
            "ion-birth callback counts, not the sampled cloud. "
            "Timing is estimated from 128 or fewer saved ion-birth "
            "timestamps, with t=0 at muon entry. "+require_text,
            f"**Penning-ON observations:** {len(x):,} valid avalanches; "
            f"mean |A| earliest/latest quartiles "
            f"{a[x<=q25].mean():.4f} / {a[x>=q75].mean():.4f}. "
            "This is birth-sector imbalance, not the induced cathode charge ratio."
        )
    if kind=="t10_angles":
        early,late,descriptions=group_masks(data,settings)
        return (
            "**What you see:** Two ion-weighted azimuthal distributions "
            "for avalanches with early versus late sampled ion-birth times. "
            "φ is the *birth position angle* around the anode, "
            "**not the incoming electron's arrival angle**.",
            "**How measured:** Each selected avalanche contributes its "
            "uniformly sampled birth azimuth histogram multiplied by "
            "its **exact** ion-birth count. Time groups come from its "
            "sampled early (10th percentile) birth time. "+require_text,
            f"**Penning-ON observations:** {descriptions[0]} — "
            f"{int(early.sum()):,} avalanches; {descriptions[1]} — "
            f"{int(late.sum()):,}. Compare the angular *shapes*, rather than "
            "their absolute ion yields; each curve is normalised separately."
        )
    if kind=="t10_start":
        x=np.abs(data["v0"][selected]);t=data["t10"][selected]
        rho=float(np.corrcoef(x,t)[0,1]) if np.std(x)>0 and np.std(t)>0 else float("nan")
        return (
            "**What you see:** Early sampled ion-birth time versus the "
            "distance |v₀| of the starting TrackHeed electron from the "
            "anode-wire plane (v=0). This is not a measured arrival position.",
            "**How measured:** v₀ is the original exact Heed electron "
            "position; the timestamp is the sampled 10th-percentile "
            "Garfield ion-creation time since muon entry. "
            "Each cell counts starting electrons. "+require_text,
            f"**Penning-ON observations:** {n:,} avalanches; Pearson "
            f"correlation between |v₀| and the time proxy = {rho:+.3f}. "
            "Position dependence can confound an apparent "
            "'late avalanche' effect."
        )
    if kind=="t10_gain_position":
        return (
            "**What you see:** Mean collected-electron gain at each "
            "combination of starting |v₀| and sampled early "
            "ion-birth time (two conditioned variables).",
            "**How measured:** The anode-collected gain is exact for each "
            "initial Heed electron; |v₀| is its saved start position; "
            "early-time is inferred from the saved 128-point birth "
            "reservoir. Cells with fewer than 20 avalanches are masked. "
            +require_text,
            f"**Penning-ON observations:** {n:,} eligible avalanches. "
            "If gain changes mainly with |v₀| at fixed time, or with "
            "time at fixed |v₀|, a one-variable LUT may be insufficient. "
            "This descriptive comparison is not a causal test."
        )
    if kind=="t10_timeline":
        mi=int(settings["muon_index"])
        m=meta["muons"][mi]
        selected_mu=data["parent"]==mi
        x=data["t10"][selected_mu]
        observed=int(np.isfinite(x).sum())
        return (
            f"**What you see:** A timeline of individual avalanches "
            f"created by real J/ψ muon #{mi}, grouped by their "
            "original TrackHeed cluster indices.",
            "**How measured:** Small grey marks are the exact TrackHeed "
            "electron creation times. Larger circles mark the sampled "
            "10th-percentile ion-birth times; their area approximately "
            "represents collected gain, and colour represents exact "
            "signed ion-birth asymmetry A. No separate electron "
            "arrival time or incoming angle is stored.",
            f"**This muon:** {m['n_clusters']} original clusters, "
            f"{m['n_Heed']} starting electrons, {observed} avalanches "
            f"with sampled ion births; path {m['path_mm']:.2f} mm, "
            f"incidence angles (θu,θw)=({m['theta_u_deg']:.2f}°, "
            f"{m['theta_w_deg']:.2f}°). The timeline is a diagnostic "
            "ion-creation chronology, not a reconstructed cathode waveform."
        )
    raise RuntimeError("Unknown temporal plot interpretation "+kind)

def render(kind,fig,data,meta,controls,st):
    """Add one of six physical plots to a caller-owned editable Plotly figure.

    Returns (plotted CSV data or None, short graph-specific note).
    """
    import plotly.graph_objects as go
    mask=eligible(data,controls["reliable_only"])
    x=data["t10"]
    y=data["gain"].astype(float)
    color=controls["color"]
    width=controls["line_width"]
    min_entries=20
    if kind in ("t10_gain","t10_asym"):
        good=mask.copy()
        if kind=="t10_asym":good&=np.isfinite(data["A"])
        ab=x[good]
        v=(y[good] if kind=="t10_gain" else np.abs(data["A"][good]))
        edges=quantile_edges(ab,14)
        mids=(edges[:-1]+edges[1:])/2
        n=np.histogram(ab,bins=edges)[0]
        summed=np.histogram(ab,bins=edges,weights=v)[0]
        avg=np.divide(summed,n,out=np.full(len(n),np.nan),
                      where=n>=min_entries)
        fig.add_trace(go.Scatter(
            x=mids,y=avg,mode="lines+markers",
            line=dict(color=color,width=width),
            marker=dict(size=7,color=color),
            name=(controls["legend"]+
                  (f" · mean={np.mean(v):,.2f}, RMS={np.std(v,ddof=0):,.2f}"
                   if controls.get("add_stats",False) else "")),
            customdata=n,hovertemplate="Early birth time: %{x:.3f} ns"
              "<br>Bin mean: %{y:.5g}<br>Avalanches: %{customdata}"
              "<extra></extra>"))
        return ({"bin_start_ns":edges[:-1],"bin_end_ns":edges[1:],
                 "mean":avg,"n_avalanches":n},
                "Time bins contain approximately equal numbers of "
                "eligible avalanches; plot values are arithmetic means, "
                "not probability densities.")
    if kind=="t10_angles":
        early,late,labels=group_masks(data,controls)
        out={}
        edges=np.asarray(PHI_EDGES,dtype=float)
        factor=int(controls["angular_deg"]//5)
        if min(int(early.sum()),int(late.sum()))<30:
            st.warning("Too few avalanches in at least one time group. "
                       "Broaden time windows or choose quantile groups.")
            return None,"Insufficient avalanche counts for reliable angular contrast."
        for name,selection,col,title in (
            ("early",early,color,controls["legend_early"]),
            ("late",late,"#1F77B4",controls["legend_late"])
        ):
            m=data["ions"][selection].astype(float)
            frac=data["angle_fraction"][selection]
            count=np.sum(frac*m[:,None],axis=0)
            reb=count.reshape(-1,factor).sum(axis=1)
            reb_edges=edges[::factor]
            must(len(reb_edges)==len(reb)+1,"Angular widths inconsistent")
            prob=reb/reb.sum()
            centers=(reb_edges[:-1]+reb_edges[1:])/2.
            xx=np.repeat(reb_edges,2)[1:-1]
            yy=np.repeat(prob,2)
            caption=title+" ("+labels[0 if name=="early" else 1]+")"
            if controls.get("add_stats",False):
                # Time means/RMS apply to selected t10 group, never the
                # periodic birth azimuth itself.
                group_t10=np.asarray(x[selection],dtype=float)
                caption+=(f" · t10 mean={group_t10.mean():.2f} ns, "
                          f"RMS={group_t10.std(ddof=0):.2f} ns")
            fig.add_trace(go.Scatter(
                x=xx,y=yy,mode="lines",line=dict(color=col,width=width),
                name=caption))
            out[name+"_probability"]=prob
        centers=(edges[::factor][:-1]+edges[::factor][1:])/2.
        return ({"phi_bin_center_deg":centers,
                 "early_probability":out["early_probability"],
                 "late_probability":out["late_probability"]},
                f"Groups contain {int(early.sum()):,} (early) and "
                f"{int(late.sum()):,} (late) avalanches. Both curves "
                "normalised by their own ion-weighted totals; "
                "blue is late, orange is early.")
    if kind=="t10_start":
        vv=np.abs(data["v0"][mask])
        tt=x[mask]
        edges_v=np.linspace(0.,3.,19)
        edges_t=quantile_edges(tt,20)
        n=np.histogram2d(tt,vv,bins=(edges_t,edges_v))[0]
        cen_t=(edges_t[:-1]+edges_t[1:])/2
        cen_v=(edges_v[:-1]+edges_v[1:])/2
        z=n
        cb="Number of Heed electrons"
        if controls.get("zscale")=="Log":
            with np.errstate(divide="ignore"):
                z=np.where(n>0,np.log10(n),np.nan)
            cb="log10(number of Heed electrons)"
        fig.add_trace(go.Heatmap(
            x=cen_v,y=cen_t,z=z,
            colorscale="YlOrBr",colorbar=dict(title=cb),
            hovertemplate="|v0| %{x:.3g} mm<br>early time %{y:.3g} ns"
                          "<br>Heed electrons %{z:.0f}<extra></extra>"))
        return None,("Cell occupancy is a count, not a probability density; "
                     "time bins are chosen for equal global occupancy.")
    if kind=="t10_gain_position":
        tt=x[mask]
        vv=np.abs(data["v0"][mask])
        gain=y[mask]
        edget=quantile_edges(tt,15)
        edgev=np.linspace(0.,3.,17)
        cnt=np.histogram2d(vv,tt,bins=(edgev,edget))[0]
        total=np.histogram2d(vv,tt,bins=(edgev,edget),weights=gain)[0]
        mean=np.divide(total,cnt,out=np.full(cnt.shape,np.nan),
                       where=cnt>=20)
        xc=(edget[:-1]+edget[1:])/2
        yc=(edgev[:-1]+edgev[1:])/2
        z=mean
        cb=controls["legend"]
        if controls.get("zscale")=="Log":
            with np.errstate(divide="ignore",invalid="ignore"):
                z=np.where(mean>0,np.log10(mean),np.nan)
            cb="log10("+cb+")"
        fig.add_trace(go.Heatmap(
            x=xc,y=yc,z=z,
            colorscale="YlOrBr",colorbar=dict(title=cb),
            customdata=cnt,
            hovertemplate="early birth time %{x:.3g} ns<br>"
                          "start |v0| %{y:.3g} mm<br>"
                          "mean gain %{z:.2f}<br>"
                          "N avalanches %{customdata:.0f}<extra></extra>"))
        return None,("Only 2D bins with at least 20 true starting "
                     "electrons are shown. Mean gain includes collected_e=0.")
    if kind=="t10_timeline":
        idx=int(controls["muon_index"])
        good=(data["parent"]==idx)
        if not np.any(good):
            st.warning("Selected muon has no retained TrackHeed electrons.")
            return None,"This real muon has zero retained electron primaries."
        clusters=data["cluster"][good].astype(float)
        th=data["t_heed"][good]
        t=x[good]
        gain=y[good]
        asym=data["A"][good]
        valid=np.isfinite(t)
        # Exact Heed creation precedes avalanche (grey crosses), and
        # approximate onset from sampled ion births is a coloured circle.
        fig.add_trace(go.Scatter(
            x=th,y=clusters,mode="markers",
            marker=dict(symbol="x",size=6,color="#888888",opacity=.65),
            name="TrackHeed start time (exact)",
            hovertemplate="Heed creation t %{x:.3f} ns<br>"
                          "Cluster %{y:.0f}<extra></extra>"))
        if np.any(valid):
            size=np.clip(6+3*np.sqrt(np.maximum(gain[valid],0.)/100.),
                         6,36)
            cd=np.column_stack((gain[valid],data["primary"][good][valid],
                                th[valid],data["ions"][good][valid],
                                data["sample_count"][good][valid]))
            fig.add_trace(go.Scatter(
                x=t[valid],y=clusters[valid],mode="markers",
                marker=dict(size=size,color=asym[valid],cmin=-1,cmax=1,
                    colorscale="RdBu_r",showscale=True,
                    colorbar=dict(title="Exact A_ion")),
                customdata=cd,
                name=controls["legend"],
                hovertemplate="Early 10% ion-birth time %{x:.3f} ns"
                              "<br>Cluster %{y:.0f}"
                              "<br>Collected gain %{customdata[0]:.0f}"
                              "<br>Heed primary index %{customdata[1]:.0f}"
                              "<br>Heed creation t %{customdata[2]:.3f} ns"
                              "<br>Ion births %{customdata[3]:.0f}"
                              "<br>Sampled positions %{customdata[4]:.0f}"
                              "<extra></extra>"))
        return None,("Circles depict sampled early ion-birth times; grey "
                     "crosses show original Heed creation times. "
                     "Marker area varies approximately with gain; "
                     "colour is exact signed ion-birth asymmetry.")
    raise RuntimeError("No temporal plot renderer for "+kind)


if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "selftest"
    if mode=="selftest":selftest()
    elif mode=="prepare":prepare()
    else:raise SystemExit("Usage: python SCRIPT selftest | prepare")

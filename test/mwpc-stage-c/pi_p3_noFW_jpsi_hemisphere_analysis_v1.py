#!/usr/bin/env python3
"""Read-only hemisphere-only analysis of the 1000 + (up to) 10000 new J/psi muons.

Physics: p3_noFW full microscopic Garfield, genuine TrackHeed primaries.
ONE and ONLY asymmetry: A_hemi = (ion_up - ion_down)/(ion_up + ion_down)
based on ALL Garfield ion births, not the 128-point-per-avalanche reservoir.
No Gaussian-cloud claim, no cathode signal, no new Garfield. All unsealed
chunks are excluded, even while the 10k campaign is actively running.

Input checks: immutable input CSV SHA, distinct Geant4 parent IDs (old/new),
chunk completion seal and the SHA of every CSV actually used, per-parent
conservation of original full ion counts and per-TrackHeed-electron counts.
All output goes into an immutable separate snapshot, never into source dirs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pi_prepare_muon_avalanche_b1 import FIELDS
from pi_run_native_muon_b1 import CROSSING_COLUMNS

REPO = Path(__file__).resolve().parents[2]
BASE = REPO/"test_runs/mwpc_stage_c"
OLD = BASE/"pi_p3_noFW_jpsi1000_native_v1"
NEW = BASE/"pi_p3_noFW_jpsi_extra10000_native_v1"
DEST = BASE/"pi_p3_noFW_jpsi_hemisphere_v1"
HIT_ID = CROSSING_COLUMNS[0]
NEEDED = ("native_muon_crossings_input.csv","muon_reference_crossing.csv",
          "muon_reference_seed.csv","muon_reference_heed_primaries.csv")
PITCH = 3.0


def insist(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path, required=None):
    with Path(path).open(newline="") as f:
        reader=csv.DictReader(f)
        insist(reader.fieldnames is not None,"Empty source CSV "+str(path))
        if required is not None:
            insist(set(required).issubset(reader.fieldnames),
                   "Wrong source CSV columns "+str(path))
        yield from reader


def read_json(path):
    return json.loads(Path(path).read_text())


def pid(row):
    return str(row["source"]),str(row["chunk"]),str(row["event"])


def crossing_id(row):
    return "|".join(str(row[x]) for x in
                    ("source","chunk","event","hit_index",
                     "track_id","station","det_id"))


def round_cpp(z):
    return math.floor(z+.5) if z>=0 else math.ceil(z-.5)


def entrance_phase_mm(u,st):
    """Exact phase used by native C++ NAnodes/EntrancePhase for 3 mm."""
    n=math.floor((307.2 if st==0 else 512.0)/PITCH+.5)
    phase=PITCH*.5 if n%2==0 else 0.0
    nearest=round_cpp((u-phase)/PITCH)*PITCH+phase
    result=u-nearest
    insist(-1.50000001<=result<=1.50000001,
           "Muon entry outside centered 3mm anode period")
    return result


def expected_A(up,down):
    return (up-down)/(up+down) if up+down else math.nan


def close_A(row,truth):
    given=row["A_ion_hemi"]
    if math.isnan(truth):
        insist(given=="","C++ undefined A_hemi is not blank")
    else:
        insist(given and abs(float(given)-truth)<2e-9,
               "Stored whole-muon A_hemi differs from exact ion counters")


def campaign(root, expected_count, required):
    input_csv=root/"input/native_muon_crossings_input.csv"
    meta_path=root/"input/input_manifest.json"
    if not input_csv.is_file() or not meta_path.is_file():
        if required:
            raise FileNotFoundError("Required frozen campaign missing: "+str(root))
        return [],None
    meta=read_json(meta_path)
    insist(meta["n_parents"]==expected_count and meta["geometry"]=="p3_noFW"
           and meta["group_size"]==20 and
           digest(input_csv)==meta["selected_input_sha256"],
           "Invalid SHA-frozen Geant4 input "+str(root))
    frozen=list(read_csv(input_csv,FIELDS))
    insist(len(frozen)==expected_count and
           len({pid(r) for r in frozen})==expected_count,
           "Frozen input parent identities are not unique")
    out=[]
    records=[]
    seal_sums=[]
    identity=set()
    for k in range(expected_count//20):
        d=root/"native_chunks"/f"chunk_{k:04d}"
        if not (d/"complete.json").is_file():
            # Never open unsealed chunk, which may be actively written.
            later=list((root/"native_chunks").glob("chunk_*/complete.json"))
            insist(all(int(p.parent.name.split("_")[1])<k for p in later),
                   "Unexpected sealed gap in campaign "+str(root))
            break
        seal=read_json(d/"complete.json")
        insist(seal["chunk"]==k and seal["n_parents"]==20 and
               seal["schema"]=="NA60_JPSI_20PARENT_P3NOFW_SHA_V1",
               "Wrong source campaign checkpoint")
        fhash=seal["file_sha256"]
        insist(all(name in fhash for name in NEEDED),"Missing SHA protection")
        for name in NEEDED:
            insist(digest(d/name)==fhash[name],
                   "SHA mismatch for native microscopic data: "+str(d/name))
        inputs=list(read_csv(d/"native_muon_crossings_input.csv",FIELDS))
        insist([crossing_id(r) for r in inputs]==
               [crossing_id(r) for r in frozen[k*20:(k+1)*20]],
               "Chunk Geant4 input does not match original SHA-frozen selection")
        source={crossing_id(r):r for r in inputs}
        stats={}
        for r in read_csv(d/"muon_reference_crossing.csv",
                          (HIT_ID,"configuration","ion_births","ion_up",
                           "ion_down","ion_boundary","A_ion_hemi",
                           "heed_primaries","failed_seeds","capped_seeds")):
            key=r[HIT_ID]
            insist(key in source and key not in stats and
                   r["configuration"]=="p3_noFW",
                   "Wrong geometry/duplicate native event identity")
            up,down,boundary=(int(r[z]) for z in
                              ("ion_up","ion_down","ion_boundary"))
            births=int(r["ion_births"])
            insist(up>=0 and down>=0 and boundary>=0 and
                   up+down+boundary==births,
                   "Nonconserved exact full Garfield ion hemisphere counts")
            a=expected_A(up,down)
            close_A(r,a)
            hit=source[key]
            event={
                "source_batch":"original1000" if expected_count==1000 else "extra10000",
                "crossing_id":key,
                "parent_key":pid(hit),
                "station":int(hit["station"]),
                "theta_u_deg":float(hit["theta_u_in_deg"]),
                "theta_w_deg":float(hit["theta_w_in_deg"]),
                "entry_u_phase_mm":entrance_phase_mm(
                    float(hit["u_in_mm"]),int(hit["station"])),
                "v_entry_mm":float(hit["v_in_mm"]),
                "n_Heed_electrons":int(r["heed_primaries"]),
                "n_ion_births":births,
                "n_up":up,"n_down":down,"n_boundary":boundary,
                "A_hemi":a,
                "failed_avalanche_calls":int(r["failed_seeds"]),
                "capped_avalanche_calls":int(r["capped_seeds"]),
            }
            stats[key]=event
        insist(len(stats)==20,"Missing genuine whole-muon response rows")

        seeds=Counter()
        counts=defaultdict(lambda: [0,0,0,0,0])
        # These are per-avalanche exact, UNSUBSAMPLED ion counters.
        # sign-flip null: each avalanche's complete up/down imbalance
        # has a randomized independent sign, magnitude and gain unchanged.
        sign_null=Counter()
        sum_t_sq=Counter()
        for r in read_csv(d/"muon_reference_seed.csv",
                          ("crossing_id","configuration","primary_index",
                           "ion_births","ion_up","ion_down","successful","capped")):
            key=r["crossing_id"]
            insist(key in stats and r["configuration"]=="p3_noFW",
                   "Unknown or wrong-geometry Garfield primary")
            ui,di=int(r["ion_up"]),int(r["ion_down"])
            bi=int(r["ion_births"])
            insist(ui>=0 and di>=0 and bi>=ui+di,
                   "Invalid complete single-avalanche ion census")
            c=counts[key]
            c[0]+=ui;c[1]+=di;c[2]+=bi
            c[3]+=int(r["successful"])==0
            c[4]+=int(r["capped"])
            seeds[key]+=1
            sign_null[key]+=(ui-di)**2
            sum_t_sq[key]+=(ui+di)**2
        primary_v=defaultdict(lambda:[0,0,0,0.0])
        for r in read_csv(d/"muon_reference_heed_primaries.csv",
                          ("crossing_id","primary_index","v_abs_cm")):
            key=r["crossing_id"]
            insist(key in stats,"Unknown original TrackHeed electron source")
            v=float(r["v_abs_cm"])
            insist(math.isfinite(v),"Nonfinite primary v")
            pv=primary_v[key]
            if v>0:pv[0]+=1
            elif v<0:pv[1]+=1
            else:pv[2]+=1
            pv[3]+=10*v
        for key,z in stats.items():
            up,down,births,failed,capped=counts[key]
            insist(seeds[key]==z["n_Heed_electrons"] and
                   up==z["n_up"] and down==z["n_down"] and
                   births==z["n_ion_births"] and
                   failed==z["failed_avalanche_calls"] and
                   capped==z["capped_avalanche_calls"],
                   "Full Garfield per-electron vs per-muon mismatch")
            pv=primary_v[key]
            insist(sum(pv[:3])==seeds[key],
                   "TrackHeed primary electron census mismatch")
            T=up+down
            z["primary_upper_fraction"]=pv[0]/seeds[key] if seeds[key] else math.nan
            z["primary_v_imbalance"]=(pv[0]-pv[1])/seeds[key] if seeds[key] else math.nan
            z["primary_v_mean_mm"]=pv[3]/seeds[key] if seeds[key] else math.nan
            z["null_sd_independent_avalanche_signs"]=(
                math.sqrt(sign_null[key])/T if T else math.nan)
            z["n_gain_effective"]=T*T/sum_t_sq[key] if sum_t_sq[key] else math.nan
            z["n_active_avalanche"]=sum(1 for _ in ()) # set below if relevant
            insist(key not in identity,"Duplicate physical crossing in campaign")
            identity.add(key)
            out.append(z)
        seal_sums.append({"chunk":k,"seal_sha256":digest(d/"complete.json")})
    if out:
        print(f"PASS {root.name}: {len(out)}/{expected_count} SHA-audited muons",flush=True)
    elif required:
        raise RuntimeError("No SHA-sealed microscopic muon chunks: "+str(root))
    receipt={"campaign":root.name,
             "selected_input_sha256":digest(input_csv),
             "physics_checkpoint_seals":seal_sums,
             "n_used_complete_muons":len(out),
             "n_expected_muons":expected_count}
    return out,receipt


def pearson(a,b):
    x=np.asarray(a,float);y=np.asarray(b,float)
    mask=np.isfinite(x)&np.isfinite(y)
    if mask.sum()<5 or np.std(x[mask])==0 or np.std(y[mask])==0:
        return math.nan
    return float(np.corrcoef(x[mask],y[mask])[0,1])


def per_bin(x,y,edges):
    x=np.asarray(x);y=np.asarray(y)
    result=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        v=y[(x>=lo)&(x<hi)&np.isfinite(y)]
        n=len(v)
        result.append((.5*(lo+hi),float(np.mean(v)) if n else math.nan,
                       float(np.std(v,ddof=1)/np.sqrt(n)) if n>1 else math.nan,n))
    return np.asarray(result,dtype=float)


def binned(ax,x,y,edges,label,min_n=12):
    z=per_bin(x,y,edges)
    sel=(z[:,3]>=min_n)
    ax.errorbar(z[sel,0],z[sel,1],yerr=z[sel,2],fmt="o-",
                capsize=3,label=label,ms=4)
    return z


def make(r, receipts, reps):
    insist(len(r)>0,"No genuine full Garfield muons available")
    pkeys=[x["parent_key"] for x in r]
    insist(len(set(pkeys))==len(r),
           "Original 1000 and extra 10000 share Geant4 Jpsi parents")
    array=lambda k:np.array([x[k] for x in r],dtype=float)
    A=array("A_hemi")
    valid=np.isfinite(A)
    insist(np.any(valid),"No defined A_hemi")
    A0=A[valid]
    angle_u=array("theta_u_deg")
    angle_w=array("theta_w_deg")
    phase=array("entry_u_phase_mm")
    number=array("n_Heed_electrons")
    prim_imb=array("primary_v_imbalance")
    sd_null=array("null_sd_independent_avalanche_signs")
    eff_n=array("n_gain_effective")
    ion=array("n_ion_births")
    # Descriptive uncertainty is clustered at the independent Geant4
    # parent event level, NEVER at the ion or electron level.
    rng=np.random.default_rng(20261105)
    boot=np.empty(reps)
    boot_ratio=np.empty(reps)
    a=A0
    b=sd_null[valid]
    insist(np.all(np.isfinite(b)),"Null moment missing for defined A")
    for i in range(reps):
        picks=rng.integers(0,len(a),size=len(a))
        boot[i]=np.mean(a[picks])
        boot_ratio[i]=np.mean(a[picks]**2)/max(1e-20,np.mean(b[picks]**2))
    ci=np.quantile(boot,[.025,.975])
    ratio_ci=np.quantile(boot_ratio,[.025,.975])
    r2=float(np.mean(a*a)/np.mean(b*b))
    total=len(r)
    n_new=sum(x["source_batch"]=="extra10000" for x in r)
    name=f"snapshot_{total:05d}_old1000_new{n_new:05d}"
    folder=DEST/name
    insist(not folder.exists(),
           "Immutable hemisphere-only snapshot already exists: "+str(folder))
    folder.mkdir(parents=True,exist_ok=False)

    # 1: Whole-muon A_hemi, exact from all ion-birth counters.
    fig,ax=plt.subplots(figsize=(9,5.4))
    ax.hist(A0,bins=np.linspace(-1,1,51),density=False,
            weights=np.ones(len(a))/len(a),histtype="step",lw=1.7)
    ax.axvline(0,ls="--",color="grey",lw=1)
    ax.axvline(float(np.mean(a)),ls="-.",lw=1.3,label=f"Mean {np.mean(a):+.4f}")
    ax.set_xlabel("Exact whole-muon ion-birth hemisphere asymmetry")
    ax.set_ylabel("Probability per bin / independent J/psi parent")
    ax.set_xlim(-1,1)
    ax.set_title(f"p3_noFW: A_hemi for {len(a)} real muons; undefined={len(r)-len(a)}")
    ax.grid(alpha=.15);ax.legend()
    fig.tight_layout();fig.savefig(folder/"01_exact_hemisphere_asymmetry.png",dpi=190);plt.close(fig)

    # 2: A sign and magnitude versus original TRUE 3D Geant4 direction.
    fig,axes=plt.subplots(2,2,figsize=(11.3,8.0),sharey="row")
    for j,(th,title) in enumerate(((angle_u,"theta_u"),(angle_w,"theta_w"))):
        ix=valid & np.isfinite(th)
        axes[0,j].scatter(th[ix],A[ix],s=6,alpha=.11,rasterized=True)
        edges=np.linspace(-27,27,13)
        binned(axes[0,j],th[ix],A[ix],edges,"Binned mean A",min_n=12)
        axes[0,j].axhline(0,ls="--",lw=.8,color="gray")
        axes[0,j].set_title("Signed "+title)
        axes[0,j].set_xlabel(title+" [deg]")
        absedge=np.array([0,2,4,6,8,10,13,16,20,30],float)
        binned(axes[1,j],np.abs(th[ix]),np.abs(A[ix]),absedge,
               "Mean |A|",min_n=12)
        axes[1,j].set_xlabel("|"+title+"| [deg]")
        axes[0,j].legend(fontsize=8);axes[1,j].legend(fontsize=8)
        for i in range(2):axes[i,j].grid(alpha=.13)
    axes[0,0].set_ylabel("Signed A_hemi")
    axes[1,0].set_ylabel("Mean |A_hemi|")
    axes[0,0].set_ylim(-1,1)
    fig.suptitle("True Geant4 incidence: hemisphere asymmetry only")
    fig.tight_layout();fig.savefig(folder/"02_hemisphere_vs_incidence.png",dpi=190);plt.close(fig)

    # 3: Position relative nearest physical anode phase, not a random cell.
    fig,axes=plt.subplots(1,2,figsize=(11.4,4.7))
    edges=np.linspace(-1.5,1.5,16)
    sel=valid
    binned(axes[0],phase[sel],A[sel],edges,"Mean signed A",min_n=10)
    axes[0].axhline(0,ls="--",color="gray",lw=.8)
    binned(axes[1],phase[sel],np.abs(A[sel]),edges,
           "Mean absolute A",min_n=10)
    for ax in axes:
        ax.set_xlabel("Muon entry u relative nearest real 3mm anode [mm]")
        ax.grid(alpha=.14);ax.legend(fontsize=8)
    axes[0].set_ylabel("A_hemi");axes[1].set_ylabel("|A_hemi|")
    fig.suptitle("Actual Geant4 entrance phase; C++ station parity preserved")
    fig.tight_layout();fig.savefig(folder/"03_hemisphere_vs_anode_phase.png",dpi=190);plt.close(fig)

    # 4: Why averaging many microscopic avalanches does NOT necessarily
    # suppress event-level A: seed multiplicity and primary cloud side.
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.9))
    mask=valid&(number>0)
    n_edges=np.array([0,15,25,35,50,70,100,150,250,600,1000],float)
    binned(axes[0],number[mask],np.abs(A[mask]),n_edges,
           "Mean |A| vs number of Heed electrons",min_n=10)
    axes[0].set_xlabel("True TrackHeed conduction electrons / muon")
    axes[0].set_ylabel("Mean |A_hemi|")
    edges_v=np.linspace(-1,1,17)
    m=valid&np.isfinite(prim_imb)
    binned(axes[1],prim_imb[m],A[m],edges_v,
           "Mean A vs primary-electron hemisphere imbalance",min_n=10)
    axes[1].axhline(0,ls="--",lw=.8,color="grey")
    axes[1].set_xlabel("(Heed N_v>0 - N_v<0) / N_Heed")
    axes[1].set_ylabel("Mean signed A_hemi")
    for ax in axes:ax.grid(alpha=.14);ax.legend(fontsize=8)
    fig.suptitle("Primary ionisation: event independence is at muon-parent level")
    fig.tight_layout();fig.savefig(folder/"04_hemisphere_vs_primaries.png",dpi=190);plt.close(fig)

    # 5: Coherence diagnostic: independent sign flips of the COMPLETE
    # per-avalanche ion difference preserve individual gain and shape, but
    # erase inter-avalanche side alignment. This is a NULL DIAGNOSTIC,
    # not a microscopic Garfield physics simulation.
    fig,axes=plt.subplots(1,2,figsize=(11.3,4.8))
    m=valid&np.isfinite(sd_null)
    axes[0].scatter(sd_null[m],np.abs(A[m]),s=7,alpha=.24,rasterized=True)
    axes[0].plot([0,1],[0,1],ls="--",lw=1,color="gray")
    axes[0].set_xlabel("Expected A rms if avalanche signs independent")
    axes[0].set_ylabel("Observed |A_hemi|")
    axes[0].set_title("Per-muon actual asymmetry vs sign-flip null scale")
    xs=np.linspace(0,1,12)
    binned(axes[1],sd_null[m],A[m]**2,xs,
           "Observed mean A²",min_n=12)
    binned(axes[1],sd_null[m],sd_null[m]**2,xs,
           "Independent-sign-flip E[A²]",min_n=12)
    axes[1].set_xlabel("Independent sign-flip null rms, by muon")
    axes[1].set_ylabel("Second moment of asymmetry")
    axes[1].set_title(f"Ratio observed/null={r2:.3f} (parent bootstrap)")
    for ax in axes:ax.grid(alpha=.13);ax.legend(fontsize=8)
    fig.suptitle("Conditional independent avalanche-side-flip diagnostic ONLY")
    fig.tight_layout();fig.savefig(folder/"05_hemisphere_avalanche_alignment_test.png",dpi=190);plt.close(fig)

    # 6: A simple true-muon 2D angle/entry-phase map, with low-stat bins
    # hidden rather than allowing small cells to mislead the conclusion.
    xx=np.linspace(-1.5,1.5,13)
    yy=np.array([-30,-15,-10,-6,-3,0,3,6,10,15,30],float)
    mask=valid
    sums,_,_=np.histogram2d(angle_u[mask],phase[mask],bins=(yy,xx),
                             weights=np.abs(A[mask]))
    num,_,_=np.histogram2d(angle_u[mask],phase[mask],bins=(yy,xx))
    z=np.divide(sums,num,out=np.full_like(sums,np.nan),where=num>=10)
    fig,ax=plt.subplots(figsize=(10,5.4))
    im=ax.pcolormesh(xx,yy,np.ma.masked_invalid(z),shading="flat",
                      cmap="viridis",vmin=0,vmax=1)
    fig.colorbar(im,ax=ax,label="Mean |A_hemi| (bins with >=10 parents)")
    ax.set_xlabel("Physical entry phase relative nearest anode [mm]")
    ax.set_ylabel("Signed true theta_u [deg]")
    ax.set_title("Muon geometry vs hemisphere asymmetry, low-count bins masked")
    fig.tight_layout();fig.savefig(folder/"06_hemisphere_angle_phase_map.png",dpi=190);plt.close(fig)

    # Exact, unsampled per-event physics and source station/angle coordinates.
    out_csv=folder/"exact_hemisphere_per_muon.csv"
    cols=[
        "source_batch","crossing_id","station",
        "theta_u_deg","theta_w_deg","entry_u_phase_mm","v_entry_mm",
        "n_Heed_electrons","n_ion_births","n_up","n_down","n_boundary",
        "A_hemi","primary_upper_fraction","primary_v_imbalance",
        "primary_v_mean_mm","null_sd_independent_avalanche_signs",
        "n_gain_effective","failed_avalanche_calls","capped_avalanche_calls"
    ]
    with out_csv.open("x",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader()
        for z in r:
            w.writerow({k:(z[k] if not isinstance(z[k],float)
                          or math.isfinite(z[k]) else "") for k in cols})
    report={
        "schema":"NA60_P3NOFW_JPSI_HEMISPHERE_ONLY_V1",
        "n_original1000":total-n_new,
        "n_extra10000":n_new,
        "n_independent_geant_parent_muons":total,
        "n_defined_hemisphere_asymmetries":int(valid.sum()),
        "n_undefined_hemisphere_asymmetries":int((~valid).sum()),
        "n_garfield_ion_births":int(np.sum(ion)),
        "n_garfield_failed_calls":int(sum(x["failed_avalanche_calls"] for x in r)),
        "n_garfield_capped_calls":int(sum(x["capped_avalanche_calls"] for x in r)),
        "mean_A_hemi":float(np.mean(a)),
        "std_A_hemi":float(np.std(a,ddof=1)),
        "mean_A_hemi_CI95_parent_bootstrap":ci.tolist(),
        "mean_abs_A_hemi":float(np.mean(np.abs(a))),
        "fraction_abs_A_greater_0p5":float(np.mean(np.abs(a)>.5)),
        "fraction_abs_A_greater_0p8":float(np.mean(np.abs(a)>.8)),
        "correlation_A_vs_signed_theta_u":pearson(A,angle_u),
        "correlation_absA_vs_abs_theta_u":pearson(np.abs(A),np.abs(angle_u)),
        "correlation_absA_vs_number_Heed":pearson(np.abs(A),number),
        "correlation_A_vs_Heed_primary_side_imbalance":pearson(A,prim_imb),
        "mean_observed_A_squared":float(np.mean(a*a)),
        "mean_independent_sign_flip_A_squared":float(np.mean(b*b)),
        "observed_to_null_A_squared_ratio":r2,
        "ratio_CI95_parent_bootstrap":ratio_ci.tolist(),
        "interpretation_warning":"Sign-flip test is a conditional null diagnostic, "
                                  "not a test of cathode charge and not proof of "
                                  "independent electron histories. Parent is sampling unit.",
        "physical_definition":"A_hemi=(N_ionbirth_v>0-N_ionbirth_v<0)/"
                              "(N_ionbirth_v>0+N_ionbirth_v<0), nearest anode; "
                              "exact unsampled Garfield ion creation census.",
        "neither_ion_drift_nor_cathode_signal":True,
        "input_checkpoint_receipts":receipts,
        "source_script_sha256":digest(Path(__file__))
    }
    report["files_sha256"]={p.name:digest(p) for p in sorted(folder.glob("*.png"))}
    report["files_sha256"][out_csv.name]=digest(out_csv)
    (folder/"summary.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("PASS p3_noFW J/psi REAL Garfield HEMISPHERE ONLY",flush=True)
    print(f"Independent physical parents {total}; old={total-n_new}; new={n_new}")
    print(f"Exact A_hemi defined={len(a)}, undefined={total-len(a)}")
    print(f"Mean {np.mean(a):+.6f}, 95% parent bootstrap CI [{ci[0]:+.5f},{ci[1]:+.5f}]")
    print(f"SD={np.std(a,ddof=1):.5f}; mean |A|={np.mean(np.abs(a)):.5f}")
    print(f"corr(|A|,abs theta_u)={pearson(np.abs(A),np.abs(angle_u)):+.5f}")
    print(f"corr(|A|,N_Heed)={pearson(np.abs(A),number):+.5f}")
    print(f"corr(A,primary v hemisphere imbalance)={pearson(A,prim_imb):+.5f}")
    print(f"Observed/null independent sign-flip A^2 ratio {r2:.4f}, "
          f"bootstrap CI [{ratio_ci[0]:.4f},{ratio_ci[1]:.4f}]")
    print("OUTPUT IMMUTABLE",folder,flush=True)


def selftest():
    insist(round_cpp(-.5)==-1 and round_cpp(.5)==1,
           "C++ period phase rounding")
    insist(abs(entrance_phase_mm(0.,0)+1.5)<1e-9,
           "Incorrect even-anode physical station parity")
    insist(abs(entrance_phase_mm(0.,1))<1e-9,
           "Incorrect odd-anode physical station parity")
    insist(abs(expected_A(30,10)-.5)<1e-12,
           "Exact ion-sector hemisphere asymmetry")
    insist(math.isnan(expected_A(0,0)),
           "Undefined A with zero ions must not be fabricated")
    d=per_bin([0,1,2,3],[1,1,3,3],[0,2,4])
    insist(len(d)==2 and d[0,3]==2 and d[1,3]==2,
           "Event-level binning sanity")
    print("PASS independent physics units, exact A_hemi, C++ anode phase parity")
    print("No Garfield jobs, source files or plotting outputs modified")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("selftest","analyze"))
    p.add_argument("--old-only",action="store_true",
                   help="Analyze only original complete 1000 even if extra10k exists")
    p.add_argument("--require-complete",action="store_true",
                   help="Require all 10000 ADDITIONAL parents SHA-sealed")
    p.add_argument("--bootstrap",type=int,default=800)
    a=p.parse_args()
    insist(a.bootstrap>=100,"Use >=100 physical parent bootstrap draws")
    if a.mode=="selftest":
        selftest()
        return
    old,proof_old=campaign(OLD,1000,required=True)
    insist(len(old)==1000,"Original 1000 reference must be complete")
    new,proof_new=([],None)
    if not a.old_only:
        new,proof_new=campaign(NEW,10000,required=a.require_complete)
    if a.require_complete:
        insist(not a.old_only and len(new)==10000,
               "Extra 10k incomplete: use regular analyze for partial snapshot")
    make(old+new,[p for p in (proof_old,proof_new) if p],a.bootstrap)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Read-only, paired-start comparison: frozen microscopic LUT500 vs fast MC.

fit/validation NEVER pooled in the primary metrics. Cell-clustered 95% CI
for paired collection differences; phi PDF is conditional on each model's
OWN surviving electrons (different survivor populations with MC attachment
disabled). No automatic fast-model acceptance or calibration here.
"""
from __future__ import annotations

import argparse,csv,hashlib,json,math,statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

GEOMS=("p3p0_anode_only","p3p0_field75p0um")
REQUIRED_MICRO={
  "geometry","split","iu","iv","sample","random_seed","u_start_mm",
  "v_start_mm","reference_collected","reference_field_wire_collected",
  "reference_status","reference_wire_index","reference_phi_hit_deg",
  "reference_drift_time_ns"}
REQUIRED_MC={
  "geometry","split","iu","iv","sample","random_seed","u_start_mm",
  "v_start_mm","mc_step_um","mc_attachment","drift_ok","status",
  "anode_collected","field_wire_collected","other_end","phi_hit_deg",
  "drift_time_ns","drift_wall_us"}


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()


def ident(row):
    return (row["geometry"],int(row["iu"]),int(row["iv"]),int(row["sample"]))


def read_micro(path):
    reference={}
    with path.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        if not REQUIRED_MICRO.issubset(reader.fieldnames or []):
            raise ValueError("Wrong selected frozen reference CSV header")
        for row in reader:
            key=ident(row)
            if key in reference:raise ValueError("Duplicate microscopic start")
            if key[0] not in GEOMS or row["split"] not in ("fit","validation"):
                raise ValueError("Unexpected reference geometry or split")
            if int(row["reference_collected"]) not in (0,1):
                raise ValueError("Invalid reference collection flag")
            reference[key]=row
    if not reference:raise ValueError("No microscopic reference starts")
    # A valid selection is produced by the SHA-pinned selector; check
    # manifest if it exists to reject truncation/accidental modifications.
    manifest=path.with_suffix(".manifest.json")
    if not manifest.exists():
        raise ValueError("Missing selector provenance manifest: "+str(manifest))
    m=json.loads(manifest.read_text(encoding="utf-8"))
    if m["selection_csv_sha256"]!=sha256(path):
        raise ValueError("Selected start table changed after creation")
    if m["selected_electrons"]!=len(reference):
        raise ValueError("Selected reference count changed")
    return reference,m


def read_mc(path,reference):
    rows={}
    attachment=set()
    step=set()
    with path.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        if not REQUIRED_MC.issubset(reader.fieldnames or []):
            raise ValueError("Wrong fast-MC CSV header")
        for row in reader:
            key=ident(row)
            if key in rows:raise ValueError("Duplicated MC trial: "+str(key))
            if key not in reference:raise ValueError("MC trial not in reference: "+str(key))
            ref=reference[key]
            if row["split"]!=ref["split"] or int(row["random_seed"])!=int(ref["random_seed"]):
                raise ValueError("Split/MC seed mismatch: "+str(key))
            for coord in ("u_start_mm","v_start_mm"):
                if not math.isclose(float(row[coord]),float(ref[coord]),rel_tol=0,abs_tol=2e-9):
                    raise ValueError("Birth coordinate changed for "+str(key))
            if sum(int(row[c]) for c in ("anode_collected",
                       "field_wire_collected","other_end"))!=1:
                raise ValueError("MC endpoint classification not exclusive")
            if int(row["drift_ok"]) not in (0,1):
                raise ValueError("Bad MC drift_ok")
            t=float(row["drift_wall_us"])
            if not math.isfinite(t) or t<0:raise ValueError("Invalid per-electron timing")
            rows[key]=row
            attachment.add(row["mc_attachment"])
            step.add(float(row["mc_step_um"]))
    if set(rows)!=set(reference):
        raise ValueError(f"MC selected-start exposure incomplete: "
                         f"{len(rows)}/{len(reference)}")
    if len(attachment)!=1 or len(step)!=1:
        raise ValueError("Mixed attachment mode or step in MC output")
    st=next(iter(step))
    if not math.isfinite(st) or st<=0:raise ValueError("Invalid MC step")
    return rows,st,next(iter(attachment))


def u_zone(u):
    a=abs(u)
    return "anode_near_u" if a<.3 else (
        "intermediate_u" if a<1.0 else "half_pitch_u")


def v_zone(v):
    a=abs(v)
    return "wire_near_v" if a<.5 else (
        "mid_gap_v" if a<2.0 else "near_cathode_v")


def phi_pdf(rows,key):
    hist=np.zeros(36,dtype=np.int64)
    times=[]
    n=0
    for ref,mc in rows:
        row=ref if key=="reference" else mc
        iscol=int(row["reference_collected"] if key=="reference"
                  else row["anode_collected"])
        if not iscol:continue
        n+=1
        phi=float(row["reference_phi_hit_deg"] if key=="reference"
                  else row["phi_hit_deg"])
        t=float(row["reference_drift_time_ns"] if key=="reference"
                else row["drift_time_ns"])
        if not(-180<=phi<=180) or not math.isfinite(phi) or not math.isfinite(t):
            raise ValueError("Collected electron has invalid phi/time")
        hist[max(0,min(35,int((phi+180.)/10.)))]+=1
        times.append(t)
    return (hist/n if n else np.zeros(36),times)


def compare(rows):
    if not rows:raise ValueError("Empty comparison region")
    n=len(rows)
    rcount=sum(int(ref["reference_collected"]) for ref,mc in rows)
    mcount=sum(int(mc["anode_collected"]) for ref,mc in rows)
    r_attach=sum(int(ref["reference_status"])==-7 for ref,mc in rows)
    m_attach=sum(int(mc["status"])==-7 for ref,mc in rows)
    differences=defaultdict(list)
    timing=[]
    for ref,mc in rows:
        k=(int(ref["iu"]),int(ref["iv"]))
        differences[k].append(int(mc["anode_collected"])-
                              int(ref["reference_collected"]))
        timing.append(float(mc["drift_wall_us"]))
    # Every box selected with fixed counts; compute cell-clustered SE for
    # the equal-cell average, NOT a pseudoreplicated electron-bin error.
    percell=np.array([np.mean(v) for v in differences.values()])
    if len(set(len(v) for v in differences.values()))!=1:
        raise ValueError("Unequal cell exposure for grouped CI")
    delta=float(percell.mean())
    stderr=float(percell.std(ddof=1)/math.sqrt(len(percell))) if len(percell)>1 else None
    phi_ref,t_ref=phi_pdf(rows,"reference")
    phi_mc,t_mc=phi_pdf(rows,"mc")
    jointly_collected=[(ref,mc) for ref,mc in rows
        if int(ref["reference_collected"]) and int(mc["anode_collected"])]
    matched_phi_ref,matched_time_ref=phi_pdf(jointly_collected,"reference")
    matched_phi_mc,matched_time_mc=phi_pdf(jointly_collected,"mc")
    return {
      "n_electrons":n,"n_independent_start_cells":len(percell),
      "reference_anode_fraction":rcount/n,
      "mc_anode_fraction":mcount/n,
      "reference_attachment_fraction":r_attach/n,
      "mc_attachment_fraction":m_attach/n,
      "delta_mc_minus_ref_anode_fraction":delta,
      "delta_95ci_cell_clustered":([delta-1.96*stderr,delta+1.96*stderr]
                                    if stderr is not None else None),
      "phi_total_variation_36bins_conditional_survival":
           float(.5*np.abs(phi_ref-phi_mc).sum()),
      "phi_ref_36bins":phi_ref.tolist(),"phi_mc_36bins":phi_mc.tolist(),
      "n_jointly_collected_at_same_birth_sites":len(jointly_collected),
      "phi_tv_36bins_jointly_collected_start_sites":
        (float(.5*np.abs(matched_phi_ref-matched_phi_mc).sum())
         if jointly_collected else None),
      "mean_drift_ref_ns_jointly_collected":
        (float(np.mean(matched_time_ref)) if matched_time_ref else None),
      "mean_drift_mc_ns_jointly_collected":
        (float(np.mean(matched_time_mc)) if matched_time_mc else None),
      "mean_drift_ref_ns":float(np.mean(t_ref)) if t_ref else None,
      "mean_drift_mc_ns":float(np.mean(t_mc)) if t_mc else None,
      "rms_drift_ref_ns":float(np.std(t_ref)) if t_ref else None,
      "rms_drift_mc_ns":float(np.std(t_mc)) if t_mc else None,
      "mean_mc_drift_us":float(np.mean(timing)),
      "median_mc_drift_us":float(np.median(timing))
    }


def report(ref_path,mc_paths,output):
    if output.exists():raise FileExistsError("Refusing to overwrite report: "+str(output))
    ref,provenance=read_micro(ref_path)
    payload={
      "schema":"PI_FAST_LUT_EXISTING_MICRO_MC_STEP_VALIDATION_V1",
      "microscopic_reference_sha256":sha256(ref_path),
      "source_micro_shas":{
        name:data["sha256"] for name,data in provenance["source_provenance"].items()},
      "description":"Only MC drifts run. Micro endpoints reused from frozen "
        "500-trial LUT. Do not infer survival accuracy with attachment off.",
      "comparisons":{}}
    plots=[]
    observed_steps=set()
    for path in mc_paths:
        runs,step,attachment=read_mc(path,ref)
        if step in observed_steps:raise ValueError("Duplicate MC step file")
        observed_steps.add(step)
        key=f"{step:g}um"
        by_group={}
        for geom in GEOMS:
            for split in ("fit","validation"):
                paired=[(ref[k],runs[k]) for k in ref
                        if k[0]==geom and ref[k]["split"]==split]
                if not paired:raise ValueError("Missing geometry/split: "+geom+"/"+split)
                base=compare(paired)
                zoned={}
                for dimension,f in (("u",u_zone),("v",v_zone)):
                    labels=sorted({f(float(row[0][f"{dimension}_start_mm"]))
                                   for row in paired})
                    zoned[dimension]={label:compare([
                      x for x in paired if f(float(x[0][f"{dimension}_start_mm"]))==label])
                         for label in labels}
                by_group[geom+"/"+split]={"global":base,"regions":zoned}
                print(f"Step {key} {geom} {split}: "
                      f"micro anode={base['reference_anode_fraction']:.4f}, "
                      f"MC={base['mc_anode_fraction']:.4f}, "
                      f"phi TV={base['phi_total_variation_36bins_conditional_survival']:.4f}, "
                      f"MC drift {base['mean_mc_drift_us']:.1f}us/e",
                      flush=True)
                if split=="validation":
                    plots.append((step,geom,base))
        payload["comparisons"][key]={
            "input":str(path),"csv_sha256":sha256(path),
            "mc_step_um":step,"mc_attachment":attachment,
            "by_geometry_split":by_group}
    payload["caveats"]=[
      "Attachment off is an intentionally unphysical diagnostic.",
      "Phi/time PDFs condition on each method's OWN anode survivors; "
      "attachment changes which electrons survive.",
      "36-bin phi TV has appreciable finite-sample fluctuation for 4 "
      "validation electrons per cell; not a precision acceptance gate.",
      "Paired collection CI clusters by spatial cell because 500 LUT "
      "microscopic replicates are the independent original trials.",
      "No eventwise Heed/Geant4 folding is done at this stage.",
      "MC CPU timing includes DriftElectron invocation only, not "
      "Magboltz table generation and serialisation."
    ]
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("NOTE: matplotlib missing; quantitative JSON created.")
    else:
        fig,axes=plt.subplots(1,2,figsize=(13,4.5),sharey=True)
        centers=np.arange(-175.,180.,10.)
        for i,geom in enumerate(GEOMS):
            items=sorted((a,b) for a,g,b in plots if g==geom)
            if not items:continue
            baseline=np.asarray(items[0][1]["phi_ref_36bins"])
            axes[i].step(centers,baseline,where="mid",
                         color="#333333",linewidth=2,label="Microscopic LUT500")
            for step,b in items:
                axes[i].step(centers,b["phi_mc_36bins"],where="mid",
                             label=f"MC no-attachment {step:g}um",alpha=.75)
            axes[i].set_title(geom+" (validation split)")
            axes[i].set_xlabel("Collected electron phi_hit [degrees]")
            axes[i].legend(frameon=False,fontsize=8)
        axes[0].set_ylabel("Probability per 10-degree bin")
        fig.suptitle("Frozen micro start sites vs un-killed MC: angular transport control")
        fig.tight_layout()
        img=output.with_suffix(".png")
        if img.exists():raise FileExistsError("Refusing PNG overwrite: "+str(img))
        fig.savefig(img,dpi=150)
        plt.close(fig)
        print("PASS plot:",img)
    print("PASS numerical paired-transport report:",output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--mc",type=Path,nargs="+",required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    report(a.reference,a.mc,a.out)

if __name__=="__main__":main()

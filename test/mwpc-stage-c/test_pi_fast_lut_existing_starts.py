#!/usr/bin/env python3
"""Small deterministic, synthetic tests; no Garfield or frozen LUT files required."""
from __future__ import annotations
import csv,hashlib,json,tempfile
from pathlib import Path
import numpy as np
from pi_fast_lut_select_reference import NCELL,GEOMS,sample_plan
from pi_fast_lut_existing_starts_report import (
  read_micro,read_mc,compare,phi_pdf)

HEADER=["geometry","split","iu","iv","sample","random_seed","u_start_mm",
    "v_start_mm","reference_collected","reference_field_wire_collected",
    "reference_status","reference_wire_index","reference_phi_hit_deg",
    "reference_drift_time_ns"]
MC_HEADER=["geometry","split","iu","iv","sample","random_seed","u_start_mm",
    "v_start_mm","mc_step_um","mc_attachment","drift_ok","status",
    "u_end_mm","v_end_mm","w_end_mm","drift_time_ns","anode_collected",
    "field_wire_collected","other_end","wire_index","phi_hit_deg","drift_wall_us"]

def main():
    for geom in GEOMS:
        p=sample_plan(geom,4,4,276014)
        assert len(p)==NCELL*8
        assert p==sample_plan(geom,4,4,276014)
        assert {v[0] for v in p.values()}=={"fit","validation"}
        assert all(v[1]>0 for v in p.values())
        for cell in range(NCELL):
            assert sum(i//500==cell for i in p)==8
    with tempfile.TemporaryDirectory(prefix="pi_fast_lut_test_") as t:
        directory=Path(t)
        refpath=directory/"reference.csv"
        mcpath=directory/"mc.csv"
        rows=[]
        sim=[]
        for geom in GEOMS:
            for split in ("fit","validation"):
                for k in range(8):
                    iu,iv=1+k,2+k
                    # Fit and validation must have disjoint frozen trial IDs.
                    sample=k+(8 if split=="validation" else 0)
                    seed=1234+sample
                    startu=-1.5+.1*(iu+.5)
                    startv=-3.+.2*(iv+.5)
                    collected=int(k!=0)
                    row=dict(geometry=geom,split=split,iu=iu,iv=iv,
                       sample=sample,random_seed=seed,u_start_mm=startu,
                       v_start_mm=startv,reference_collected=collected,
                       reference_field_wire_collected=0,
                       reference_status=-5 if collected else -7,
                       reference_wire_index=0,
                       reference_phi_hit_deg=90. if collected else "nan",
                       reference_drift_time_ns=17. if collected else "nan")
                    rows.append(row)
                    sim.append(dict(
                       geometry=geom,split=split,iu=iu,iv=iv,
                       sample=sample,random_seed=seed,u_start_mm=startu,
                       v_start_mm=startv,mc_step_um=5,mc_attachment="off",
                       drift_ok=1,status=-5,u_end_mm=0.,v_end_mm=.015,
                       w_end_mm=0.,drift_time_ns=18.,
                       anode_collected=1,field_wire_collected=0,other_end=0,
                       wire_index=0,phi_hit_deg=90.,drift_wall_us=180.))
        with refpath.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=HEADER)
            w.writeheader()
            w.writerows(rows)
        sha=hashlib.sha256(refpath.read_bytes()).hexdigest()
        refpath.with_suffix(".manifest.json").write_text(json.dumps({
            "selection_csv_sha256":sha,"selected_electrons":len(rows),
            "source_provenance":{geom:{"sha256":"test"} for geom in GEOMS}}))
        with mcpath.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=MC_HEADER)
            w.writeheader()
            w.writerows(sim)
        # The anti-duplicate check must reject a deliberately repeated
        # microscopic trial ID, including duplicate copies across splits.
        bad=rows.copy()
        bad[-1]={**bad[-1],"sample":bad[-1]["sample"]-8,
                 "split":"validation"}
        with refpath.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=HEADER)
            writer.writeheader()
            writer.writerows(bad)
        try:
            read_micro(refpath)
        except ValueError as e:
            assert "Duplicate microscopic start" in str(e)
        else:
            raise AssertionError("Failed to reject duplicate reference ID")
        with refpath.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=HEADER)
            writer.writeheader()
            writer.writerows(rows)
        # The regenerated content must match the original frozen hash.
        assert hashlib.sha256(refpath.read_bytes()).hexdigest()==sha
        ref,_=read_micro(refpath)
        mc,step,attach=read_mc(mcpath,ref)
        assert len(mc)==len(ref)==len(rows)
        assert step==5. and attach=="off"
        for geom in GEOMS:
            pairs=[(ref[k],mc[k]) for k in ref
                   if k[0]==geom and ref[k]["split"]=="validation"]
            stats=compare(pairs)
            assert np.isclose(stats["reference_anode_fraction"],.875)
            assert stats["mc_anode_fraction"]==1.
            assert np.isclose(stats["delta_mc_minus_ref_anode_fraction"],.125)
            assert np.isclose(stats["phi_total_variation_36bins_conditional_survival"],0.)
            assert stats["n_independent_start_cells"]==8
            assert stats["delta_95ci_cell_clustered"] is not None
        broken=sim.copy()
        broken[0]={**broken[0],"u_start_mm":999}
        with mcpath.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=MC_HEADER)
            w.writeheader();w.writerows(broken)
        try:
            read_mc(mcpath,ref)
        except ValueError as e:
            assert "Birth coordinate changed" in str(e)
        else:
            raise AssertionError("Fail-open paired-start coordinates")
    print("PASS fast-LUT: 900-cell reproducible selection, fit/validation,"
          " complete paired IDs, spatial cell CI, phi PDFs and mismatch rejection")

if __name__=="__main__":main()

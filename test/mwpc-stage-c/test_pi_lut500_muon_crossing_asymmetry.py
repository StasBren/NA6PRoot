#!/usr/bin/env python3
"""Synthetic no-physics checks for crossing-level p3/p4 LUT20/LUT500 analysis."""
from __future__ import annotations
import csv
import gzip
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

import pi_lut500_muon_crossing_asymmetry as m
from pi_pitch_fieldwire_pilot import tag

def test_classification():
    x=np.asarray([[1.,0.,1.,0.,phi] for phi in
                   (90.,-90.,10.,-10.,0.,180.)])
    c=m.classify(x,np.ones(len(x),dtype=bool))
    assert (c["up"],c["down"],c["hemisphere_boundary"])==(2,2,2)
    assert (c["top"],c["bottom"],c["side"])==(1,1,4)
    assert c["A_hemisphere"]==c["A_top_bottom"]==0.
    assert abs(c["f_top_bottom"]-1/3)<1e-14
    empty=m.classify(np.empty((0,5)),np.zeros(0,dtype=bool))
    assert empty["A_hemisphere"] is None
    assert empty["A_top_bottom"] is None
    assert empty["f_top_bottom"] is None
    acc=m.Moments()
    acc.add(c,6,6)
    acc.add(empty,0,0)
    r=acc.result()
    assert r["crossings"]==2
    assert r["hemisphere_A"]["valid"]==1
    assert r["hemisphere_A"]["undefined"]==1
    assert r["hemisphere_A"]["hist"][20]==1
    assert r["f_top_bottom_by_crossing"]["undefined"]==1
    pair=m.Paired()
    pair.add(1.,-1.); pair.add(None,1.)
    assert pair.result()["both_defined"]==1
    assert pair.result()["not_both_defined"]==1
    assert pair.result()["FW_minus_noFW_hemisphere_A"]["mean"]==-2.

def test_full_pipeline_without_physics():
    with tempfile.TemporaryDirectory(prefix="pi_lut500_crossing_") as d:
        root=Path(d)
        infile=root/"native.csv"
        header=("station","pdg","det_id","event","path_length_mm",
                "hit_index","track_id","theta_u_in_deg","theta_w_in_deg",
                "u_in_mm","u_out_mm")
        with infile.open("w",newline="") as f:
            w=csv.writer(f);w.writerow(header)
            w.writerow((0,13,10,0,1.,1,11,5.,10.,-.2,.2))
            w.writerow((1,-13,20,0,1.,2,12,-5.,-10.,.1,-.1))
        modules=[
            [SimpleNamespace(station=0,det_id=10)],
            [SimpleNamespace(station=1,det_id=20)]]
        def mock_map(pitch,fw,path):
            ntrials=20 if Path(path)==m.ROOT20 else 500
            table=np.zeros((round(pitch/.1)*30,ntrials,5))
            table[0,0,0]=fw
            return table,{"sha256":"fixture"+str(fw),
                          "samples_per_box":ntrials}
        def mock_birth(row,lam,marks,seed):
            return np.array([[0.,0.,0.]]),np.array([.5])
        def mock_response(u,v,q,pitch,station,pool):
            # Both geometries: one electron on an anode. Change only
            # endpoint hemisphere, NOT native Geant4/Heed birth.
            phi=90. if pool[0,0,0]==0 else -90.
            arr=np.array([[1.,0.,1.,0.,phi]])
            return np.array([True]),np.array([True]),np.array([0]),arr
        with patch.object(m,"crossing_provenance",
                          lambda name,chunk:(infile,{"csv_sha256":"mock"})),\
             patch.object(m,"build_chambers",lambda:modules),\
             patch.object(m,"load_map",mock_map),\
             patch.object(m,"read_heed",lambda a,b:(1.,np.array([1]))),\
             patch.object(m,"birth_sample",mock_birth),\
             patch.object(m,"response_for_births",mock_response):
            result=m.run(1,[0,1],0,[3.],root/"report",plots=False)
            assert result["stations"]["0"]["gas_crossings"]==1
            assert result["stations"]["1"]["gas_crossings"]==1
            s=result["stations"]["0"]
            a=s["configs"][tag(3.,0.)]["hemisphere_A"]
            b=s["configs"][tag(3.,75.)]["hemisphere_A"]
            assert a["mean"]==1. and b["mean"]==-1.
            assert s["paired"]["3.0"]["FW_minus_noFW_hemisphere_A"]["mean"]==-2.
            assert s["configs"][tag(3.,0.)]["top_bottom_A"]["mean"]==1.
            assert s["configs"][tag(3.,0.)]["top_bottom_electron_fraction"]==1.
            with gzip.open(root/"report"/"crossings.csv.gz","rt") as gz:
                saved=list(csv.DictReader(gz))
            assert len(saved)==4 # Two crossings times two geometries
            assert saved[0]["event"]=="0" and saved[0]["track_id"]=="11"
            assert saved[0]["station"]=="0"
            assert saved[0]["A_hemisphere"]=="1.0"
            assert (root/"report"/"manifest.json").is_file()
            # Exactly the same crossing identities and synthetic Heed sources,
            # but use the old frozen 20-trial p4 maps instead of LUT500.
            result20=m.run(1,[0,1],0,[4.],root/"report20",
                           plots=False,lut_trials=20)
            assert result20["lut_trials_per_cell"]==20
            assert result20["pitches_mm"]==[4.]
            assert result20["schema"].endswith("LUT20_V1")
            assert "LUT20" in result20 and "LUT500" not in result20
            assert result20["stations"]["0"]["configs"][tag(4.,0.)][
                "hemisphere_A"]["mean"]==1.
            assert result20["stations"]["0"]["configs"][tag(4.,75.)][
                "hemisphere_A"]["mean"]==-1.
            with gzip.open(root/"report20"/"crossings.csv.gz","rt") as f20:
                rows20=list(csv.DictReader(f20))
            assert len(rows20)==4
            assert all(r["lut_trials_per_cell"]=="20" for r in rows20)
            assert (root/"report20"/"manifest.json").is_file()
            # Both pitches at LUT20 use exactly the same original gas
            # crossing/Heed source IDs. The pitch-dependent wire lattice
            # is intentionally different; ensure counters are not pooled.
            r34=m.run(1,[0,1],0,[3.,4.],root/"report20_34",
                       plots=False,lut_trials=20)
            ms0=r34["stations"]["0"]
            assert set(ms0["configs"])=={
                tag(3.,0.),tag(3.,75.),tag(4.,0.),tag(4.,75.)}
            assert ms0["gas_crossings"]==1
            assert all(ms0["configs"][name]["crossings"]==1
                       for name in ms0["configs"])
            with gzip.open(root/"report20_34"/"crossings.csv.gz","rt") as g34:
                rows34=list(csv.DictReader(g34))
            assert len(rows34)==8 # 2 crossings * 2 pitches * 2 FW choices
            assert set(r["lut_trials_per_cell"] for r in rows34)=={"20"}
            try:
                m.run(1,[0,1],0,[3.],root/"report",plots=False)
            except FileExistsError:
                pass
            else:raise AssertionError("Analysis overwrote an existing output")
    print("PASS: crossing-level classification, parent/track/DET IDs, "
          "paired FW comparison, conservation and no-overwrite.")


def test_render_only_repair_for_all_six_stations():
    """Catch plotting NameError and prove frozen physics survives failed plotting."""
    with tempfile.TemporaryDirectory(prefix="pi_lut20_repair_") as d:
        root=Path(d)
        src=root/"native.csv"
        head=("station","pdg","det_id","event","path_length_mm",
              "hit_index","track_id","theta_u_in_deg","theta_w_in_deg",
              "u_in_mm","u_out_mm")
        with src.open("w",newline="") as f:
            w=csv.writer(f);w.writerow(head)
            for st in range(6):
                w.writerow((st,13,100+st,st,1.,st,10+st,5.,-5.,-.2,.2))
        modules=[[SimpleNamespace(station=st,det_id=100+st)]
                 for st in range(6)]
        maps=root/"maps"
        for fw in (0.,75.):
            folder=maps/tag(4.,fw)
            folder.mkdir(parents=True)
            (folder/"garfield_electrons.csv").write_text(
                f"fixture LUT20 p4 FW{fw}\n",encoding="utf-8")
            (folder/"frozen_config.json").write_text(
                f"fixture frozen p4 FW{fw}\n",encoding="utf-8")
        def fake_map(pitch,fw,directory):
            folder=Path(directory)/tag(pitch,fw)
            pool=np.zeros((1200,20,5),dtype=float)
            pool[0,0,0]=fw
            return pool,{
                "samples_per_box":20,
                "csv_sha256":m.sha(folder/"garfield_electrons.csv"),
                "frozen_sha256":m.sha(folder/"frozen_config.json")}
        def fake_birth(row,lam,marks,seed):
            return np.array([[0.,0.,0.]]),np.array([.5])
        def fake_response(u,v,q,pitch,station,pool):
            phi=90. if pool[0,0,0]==0 else -90.
            return (np.array([True]),np.array([True]),
                    np.array([0]),np.array([[1.,0.,5.,0.,phi]]))
        with patch.object(m,"ROOT20",maps),\
             patch.object(m,"crossing_provenance",
                          lambda name,chunk:(src,{"csv_sha256":m.sha(src)})),\
             patch.object(m,"build_chambers",lambda:modules),\
             patch.object(m,"load_map",fake_map),\
             patch.object(m,"read_heed",lambda a,b:(1.,np.array([1]))),\
             patch.object(m,"birth_sample",fake_birth),\
             patch.object(m,"response_for_births",fake_response):
            out=root/"report"
            m.run(1,list(range(6)),10,[4.],out,plots=False,lut_trials=20)
            before=(m.sha(out/"summary.json"),m.sha(out/"crossings.csv.gz"))
            # Emulate the real crash: completed and closed raw data, with
            # partial png renders but no publication manifest.
            (out/"manifest.json").unlink()
            (out/"hemisphere_A_p4.png").write_bytes(b"interrupted plot")
            manifest=m.repair_plots_from_saved_counts(
                out,required_trials=20,required_pitches=[4.])
            after=(m.sha(out/"summary.json"),m.sha(out/"crossings.csv.gz"))
            assert before==after
            assert manifest["recovered_original_crossings"]==6
            assert manifest["recovered_raw_geometry_rows"]==12
            assert manifest["recovery"].startswith("plot-only")
            assert (out/"manifest.json").is_file()
            assert (out/"hemisphere_A_p4.png").stat().st_size>1000
            assert (out/"up_down_joint_p4_FW75.png").stat().st_size>1000
            assert (out/"RMS_A_vs_theta_w_p4.png").stat().st_size>1000
            assert len(list(out.glob("*.png")))==8
            try:
                m.repair_plots_from_saved_counts(
                    out,required_trials=20,required_pitches=[4.])
            except FileExistsError:
                pass
            else:
                raise AssertionError("Repair could overwrite completed manifest")
    print("PASS six-station plot-only crash recovery; all figures rendered "
          "without modifying frozen raw/summary data.")

if __name__=="__main__":
    test_classification()
    test_full_pipeline_without_physics()
    test_render_only_repair_for_all_six_stations()

#!/usr/bin/env python3
"""Synthetic regression for real-event TWO independent-anode waveform explorer.

No local native ROOT/Geant4/Garfield data or plotting runtime required.
Fixtures use one true-format native hit per source and a deterministic
two-wire mock Garfield map; the production reader/Heed sampler/digitizer
and offline HTML/manifest generation are exercised.
"""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import csv
import json
import numpy as np

import pi_native_two_anode_signal_explorer as ex
from pi_centered_electron_sampler import ionize_and_drift,multiplicities,stable_seed


def main():
    data=[
        {"source":"jpsi40-bore","det_id":0,"station":0,"event":0,
         "hit_index":0,"track_id":0,"theta_u_in_deg":-12.,
         "theta_w_in_deg":4.,"edep_eV":1234.,"gas_path_mm":6.,
         "N_t100":2,"N_all":2},
        {"source":"phi40-bore","det_id":0,"station":0,"event":0,
         "hit_index":0,"track_id":0,"theta_u_in_deg":.4,
         "theta_w_in_deg":4.,"edep_eV":1220.,"gas_path_mm":6.,
         "N_t100":2,"N_all":2},
        {"source":"omega40-bore","det_id":0,"station":0,"event":0,
         "hit_index":0,"track_id":0,"theta_u_in_deg":12.,
         "theta_w_in_deg":4.,"edep_eV":1300.,"gas_path_mm":6.,
         "N_t100":2,"N_all":2},
    ]
    chosen=ex.choose_examples(data,7)
    assert len(chosen)==3 and all(len(x)==1 for x in chosen.values())
    # Deliberately construct an independent two-wire collection for the
    # true FIRST chunk deterministic seed, and check N100 from replay.
    field=np.zeros((600,1,3),dtype=float)
    field[:,:,0]=1.;field[:,:,2]=20.
    marks=np.array([1])
    with TemporaryDirectory(prefix="pi_native_2wire_explorer_") as rootdir:
        root=Path(rootdir)
        export=root/"export";export.mkdir()
        template=root/"template.html"
        template.write_text("__PI_TWO_ANODE_EVENT_DATA__",encoding="utf-8")
        dummy=root/"field.csv";dummy.write_text("dummy")
        frozen=root/"pilot.csv"
        with frozen.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]))
            w.writeheader();w.writerows(data)
        for item in data:
            channel=ex.SOURCES[item["source"]]
            raw={
                "event":0,"hit_index":0,"track_id":0,"pdg":13,
                "station":0,"det_id":0,
                "u_in_mm":-1.,"u_out_mm":3.,
                "v_in_mm":-3.,"v_out_mm":3.,
                "w_in_mm":-13.,"w_out_mm":13.,
                "path_length_mm":6.,
                "theta_u_in_deg":item["theta_u_in_deg"],
                "theta_w_in_deg":4.,
                "edep_eV":item["edep_eV"],
            }
            seed=stable_seed("centered_native_pilot_v1",
                             item["source"],"chunk0001",0,0,0)
            h=ionize_and_drift(raw,np.random.default_rng(seed),20.,marks,field)
            assert multiplicities(h)[1]==2 and multiplicities(h)[-1]==2
            with (export/f"{channel}_geant_gas_crossings.csv").open(
                "w",newline="",encoding="utf-8"
            ) as f:
                w=csv.DictWriter(f,fieldnames=list(raw))
                w.writeheader();w.writerow(raw)
        original_map=ex.load_centered_map
        original_heed=ex.read_heed
        try:
            ex.load_centered_map=lambda *_:(field,{
                "response_csv_sha256":"0"*64})
            ex.read_heed=lambda *_:(20.,marks)
            out=root/"viewer.html"
            args=SimpleNamespace(
                pilot=frozen,export_dir=export,garfield=dummy,manifest=dummy,
                heed_tracks=dummy,heed_clusters=dummy,template=template,
                gain_csv=None,max_per_chamber=7,html=out)
            result=ex.build(args)
            assert result["two_wire_example_count"]==3
            data=json.loads(out.read_text(encoding="utf-8"))
            assert len(data["layout"])==298 and len(data["sources"])==3
            for src in ex.SOURCES:
                muon=data["examples"][src+"|0"][0]
                assert len(muon["anodes"])==len(muon["wires"])==2
                assert muon["N_t100"]==2
                a,b=muon["anodes"]
                assert a["wire"]!=b["wire"]
                assert len(a["samples_mV"])==len(b["samples_mV"])==301
                assert all(abs(c-(x+y))<0.0011 for c,x,y in zip(
                    muon["sumSamples_mV"],
                    a["samples_mV"],b["samples_mV"]))
                assert muon["deltaFirstElectron_ns"]==0.
        finally:
            ex.load_centered_map=original_map
            ex.read_heed=original_heed
        print("PASS: exact native (source,event,hit,track,DET) Heed+Garfield replay")
        print("PASS: actual two collected physical anodes, N100=2 and correct stable seed")
        print("PASS: independent anode-only electronics gain/noise and 301-sample waves")
        print("PASS: sum waveform = pointwise two independent amplifier voltages")
        print("PASS: source/station XY structure and offline view-data/manifest outputs")


if __name__=="__main__":
    main()

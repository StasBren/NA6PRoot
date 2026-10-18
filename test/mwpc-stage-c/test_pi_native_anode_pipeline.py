#!/usr/bin/env python3
"""Synthetic new-Geant4 -> true local gas hit -> Heed/Garfield endpoint
-> anode-only waveform check. No ROOT files, Geant4 or Garfield executable.
"""
import csv
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from prepare_pi_native_hit_explorer import convert
from build_pi_muon_track_explorer import describe_row
from pi_anode_digitizer import threshold_features, AnodeConfig


def main():
    with TemporaryDirectory(prefix="pi_native_electronics_") as tmp:
        root=Path(tmp)
        inp=root/"native_hits.csv"
        output=root/"native_viewer.csv"
        data={
            "event":"3","hit_index":"5","track_id":"19",
            "pdg":"13","station":"1","det_id":"32",
            "x_in_cm":"0","y_in_cm":"0","z_in_cm":"340.1",
            "x_out_cm":"0.1","y_out_cm":"0.01","z_out_cm":"340.7",
            "u_in_mm":"-.2","v_in_mm":"-3",
            "u_out_mm":"-.1","v_out_mm":"3",
            "path_length_mm":"6.1","edep_eV":"1200.",
            "time_ns":"21.","px_in_GeV":"0.3","py_in_GeV":"0.2",
            "pz_in_GeV":"20.",
            "theta_u_in_deg":"0.573","theta_w_in_deg":"0.859",
        }
        with inp.open("w",newline="",encoding="utf-8") as stream:
            w=csv.DictWriter(stream,fieldnames=list(data))
            w.writeheader();w.writerow(data)
            w.writerow({**data,"event":"4","pdg":"11"})
        assert convert(inp,output,"jpsi40-bore","native_jpsi_chunk0001")==1
        rows=list(csv.DictReader(output.open(newline="",encoding="utf-8")))
        assert len(rows)==1
        row=rows[0]
        assert row["native_hit"]=="1" and row["pdg"]=="13"
        assert math.isclose(float(row["edep_eV"]),1200.)
        # Synthetic saved Garfield map: each of 600 0.2-mm boxes collects
        # onto its nearest wire, and every drift time is 50ns.
        offsets=np.arange(600,dtype=np.int64)
        sizes=np.ones(600,dtype=np.int64)
        response=np.zeros((600,3),dtype=float)
        response[:,0]=1.
        response[:,2]=50.
        models=(3.0,np.array([1,2,3],dtype=int),
                (offsets,sizes,response))
        result=describe_row(row,models,True,None)
        assert result["isRealGeantGasHit"] is True
        assert result["edep_eV"]==1200.
        assert result["time_ns"]==21.
        assert result["length"]==6.1
        assert len(result["geom"])==1
        assert result["charges"] is not None
        assert result["nClusters"]>=0
        charged=sum(e[1] for e in result["charges"])
        assert result["lostElectrons"]==0
        assert sum(v["nPrimaryCollected"]
                   for v in result["anodeSignals"].values())==charged
        assert all(len(v["samples_mV"])==301
                   for v in result["anodeSignals"].values())
        assert all(v["nUnknownOrOutOfWindow"]==0
                   for v in result["anodeSignals"].values())
        result2=describe_row(row,models,True,None)
        assert result2==result, "Native event/track/hit reproducibility failed"
        cfg=AnodeConfig()
        for v in result["anodeSignals"].values():
            assert not threshold_features(v["samples_mV"],cfg,100000.)["fired"]
    print("PASS native Geant4 TGeo-local true gas-hit Edep and muon-PDG converter")
    print("PASS native Heed/Garfield individual electron endpoints -> anode-only analog")
    print("PASS reproducible gain/shaping and no false conversion of Geant4 Edep")


if __name__=="__main__":
    main()

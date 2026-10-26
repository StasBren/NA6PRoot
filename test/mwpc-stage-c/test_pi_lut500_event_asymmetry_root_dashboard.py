#!/usr/bin/env python3
"""Synthetic parent/DET conservation and ROOT dashboard smoke, no large data."""
import csv
import gzip
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pi_lut500_event_asymmetry_root_dashboard as d


def full_test():
    with tempfile.TemporaryDirectory(prefix="mwpc_asym_root_") as tmp:
        root=Path(tmp)
        source=root/"fake_crossings.csv"
        header=["station","pdg","det_id","event","path_length_mm",
                "hit_index","track_id"]
        rows=[
            [0,13,10,0,1,0,0],  # upper in DET 10
            [0,13,11,0,1,1,1],  # lower in DET 11: station A=0
            [0,13,10,1,1,2,2],  # upper in DET 10
            [0,13,11,2,0,3,3],  # zero-length crossing -> undefined A
            [1,13,20,0,1,4,4],  # station MS1 independent
        ]
        with source.open("w",encoding="utf-8",newline="") as f:
            w=csv.writer(f);w.writerow(header);w.writerows(rows)
        modules=[[SimpleNamespace(station=0,det_id=10,row=0,col=0,
                        layer="A",x=-1.,y=0.,half_x=.3,half_y=.3),
                  SimpleNamespace(station=0,det_id=11,row=0,col=1,
                        layer="A",x=+1.,y=0.,half_x=.3,half_y=.3)],
                 [SimpleNamespace(station=1,det_id=20,row=0,col=0,
                        layer="B",x=0.,y=0.,half_x=.3,half_y=.3)]]
        # Same generated birth positions for both Garfield geometries;
        # endpoint phi chosen from original hit_index for fixture only.
        def synthetic_birth(row,lam,marks,seed):
            return np.array([[0.,0.,float(row["hit_index"])]]),np.array([.3])
        def synthetic_response(u,v,q,pitch,st,table):
            index=int(table[0,0,0]) # marker gives noFW=0, FW=1
            hitnum=int(v[0]) if len(v) else 0
            phi=+90. if hitnum in (0,2,4) else -90.
            samples=np.array([[1.,0.,8.,0.,phi]]) if len(v) else np.empty((0,5))
            return np.ones(len(v),dtype=bool),np.ones(len(v),dtype=bool),np.zeros(len(v)),samples
        def synthetic_map(pitch,diameter,root):
            arr=np.zeros((900,500,5),dtype=float)
            arr[0,0,0]=1 if diameter>0 else 0
            return arr,{"sha256":"fixture-"+str(diameter)}
        raw=root/"view.html"
        raw.write_text("<html><body><script>const D=__PAYLOAD__;</script></body></html>")
        with patch.object(d,"TEMPLATE",raw),\
             patch.object(d,"crossing_provenance",lambda source_key,c:(source,{"csv_sha256":"fixture"})),\
             patch.object(d,"build_chambers",lambda:modules),\
             patch.object(d,"load_map",synthetic_map),\
             patch.object(d,"read_heed",lambda a,b:(.1,np.array([1.]))),\
             patch.object(d,"birth_sample",synthetic_birth),\
             patch.object(d,"response_for_births",synthetic_response),\
             patch.object(d,"coincident_zero_chord",lambda row:True):
            out=root/"dash"
            result=d.summary_for_chunks(1,[0,1],0,out,True)
            assert result["number_of_source_chunks"]==1
            assert result["partial_exposure"] is False
            assert result["n_stations_with_computed_results"]==2
            st=result["stations"]["0"]
            parent=st["station_parent"]["configurations"][d.CONFIGS[0]]
            assert parent["parents_with_gas_crossing"]==3
            assert parent["valid_A"]==2 and parent["undefined_A"]==1
            assert parent["hist"][20]==1 # A=0 in station parent 0
            # DO NOT sum per-chamber histograms to get station A!
            c10=st["chambers"]["10"]["configurations"][d.CONFIGS[0]]
            c11=st["chambers"]["11"]["configurations"][d.CONFIGS[0]]
            assert c10["valid_A"]==2 and c11["valid_A"]==1
            assert c10["hist"][40]==2 and c11["hist"][0]==1
            assert parent["mean_A"]==.5 and parent["rms_A"]==.5
            assert st["station_parent"]["paired"]["both_valid_parents"]==2
            assert result["stations"]["1"]["station_parent"]["configurations"][d.CONFIGS[0]]["valid_A"]==1
            assert (out/"dashboard.html").is_file()
            assert (out/"manifest.json").is_file()
            assert "__PAYLOAD__" not in (out/"dashboard.html").read_text()
            with gzip.open(out/"event_counts.csv.gz","rt") as f:
                data=list(csv.DictReader(f))
            assert len(data)==8 # 3 MS0 station, 3 MS0 DET, 1 MS1 station, 1 MS1 DET
            try:
                d.summary_for_chunks(1,[0,1],0,out,True)
            except FileExistsError:pass
            else:raise AssertionError("Could overwrite dashboard without warning")
    print("PASS ROOT event-asymmetry fixture: station and DET parent pooling, "
          "valid/undefined conservation, raw event gzip, immutable output")


if __name__=="__main__":
    full_test()

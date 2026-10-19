#!/usr/bin/env python3
"""Tiny end-to-end first-chunk native timing pilot fixture (no real Garfield).

Checks all three source files, TGeo local u/v/w, all six stations including
zero exposure, angle grouping, histograms, per-hit PDF, HTML, report and ZIP
prerequisites; NEVER claims this toy data is physics or modifies frozen runs.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

import pi_native_centered_timing_pilot as pilot


def main():
    old_loader=pilot.load_centered_map
    old_heed=pilot.read_heed
    old_expected=pilot.FIRST_CHUNK_EXPECTED
    try:
        with TemporaryDirectory(prefix="pi_centered_timing_fixture_") as tmp:
            root=Path(tmp)
            exp=root/"exports";exp.mkdir()
            out=root/"outputs"
            fake=root/"dummy.csv";fake.write_text("unused\n")
            temp=root/"dummy.html"
            temp.write_text("__CENTERED_PILOT_PAYLOAD__",encoding="utf-8")
            # One valid native muon gas hit and one NON-muon detector hit per
            # channel. PDG filter must not count nonmuons in PI denominators.
            row={
                "event":0,"hit_index":1,"track_id":0,"pdg":13,
                "station":0,"det_id":0,"path_length_mm":6.,
                "edep_eV":1200.,
                "u_in_mm":0.,"u_out_mm":0.,
                "v_in_mm":-3.,"v_out_mm":3.,
                "w_in_mm":-15.,"w_out_mm":15.,
                "theta_u_in_deg":1.,"theta_w_in_deg":0.,
            }
            for src,channel in pilot.SOURCES.items():
                with (exp/f"{channel}_geant_gas_crossings.csv").open(
                    "w",newline="",encoding="utf-8"
                ) as stream:
                    w=csv.DictWriter(stream,fieldnames=list(row))
                    w.writeheader();w.writerow(row)
                    w.writerow({**row,"pdg":211,"hit_index":2})
            pool=np.zeros((600,1,3),dtype=float)
            pool[:,:,0]=1.;pool[:,:,2]=20.
            pilot.load_centered_map=lambda *_: (
                pool,{"response_csv_sha256":"3"*64})
            pilot.read_heed=lambda *_: (15.,np.array([1],dtype=np.int32))
            pilot.FIRST_CHUNK_EXPECTED={src:1 for src in pilot.SOURCES}
            args=SimpleNamespace(export_dir=exp,garfield=fake,manifest=fake,
                heed_tracks=fake,heed_clusters=fake,template=temp,
                html=root/"built.html",out=out)
            pilot.run(args)
            with (out/"pi_centered_native_hit_multiplicities.csv").open() as stream:
                hits=list(csv.DictReader(stream))
            assert len(hits)==3
            assert all(int(r["N_all"])==1 for r in hits)
            assert all(int(r["n_missing_time"])==0 for r in hits)
            with (out/"pi_centered_time_angle_multiplicity.csv").open() as stream:
                table=list(csv.DictReader(stream))
            assert len(table)>0
            for src in pilot.SOURCES:
                for st in range(6):
                    baseline=[r for r in table if r["source"]==src
                         and int(r["station"])==st and int(r["det_id"])==-1
                         and r["angle_axis"]=="u" and int(r["angle_bin"])==-1
                         and r["window_ns"]=="all"]
                    assert len(baseline)==1
                    assert int(baseline[0]["gas_hit_segments"])==(1 if st==0 else 0)
            payload=json.loads(args.html.read_text(encoding="utf-8"))
            assert payload["meta"]["geometry_id"]==pilot.GEOMETRY_ID
            assert payload["groups"]["jpsi40-bore|0|-1|u|-1"]["n"]==1
            assert payload["groups"]["omega40-bore|5|-1|u|-1"]["n"]==0
            assert sum(payload["groups"]["phi40-bore|0|-1|u|-1"]["time"])>0
            print("PASS: full tiny native Geant4 first-chunk aggregation end-to-end")
            print("PASS: 18 source×station totals incl empty stations, nonmuon filter")
            print("PASS: signed angle/time PDF, per-hit physical-wire CSV and HTML payload")
    finally:
        pilot.load_centered_map=old_loader
        pilot.read_heed=old_heed
        pilot.FIRST_CHUNK_EXPECTED=old_expected


if __name__=="__main__":
    main()

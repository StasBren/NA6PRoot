#!/usr/bin/env python3
"""Six-station synthetic regression: crossing P(A) dashboard, SHA + conservation."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

from pi_crossing_asymmetry_dashboard import assemble, digest, render
from pi_lut500_muon_crossing_asymmetry import Moments
from pi_pitch_fieldwire_pilot import tag


def fixture(root):
    configs=[tag(p,fw) for p in (3.,4.) for fw in (0.,75.)]
    stations={}
    layout=[]
    for st in range(6):
        obj=SimpleNamespace(
            station=st,det_id=st+100,x=float(st),y=0.,
            half_x=.45,half_y=.75,layer="A")
        layout.append([obj])
        stat={}
        for cfg in configs:
            m=Moments()
            m.add({"up":6,"down":4,"hemisphere_boundary":0,
                   "top":3,"bottom":3,"side":4,
                   "A_hemisphere":.2,"A_top_bottom":0.,
                   "f_top_bottom":.6},10,10)
            m.add({"up":0,"down":0,"hemisphere_boundary":0,
                   "top":0,"bottom":0,"side":0,
                   "A_hemisphere":None,"A_top_bottom":None,
                   "f_top_bottom":None},0,0)
            stat[cfg]=m.result()
        stations[str(st)]={
            "gas_crossings":2, "configs":stat,
            "det":{str(st+100):stat},
            "paired":{},"angle":{}}
    original={
      "schema":"PI_JPSI_NATIVE_MUON_GAS_CROSSING_ASYM_P3_P4_LUT20_V1",
      "source":"jpsi40-bore","lut_trials_per_cell":20,
      "pitches_mm":[3.,4.],"stations_selected":list(range(6)),
      "source_chunks":1,"generated_parent_slots":50000,
      "partial_exposure":True,
      "LUT20":{cfg:{"csv_sha256":"synthetic"} for cfg in configs},
      "definitions":{},"stations":stations
    }
    root.mkdir()
    (root/"summary.json").write_text(json.dumps(original)+"\n")
    (root/"crossings.csv.gz").write_bytes(b"synthetic raw placeholder")
    manifest={
        "schema":original["schema"],"complete_six_station":False,
        "generated_parents":50000,
        "garfield_lut_sha256":original["LUT20"],
        "sha256_outputs":{
            name:digest(root/name) for name in ("summary.json","crossings.csv.gz")
        }
    }
    (root/"manifest.json").write_text(json.dumps(manifest)+"\n")
    return layout


def run_test():
    with tempfile.TemporaryDirectory(prefix="crossing_dash_test_") as tmp:
        root=Path(tmp)/"physics"
        geom=fixture(root)
        frozen={f:digest(root/f) for f in
                ("summary.json","crossings.csv.gz","manifest.json")}
        payload,checks=assemble(root,verify_raw=True,geometry=geom)
        assert payload["all_crossings"]==12
        assert payload["generated_parent_slots"]==50000
        assert payload["stations"]["0"]["configs"][tag(3.,0.)][
            "hemisphere_A"]["valid"]==1
        assert checks["summary.json"]==frozen["summary.json"]
        path=render(root,geometry=geom)
        page=path.read_text()
        assert "ONE physical muon gas crossing" in page
        assert "parent-event dashboard" in page
        assert "__PAYLOAD__" not in page
        assert "LUT20" in page and "NA60+ Simulation" in page
        assert (root/"crossing_dashboard.provenance.json").is_file()
        assert render(root,geometry=geom)==path # idempotent
        assert frozen=={f:digest(root/f) for f in frozen}
        try:
            path.write_text("accidental unrelated preexisting dashboard")
            render(root,geometry=geom)
        except FileExistsError:
            pass
        else:
            raise AssertionError("Unapproved dashboard replacement")
        render(root,replace=True,geometry=geom)
        assert "ONE physical muon gas crossing" in path.read_text()
        (root/"summary.json").write_text((root/"summary.json").read_text()+" ")
        try:
            render(root,geometry=geom)
        except ValueError as e:
            assert "SHA" in str(e)
        else:
            raise AssertionError("Physics checksum corruption went unnoticed")
    print("PASS: six-station crossing-only dashboard; hist and DET "
          "conservation, SHA protection, no parent pooling, no physics overwrite")


if __name__=="__main__":
    run_test()

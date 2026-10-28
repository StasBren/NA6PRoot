#!/usr/bin/env python3
"""Offline B1 accounting regression; NO Garfield or Geant4 executed."""
from __future__ import annotations
import csv
import math
import tempfile
from pathlib import Path

from pi_prepare_muon_avalanche_b1 import validate, FIELDS
from pi_run_native_muon_b1 import GEOMS, CROSSING_COLUMNS, check_crossings
from pi_benchmark_native_muon_b1 import summarize, percent

def crossing():
    r={k:"0" for k in FIELDS}
    r.update(source="jpsi40-bore",chunk="1",event="4",hit_index="9",
             track_id="2",station="0",det_id="17",pdg="13",
             u_in_mm="0",v_in_mm="3",w_in_mm="0",
             u_out_mm="0",v_out_mm="-3",w_out_mm="0",
             theta_u_in_deg="180",theta_w_in_deg="180",
             path_length_mm="6",edep_eV="800",
             pu_in_GeV="0",pv_in_GeV="-12",pw_in_GeV="0",
             momentum_GeV="12",heed_seed="12345")
    return r

def row(r,name):
    h={c:"" for c in CROSSING_COLUMNS}
    h.update({
     CROSSING_COLUMNS[0]:"|".join(r[k] for k in
       ("source","chunk","event","hit_index","track_id","station","det_id")),
     "configuration":name,"station":r["station"],"det_id":r["det_id"],
     "theta_u_in_deg":r["theta_u_in_deg"],
     "theta_w_in_deg":r["theta_w_in_deg"],
     "path_mm":"6","momentum_GeV":"12","heed_clusters":"3",
     "heed_primaries":"5","heed_discarded":"0","heed_failed":"0",
     "collected_avalanche_electrons":"9","ion_births":"9",
     "garfield_ion_count":"9","zero_seeds":"1","capped_seeds":"0",
     "failed_seeds":"0","ion_up":"5","ion_down":"4","ion_boundary":"0",
     "ion_top":"3","ion_bottom":"2","ion_side":"4",
     "A_ion_hemi":str(1/9),"A_ion_TB":str(1/5),
     "e_up":"4","e_down":"5","e_top":"2","e_bottom":"3",
     "A_e_hemi":str(-1/9),"A_e_TB":str(-1/5),
     "G_eff_per_Heed_primary":"1.8"})
    return h

def run():
    r=crossing()
    res=validate(r)
    assert res is not None and math.isclose(res[0],12)
    bad=dict(r,pv_in_GeV="12")
    assert validate(bad) is None, "Reversed incoming muon was accepted"
    bad=dict(r,u_out_mm="300")
    assert validate(bad) is None, "Muon near finite u-edge was accepted"
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/"crossing.csv"
        # Python csv.DictWriter with newline='' writes CRLF on Linux.
        # C++ ReadCrossings historically failed at exact header comparison
        # because std::getline strips LF but not CR.
        raw=Path(d)/"native_input.csv"
        with raw.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=FIELDS)
            writer.writeheader()
            writer.writerow(r)
        assert raw.read_bytes().split(b"\n")[0].endswith(b"\r")
        good=[row(r,k) for k in GEOMS]
        def write(rows):
            with path.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=CROSSING_COLUMNS)
                w.writeheader();w.writerows(rows)
        write(good)
        stats=check_crossings(path,[r])
        assert all(d["n_crossings"]==1 for d in stats.values())
        assert all(d["n_avalanche_ion_births"]==9 for d in stats.values())
        badrows=[dict(x) for x in good]
        badrows[0]["ion_top"]="4"
        write(badrows)
        try:check_crossings(path,[r])
        except ValueError:pass
        else:raise AssertionError("Failed to reject inconsistent ion-sector sum")
        badrows=[dict(x) for x in good]
        badrows[0]["heed_primaries"]="4"
        write(badrows)
        try:check_crossings(path,[r])
        except ValueError:pass
        else:raise AssertionError("Failed to reject non-matching Heed primaries")
    # Tiny synthetic BENCHMARK_ONLY regression; four geometry samples with
    # identical real Heed primary counts and known wall-clock arithmetic.
    timed=[]
    for geom in GEOMS:
        for j in (0,2):
            timed.append({"crossing_id":"|".join(
                r[k] for k in ("source","chunk","event","hit_index","track_id","station","det_id")),
                "configuration":geom,
                "primary_index":str(j),
                "n_true_heed_primaries":"5",
                "wall_seconds":"0.2" if j==0 else "0.4",
                "avalanche_electrons":"30","avalanche_ions":"20",
                "ion_births":"20","collected_e":"20"})
    stats=summarize(timed,[r])
    assert len(stats)==4
    assert all(abs(d["mean_seconds_per_seed"]-0.3)<1e-10 for d in stats.values())
    assert all(abs(d["rough_full_native_muon_2event_time_seconds"]-1.5)<1e-10
               for d in stats.values())
    assert percent([1.,2.,3.],.5)==2.
    repeated=[dict(x) for x in timed]
    repeated[0]["primary_index"]=repeated[1]["primary_index"]
    try:summarize(repeated,[r])
    except ValueError:pass
    else:raise AssertionError("Timing benchmark duplicate primary was not rejected")
    print("PASS real muon momentum/chord, finite-edge rejection, B1 accounting and 4-geometry timing benchmark regression tests")

if __name__=="__main__":run()

#!/usr/bin/env python3
"""Synthetic unit tests, no Garfield/Geant4/native files needed."""
import csv
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from pi_pitch_fieldwire_paired_analysis import (
    analyse, match_maps, metrics, spatial_profiles,
)
from pi_pitch_fieldwire_paired_html import render
from pi_pitch_fieldwire_pilot import tag, NV, DU


def test_paired_metrics():
    bare=np.array([[1,1,0],[1,0,1]],dtype=bool)
    field=np.array([[1,0,1],[1,0,0]],dtype=bool)
    ta=np.array([[10.,20.,np.nan],[30.,np.nan,50.]])
    tb=np.array([[15.,np.nan,40.],[35.,np.nan,np.nan]])
    aa=np.where(bare,20.,np.nan)
    ab=np.where(field,25.,np.nan)
    m=metrics(bare,field,ta,tb,aa,ab)
    assert m["n_both_collected"] == 2
    assert m["n_paired_valid_times"] == 2
    assert abs(m["paired_mean_time_delta_both_ns"]-5.)<1.e-12
    assert abs(m["mean_time_no_ns"]-27.5)<1.e-12
    assert abs(m["mean_time_field_ns"]-30.)<1.e-12
    assert abs(m["delta_conditional_mean_time_ns"]-2.5)<1.e-12
    assert abs(m["delta_anode_fraction_pp"]+100/6)<1.e-12
    pair={
        "pitch":3., "samples":3, "bare":bare, "field":field,
        "bare_time":ta, "field_time":tb,
        "bare_angle":aa, "field_angle":ab,
        "transition_counts":{"0->0":1,"0->1":1,"1->0":2,"1->1":2},
        "status_bare":{"0":6}, "status_field":{"0":6},
        "field_wire_ends":0,"input_sha256":{},
    }
    report=analyse(pair,replicas=60)
    assert report["n_electron_pairs"] == 6
    assert math.isfinite(
        report["bootstrap_standard_errors"]["paired_mean_time_delta_both_ns"])
    assert report["point_estimate"]["paired_mean_time_delta_both_ns"] == 5.


def test_full_map_matching(tmp):
    config={
        "pitch_mm":3., "field_wire_diam_um":0.,
        "gap_minus_mm":3., "gap_plus_mm":3.,
        "anode_voltage_v":1800., "anode_diameter_um":30.,
        "gas":"ArCO2","magnetic_field_t":0., "n_u":30,"n_v":30,
        "du_mm":DU,"dv_mm":.2, "samples_per_box":2,
        "seed":260020, "cpp_sha256":"same", "binary_sha256":"same",
        "mode":"pilot",
    }
    schema=("iu","iv","box_id","sample","collected",
            "field_wire_collected","drift_time_ns","incidence_normal_deg",
            "status","u_start_mm","v_start_mm","random_seed","wire_index")
    for fw in (0.,75.):
        dest=tmp/tag(3.,fw)
        dest.mkdir()
        (dest/"frozen_config.json").write_text(json.dumps(
            {**config,"field_wire_diam_um":fw}))
        with (dest/"garfield_electrons.csv").open(
                "w",newline="",encoding="utf-8") as file:
            w=csv.DictWriter(file,fieldnames=schema)
            w.writeheader()
            for iu in range(30):
                for iv in range(NV):
                    for k in range(2):
                        ib=iu*NV+iv
                        # Identical physical starts and seeds in both setups.
                        c=(k==0 or (fw==0 and iv%2==0))
                        # The original C++ safeguards electron births inside
                        # a grounded field wire: only that geometry moves the
                        # point to the box center. Exercise this actual case.
                        relocated = (iu == 0 and iv == 14 and k == 0)
                        start_u = (-1.49 if relocated and fw == 0.
                                   else -1.5+(iu+.5)*DU)
                        start_v = (-.01 if relocated and fw == 0.
                                   else -3+(iv+.5)*.2)
                        w.writerow({
                            "iu":iu, "iv":iv, "box_id":ib,
                            "sample":k, "collected":int(c),
                            "field_wire_collected":0,
                            "drift_time_ns":15+iu if c else "nan",
                            "incidence_normal_deg":25 if c else "nan",
                            "status":-1 if not c else -5,
                            "u_start_mm":start_u,
                            "v_start_mm":start_v,
                            "random_seed":260020+ib*2+k,
                            "wire_index":0,
                        })
    pair=match_maps(3.,tmp)
    assert pair["bare"].shape == (900,2)
    assert pair["field"].shape == (900,2)
    assert pair["transition_counts"]["1->0"] == 450
    assert pair["transition_counts"]["1->1"] == 899
    assert pair["transition_counts"]["0->0"] == 450
    assert pair["transition_counts"]["0->1"] == 0
    assert len(pair["start_relocations"]) == 1
    assert pair["start_relocations"][0]["box_id"] == 14
    assert pair["start_relocations"][0]["sample"] == 0
    assert np.count_nonzero(pair["paired_start_valid"]) == 1799
    profile=spatial_profiles(pair)
    assert len(profile)==30
    assert profile[0]["u_bin"]==0
    assert abs(profile[0]["delta_anode_fraction_pp"]+1500/59)<1.e-12
    # Offline HTML is pure presentation; original simulation arrays untouched.
    analysis=analyse(pair,replicas=30)
    mock_report={"field_wire_diameter_um":75.,
                 "groups":[{**analysis,"pitch_mm":p}
                           for p in (3.,3.5,4.)]}
    sanitized=[{key:("nan" if isinstance(value,float) and
                    not math.isfinite(value) else value)
                for key,value in row.items()} for row in profile]
    page=render(mock_report,sanitized,tmp/"paired.html")
    contents=page.read_text(encoding="utf-8")
    assert "NA60+ MWPC" in contents
    assert "application/json" in contents
    assert "Same initial electron" in contents
    assert "75.0 µm" in contents
    # Unexplained drifts in the birth coordinate must be rejected, not
    # casually hidden as if they were valid matched electron trajectories.
    f=tmp/tag(3.,75.)/"garfield_electrons.csv"
    original=f.read_text(encoding="utf-8")
    lines=original.splitlines()
    header=lines[0].split(",")
    index=header.index("u_start_mm")
    invalid=lines[1].split(",")
    invalid[index]=str(float(invalid[index])+.04)
    lines[1]=",".join(invalid)
    f.write_text("\n".join(lines)+"\n",encoding="utf-8")
    try:
        match_maps(3.,tmp)
    except ValueError as exc:
        assert "Unexplained Garfield start mismatch" in str(exc)
    else:
        raise AssertionError("Unexplained start-coordinate mismatch accepted")


def main():
    test_paired_metrics()
    with TemporaryDirectory(prefix="mwpc_paired_") as t:
        test_full_map_matching(Path(t))
    print("PASS: true paired electron starts, collection transitions, "
          "common-collected time delta, stratified errors and u profiles")


if __name__=="__main__":
    main()

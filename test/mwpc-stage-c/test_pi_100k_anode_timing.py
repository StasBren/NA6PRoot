#!/usr/bin/env python3
"""Simulation-free regression: 18 small, conserved 100k station JSON files."""
from __future__ import annotations
import json
import math
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from analyze_pi_100k_anode_timing import (
    CONFIGS,SOURCES,analyze_file,sum_across_stations,write_reports,stats
)
from pi_six_geometry_native_timing import n_anodes
from pi_pitch_fieldwire_pilot import tag


def group(source,st,cfg,det,wire,times,primaries,segments):
    hist=[0]*31
    for t in times:hist[min(30,int(t//10))]+=1
    return dict(source=source,station=st,config=cfg,
                det_id=det,anode=wire,axis="all",anglebin=-1,
                segments=segments,primaries=primaries,collected=len(times),
                known=len(times),unknown=0,edge=0,
                hist=hist,angle_hist=[0]*18,phi_hist=[0]*36,
                angle_n=0,phi_n=0,
                sum_drift_ns=float(sum(times)),
                sum_drift_ns2=float(sum(t*t for t in times)))


def write_mock(root,channel,st):
    # Production merged JSON directories use "jpsi"/"phi"/"omega",
    # while the groups and meta retain "jpsi40-bore" etc.
    source=SOURCES[channel]
    rows=[]
    for pitch in (3.,3.5,4.):
        for field in (0.,75.):
            cfg=tag(pitch,field)
            n=n_anodes(st,pitch)
            # A complete station with one chamber, same original exposure
            # but deliberately distinct central and outer timing.
            center=[25.,25.,25.]
            outer=[80.,80.]
            all_times=center+outer
            rows.extend([
                group(source,st,cfg,-1,-1,all_times,8,7),
                group(source,st,cfg,0,-1,all_times,8,7),
                group(source,st,cfg,0,n//2,center,0,1),
                group(source,st,cfg,0,0,outer,0,1),
            ])
    file=root/channel.lower()/f"MS{st}.json"
    file.parent.mkdir(parents=True,exist_ok=True)
    meta={"n_generated_parent_events_per_source":100000,
          "geometries":6,"garfield_luts":{cfg:{} for cfg in CONFIGS},
          "source_station":f"{source} MS{st}",
          "n_muon_gas_crossings_station":7}
    file.write_text(json.dumps(
        {"meta":meta,"group_count":len(rows),"groups":rows},
        ensure_ascii=False,separators=(",",":")),encoding="utf-8")
    return file


def main():
    with tempfile.TemporaryDirectory(prefix="pi_100k_anode_synthetic_") as folder:
        root=Path(folder)/"merged"
        records={}
        for key,physical in SOURCES.items():
            for st in range(6):
                actual=write_mock(root,key,st)
                assert actual==root/key.lower()/f"MS{st}.json"
                records[(key,st)]=analyze_file(key,st,root)
                c=records[(key,st)]["p3p0_anode_only"]
                assert c["segments"]==7
                assert c["primary_electrons"]==8
                assert c["collected_electrons"]==5
                assert c["known_times"]==5
                assert math.isclose(c["mean_ns"],47,abs_tol=1e-12)
                assert math.isclose(c["collection_fraction"],5/8)
                assert math.isclose(c["fast_collected_fraction_all_primaries"],5/8)
                assert math.isclose(c["outer_minus_center_mean_ns"],55.)
                assert c["center_known_electrons"]==3
                assert c["outer_known_electrons"]==2
                assert c["per_position_bin"][0]["known_times"]==2
        totals=sum_across_stations(records)
        a=totals[("Jpsi","p3p0_anode_only")]
        assert a["segments"]==42 and a["known_times"]==30
        assert math.isclose(a["mean_ns"],47)
        assert math.isclose(a["centered_rms_ns"],math.sqrt(726),abs_tol=1e-10)
        out=Path(folder)/"analysis"
        with redirect_stdout(StringIO()):
            write_reports(records,totals,out)
        page=(out/"anode_timing_report_en.html").read_text()
        assert "id=\"anode\"" in page and "id=\"fast\"" in page
        assert "id=\"time\"" in page and "100,000 generated parents" in page
        assert (out/"station_summary.csv").is_file()
        assert (out/"summary.txt").is_file()
        # Corrupt one wire's time histogram: strict station-vs-anode
        # conservation must reject it, not silently smooth missing hits.
        bad=root/"jpsi"/"MS0.json"
        data=json.loads(bad.read_text())
        wire=next(r for r in data["groups"] if r["anode"]>=0)
        wire["hist"][2]+=1
        bad.write_text(json.dumps(data,separators=(",",":")))
        try:
            analyze_file("Jpsi",0,root)
        except (AssertionError,ValueError):
            pass
        else:
            raise AssertionError("Corrupted per-anode time spectrum passed QA")
    print("PASS: streamed 18 synthetic source/station native 100k manifests; "
          "correct 31-bin conservation, true electron-weighted mean/RMS, "
          "outer-minus-center anode timing, six-field English HTML and "
          "corruption rejection; NO real Geant4/Garfield jobs.")


if __name__=="__main__":
    main()

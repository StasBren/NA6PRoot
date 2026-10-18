#!/usr/bin/env python3
"""Synthetic fixture for the event-aligned, geometry-only Muon Track Explorer.

No historical input/Geant4/Garfield or charge-map CPU required.
"""
import csv
import json
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory

from mwpc_checkerboard_projection import SOURCES, build_chambers
from build_pi_muon_track_explorer import build, nominal_geometry_anodes


def main():
    groups=build_chambers()
    with TemporaryDirectory(prefix="pi_track_fixture_") as tmp:
        root=Path(tmp)
        inp=root/"crossings"
        inp.mkdir()
        fields=[
            "source","job","event","track","pdg","station","det_id",
            "layer","x_in_cm","y_in_cm","z_in_cm",
            "x_out_cm","y_out_cm","z_out_cm","chamber_center_y_cm",
            "gas_center_z_cm","theta_u_deg","theta_w_deg",
            "gas_path_mm","p_GeV",
        ]
        for src in SOURCES:
            with (inp/f"crossings_{src}_ms2_530.csv").open(
                "w",newline="",encoding="utf-8"
            ) as stream:
                writer=csv.DictWriter(stream,fieldnames=fields)
                writer.writeheader()
                for station in range(6):
                    ch=groups[station][0]
                    z=(300,340,530,590,810,850)[station]+ch.z_offset
                    for track,pdg in ((11,13),(12,-13)):
                        writer.writerow({
                            "source":src,"job":"subrun01","event":"42",
                            "track":track,"pdg":pdg,"station":station,
                            "det_id":ch.det_id,"layer":ch.layer,
                            "x_in_cm":ch.x,"y_in_cm":ch.y,
                            "z_in_cm":z-.3,"x_out_cm":ch.x,
                            "y_out_cm":ch.y,"z_out_cm":z+.3,
                            "chamber_center_y_cm":ch.y,
                            "gas_center_z_cm":z,
                            "theta_u_deg":0.,"theta_w_deg":0.,
                            "gas_path_mm":6.,"p_GeV":20.,
                        })
        out=root/"track.html"
        build(Namespace(
            crossings=inp,geometry_only=True,sample_modulo=1,
            max_events_per_source=5,
            template=Path(__file__).with_name(
                "pi_muon_track_explorer_template.html"
            ),
            out=out,
        ))
        html=out.read_text(encoding="utf-8")
        assert "__PI_TRACK_DATA__" not in html
        marker='<script id="trackData" type="application/json">'
        data=json.loads(html.split(marker,1)[1].split("</script>",1)[0])
        assert len(data["sources"])==3
        assert not data["meta"]["chargeModel"]
        for src in SOURCES:
            events=data["sources"][src]
            assert len(events)==1
            event=events[0]
            assert event["job"]=="subrun01" and event["event"]=="42"
            assert len(event["tracks"])==2
            assert {x["pdg"] for x in event["tracks"]}=={13,-13}
            for track in event["tracks"]:
                assert track["stations"]==[0,1,2,3,4,5]
                assert len(track["points"])==6
                for point in track["points"]:
                    assert len(point["geom"])==1
                    assert point["charges"] is None
                    assert point["length"]==6.
        sample={
            "station":0,"y_in_cm":0.,"y_out_cm":0.,
            "chamber_center_y_cm":0.,
        }
        assert nominal_geometry_anodes(sample)==[38]
        sample["station"]=1
        assert nominal_geometry_anodes(sample)==[64]
        print("PASS: complete six-station same-event mu+/mu- association")
        print("PASS: finite 77/128-wire central cell selection and source grouping")
        print("PASS: offline geometry-only HTML JSON embedded without Geant4")
if __name__=="__main__":
    main()

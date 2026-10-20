#!/usr/bin/env python3
"""Synthetic full-coverage histogram fixture for signed-angle map builder."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_native_full_campaign_qa import ANGLE
from pi_native_angle_position import FIELDS, build


def main():
    with TemporaryDirectory(prefix="native_pi_angle_map_") as t:
        root=Path(t)
        (root/"pi_native_all_QA_report.json").write_text(
            json.dumps({"full_campaign_analyzed":True,
                        "processed_chunks_per_source":67}),
            encoding="utf-8")
        chambers=[c for group in build_chambers() for c in group]
        with (root/"pi_native_all_chamber_quality.csv").open("w",newline="") as stream:
            w=csv.DictWriter(stream,fieldnames=("source","station","det_id","gas_hit_rows"))
            w.writeheader()
            for src in SOURCES:
                for c in chambers:
                    w.writerow(dict(source=src,station=c.station,
                                    det_id=c.det_id,gas_hit_rows=3))
        with (root/"pi_native_all_station_quality.csv").open("w",newline="") as stream:
            w=csv.DictWriter(stream,fieldnames=("source","station","gas_hit_rows"))
            w.writeheader()
            for src in SOURCES:
                for st,group in enumerate(build_chambers()):
                    w.writerow(dict(source=src,station=st,gas_hit_rows=3*len(group)))
        fields=("source","det_id","variable","bin_lower","bin_upper",
                "hit_count","gas_hit_exposure")
        with (root/"pi_native_all_histograms.csv").open("w",newline="") as stream:
            w=csv.DictWriter(stream,fieldnames=fields)
            w.writeheader()
            for src in SOURCES:
                for c in chambers:
                    for field in FIELDS:
                        for i in range(len(ANGLE)-1):
                            w.writerow(dict(
                                source=src,det_id=c.det_id,variable=field,
                                bin_lower=ANGLE[i],bin_upper=ANGLE[i+1],
                                hit_count=3 if ANGLE[i]==0 else 0,
                                gas_hit_exposure=3,
                            ))
        out=root/"angle.html"
        template=Path(__file__).with_name("pi_native_angle_position_template.html")
        build(root,out,template)
        html=out.read_text(encoding="utf-8")
        marker='<script id="angleData" type="application/json">'
        payload=json.loads(html.split(marker,1)[1].split("</script>",1)[0])
        assert len(payload["geometry"])==298
        assert len(payload["module"])==3
        assert len(payload["module"][SOURCES[0]])==298
        assert payload["module"][SOURCES[0]]["0"]["n"]==3
        assert sum(payload["module"][SOURCES[1]]["278"]["angles"]["theta_w_in_deg"])==3
        assert len(payload["module"][SOURCES[2]]["34"]["angles"]["theta_u_out_deg"])==len(ANGLE)-1
        # Ensure the same page truly embeds local u-w data while its
        # existing global angle histograms retain full QA exposure.
        local_path=root/"synthetic_local_uw.json"
        local_data={
            "meta":{"scope":"FIRST_CHUNK_3000_EVENTS_PER_SOURCE",
                    "nu":16,"nw":24,"full_independent_qa_passed":False},
            "modules":{
                src:{str(c.det_id):{
                    "n":3,"outside":0,
                    "all":[[0,1],[1,2]],
                    "selected":{field:{
                        "positive":[[0,1]],"far10":[],
                        "near2":[[0,1],[1,2]]
                        } for field in FIELDS}
                } for c in chambers}
                for src in SOURCES},
        }
        local_path.write_text(json.dumps(local_data),encoding="utf-8")
        build(root,out,template,local_path)
        joined=json.loads(out.read_text(encoding="utf-8").split(
            marker,1)[1].split("</script>",1)[0])
        assert joined["local"]["meta"]["scope"]==(
            "FIRST_CHUNK_3000_EVENTS_PER_SOURCE")
        assert joined["local"]["modules"][SOURCES[0]]["0"]["all"]==[
            [0,1],[1,2]]
        assert joined["module"][SOURCES[0]]["0"]["n"]==3
        assert joined["local"]["modules"][SOURCES[2]]["278"][
            "selected"]["theta_u_in_deg"]["near2"]==[[0,1],[1,2]]
        local_data["meta"]["scope"]="FULL_201_CHUNKS_200000_EVENTS_PER_SOURCE"
        local_path.write_text(json.dumps(local_data),encoding="utf-8")
        try:
            build(root,out,template,local_path)
        except ValueError as exc:
            assert "not independently verified" in str(exc)
        else:
            raise AssertionError("Unvalidated full local-u-w source was accepted")
        local_data["meta"]["full_independent_qa_passed"]=True
        local_data["modules"][SOURCES[0]]["0"]["n"]=4
        local_data["modules"][SOURCES[0]]["0"]["outside"]=1
        local_path.write_text(json.dumps(local_data),encoding="utf-8")
        try:
            build(root,out,template,local_path)
        except ValueError as exc:
            assert "exposure differ" in str(exc)
        else:
            raise AssertionError("Full local-u-w mismatch with full Geant4 QA accepted")
        print("PASS: embedded pilot native u-w maps alongside full native angle QA")
        print("PASS: refuse unvalidated or mismatched full 201-chunk local map")
        print("PASS: all 894 physical module/source PDFs and 4 signed angle fields")
        print("PASS: native complete-campaign provenance, histogram conservation and offline HTML JSON")


if __name__=="__main__":
    main()

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
        print("PASS: all 894 physical module/source PDFs and 4 signed angle fields")
        print("PASS: native complete-campaign provenance, histogram conservation and offline HTML JSON")


if __name__=="__main__":
    main()

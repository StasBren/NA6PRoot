#!/usr/bin/env python3
"""Bounded-memory full-JSON stream parser and offline source-station shard smoke."""
import json
import tempfile
from pathlib import Path
import re

from pi_six_geometry_native_dashboard_pages import (
    read_groups_stream, generate,
)
from pi_six_geometry_native_timing import SOURCES, tag, PITCHES
from pi_six_geometry_native_timing_html import HTML
from mwpc_checkerboard_projection import build_chambers
from pi_six_geometry_layout_v3 import VERSION


def main():
    # UI regression: the first station page used to render empty because
    # draw() referred to el.stats and el.coverage without binding either.
    selectors=json.loads(re.search(r'const ids=(\[[^;]+\]);', HTML).group(1))
    element_refs=set(re.findall(r'\bel\.([A-Za-z_]\w*)', HTML))
    extra={"stats","coverage"}
    assert element_refs.issubset(set(selectors)|extra), (
        "Unregistered native HTML DOM element references",
        sorted(element_refs-set(selectors)-extra))
    for name in extra:
        assert (name+':document.getElementById("'+name+'")') in HTML, (
            f"Missing real DOM binding for {name}")
        assert 'id="'+name+'"' in HTML
    configs=[tag(p,fw) for p in PITCHES for fw in (0.,75.)]
    assert len(configs)==6
    meta={
        "garfield_luts":{k:{"test_only":True} for k in configs},
        "gas_hit_counts":{s:1 for s in SOURCES},
        "angle_bins_deg":[-180,-5,0,5,180],
        "input_model":"synthetic test only",
        "statistical_errors":"unavailable",
        "data_scope":"SYNTHETIC",
        "finite_anodes_by_pitch":{},
    }
    rows=[]
    det_by_st={st:next(c.det_id for c in build_chambers()[st])
               for st in range(6)}
    for src in SOURCES:
        for st in range(6):
            for cfg in configs:
                row={
                    "source":src,"config":cfg,"station":st,
                    "det_id":-1,"anode":-1,"axis":"all","anglebin":-1,
                    "segments":1,"primaries":3,"collected":2,"known":2,
                    "unknown":0,"edge":0,
                    "hist":[2]+[0]*30,
                    "angle_hist":[1]+[0]*17,"phi_hist":[2]+[0]*35,
                    "angle_n":1,"phi_n":2,
                    "mean_ns":5.,"rms_ns":1.,
                    "mean_angle_deg":10.,"rms_angle_deg":0.,
                }
                rows.append(row)
                rows.append({**row,"det_id":det_by_st[st]})
    with tempfile.TemporaryDirectory(prefix="six_native_pages_") as t:
        tmp=Path(t)
        original=tmp/"full.json"
        original.write_text(json.dumps({"meta":meta,"group_count":len(rows),
                                        "groups":rows},
                                       separators=(",",":")),encoding="utf-8")
        full=list(read_groups_stream(original))
        assert full[0]["group_count"]==len(rows)
        assert full[1:]==rows
        idx=generate(original,tmp/"pages")
        assert idx.exists() and "MS5.html" in idx.read_text()
        station_page=tmp/"pages"/list(SOURCES)[0]/"MS0.html"
        assert station_page.is_file()
        page_text=station_page.read_text(encoding="utf-8")
        assert "MWPC-NATIVE-STATIC-FALLBACK-20261021-v2" in page_text
        assert "Резервная гистограмма" in page_text
        assert 'id="staticFallback"' in page_text
        assert '<rect x="' in page_text  # Python-generated physical histogram
        assert "document.getElementById(\"staticFallback\").hidden=true" in page_text
        compact=json.loads((tmp/"pages"/
                            "pi_six_geometry_station_reference.json").read_text())
        assert len(compact["groups"])==len(rows)//2
        assert compact["groups"][0]["mean_ns"]==5.
        assert (tmp/"pages"/list(SOURCES)[-1]/"MS5.html").is_file()
        new_portal=tmp/"pages"/"pi_pitch_fieldwire_timing_dashboard_v3.html"
        assert new_portal.is_file()
        assert "XY-карта" in new_portal.read_text(encoding="utf-8")
        v3page=tmp/"pages"/list(SOURCES)[0]/"MS0_layout_v3.html"
        assert v3page.is_file()
        txt=v3page.read_text(encoding="utf-8")
        assert VERSION in txt
        assert 'id="map"' in txt and 'id="time"' in txt
        assert 'id="arrival"' in txt and 'id="comparison"' in txt
        assert "P(N≥2)" in txt and "Это не P(N≥2)" in txt
        assert "Резервная гистограмма" in txt
        assert "__STATIC_FALLBACK__" not in txt
    print("PASS: stream native JSON to 18 XY/pitch/field-wire PI V3 "
          "plus original 18 reports and exact bootstrap station reference")


if __name__=="__main__":
    main()

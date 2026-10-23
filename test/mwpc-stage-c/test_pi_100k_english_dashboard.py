#!/usr/bin/env python3
"""Read-only smoke test: English localisation is lossless for physics JSON."""
from __future__ import annotations
import json
import tempfile
from pathlib import Path

from build_pi_100k_english_dashboard import translate_page, CYRILLIC
from pi_six_geometry_layout_v3 import TEMPLATE, make_portal
from pi_six_geometry_native_timing_html import build as build_detail
from pi_six_geometry_native_timing import PITCHES
from pi_pitch_fieldwire_pilot import tag

CONFIGS={tag(p,d):{} for p in PITCHES for d in (0.,75.)}
REFERENCE={
    "source":"jpsi40-bore","config":"p3p0_anode_only","station":0,
    "det_id":-1,"anode":-1,"axis":"all","anglebin":-1,
    "segments":33,"primaries":104,"collected":80,
    "known":80,"unknown":0,"edge":0,
    "mean_ns":25.,"rms_ns":11.,"hist":[0,0,80]+[0]*28,
    "angle_n":80,"angle_hist":[80]+[0]*17,
    "phi_n":80,"phi_hist":[80]+[0]*35,
}
META={
    "model_version":"synthetic-test-only",
    "n_generated_parent_events_per_source":100000,
    "data_scope":"SYNTHETIC TEST ONLY; not real data",
    "garfield_luts":CONFIGS,
    "gas_hit_counts":{"jpsi40-bore":33},
    "input_model":"synthetic",
    "statistical_errors":"none",
    "finite_anodes_by_pitch":{},
    "angle_bins_deg":[-180,0,180],
}


def assert_payload_unchanged(old,new):
    from build_pi_100k_english_dashboard import DATA_TAG
    before=DATA_TAG.search(old).group(2)
    after=DATA_TAG.search(new).group(2)
    assert before==after
    assert json.loads(after)["meta"]["n_generated_parent_events_per_source"]==100000


def main():
    with tempfile.TemporaryDirectory(prefix="mwpc_v3_english_synthetic_") as t:
        base=Path(t)
        directory=base/"jpsi40-bore"
        directory.mkdir()
        payload=dict(meta=META,groups=[REFERENCE],group_count=1)
        # Real V3 template with real dynamic JS; use a synthetic full-fold
        # JSON placeholder and the real first-chunk -> 100k displayed edits.
        v3=TEMPLATE.read_text(encoding="utf-8")
        v3=v3.replace("__PAYLOAD__",json.dumps(payload,ensure_ascii=False))
        v3=v3.replace("__STATIC_FALLBACK__","<p>Fallback test</p>")
        v3=v3.replace(
            "Трёхтысячный первый чанк J/ψ / ϕ / ω на источник.",
            "100,000 сгенерированных родительских событий на источник; "
            "объединены проверенные независимые Geant4-чанки.")
        v3=v3.replace(
            "истинные native Geant4 газовые пересечения первого чанка "
            "(3 000 сгенерированных событий на канал)",
            "истинные native Geant4 газовые пересечения из независимых "
            "чанков (по 100,000 сгенерированных событий на канал)")
        v3=v3.replace(
            "18 детальных страниц и исходные 496 809 групп не изменяются.",
            "Оригинальные отчёты на 3 000 событий сохранены отдельно.")
        # Add exactly the same source navigation expected from production.
        ru=directory/"MS0_layout_v3.html"
        en=directory/"MS0_layout_v3_en.html"
        ru.write_text(v3,encoding="utf-8")
        translate_page(ru,en,"layout")
        text=en.read_text(encoding="utf-8")
        assert not CYRILLIC.search(text)
        assert_payload_unchanged(v3,text)
        assert "../pi_pitch_fieldwire_timing_dashboard_v3_en.html" in text
        assert 'MS"+S.st+"_layout_v3_en.html"' in text
        assert "id=\"timeScale\"" in text
        assert "id=\"arrival\"" in text
        # Render real detailed page HTML + real Python static fallback,
        # using small but physically conserved synthetic histogram rows.
        detail_json=base/"mock_station.json"
        detail_json.write_text(json.dumps(payload),encoding="utf-8")
        detail_ru=directory/"MS0.html"
        detail_en=directory/"MS0_en.html"
        build_detail(detail_json,detail_ru)
        translate_page(detail_ru,detail_en,"detail")
        text2=detail_en.read_text(encoding="utf-8")
        assert not CYRILLIC.search(text2)
        assert_payload_unchanged(detail_ru.read_text(encoding="utf-8"),text2)
        assert 'id="yscale"' in text2
        # Original portal source with genuine 100k labels/relative links.
        make_portal(base,[("jpsi40-bore",0)],
                    n_generated_parent_events_per_source=100000)
        portal_ru=base/"pi_pitch_fieldwire_timing_dashboard_v3.html"
        portal_en=base/"pi_pitch_fieldwire_timing_dashboard_v3_en.html"
        translate_page(portal_ru,portal_en,"portal")
        portal=portal_en.read_text(encoding="utf-8")
        assert not CYRILLIC.search(portal)
        assert "jpsi40-bore/MS0_layout_v3_en.html" in portal
        assert 'lang="en"' in portal
        # Wrong exposure MUST be rejected, never silently labelled 100k.
        bad=v3.replace(
            '"n_generated_parent_events_per_source": 100000',
            '"n_generated_parent_events_per_source": 3000')
        ru.write_text(bad,encoding="utf-8")
        try:
            translate_page(ru,en,"layout",write=False)
        except ValueError as error:
            assert "parent exposure" in str(error)
        else:
            raise AssertionError("English 100k builder accepted 3000 parents")
    print("PASS: English MWPC V3 interactive map, full per-anode timing "
          "page, static Python fallback and navigation portal; "
          "numerical JSON preserved byte-for-byte; no Cyrillic "
          "in all three HTML types; false exposure rejected")


if __name__=="__main__":
    main()

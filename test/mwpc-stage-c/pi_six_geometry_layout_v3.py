#!/usr/bin/env python3
"""Make a PI V2-style clickable XY-map/timing dashboard from existing six LUTs.

Reads only the source/station JSON slice already created by the v1 pages
streamer. Does NOT re-run Geant4, Heed, Garfield or alter original results.
The data DO NOT contain V2 P(N>=2) gas-hit multiplicities or 2D electron
birth coordinate/time response for SIX configs. The V3 XY heatmap is
explicitly labelled P(primary electron collected with drift < 100ns)
and the lower panels show the TWO *actually stored* arrival-angle PDFs
and six-configuration mean/RMS comparisons.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

from mwpc_checkerboard_projection import build_chambers
from pi_six_geometry_native_timing_html import static_native_spectrum

HERE=Path(__file__).resolve().parent
TEMPLATE=HERE/"pi_six_geometry_layout_v3_template.html"
VERSION="NA60-SIX-GEOMETRY-XY-V3-MULTIDET-20261025"


def construct_payload(blob,source,station):
    """Validate frozen full-source/station groups and add physical 298-det map."""
    if not 0<=station<6:
        raise ValueError("Bad native station")
    meta=blob["meta"]
    confs=list(meta.get("garfield_luts",{}))
    if len(confs)!=6:
        raise ValueError("V3 expects exactly six frozen Garfield configurations")
    for pitch,wire in ((3.,0.),(3.,75.),(3.5,0.),(3.5,75.),(4.,0.),(4.,75.)):
        name=("p"+str(pitch).replace(".","p")+"_"+(
            "anode_only" if wire==0 else "field75p0um"))
        if name not in confs:
            raise ValueError(f"V3 missing pinned pitch/field-wire map {name}")
    data=blob["groups"]
    if not data or any(g["source"]!=source or g["station"]!=station for g in data):
        raise ValueError("Non-station, empty or mixed source data")
    group_by_key={}
    for row in data:
        key=(row["config"],row["det_id"],row["anode"],
             row["axis"],row["anglebin"])
        if key in group_by_key:
            raise ValueError("Duplicate full-fold group in V3")
        group_by_key[key]=row
        if (len(row["hist"])!=31 or
            sum(row["hist"])!=row["known"] or
            row["collected"]!=row["known"]+row["unknown"] or
            len(row["angle_hist"])!=18 or
            sum(row["angle_hist"])!=row["angle_n"] or
            len(row["phi_hist"])!=36 or
            sum(row["phi_hist"])!=row["phi_n"]):
            raise ValueError("Native time/electron angle conservation failed")
    for cfg in confs:
        summary=group_by_key.get((cfg,-1,-1,"all",-1))
        if summary is None:
            raise ValueError("Missing V3 source/station/config full summary")
        # For time PDF on a station, distinct module rows must reproduce
        # the station's 31-bin native microscopic electron histogram exactly.
        modules=[g for g in data if g["config"]==cfg and
                 g["det_id"]>=0 and g["anode"]==-1 and
                 g["axis"]=="all" and g["anglebin"]==-1]
        if sum(g["known"] for g in modules)!=summary["known"]:
            raise ValueError(f"V3 chamber/station electron mismatch: {cfg}")
        for i in range(31):
            if sum(g["hist"][i] for g in modules)!=summary["hist"][i]:
                raise ValueError(f"V3 time bin {i} module mismatch: {cfg}")
    layout=[
        {"id":c.det_id,"st":c.station,"x":c.x,"y":c.y,
         "hx":c.half_x,"hy":c.half_y,"row":c.row,
         "col":c.col,"layer":c.layer}
        for group in build_chambers() for c in group if c.station==station
    ]
    if len(layout)!=len({c["id"] for c in layout}) or len(layout)==0:
        raise AssertionError("Duplicate/empty physical XY DET coordinates")
    payload={
        "meta":meta,"sourceNames":list(meta["gas_hit_counts"]),
        "selected_source":source,"selected_station":station,
        "layout":layout,"groups":data,"dashboard_version":VERSION,
        "scope":(
            f"{int(meta.get('n_generated_parent_events_per_source',3000)):,} "
            "GENERATED Geant4 dimuon parents/source; B=0 transport lookup, "
            "six p/field-wire configurations."
        ),
    }
    return payload


def render(source_file,out,source,station,template=TEMPLATE):
    blob=json.loads(Path(source_file).read_text(encoding="utf-8"))
    data=construct_payload(blob,source,station)
    template_text=Path(template).read_text(encoding="utf-8")
    if template_text.count("__PAYLOAD__")!=1 or \
            template_text.count("__STATIC_FALLBACK__")!=1:
        raise ValueError("Incomplete V3 HTML template placeholders")
    fallback=static_native_spectrum({
        "meta":data["meta"],"groups":data["groups"]})
    page=template_text.replace("__PAYLOAD__",json.dumps(
        data,ensure_ascii=False,separators=(",",":"),allow_nan=False
    ).replace("</","<\\/")).replace("__STATIC_FALLBACK__",fallback)
    n_parents=int(data["meta"].get("n_generated_parent_events_per_source",3000))
    if n_parents!=3000:
        if n_parents<3000:
            raise ValueError("Unsupported native V3 parent exposure")
        page=page.replace(
            "Трёхтысячный первый чанк J/ψ / ϕ / ω на источник.",
            f"{n_parents:,} сгенерированных родительских событий на "
            "источник; объединены проверенные независимые Geant4-чанки."
        ).replace(
            "истинные native Geant4 газовые пересечения первого чанка "
            "(3 000 сгенерированных событий на канал)",
            "истинные native Geant4 газовые пересечения из независимых "
            f"чанков (по {n_parents:,} сгенерированных событий на канал)"
        ).replace(
            "18 детальных страниц и исходные 496 809 групп не изменяются.",
            "Оригинальные отчёты на 3 000 событий сохранены отдельно."
        )
        if "Трёхтысячный первый чанк" in page:
            raise AssertionError("High-stat V3 page still labelled 3000")
    
    dest=Path(out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(page,encoding="utf-8")
    print(f"V3 clickable XY + 6 fields {source} MS{station}: "
          f"{dest.stat().st_size/1048576:.2f} MiB {dest}",flush=True)
    return dest


def make_portal(output,refs,n_generated_parent_events_per_source=None):
    """Browsable local launcher, source-station URLs are all relative."""
    cards=[]
    for src,st in refs:
        href=f"{src}/MS{st}_layout_v3.html"
        cards.append(
            '<a class="card" href="'+html.escape(href,quote=True)+'">'
            '<b>'+html.escape({"jpsi40-bore":"J/ψ","phi40-bore":"ϕ",
            "omega40-bore":"ω"}.get(src,src))+'</b>'
            ' <span>MS'+str(st)+' · интерактивная XY-карта камер</span></a>')
    source="""<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ · PI · шесть конфигураций проволок</title>
<style>body{background:#091626;color:#ecf7ff;font:15px/1.55 system-ui,
sans-serif;max-width:1110px;margin:35px auto;padding:0 20px}
h1{font-size:26px}p{color:#b3cadc}
.grid{display:grid;grid-template-columns:repeat(auto-fill,
minmax(250px,1fr));gap:12px;margin-top:24px}
.card{color:#ecf7ff;text-decoration:none;display:flex;flex-direction:column;
gap:4px;border:1px solid #45617a;background:#152d42;border-radius:9px;
padding:14px}.card:hover{border-color:#59d9c5;background:#24435a}
.card span{color:#b3cadc;font-size:13px}a{color:#8de9e6}</style>
</head><body><h1>NA60+ · время дрейфа · 6 геометрий проволок</h1>
<p>Интерфейс как в исходном pi_centered_timing_dashboard_v2:
кликабельная XY-карта физических камер слева и время дрейфа справа.
Шаг p=3, 3.5, 4 мм; без и с заземлённой полевой проволокой
диаметром 75 мкм. Ниже — углы прихода электронов и сравнительная таблица.
Источник и станцию можно менять на каждой странице.</p>
<p><a href="index.html">Упрощённый обзор</a> ·
<a href="../pi_six_geometry_native_timing.html">
Первоначальный подробный дашборд</a></p>
<div class="grid">__CARDS__</div>
<p>Данные: реальный нативный первый чанк Geant4/Heed + шесть Garfield
LUT; не пороговая электроника. Базовые результаты не изменяются.</p>
</body></html>"""
    source=source.replace("__CARDS__","\n".join(cards))
    if n_generated_parent_events_per_source is None:
        # High-stat merge is a different, independent output directory.
        # Read its already-written VERIFIED native exposure metadata;
        # never infer an exposure from the directory name.
        n_generated_parent_events_per_source=3000
        if refs:
            src,st=refs[0]
            channel={"jpsi40-bore":"jpsi","phi40-bore":"phi",
                     "omega40-bore":"omega"}.get(src)
            if channel:
                merged=Path(output).parent/"merged"/channel/f"MS{st}.json"
                if merged.is_file():
                    info=json.loads(merged.read_text(encoding="utf-8"))
                    n_generated_parent_events_per_source=int(
                        info["meta"].get("n_generated_parent_events_per_source",3000))
    nparents=int(n_generated_parent_events_per_source)
    if nparents!=3000:
        if nparents<3000:raise ValueError("Unsupported high-stat portal exposure")
        source=source.replace(
            "реальный нативный первый чанк Geant4/Heed + шесть Garfield",
            f"{nparents:,} сгенерированных родительских событий на источник; "
            "Geant4/Heed + шесть Garfield"
        )
        source=source.replace(
            '<p><a href="index.html">Упрощённый обзор</a> ·\n'
            '<a href="../pi_six_geometry_native_timing.html">\n'
            'Первоначальный подробный дашборд</a></p>',"")
        if "реальный нативный первый чанк" in source:
            raise AssertionError("Old first-chunk label in high-stat portal")
    portal=Path(output)/"pi_pitch_fieldwire_timing_dashboard_v3.html"
    portal.write_text(source,encoding="utf-8")
    return portal

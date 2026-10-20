#!/usr/bin/env python3
"""Link only newest generated PI dashboards; preserve all physical payloads.

Refreshes navigation in already-existing HTML in place (not ROOT/Heed/
Garfield). Missing new plots are omitted from menus; older projected
and 2+4mm pages remain on disk but do not appear in the new navigation.
"""
from __future__ import annotations
import argparse
import html
from pathlib import Path
import re

REPO=Path(__file__).resolve().parents[2]
PAGES=REPO/"test_runs/mwpc_stage_c/pi_fixed_wire_maps"
HUB="pi_latest_dashboards.html"
DASH=[
 ("pi_native_angle_position.html","Геометрия и углы · XY / u–w",
  "Полный Geant4: 200 000 событий/канал; локальные u–w проекции при наличии полного расчёта",
  "run_pi_native_local_uw.sh all"),
 ("pi_centered_timing_dashboard_v2.html","Время дрейфа и кратность",
  "Пилот: первые 3 000 событий/канал, Heed + симметричная Garfield 3+3 мм",
  "run_pi_centered_timing_dashboard_v2.sh"),
 ("pi_two_anode_electron_explorer.html","Один мюон · два анода",
  "Первичные электроны по времени, два анода и настоящая Geant4 траектория u–v",
  "run_pi_two_anode_electron_explorer.sh"),
 ("pi_centered_true_1d_cloud.html","Облако u(t), v(t), w(t)",
  "Микроскопические Garfield траектории, выбранные углы мюонов; отдельный пилот",
  "run_pi_centered_true_1d_cloud.sh"),
]
ADDITIONAL=["pi_station_interactive.html","pi_muon_track_explorer_native.html"]
HEAD=re.compile(r"(<header\b[^>]*>)(.*?)(</header>)",re.I|re.S)
DATA=re.compile(r'(<script\b[^>]*\btype=["\']application/json["\'][^>]*>)(.*?)(</script>)',re.I|re.S)
STYLE="""<style id="pi-latest-nav-style">
.pi-latest-nav{display:flex;flex-wrap:wrap;align-items:center;gap:7px 9px;
margin:14px 0 3px;padding:10px;border:1px solid #536e80;border-radius:9px;
background:rgba(9,31,47,.75);font:12px/1.4 system-ui,"Segoe UI",sans-serif}
.pi-latest-nav .pi-latest-name{font-weight:750;color:#d4eaf3;margin-right:9px}
.pi-latest-nav a,.pi-latest-nav .pi-latest-active{display:inline-block;
padding:7px 10px;border:1px solid #547a91;border-radius:6px;
text-decoration:none;color:#d1f7f2;background:rgba(35,66,81,.7)}
.pi-latest-nav a:hover{border-color:#84eee3;color:#fff}
.pi-latest-nav .pi-latest-active{font-weight:750;color:#071f2c;
background:#85e5d6;border-color:#85e5d6}
</style>"""


def navbar(current,present):
    menu=[(HUB,"Обзор")]+[(f,title) for f,title,*_ in DASH if f in present]
    pieces=['<nav class="pi-latest-nav" aria-label="Актуальные dashboard NA60+">',
            '<span class="pi-latest-name">NA60+ / PI</span>']
    for filename,label in menu:
        if filename==current:
            pieces.append('<span class="pi-latest-active" aria-current="page">'+
                          html.escape(label)+'</span>')
        else:
            pieces.append('<a href="'+filename+'">'+html.escape(label)+'</a>')
    pieces.append("</nav>")
    return "\n".join(pieces)


def update_page(source,current,present):
    m=HEAD.search(source)
    if m is None:
        raise ValueError("No header in "+current)
    contents=m.group(2)
    # All current templates put prior dashboard links at END of header.
    # Historic first station page wraps them in an independent paragraph.
    cut=re.search(r'(?is)<nav\b[^>]*>|<p\s+style="margin-top:9px"[^>]*>',contents)
    if cut is None:
        cut=re.search(r'(?is)<a\b[^>]*\bhref\s*=\s*["\'][^"\']+\.html',contents)
    clean=contents[:cut.start()] if cut else contents
    revised=(source[:m.start()]+m.group(1)+"\n"+clean.rstrip()+"\n"+
             navbar(current,present)+"\n"+m.group(3)+source[m.end():])
    if '<style id="pi-latest-nav-style">' not in revised:
        close=re.search(r"</head\s*>",revised,re.I)
        if close is None:
            raise ValueError("No </head> in "+current)
        revised=revised[:close.start()]+STYLE+"\n"+revised[close.start():]
    original_data=[part.group(2) for part in DATA.finditer(source)]
    after_data=[part.group(2) for part in DATA.finditer(revised)]
    if original_data!=after_data:
        raise AssertionError("Changed embedded physical data payload")
    if revised.count('class="pi-latest-nav"')!=1:
        raise AssertionError("Navigation not unique")
    header=HEAD.search(revised).group(2)
    for link in re.findall(r'<a\s+href="([^"]+\.html)"',header):
        if link!=HUB and link not in present:
            raise AssertionError("Nonexistent dashboard link: "+link)
    return revised


def index_page(present):
    cards=[]
    for filename,title,scope,runner in DASH:
        available=filename in present
        headline=('<a href="'+filename+'">'+html.escape(title)+'</a>'
                  if available else html.escape(title))
        status=("Готово · открыть" if available else
                "Пока не создан · "+runner)
        cards.append('<article><h2>'+headline+'</h2><p>'+
                     html.escape(scope)+'</p><p class="status">'+
                     html.escape(status)+'</p></article>')
    return """<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ · Актуальные dashboard PI</title>
<style>:root{color-scheme:dark}
body{background:#0b1928;color:#eaf4ff;margin:0;font:14px/1.55 system-ui,"Segoe UI",sans-serif}
main{max-width:1120px;margin:auto;padding:25px}h1{font-size:26px}
h2{font-size:18px}a{color:#87e5d8}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:14px}
article{border:1px solid #45657a;background:#162e44;border-radius:11px;padding:17px}
p{color:#c2d6e4}.status{color:#95e2d2;font-size:12px}
</style></head><body><main><header>
<p>NA60+ · MWPC · PI</p><h1>Актуальные интерактивные результаты</h1>
<p>Открываются только страницы, которые уже существуют на этом компьютере.
Полный Geant4 и пилоты Garfield имеют разные объёмы статистики.
Исторические модели 2+4 мм и прежние условные напряжения не входят в меню.</p>
</header><section class="grid">__CARDS__</section>
<p>Переход между страницами не переносит выбор конкретного события.
Обновление ссылок не пересчитывает ROOT, Garfield или Heed.</p>
</main></body></html>""".replace("__CARDS__","\n".join(cards))


def run(folder):
    folder=Path(folder)
    if not folder.is_dir():
        raise FileNotFoundError(folder)
    present={f for f,*_ in DASH if (folder/f).is_file()}
    result=[]
    for filename in [f for f,*_ in DASH]+ADDITIONAL:
        file=folder/filename
        if not file.is_file():
            continue
        old=file.read_text(encoding="utf-8")
        new=update_page(old,filename,present)
        if update_page(new,filename,present)!=new:
            raise AssertionError("Not idempotent: "+filename)
        if old!=new:
            file.write_text(new,encoding="utf-8")
        result.append(filename)
    (folder/HUB).write_text(index_page(present),encoding="utf-8")
    print("PASS: latest-only inter-dashboard navigation (no physics data changed)")
    print("Active:",", ".join(sorted(present)))
    print("Refreshed:",", ".join(result))
    print("Start page:",folder/HUB)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder",type=Path,default=PAGES)
    run(p.parse_args().folder)


if __name__=="__main__":
    main()

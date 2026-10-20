#!/usr/bin/env python3
"""No physics dependency regression for offline latest-only dashboard links."""
from pathlib import Path
from tempfile import TemporaryDirectory
from refresh_pi_dashboard_navigation import (
    DASH,HUB,HEAD,update_page,index_page,run,
)


def document(title,payload):
    return (
      '<!doctype html><html lang="ru"><head><title>'+title+
      '</title></head><body><header>\n<h1>'+title+'</h1>\n'+
      '<p>Native Geant4 + correct Garfield 3+3 physics</p>\n'+
      '<a href="pi_station_interactive.html">Legacy projection</a> · '+
      '<a href="pi_native_drift_times.html">Legacy 2+4mm</a>\n'+
      '</header><main><script id="data" type="application/json">'+
      payload+'</script></main></body></html>'
    )


def test():
    physics='{"source":"jpsi40-bore","native_events":200000}'
    current=document("Native full angle and local u-w",physics)
    present={DASH[0][0],DASH[1][0]}
    refreshed=update_page(current,DASH[0][0],present)
    assert update_page(refreshed,DASH[0][0],present)==refreshed
    assert refreshed.count('class="pi-latest-nav"')==1
    h=HEAD.search(refreshed).group(2)
    assert 'href="pi_native_drift_times.html"' not in h
    assert 'href="pi_station_interactive.html"' not in h
    assert f'href="{DASH[1][0]}"' in h
    assert f'href="{DASH[2][0]}"' not in h
    assert f'href="{DASH[3][0]}"' not in h
    assert physics in refreshed
    assert 'aria-current="page"' in h

    with TemporaryDirectory(prefix="pi_latest_dashboard_nav_") as d:
        folder=Path(d)
        (folder/DASH[0][0]).write_text(current,encoding="utf-8")
        (folder/DASH[1][0]).write_text(
            document("Centered timing",physics),encoding="utf-8")
        (folder/"pi_station_interactive.html").write_text(
            document("Previous dashboard",physics),encoding="utf-8")
        run(folder)
        assert (folder/HUB).is_file()
        main=(folder/HUB).read_text(encoding="utf-8")
        assert DASH[0][0] in main and DASH[1][0] in main
        assert 'href="'+DASH[3][0]+'"' not in main
        original={
            p.name:p.read_bytes()
            for p in folder.glob("*.html")
        }
        run(folder)
        assert original=={p.name:p.read_bytes() for p in folder.glob("*.html")}
        # Once an optional microscopic cloud file actually exists,
        # every current dashboard must gain its working link.
        (folder/DASH[3][0]).write_text(
            document("Real Garfield microscopic trajectory",physics),
            encoding="utf-8")
        run(folder)
        for filename in (DASH[0][0],DASH[1][0],DASH[3][0]):
            h=HEAD.search((folder/filename).read_text()).group(2)
            if filename!=DASH[3][0]:
                assert 'href="'+DASH[3][0]+'"' in h
            else:
                assert "Облако u(t), v(t), w(t)" in h
        assert physics in (folder/DASH[3][0]).read_text()

    print("PASS: current pages crosslinked; archived 2+4mm and projection pages absent")
    print("PASS: unavailable microscopic dashboard omitted until its HTML exists")
    print("PASS: embedded Geant4/Garfield JSON byte-identical, output refresh idempotent")
    print("PASS: current hub with runnable/correct present or absent status")


if __name__=="__main__":
    test()

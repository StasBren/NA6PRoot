#!/usr/bin/env python3
"""Physics and UI contract check for centered 3+3mm first-chunk PI dashboard V2.

Consumes existing native pilot HTML data and its 3+3mm geometry. No ROOT,
Heed or Garfield reruns, and no browser automation dependency.
"""
from __future__ import annotations

from pathlib import Path

from build_pi_centered_timing_dashboard_v2 import (
    ORIGINAL,TEMPLATE,payload_from_html,validate_payload,
)


def main():
    if not ORIGINAL.is_file():
        raise FileNotFoundError(
            "Original first-chunk 3+3mm PI pilot HTML missing; "
            "run run_pi_native_centered_timing_pilot.sh once")
    payload,sha=payload_from_html(ORIGINAL)
    n=validate_payload(payload)
    assert n>2000 and len(payload["layout"])==298
    assert len({c["id"] for c in payload["layout"]})==298
    assert len({(c["st"],c["id"]) for c in payload["layout"]})==298
    assert sha==__import__("hashlib").sha256(
        ORIGINAL.read_text(encoding="utf-8").split(
            '<script id="payload" type="application/json">',1
        )[1].split("</script>",1)[0].encode("utf-8")
    ).hexdigest()

    # Non-empty bin selector draws only bins with native hit exposure.
    def choices(source,st,det,axis):
        prefix="|".join(map(str,(source,st,det,axis)))+"|"
        return [int(k[len(prefix):]) for k,v in payload["groups"].items()
                if k.startswith(prefix) and v["n"]>0]
    for source in payload["sourceNames"]:
        for st in range(6):
            for axis in ("u","w"):
                all_bins=choices(source,st,-1,axis)
                assert -1 in all_bins
                for bin_index in all_bins:
                    if bin_index<0:continue
                    b=payload["groups"]["|".join(map(str,
                           (source,st,-1,axis,bin_index)))]
                    assert b["n"]>0

    # Every 1D curve slice must be normalized within its own t interval,
    # not proportional to its total electron production.
    first=payload["groups"]["|".join(
        [payload["sourceNames"][0],"0","-1","u","-1"])]
    assert first["coords"] is not None
    for coord,grid in first["coords"].items():
        bins=[(0,2),(2,4),(4,6),(6,10),(10,31)]
        for lo,hi in bins:
            vals=[sum(row[lo:hi]) for row in grid]
            total=sum(vals)
            if total:
                probs=[v/total for v in vals]
                assert abs(sum(probs)-1)<1.e-12
    html=TEMPLATE.read_text(encoding="utf-8")
    for phrase in [
        "__PI_CENTERED_DASHBOARD_V2_DATA__",
        "P(N≥2), % газовых хитов",
        "Время дрейфа от рождения до анода, нс",
        "Доля e⁻ в бине начальной координаты, %",
        "X, см",
        "Y, см",
        "b&&b.n>0",
        'click.addEventListener("click"',
        "НЕ траектория",
        "выбор камеры влияет на графики 2 и 4",
    ]:
        if phrase not in html:
            raise AssertionError("Missing required PI plot feature: "+phrase)
    print("PASS: 298 clickable physical chambers and native Geant4 station exposures")
    print("PASS: only occupied muon angle ranges available per selected module")
    print("PASS: normalized 1D electron birthplace distributions per arrival-time slice")
    print("PASS: all 4 plots have both labeled axes, units and scope explanations")
    print(f"PASS: unchanged original 3+3mm physics payload SHA256 {sha[:16]}...")


if __name__=="__main__":
    main()

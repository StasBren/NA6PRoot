#!/usr/bin/env python3
"""Render an ATLAS-style NA60+ CROSSING-level P(A) dashboard without redoing physics.

Reads already produced pi_lut500_muon_crossing_asymmetry.py products:
summary.json, crossings.csv.gz and manifest.json. Validates SHA256 and
station/physical-DET crossing/histogram conservation. It NEVER uses the
old parent-event dashboard output and NEVER modifies the original summary,
event CSV, provenance manifest, Garfield LUT, Heed or Geant4 data.

Usage:
  python3 pi_crossing_asymmetry_dashboard.py --input-dir PATH
  python3 pi_crossing_asymmetry_dashboard.py --input-dir PATH --replace-dashboard
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

from mwpc_checkerboard_projection import build_chambers
from pi_pitch_fieldwire_pilot import tag

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "pi_crossing_asymmetry_dashboard_template.html"
EXPECTED_STATS = ("hemisphere_A", "top_bottom_A")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def valid_distribution(dist, expected_crossings):
    n = int(expected_crossings)
    valid, undefined = int(dist["valid"]), int(dist["undefined"])
    hist = dist["hist"]
    if (len(hist) != 41 or valid < 0 or undefined < 0 or
            valid + undefined != n or
            any(not isinstance(x, int) or x < 0 for x in hist) or
            sum(hist) != valid):
        raise ValueError("Non-conserved crossing histogram or undefined A count")
    if valid == 0:
        if dist["mean"] is not None or dist["rms"] is not None:
            raise ValueError("No valid crossings but finite A moments")
    else:
        if not all(isinstance(dist[k], (int, float)) and
                   math.isfinite(dist[k]) for k in ("mean", "rms")):
            raise ValueError("Nonfinite crossing mean/RMS")
        if not -1.00001 <= dist["mean"] <= 1.00001 or not 0 <= dist["rms"] <= 1.00001:
            raise ValueError("Unphysical asymmetry mean or RMS")


def assemble(input_dir, verify_raw=True, geometry=None):
    inp = Path(input_dir)
    summary_file = inp / "summary.json"
    raw_file = inp / "crossings.csv.gz"
    manifest_file = inp / "manifest.json"
    summary = json.loads(summary_file.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    trials = summary.get("lut_trials_per_cell")
    pitches = summary.get("pitches_mm")
    stations = summary.get("stations_selected")
    if (trials not in (20, 500) or not isinstance(pitches, list) or
            not pitches or pitches != sorted(set(pitches)) or
            any(p not in (3., 4.) for p in pitches) or
            stations != list(range(6)) or summary.get("source") != "jpsi40-bore" or
            summary.get("schema") !=
            "PI_JPSI_NATIVE_MUON_GAS_CROSSING_ASYM_P3_P4_LUT"+str(trials)+"_V1"):
        raise ValueError("Not a validated six-station J/psi crossing-level source")
    if manifest.get("schema") != summary["schema"]:
        raise ValueError("Wrong manifest for crossing-level physics")
    hashes = manifest.get("sha256_outputs", {})
    if hashes.get("summary.json") != digest(summary_file):
        raise ValueError("Frozen summary.json SHA mismatch")
    if verify_raw:
        if hashes.get("crossings.csv.gz") != digest(raw_file):
            raise ValueError("Frozen per-crossing CSV SHA mismatch")
    elif not raw_file.is_file() or "crossings.csv.gz" not in hashes:
        raise ValueError("Missing frozen per-crossing raw file")
    if summary.get("partial_exposure") != (not manifest.get("complete_six_station")):
        raise ValueError("Crossing exposure flags disagree")
    if manifest.get("generated_parents") != summary.get("generated_parent_slots"):
        raise ValueError("Parent exposure mismatch")
    if manifest.get("garfield_lut_sha256") != summary.get("LUT"+str(trials)):
        raise ValueError("Frozen Garfield provenance mismatch")

    configs = {tag(p, fw) for p in pitches for fw in (0., 75.)}
    if set(summary.get("LUT"+str(trials), {})) != configs:
        raise ValueError("Missing LUT geometry in summary")

    model = geometry if geometry is not None else build_chambers()
    layout = {}
    for layer in model:
        for c in layer:
            st = int(c.station)
            if st not in range(6) or int(c.det_id) in layout.get(st, {}):
                raise ValueError("Invalid or duplicate physical chamber")
            layout.setdefault(st, {})[int(c.det_id)] = {
                "id": int(c.det_id),
                "x": float(c.x), "y": float(c.y),
                "hx": float(c.half_x), "hy": float(c.half_y),
                "layer": str(c.layer),
            }
    all_crossings = 0
    stations_out = {}
    for st in stations:
        station = summary["stations"][str(st)]
        exposure = int(station["gas_crossings"])
        if exposure <= 0:
            raise ValueError("No physical muon crossings for MS"+str(st))
        all_crossings += exposure
        complete = station["configs"]
        det_data = station["det"]
        if set(complete) != configs:
            raise ValueError("Station missing a geometry response")
        if not set(int(k) for k in det_data).issubset(layout[st]):
            raise ValueError("Unknown physical DET in crossing summary")
        for cfg in configs:
            main = complete[cfg]
            if int(main["crossings"]) != exposure:
                raise ValueError("Station crossing exposure inconsistent")
            for obs in EXPECTED_STATS:
                valid_distribution(main[obs], exposure)
            for det, data in det_data.items():
                if cfg not in data:
                    raise ValueError("Physical DET missing a configured geometry")
                sub = data[cfg]
                for obs in EXPECTED_STATS:
                    valid_distribution(sub[obs], int(sub["crossings"]))
            for field in ("crossings", "n_anode_collected"):
                if sum(int(d[cfg][field]) for d in det_data.values()) != int(main[field]):
                    raise ValueError("Station/DET "+field+" conservation failure")
            for obs in EXPECTED_STATS:
                if sum(int(d[cfg][obs]["valid"]) for d in det_data.values()) != main[obs]["valid"]:
                    raise ValueError("Station/DET valid-A crossing conservation failure")
                bins = [sum(int(d[cfg][obs]["hist"][i]) for d in det_data.values())
                        for i in range(41)]
                if bins != main[obs]["hist"]:
                    raise ValueError("Station/DET crossing P(A) conservation failure")
        stations_out[str(st)] = {
            "gas_crossings":exposure,
            "configs":{cfg:{
                "crossings":int(complete[cfg]["crossings"]),
                **{obs:complete[cfg][obs] for obs in EXPECTED_STATS},
            } for cfg in sorted(configs)},
            "det":{det:{cfg:{
                "crossings":int(values[cfg]["crossings"]),
                **{obs:values[cfg][obs] for obs in EXPECTED_STATS}
            } for cfg in sorted(configs)} for det, values in det_data.items()},
            "layout":sorted(layout[st].values(), key=lambda c:c["id"]),
        }

    payload = {
        "schema":summary["schema"],
        "source":summary["source"],
        "lut_trials_per_cell":trials,
        "source_chunks":summary["source_chunks"],
        "generated_parent_slots":summary["generated_parent_slots"],
        "partial_exposure":summary["partial_exposure"],
        "all_crossings":all_crossings,
        "stations_selected":stations,
        "pitches_mm":pitches,
        "definition":summary["definitions"],
        "stations":stations_out,
    }
    return payload, {p.name:digest(p) for p in (summary_file, raw_file)
                     if p.name=="summary.json" or verify_raw}


def render(input_dir, out=None, replace=False, verify_raw=True, geometry=None,
           template=TEMPLATE):
    inp = Path(input_dir)
    out = Path(out) if out is not None else inp / "crossing_dashboard.html"
    source_hashes = {name:digest(inp/name) for name in ("summary.json",)}
    payload, checks = assemble(inp,verify_raw=verify_raw,geometry=geometry)
    if any(digest(inp/name) != value for name,value in source_hashes.items()):
        raise ValueError("Physics summary changed while preparing dashboard")
    if not template.is_file():
        raise FileNotFoundError(template)
    html = template.read_text(encoding="utf-8")
    if html.count("__PAYLOAD__") != 1:
        raise ValueError("Dashboard template must have exactly one data marker")
    safe = json.dumps(payload,ensure_ascii=False,separators=(",",":"),
                      allow_nan=False).replace("</","<\\/")
    html = html.replace("__PAYLOAD__",safe)
    if "parent event" not in html or "physical muon gas crossing" not in html:
        raise ValueError("Missing crossing-level semantics in dashboard")
    if out.exists() and not replace:
        if out.read_text(encoding="utf-8") == html:
            print("PASS: existing dashboard is byte-identical, unchanged:",out)
            return out
        raise FileExistsError("Dashboard exists; use --replace-dashboard for HTML-only refresh: "+str(out))
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp = out.with_name(out.name+".tmp")
    if tmp.exists():
        raise FileExistsError("Dashboard temp file exists: "+str(tmp))
    try:
        tmp.write_text(html,encoding="utf-8")
        if digest(inp/"summary.json") != source_hashes["summary.json"]:
            raise ValueError("Summary changed while writing dashboard")
        os.replace(tmp,out)
    finally:
        if tmp.exists():tmp.unlink()
    receipt = {
        "dashboard_file":str(out),
        "dashboard_sha256":digest(out),
        "source_input_dir":str(inp),
        "physics_sha256":checks,
        "generated_parent_slots":payload["generated_parent_slots"],
        "physical_muon_gas_crossings":payload["all_crossings"],
        "statistical_unit":"one real gas crossing; no parent pooling",
    }
    receipt_path=out.with_suffix(".provenance.json")
    receipt_tmp=receipt_path.with_suffix(".provenance.json.tmp")
    receipt_tmp.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    os.replace(receipt_tmp,receipt_path)
    print("PASS: standalone crossing-level dashboard:",out)
    print("Validated muon gas crossings:",payload["all_crossings"],
          "LUT trials per cell:",payload["lut_trials_per_cell"])
    print("Original physics summary/raw CSV/manifest were not changed.")
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--replace-dashboard",action="store_true",
                    help="Replace HTML-only dashboard, never physical data or manifest")
    ap.add_argument("--skip-raw-sha",action="store_true",
                    help="Only check manifest and summary checksum; quicker, less exhaustive")
    args=ap.parse_args()
    render(args.input_dir,args.out,args.replace_dashboard,
           verify_raw=not args.skip_raw_sha)


if __name__=="__main__":
    main()

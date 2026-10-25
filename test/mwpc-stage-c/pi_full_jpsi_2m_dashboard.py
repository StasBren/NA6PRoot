#!/usr/bin/env python3
"""FULL 2,000,000-parent J/psi interactive six-geometry MWPC dashboard.

Read ONLY the 40 *completed* Geant4/Heed/6-LUT stage-C gzip folds
for each of the six stations (240 files). Validate SHA256 and immutable
HEED/GARFIELD/fold provenance, aggregate unrounded electron moments,
validate station/chamber/physical-anode electron conservation, generate
the original physical XY V3 and detailed per-anode interactive dashboards,
translate pages to English without modifying numeric JSON payloads.

No new Geant4, Heed, Garfield, parent simulation or fold. The still-running
6M suite, original 100k report and 20-trial frozen LUTs are NEVER changed.

Usage:
 python3 test/mwpc-stage-c/pi_full_jpsi_2m_dashboard.py check
 python3 test/mwpc-stage-c/pi_full_jpsi_2m_dashboard.py build
 python3 test/mwpc-stage-c/pi_full_jpsi_2m_dashboard.py verify
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import html
import json
import os
from pathlib import Path

from pi_2m_native_six_geometry_fold import (
    OUTPUT, VERSION, SOURCE, FIELDS, shard_path, sha_file, restore,
    key, as_row,
)
from pi_six_geometry_native_timing import (
    Count, ANGLES, n_anodes, validate_conservation, WIDTHS, PITCHES, tag,
)
from pi_six_geometry_layout_v3 import render as render_xy, make_portal
from pi_six_geometry_native_timing_html import build as render_detail
from build_pi_100k_english_dashboard import translate_page

HERE=Path(__file__).resolve().parent
BASE=OUTPUT.parent/"pi_full_jpsi_2m_dashboard_v1"
MERGED=BASE/"merged"/"jpsi"
PAGES=BASE/"pages"
NAME=SOURCE["Jpsi"]
CHANNEL="Jpsi"
N_CHUNKS=40
N_PARENTS=2_000_000
PLOTS=("p3p0_anode_only","p3p0_field75p0um",
       "p3p5_anode_only","p3p5_field75p0um",
       "p4p0_anode_only","p4p0_field75p0um")
FINGERPRINT_FIELDS=(
    "garfield_map_provenance","heed_tracks_sha256",
    "heed_clusters_sha256","fold_source_sha256",
    "fold_driver_sha256","rng_semantics",
)
SCHEMA="NA60_STAGEC_JPSI_2M_FULL_XY_TIMING_ANODE_ANGLE_V1"


def atomic_json(path,payload):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".incomplete")
    try:
        tmp.write_text(json.dumps(payload,separators=(",",":"),
                            ensure_ascii=False,allow_nan=False)+"\n",
                       encoding="utf-8")
        os.replace(tmp,path)
    finally:
        if tmp.exists():tmp.unlink()


def checked_manifest(chunk,station,*,verify_sha):
    location=shard_path(OUTPUT/"full",CHANNEL,chunk,station)
    manifest=location.with_suffix(".manifest.json")
    if not location.is_file() or not manifest.is_file():
        raise FileNotFoundError(f"Missing completed immutable Heed/LUT shard: {location}")
    m=json.loads(manifest.read_text(encoding="utf-8"))
    frozen=m.get("immutable_inputs",{})
    expected={
        "model":VERSION,"source":NAME,"chunk":chunk,"station":station,
        "generated_parent_events_in_input_chunk":50_000,
        "pilot_muon_gas_segment_limit":0,
    }
    if any(frozen.get(k)!=v for k,v in expected.items()) or \
       m.get("n_generated_parent_events_in_chunk")!=50_000 or \
       m.get("n_groups",0)<6 or m.get("n_muon_gas_hits",0)<1 or \
       not m.get("gzip_sha256"):
        raise ValueError(f"Invalid frozen full-channel provenance: {manifest}")
    fingerprints={k:frozen.get(k) for k in FINGERPRINT_FIELDS}
    luts=frozen.get("garfield_map_provenance",{})
    if set(luts)!=set(PLOTS) or any(v.get("samples_per_box")!=20
                                      for v in luts.values()):
        raise ValueError(f"Wrong original 20-trial six Garfield LUTs: {manifest}")
    if verify_sha and sha_file(location)!=m["gzip_sha256"]:
        raise ValueError(f"SHA256 failure: {location}")
    return location,manifest,m,fingerprints


def precheck(verbose=True):
    """Cheap complete 240-manifest preflight before any big output is touched."""
    stats=[]
    universal=None
    for st in range(6):
        row={"station":st,"n_completed_chunks":0,"n_gas_crossings":0,
             "n_generated_parent_events":N_PARENTS}
        for i in range(1,N_CHUNKS+1):
            _,_,m,f=checked_manifest(i,st,verify_sha=False)
            if universal is None:universal=f
            if universal!=f:
                raise ValueError(f"Physics mismatch: MS{st} chunk{i:04d}")
            row["n_completed_chunks"]+=1
            row["n_gas_crossings"]+=m["n_muon_gas_hits"]
        stats.append(row)
        if verbose:print(f"READY J/psi MS{st}: 40/40 SHA-manifested shards; "
                         f"{row['n_gas_crossings']:,} muon gas segments",flush=True)
    print(f"PASS: 240/240 J/psi fold checkpoints, same frozen Garfield/Heed "
          "physics, 2M generated parents/channel; SHA check occurs on build.")
    return stats


def group_quality(row):
    if (len(row.get("hist",[]))!=31 or
        sum(row["hist"])!=row["known"] or
        len(row.get("angle_hist",[]))!=18 or
        sum(row["angle_hist"])!=row["angle_n"] or
        len(row.get("phi_hist",[]))!=36 or
        sum(row["phi_hist"])!=row["phi_n"] or
        row["known"]+row["unknown"]!=row["collected"]):
        raise ValueError(f"Corrupt Garfield histogram conservation: {row.get('config')}")


def sum_bins(groups,cfg,det,axis):
    return [r for r in groups.values()
            if r[0]==cfg and r[2]==det and r[4]==axis]


def make_station(st):
    merged=defaultdict(Count)
    expected_physics=None
    exposure=0
    witnesses=[]
    for chunk in range(1,N_CHUNKS+1):
        path,manifest,m,f=checked_manifest(chunk,st,verify_sha=True)
        if expected_physics is None:expected_physics=f
        if expected_physics!=f:
            raise ValueError(f"Mixed frozen physics MS{st}, chunk{chunk}")
        # Reading an existing immutable .gz is safe while a different source
        # is actively being folded by the ongoing production launcher.
        with __import__("gzip").open(path,"rt",encoding="utf-8") as inp:
            doc=json.load(inp)
        if doc.get("meta")!=m["immutable_inputs"] or \
           doc.get("group_count")!=m["n_groups"] or \
           len(doc.get("groups",[]))!=m["n_groups"]:
            raise ValueError(f"Truncated or mismatched gzip fold: {path}")
        per_cfg={}
        seen=set()
        for r in doc["groups"]:
            if r.get("source")!=NAME or r.get("station")!=st or \
               r.get("config") not in PLOTS:
                raise ValueError(f"Foreign source/station/geometry: {path}")
            k=key(r)
            if k in seen:raise ValueError(f"Duplicate {k} in {path}")
            seen.add(k)
            group_quality(r)
            restore(merged[k],r)
            if r["det_id"]==-1 and r["anode"]==-1 and \
               r["axis"]=="all" and r["anglebin"]==-1:
                per_cfg[r["config"]]=r["segments"]
        if set(per_cfg)!=set(PLOTS) or \
           any(v!=m["n_muon_gas_hits"] for v in per_cfg.values()):
            raise ValueError(f"Station source exposure mismatch: {path}")
        exposure+=m["n_muon_gas_hits"]
        witnesses.append({
            "chunk":chunk,"gzip_sha256":m["gzip_sha256"],
            "immutable_manifest_sha256":sha_file(manifest),
            "n_muon_gas_hits":m["n_muon_gas_hits"],
            "n_groups":m["n_groups"],
        })
        del doc
        print(f"MS{st} chunk{chunk:04d}/0040: "
              f"{m['n_muon_gas_hits']:,} actual muon gas crossings; "
              f"unique accumulated groups {len(merged):,}",flush=True)

    validate_conservation(merged)
    collected={}
    primaries={}
    heads={}
    all_groups=[]
    # Histogram conservation by angle projection in every station and chamber.
    for k,bucket in sorted(merged.items()):
        row=as_row(k,bucket)
        group_quality(row)
        all_groups.append(row)
    # Only 6 large station summaries need to enter the downstream comparative table.
    all_rows={(g["config"],g["det_id"],g["anode"],
               g["axis"],g["anglebin"]):g for g in all_groups}
    station_heed=None
    for cfg in PLOTS:
        head=all_rows.get((cfg,-1,-1,"all",-1))
        if head is None or head["segments"]!=exposure:
            raise ValueError(f"Missing or inconsistent station MS{st} {cfg}")
        if station_heed is None:station_heed=head["primaries"]
        elif head["primaries"]!=station_heed:
            raise ValueError("Shared Heed sampling broken between six LUTs")
        for axis in ("u","w"):
            binned=[all_rows.get((cfg,-1,-1,axis,bi))
                    for bi in range(len(ANGLES)-1)]
            for dimension in ("segments","primaries","known",
                              "collected","angle_n","phi_n"):
                if sum(0 if r is None else r[dimension] for r in binned)!= \
                   head[dimension]:
                    raise ValueError(f"Muon θ_{axis} bin sum != station for "
                                     f"{cfg} {dimension} MS{st}")
        collected[cfg]=head["collected"]
        primaries[cfg]=head["primaries"]
        heads[cfg]=head

    meta={
        "model_version":VERSION,"dashboard_version":SCHEMA,
        "data_scope":"COMPLETE 40 Geant4 50,000-parent chunks: "
                     "2,000,000 GENERATED J/psi parent events; "
                     "all station muon gas crossings; common Heed births "
                     "and six original Garfield 20-trial/cell LUT responses",
        "n_generated_parent_events_per_source":N_PARENTS,
        "gas_hit_counts":{NAME:exposure},
        "n_muon_gas_crossings_station":exposure,
        "n_heed_primary_electrons_station":station_heed,
        "source_station":f"{NAME} MS{st}",
        "source_channel":"Jpsi","geometries":6,
        "garfield_luts":expected_physics["garfield_map_provenance"],
        "angle_bins_deg":ANGLES.tolist(),
        "time_bin_width_ns":10,
        "time_edges_ns":list(range(0,301,10))+["overflow"],
        "electron_angle_edges_deg":list(range(0,181,10)),
        "electron_phi_hit_edges_deg":list(range(-180,181,10)),
        "finite_anodes_by_pitch":{str(p):[n_anodes(s,p) for s in range(6)]
                                  for p in PITCHES},
        "physical_gas_width_mm":list(WIDTHS),
        "input_model":"True TGeo-local Geant4 muon gas crossings, "
                      "Poisson/Heed-based primary ionization, and common "
                      "birth-coordinate/quantile ensemble fed to six "
                      "separately frozen sampled Garfield response LUTs",
        "angle_semantics":"Geant4 signed θ_u/θ_w entrance angle; "
                          "electron incidence_normal_deg is from LAST "
                          "microscopic Garfield step, φ_hit is circular",
        "simulation_limitations":"B=0 ideal field; 20 independent Garfield "
                          "drift trials per u/v cell, reused via lookup; "
                          "no avalanche, ionic pulse, cathode strips or "
                          "front-end electronics; field wire 75um is an "
                          "unverified engineering study assumption",
        "statistical_errors":"Parent-event/cluster bootstrap not yet reported; "
                             "reused Garfield LUT electrons are not "
                             "independent microscopic trajectories",
        "chunk_manifest_sha256":[r["immutable_manifest_sha256"] for r in witnesses],
        "chunk_gzip_sha256":[r["gzip_sha256"] for r in witnesses],
    }
    target=MERGED/f"MS{st}.json"
    atomic_json(target,{"meta":meta,"group_count":len(all_groups),
                        "groups":all_groups})
    del all_groups, merged, all_rows
    dest=PAGES/NAME
    render_xy(target,dest/f"MS{st}_layout_v3.html",NAME,st)
    render_detail(target,dest/f"MS{st}.html")
    translate_page(dest/f"MS{st}_layout_v3.html",
                   dest/f"MS{st}_layout_v3_en.html",
                   "layout",expected_parents=N_PARENTS)
    translate_page(dest/f"MS{st}.html",
                   dest/f"MS{st}_en.html",
                   "detail",expected_parents=N_PARENTS)
    # Fixed historical translator string (not the scientific data payload).
    detail_path=dest/f"MS{st}_en.html"
    dt=detail_path.read_text(encoding="utf-8")
    if "MWPC-NATIVE-CHUNK-MERGE-100K-V1" in dt:
        dt=dt.replace("MWPC-NATIVE-CHUNK-MERGE-100K-V1",
                      "MWPC-NATIVE-CHUNK-MERGE-2M-V1")
        t=detail_path.with_name(detail_path.name+".incomplete")
        t.write_text(dt,encoding="utf-8")
        os.replace(t,detail_path)
    summary={
        "schema":SCHEMA,"station":st,"source":"Jpsi",
        "generated_parent_events":N_PARENTS,
        "gas_crossings":exposure,"heed_primaries":station_heed,
        "n_independent_garfield_trials_per_box":20,
        "six_geometries":{},
        "sha256_input_files":witnesses,
        "files":{
          name:sha_file(dest/name)
          for name in (f"MS{st}_layout_v3.html",f"MS{st}.html",
                       f"MS{st}_layout_v3_en.html",f"MS{st}_en.html")
        },
        "merged_json_sha256":sha_file(target)
    }
    for cfg in PLOTS:
        head=heads[cfg]
        summary["six_geometries"][cfg]={
           "n_anode_collected":head["collected"],
           "n_heed_primaries":head["primaries"],
           "collection_fraction":head["collected"]/head["primaries"],
           "known_time_electrons":head["known"],
           "mean_drift_ns":head["mean_ns"],
           "rms_drift_ns":head["rms_ns"],
        }
    atomic_json(BASE/"station_summaries"/f"MS{st}.json",summary)
    print(f"PASS FULL J/psi 2M MS{st}: {exposure:,} muon gas "
          f"crossings; {station_heed:,} primary electrons; "
          f"all six LUT geometries, physical chambers, anodes, "
          f"signed muon entrance and electron arrival-angle histograms",
          flush=True)
    return summary


def reusable_station(st):
    """Re-enter an interrupted run without rebuilding finished station pages.

    Resuming never trusts only filenames: scientific source-manifest SHAs
    and local output SHA256 are compared against the station publication.
    """
    publication=BASE/"station_summaries"/f"MS{st}.json"
    if not publication.is_file():return False
    try:
        info=json.loads(publication.read_text(encoding="utf-8"))
        if info.get("schema")!=SCHEMA or info.get("station")!=st or \
           info.get("generated_parent_events")!=N_PARENTS or \
           set(info.get("six_geometries",{}))!=set(PLOTS):
            return False
        if sha_file(MERGED/f"MS{st}.json")!=info["merged_json_sha256"]:
            return False
        dest=PAGES/NAME
        for name,digest in info["files"].items():
            if sha_file(dest/name)!=digest:return False
        sources=info.get("sha256_input_files",[])
        if len(sources)!=N_CHUNKS:return False
        for item in sources:
            location,manifest,m,_=checked_manifest(
                item["chunk"],st,verify_sha=False)
            if (m["gzip_sha256"]!=item["gzip_sha256"] or
                sha_file(manifest)!=item["immutable_manifest_sha256"]):
                return False
        return True
    except (KeyError,ValueError,OSError,TypeError):
        return False


def html_index(summaries):
    # Standalone landing page, navigates all six independent offline pages.
    blocks=[]
    total=sum(x["gas_crossings"] for x in summaries)
    for x in summaries:
        n=x["station"]
        nums=x["six_geometries"]["p3p0_anode_only"]
        blocks.append(
          f'<a class="card" href="{NAME}/MS{n}_layout_v3_en.html">'
          f'<strong>MS{n}</strong><span>{x["gas_crossings"]:,} muon gas crossings'
          f' · {x["heed_primaries"]:,} Heed primary electrons</span>'
          f'<span>p=3 mm, no field wire: mean drift {nums["mean_drift_ns"]:.2f} ns</span>'
          f'<small>Open clickable XY chamber map →</small></a>')
    page="""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ · J/ψ 2M parents · Six-geometry complete MWPC dashboard</title>
<style>
:root{color-scheme:dark}*{box-sizing:border-box}
body{margin:0;background:#091626;color:#edf6fb;font:15px/1.55 system-ui,Arial,sans-serif}
main{max-width:1100px;margin:auto;padding:35px 25px}h1{font-size:31px;letter-spacing:-.03em}
.intro{color:#b4cbd9}.meta{display:flex;gap:10px;flex-wrap:wrap;margin:20px 0}
.meta b{background:#1c4053;color:#7be3d5;border-radius:10px;padding:11px 15px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));
 gap:12px;margin:22px 0}.card{display:flex;flex-direction:column;gap:8px;
border:1px solid #335b75;border-radius:12px;padding:23px;
background:#153047;color:#eaf5fc;text-decoration:none}.card:hover{border-color:#4bd5c3}
.card strong{font-size:25px}.card span{font-size:13px;color:#b5ccdb}
.card small{font-size:12px;color:#69e5d4}.notice{margin-top:27px;border-left:3px solid
#f1bf83;padding:12px 16px;color:#c5d6df;background:#173047}
a{color:#83e0e2}
</style></head><body><main>
<p style="color:#5adccf;font-weight:700;letter-spacing:.12em">NA60+ / STAGE C / FULL PARENT SAMPLE</p>
<h1>J/ψ · 2,000,000 generated parents</h1>
<p class="intro">A complete native muon-gas-crossing, common Heed birth,
six-geometry sampled Garfield LUT dashboard. Choose one of six muon
spectrometer stations. On each page: clickable physical XY chamber map,
individual DET modules and numbered anodes, 10 ns electron drift time,
signed muon entrance angles (θu/θw), polar electron approach angles,
signed φ, and six-way pitch / field-wire comparisons.</p>
<div class="meta">
<b>2,000,000 generated parent events</b>
<b>__TOTAL__ gas-crossing segments across six stations</b>
<b>Six geometries · 20 Garfield independent trials/cell</b></div>
<div class="grid">__CARDS__</div>
<p><a href="pi_pitch_fieldwire_timing_dashboard_v3_en.html">
Alternative compact source/station navigation</a> ·
<a href="../summary.csv">Numerical statistics (CSV)</a></p>
<div class="notice"><strong>Definitions:</strong> gas crossings are NOT unique muons,
and repeated Heed electron responses are NOT independent Garfield drift
trials. Measured time is birth-to-anode drift, not threshold or front-end
readout time. Frozen B=0 Garfield geometry; 75 μm field-wire diameter is
provisional. No avalanche or detector-efficiency inference; parent-event
bootstrap uncertainty not yet reported. The original J/ψ/ϕ/ω 6M production
is left unchanged.</div>
</main></body></html>"""
    page=page.replace("__TOTAL__",f"{total:,}").replace("__CARDS__",
                                                   "\n".join(blocks))
    location=PAGES/"index.html"
    tmp=location.with_name("index.html.incomplete")
    tmp.write_text(page,encoding="utf-8")
    os.replace(tmp,location)
    return location


def finalize():
    summaries=[]
    for st in range(6):
        file=BASE/"station_summaries"/f"MS{st}.json"
        if not file.is_file():raise FileNotFoundError(f"Missing station {st}: {file}")
        doc=json.loads(file.read_text(encoding="utf-8"))
        if doc.get("generated_parent_events")!=N_PARENTS or \
           doc.get("station")!=st or set(doc.get("six_geometries",{}))!=set(PLOTS):
            raise ValueError(f"Invalid 2M Jpsi summary: {file}")
        summaries.append(doc)
    csvfile=BASE/"summary.csv"
    with csvfile.open("w",newline="",encoding="utf-8") as output:
        w=csv.writer(output)
        w.writerow(("source","generated_parent_events","station",
                    "muon_gas_crossings","heed_primary_electrons",
                    "config","collected_electrons","collection_fraction",
                    "known_time_electrons","mean_drift_ns","rms_drift_ns",
                    "independent_garfield_trials_per_spatial_cell"))
        for row in summaries:
            for cfg,vals in row["six_geometries"].items():
                w.writerow(("Jpsi",N_PARENTS,row["station"],
                   row["gas_crossings"],row["heed_primaries"],cfg,
                   vals["n_anode_collected"],vals["collection_fraction"],
                   vals["known_time_electrons"],vals["mean_drift_ns"],
                   vals["rms_drift_ns"],20))
    refs=[(NAME,st) for st in range(6)]
    portal=make_portal(PAGES,refs,n_generated_parent_events_per_source=N_PARENTS)
    translate_page(portal,PAGES/"pi_pitch_fieldwire_timing_dashboard_v3_en.html",
                   "portal",expected_parents=N_PARENTS)
    index=html_index(summaries)
    out=BASE/"complete_dashboard_manifest.json"
    payload={
       "schema":SCHEMA,"source":"Jpsi",
       "generated_parent_events":N_PARENTS,
       "n_stations":6,"verified_chunk_station_shards":240,
       "gas_crossings_all_station_exposures":sum(x["gas_crossings"] for x in summaries),
       "physics":"Geant4/Heed original six 20-trial Garfield LUTs; "
                 "separate higher-stat 50-trial LUT campaign not substituted",
       "scope_limitations":"No new Garfield, Heed or Geant4; no parent bootstrap; "
                           "individual electrons are not statistically independent "
                           "Geant4 parent events.",
       "station_reports":{f"MS{v['station']}":{
           "sha256":sha_file(BASE/"station_summaries"/f"MS{v['station']}.json"),
           "muon_gas_crossings":v["gas_crossings"],
           "heed_primaries":v["heed_primaries"],
       } for v in summaries},
       "files_sha256":{
           "index.html":sha_file(index),
           "portal_en":sha_file(PAGES/"pi_pitch_fieldwire_timing_dashboard_v3_en.html"),
           "summary_csv":sha_file(csvfile),
       },
    }
    atomic_json(out,payload)
    print("PASS: Full J/psi 2M complete-six-station dashboard:",index,flush=True)
    print("Numerical summary:",csvfile,flush=True)
    return payload


def verify():
    done=BASE/"complete_dashboard_manifest.json"
    if not done.is_file():raise FileNotFoundError(done)
    complete=json.loads(done.read_text(encoding="utf-8"))
    if complete.get("verified_chunk_station_shards")!=240 or \
       complete.get("generated_parent_events")!=N_PARENTS:
        raise ValueError("Full Jpsi dashboard provenance incorrect")
    for st in range(6):
        info=json.loads((BASE/"station_summaries"/f"MS{st}.json").read_text())
        merged=MERGED/f"MS{st}.json"
        if sha_file(merged)!=info["merged_json_sha256"]:
            raise ValueError("Merged Jpsi station JSON changed: "+str(merged))
        path=PAGES/NAME
        for relative,expected in info["files"].items():
            if sha_file(path/relative)!=expected:
                raise ValueError("Station HTML checksum changed: "+relative)
        for item in info["sha256_input_files"]:
            location,manifest,m,_=checked_manifest(
                item["chunk"],st,verify_sha=False)
            if m["gzip_sha256"]!=item["gzip_sha256"] or \
               sha_file(manifest)!=item["immutable_manifest_sha256"]:
                raise ValueError("Original frozen source manifest changed: "+str(manifest))
        print(f"PASS verified MS{st}: 40 immutable file SHAs and four "
              "HTML variants; station merged SHA",flush=True)
    for rel,digest in complete["files_sha256"].items():
        file={"index.html":PAGES/"index.html",
              "portal_en":PAGES/"pi_pitch_fieldwire_timing_dashboard_v3_en.html",
              "summary_csv":BASE/"summary.csv"}[rel]
        if sha_file(file)!=digest:raise ValueError("Dashboard publication changed: "+str(file))
    print("PASS: complete 2M Jpsi dashboard files for all SIX stations validated")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("check","build","verify"))
    args=p.parse_args()
    if args.mode=="verify":verify();return
    precheck()
    if args.mode=="check":return
    for st in range(6):
        if reusable_station(st):
            print(f'SKIP already SHA-verified complete J/psi MS{st} dashboard',flush=True)
        else:
            make_station(st)
    finalize()
    verify()

if __name__=="__main__":main()

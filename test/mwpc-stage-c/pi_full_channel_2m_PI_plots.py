#!/usr/bin/env python3
"""READ-ONLY PI plots using ALL 40 completed 50k-parent six-LUT fold shards.

No Geant4, Heed, Garfield, or recomputation of frozen 6M production.
Only one source/station is read at a time to limit memory; the 240 immutable
Jpsi shards are NOT modified. Full source exposure is 2,000,000 generated
parent events, NOT 2,000,000 muon hits or 2,000,000 Garfield trajectories.

Example:
 python3 test/mwpc-stage-c/pi_full_channel_2m_PI_plots.py check --channel Jpsi --station 4
 python3 test/mwpc-stage-c/pi_full_channel_2m_PI_plots.py plots --channel Jpsi --station 4 --det 143

The full original 6M campaign is free to continue running on other
source/station shards while this program reads completed Jpsi checkpoints.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
import os
from pathlib import Path

import numpy as np

from pi_2m_native_six_geometry_fold import (
    OUTPUT, SOURCE, CHANNELS, VERSION, sha_file, shard_path,
    restore, key, as_row,
)
from pi_six_geometry_native_timing import Count, ANGLES, n_anodes
import pi_pi_meeting_preliminary_20261025 as plotter

CONFIGS=tuple(plotter.CONFIGS)
SAME_PHYSICS_KEYS=(
    "garfield_map_provenance","heed_tracks_sha256",
    "heed_clusters_sha256","fold_source_sha256",
    "fold_driver_sha256","rng_semantics",
)
BASE=OUTPUT.parent/"pi_full_channel_2m_PI_v1"
GROUP_TYPES=("station","chamber","anodes",
             "theta_u_station","theta_w_station",
             "theta_u_chamber","theta_w_chamber")

def file_manifest(ch,index,station,expected_physics=None,checksum=True):
    out=OUTPUT/"full"
    path=shard_path(out,ch,index,station)
    mfile=path.with_suffix(".manifest.json")
    if not mfile.is_file() or not path.is_file():
        raise FileNotFoundError(
            f"Missing completed 2M input shard: {path} and {mfile}; "
            "never substitute unfinished .incomplete data"
        )
    m=json.loads(mfile.read_text(encoding="utf-8"))
    frozen=m.get("immutable_inputs",{})
    expected={
        "model":VERSION, "source":SOURCE[ch],"chunk":index,
        "station":station,"pilot_muon_gas_segment_limit":0,
        "generated_parent_events_in_input_chunk":50000,
    }
    if (any(frozen.get(k)!=v for k,v in expected.items()) or
        m.get("n_generated_parent_events_in_chunk")!=50000 or
        m.get("n_muon_gas_hits",0)<=0 or
        m.get("n_groups",0)<=0):
        raise ValueError(f"Inconsistent immutable 2M physics manifest: {mfile}")
    fingerprints={k:frozen.get(k) for k in SAME_PHYSICS_KEYS}
    if expected_physics is not None and fingerprints!=expected_physics:
        raise ValueError(f"Mixed Heed/Garfield/fold-source physics: {mfile}")
    lut=frozen.get("garfield_map_provenance",{})
    if set(lut)!=set(CONFIGS) or any(v.get("samples_per_box")!=20
                                      for v in lut.values()):
        raise ValueError(f"Not the intended original SIX × 20-trial LUTs: {mfile}")
    if checksum and sha_file(path)!=m.get("gzip_sha256"):
        raise ValueError(f"SHA256 corrupted/incomplete folded gzip: {path}")
    return path,mfile,m,fingerprints

def which_group(row,station,det):
    did=int(row["det_id"])
    wire=int(row["anode"])
    axis=row["axis"]
    bi=int(row["anglebin"])
    if did==-1 and wire==-1:
        if axis=="all" and bi==-1:return "station"
        if axis in ("u","w") and 0<=bi<len(ANGLES)-1:
            return f"theta_{axis}_station"
    if did==det:
        if wire==-1:
            if axis=="all" and bi==-1:return "chamber"
            if axis in ("u","w") and 0<=bi<len(ANGLES)-1:
                return f"theta_{axis}_chamber"
        if wire>=0 and axis=="all" and bi==-1:
            return "anodes"
    return None

def reduced_key(row,group_type):
    cfg=row["config"]
    if group_type=="anodes":return (cfg,group_type,int(row["anode"]))
    if group_type.startswith("theta_"):
        return (cfg,group_type,int(row["anglebin"]))
    return (cfg,group_type,-1)

def convert(bucket):
    return bucket.done()

def build_report(ch,station,det,out,*,precheck=False):
    if ch not in CHANNELS or station not in range(6) or det<0:
        raise ValueError("Invalid channel/station/det")
    if precheck:
        previous=None
        for i in range(1,41):
            _,_,_,previous=file_manifest(ch,i,station,previous,checksum=False)
        print(f"PASS: all 40 manifests+gzip present for {ch} MS{station}; "
              "fingerprints match. SHA checks and gzip decode are in plots mode.")
        return None
    output={
        "schema":"NA60_PI_FULL_CHANNEL_2M_SIX_GEOM_READONLY_V1",
        "source":ch,"station":station,"chosen_detector":det,
        "generated_parent_events":2000000,
        "exposure":"2,000,000 generated parents from 40 SHA256-verified "
                   "50k Geant4 input chunks, one source; station gas-crossing "
                   "exposure and electron counts in data.station",
        "angle_semantics":"Signed muon theta_u,theta_w in detector-local axes; "
                  "Garfield electron incidence_normal_deg is last microscopic "
                  "step w.r.t. local anode surface normal; phi_hit_deg is "
                  "last-step azimuth, NOT line direction along the anode",
        "muon_angle_bin_edges_deg":ANGLES.tolist(),
        "electron_incidence_bin_edges_deg":list(range(0,181,10)),
        "electron_phi_bin_edges_deg":list(range(-180,181,10)),
        "time_bin_width_ns":10,
        "time_overflow":"hist[30] = t>=300ns, included in known-time denominator",
        "normalization":"Known collected electrons for drift time; finite "
            "arrival direction among anode-collected electrons for angles; "
            "electron counts are not independent Garfield LUT trials",
        "limitations":"Single-electron Garfield LUT 20 independent drift trials/box; "
            "B=0; no avalanche, electronics response or parent-event bootstrap; "
            "field-wire 75um diameter a provisional study parameter",
        "chunks_used":list(range(1,41)),
        "input_shard_sha256":{},
        "data":{cfg:{"station":None,"chamber":None,"anodes":[],
                     "theta_u_station":{},"theta_w_station":{},
                     "theta_u_chamber":{},"theta_w_chamber":{}} for cfg in CONFIGS}
    }
    agg=defaultdict(Count)
    previous=None
    total_gas_segments=0
    for i in range(1,41):
        path,manifest_file,m,previous=file_manifest(ch,i,station,previous)
        print(f"MERGE (READ ONLY) {ch} MS{station} chunk{i:04d}/0040 "
              f"gas segments={m['n_muon_gas_hits']:,}",flush=True)
        with gzip.open(path,"rt",encoding="utf-8") as stream:
            doc=json.load(stream)
        if doc.get("meta")!=m["immutable_inputs"]:
            raise ValueError(f"Frozen payload/manifest mismatch: {path}")
        if doc.get("group_count")!=m["n_groups"] or \
           len(doc.get("groups",[]))!=m["n_groups"]:
            raise ValueError(f"Incomplete gzip grouped response: {path}")
        shard_station={}
        seen=set()
        for row in doc["groups"]:
            if row.get("source")!=SOURCE[ch] or row.get("station")!=station:
                raise ValueError(f"Unexpected source/station mixed in {path}")
            cfg=row["config"]
            if cfg not in CONFIGS:
                raise ValueError(f"Unexpected Garfield geometry in {path}")
            kind=which_group(row,station,det)
            if kind is None:continue
            group_key=reduced_key(row,kind)
            if group_key in seen:raise ValueError(f"Duplicate {group_key} in {path}")
            seen.add(group_key)
            plotter.validate_group(row)
            restore(agg[group_key],row)
            if kind=="station":
                shard_station[cfg]=int(row["segments"])
        if set(shard_station)!=set(CONFIGS) or \
           any(n!=m["n_muon_gas_hits"] for n in shard_station.values()):
            raise ValueError(f"Muon exposure mismatch across six LUTs: {path}")
        total_gas_segments+=m["n_muon_gas_hits"]
        output["input_shard_sha256"][f"chunk{i:04d}"]=m["gzip_sha256"]
    # Publication-quality internal conservation, with no fake independent-electron errors.
    station_muon_counts=[]
    primary_counts=[]
    for cfg in CONFIGS:
        target=output["data"][cfg]
        for kind in ("station","chamber"):
            b=agg.get((cfg,kind,-1))
            if b is None:raise ValueError(f"Missing {kind}: {cfg}")
            target[kind]=convert(b)
            plotter.validate_group(target[kind])
        pitch=float(cfg.split("_")[0][1:].replace("p","."))
        expected_wires=n_anodes(station,pitch)
        for wire in range(expected_wires):
            b=agg.get((cfg,"anodes",wire))
            if b is None:continue
            r=convert(b)
            if r["known"]<0:raise ValueError("Negative anode known-time count")
            target["anodes"].append({
                 "anode":wire,"known":r["known"],"hist":r["hist"],
                 "mean_ns":r["mean_ns"],"rms_ns":r["rms_ns"]
            })
        for axis in ("u","w"):
            for level in ("station","chamber"):
                kind=f"theta_{axis}_{level}"
                for abin in range(len(ANGLES)-1):
                    b=agg.get((cfg,kind,abin))
                    if b is not None:
                        r=convert(b)
                        values={k:r[k] for k in (
                           "segments","angle_n","angle_hist","phi_n","phi_hist"
                        )}
                    else:
                        values={"segments":0,"angle_n":0,"angle_hist":[0]*18,
                                "phi_n":0,"phi_hist":[0]*36}
                    target[kind][str(abin)]=values
                rows=target[kind].values()
                for key,full_key in (("segments","segments"),("angle_n","angle_n"),
                                     ("phi_n","phi_n")):
                    if sum(x[key] for x in rows)!=target[level][full_key]:
                        raise ValueError(f"Wrong {kind} conservation {key}: {cfg}")
        chamber=target["chamber"]
        if (sum(x["known"] for x in target["anodes"])!=chamber["known"] or
            [int(sum(x["hist"][j] for x in target["anodes"]))
                 for j in range(31)]!=chamber["hist"]):
            raise ValueError(f"Physical anode drift-hist conservation failed: {cfg}")
        station_row=target["station"]
        if station_row["segments"]!=total_gas_segments:
            raise ValueError(f"Expected {total_gas_segments} muon gas segments "
                             f"but {cfg} has {station_row['segments']}")
        station_muon_counts.append(station_row["segments"])
        primary_counts.append(station_row["primaries"])
    if len(set(station_muon_counts))!=1 or len(set(primary_counts))!=1:
        raise ValueError("Six configurations use different parent/exposure")
    output["n_muon_gas_segments_station"]=total_gas_segments
    output["n_heed_primary_electrons_station"]=primary_counts[0]
    output["garfield_trials_per_box"]=20
    output["garfield_lut_variant_count"]=6

    out.mkdir(parents=True,exist_ok=True)
    jsonpath=out/"selected_data.json"
    temporary=jsonpath.with_name(jsonpath.name+".incomplete")
    temporary.write_text(json.dumps(output,separators=(",",":"),
                 ensure_ascii=False,allow_nan=False),encoding="utf-8")
    os.replace(temporary,jsonpath)

    stats=out/"simulation_statistics.csv"
    with stats.open("w",newline="",encoding="utf-8") as stream:
        w=csv.writer(stream)
        w.writerow(("source","station","det_id","pitch_mm","field_wire_um",
            "generated_parent_events","gas_crossings_station",
            "gas_crossings_chamber","heed_primary_electrons_station",
            "heed_primary_electrons_chamber","anode_collected_station",
            "anode_collected_chamber","collection_fraction_station",
            "drift_known_station","mean_drift_ns_station",
            "drift_rms_ns_station","garfield_trials_per_cell"))
        for cfg in CONFIGS:
            row=output["data"][cfg]
            pitch=float(cfg.split("_")[0][1:].replace("p","."))
            fw=75 if "field" in cfg else 0
            st=row["station"];cham=row["chamber"]
            w.writerow((ch,output["station"],det,pitch,fw,2000000,
                st["segments"],cham["segments"],st["primaries"],
                cham["primaries"],st["collected"],cham["collected"],
                st["collected"]/st["primaries"],st["known"],
                st["mean_ns"],st["rms_ns"],20))
    # Reuse all prior graph definitions, but dynamic 2M-parent text.
    for level in ("station","chamber"):
        plotter.plot_six_times(output,level,out)
    plotter.plot_overlay_times(output,out)
    plotter.plot_numbered_anode_means(output,out)
    plotter.plot_anode_time_heatmaps(output,out)
    for axis in ("u","w"):
        plotter.plot_angle_heatmaps(output,axis,out)
        plotter.plot_angle_overlays(output,axis,out)
        plotter.plot_azimuth_heatmaps(output,axis,out)
    plotter.plot_station_chamber_comparison(output,out)
    (out/"README.txt").write_text(
        f"Source: {ch}, station MS{output['station']}, chamber DET{det}\n"
        f"Exactly 40 completed SHA256-verified 50k-parent chunks = 2M generated "
        f"{ch} parents.\n"
        f"Actual station muon gas crossings: {total_gas_segments:,}; "
        f"Heed primary electrons: {primary_counts[0]:,}.\n"
        "Each of six configurations processes the SAME Heed births and "
        "Geant4 muon gas segments using 20-trial-per-cell Garfield LUTs.\n"
        "Angles describe the final microscopic electron step near a wire; "
        "azimuth and local-normal incidence are different observables.\n"
        "Gas crossings are not unique muons; electron samples from reused "
        "finite LUTs are not independent Garfield microscopic trials.\n"
        "This is read-only analysis. No new Geant4/Heed/Garfield physics was run.\n"
        "No parent-event bootstrap uncertainties, avalanches or ASIC responses.\n",
        encoding="utf-8")
    print(f"PASS: ALL 40 exact {ch} chunks, MS{output['station']}, 2M parents, "
          f"{total_gas_segments:,} gas crossings, "
          f"{primary_counts[0]:,} Heed primary electrons")
    print("REPORT (ZIP-ready images/JSON/CSV):",out,flush=True)
    return output

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("check","plots"))
    parser.add_argument("--channel",choices=CHANNELS,default="Jpsi")
    parser.add_argument("--station",type=int,choices=range(6),default=4)
    parser.add_argument("--det",type=int,default=143)
    parser.add_argument("--out",type=Path,default=None)
    a=parser.parse_args()
    out=a.out or (BASE/f"{a.channel}_MS{a.station}_DET{a.det}")
    build_report(a.channel,a.station,a.det,out,precheck=a.mode=="check")

if __name__=="__main__":
    main()

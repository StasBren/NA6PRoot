#!/usr/bin/env python3
"""Build a real-data offline single-muon explorer from historical Stage-C crossings.

Each gas-chamber crossing is a straight-ray projection of a stored Geant4
reference-plane state, NOT a new native gas hit. Group by the exact
(source,job,event,track,pdg), never just event or spatial proximity.

Optionally realize Heed cluster ionization + SAVED Garfield single-electron
transport map on a FIXED finite nominal 4-mm grid independently for the
selected crossing rows. This is one REPRODUCIBLE stochastic realization,
NOT experimentally fired electronics or a unique prediction.

A deterministic source/job/event SHA256 sampler selects a tiny fraction of
the historical projected crossing CSVs; this does not rerun Geant4 or Garfield.
There is no need to resimulate all ~2.5m chamber crossings to inspect tracks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict
from collections import defaultdict
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_fixed_wire_occupancy import (
    cluster_segment, cluster_native_gas_segment,
    collected_electron_endpoints, n_wires, nominal_wire_index,
)
from pi_anode_digitizer import (
    AnodeConfig, digitize_hit, empirical_gains,
)
from run_pi_per_chamber_response import read_garfield, read_heed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
GEOMETRY = build_chambers()
LOOKUP = {c.det_id: c for station in GEOMETRY for c in station}


def stable_u64(*parts):
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "big")


def nominal_geometry_anodes(row):
    station = int(row["station"])
    n = n_wires(station)
    if row.get("native_hit") == "1":
        local=np.array([float(row["u_in_mm"]),float(row["u_out_mm"])])
    else:
        local = np.array([
            (float(row["y_in_cm"]) - float(row["chamber_center_y_cm"])) * 10.,
            (float(row["y_out_cm"]) - float(row["chamber_center_y_cm"])) * 10.,
        ])
    a, b = nominal_wire_index(local, station, 0.0)
    return list(range(max(0, min(int(a), int(b))),
                      min(n - 1, max(int(a), int(b))) + 1))


def get_models(args):
    if args.geometry_only:
        return None
    for p in (args.heed_tracks, args.heed_clusters, args.garfield):
        if not p.is_file():
            raise FileNotFoundError(
                f"Microscopic response library missing: {p}. "
                "Use --geometry-only to create a track-only browser."
            )
    lam, marks = read_heed(args.heed_tracks, args.heed_clusters)
    transport = read_garfield(args.garfield)
    return lam, marks, transport


def describe_row(row, models, analog=False, gain_bank=None):
    station, det = int(row["station"]), int(row["det_id"])
    if LOOKUP[det].station != station:
        raise ValueError("Inconsistent detector/station in projected crossings")
    x0, x1 = float(row["x_in_cm"]), float(row["x_out_cm"])
    y0, y1 = float(row["y_in_cm"]), float(row["y_out_cm"])
    z0, z1 = float(row["z_in_cm"]), float(row["z_out_cm"])
    x = (x0+x1)/2
    y = (y0+y1)/2
    z = (z0+z1)/2
    p = float(row["p_GeV"])
    if not all(math.isfinite(v) for v in (x,y,z,p)) or p<=0:
        raise ValueError("Nonfinite projected crossing or invalid momentum")
    result = {
        "st":station, "det":det, "layer":row["layer"], "x":round(x,5),
        "y":round(y,5), "z":round(z,5),
        "xIn":round(x0,5), "yIn":round(y0,5), "zIn":round(z0,5),
        "xOut":round(x1,5), "yOut":round(y1,5), "zOut":round(z1,5),
        "p":round(p,5),
        "thetaU":round(float(row["theta_u_deg"]),4),
        "thetaW":round(float(row["theta_w_deg"]),4),
        "length":round(float(row["gas_path_mm"]),4),
        "geom":nominal_geometry_anodes(row),
        "charges":None, "lostElectrons":None, "nClusters":None,
        "anodeSignals":None,
        "isRealGeantGasHit":row.get("native_hit")=="1",
        "edep_eV":float(row["edep_eV"]) if row.get("native_hit")=="1" else None,
        "time_ns":float(row["time_ns"]) if row.get("native_hit")=="1" else None,
    }
    if models is not None:
        lam, marks, transport = models
        seed = stable_u64("track-explorer-Heed-Garfield-v1",
                          row["source"], row["job"], row["event"],
                          row["track"], station, det)
        rng = np.random.default_rng(seed)
        sampler=(cluster_native_gas_segment
                 if row.get("native_hit")=="1" else cluster_segment)
        u, v, ncl, _ = sampler(row, rng, lam, marks)
        wire_ids, arrival_times, lost = collected_electron_endpoints(
            u,v,station,0.,transport,rng.random(u.size)
        )
        nwire=n_wires(station)
        counts=np.bincount(wire_ids,minlength=nwire)
        early_mask=(np.isfinite(arrival_times) &
                    (arrival_times >= 0.) & (arrival_times <= 100.))
        early=np.bincount(wire_ids[early_mask],minlength=nwire)
        if np.any(early > counts):
            raise ValueError("Drift-time count cannot exceed collected charge")
        if analog:
            anode_rng=np.random.default_rng(stable_u64(
                "PI_ANODE_SURROGATE_V1",row["source"],row["job"],
                row["event"],row["track"],station,det,
            ))
            endpoints=zip(wire_ids.tolist(),arrival_times.tolist())
            result["anodeSignals"]=digitize_hit(
                endpoints,anode_rng,AnodeConfig(),gain_bank
            )
        result["charges"] = [
            [int(j),int(counts[j]),int(early[j])]
            for j in np.flatnonzero(counts)
        ]
        result["lostElectrons"] = int(lost)
        result["nClusters"] = int(ncl)
    return result


def select_rows(args):
    # A deterministic EVENT-level sample (not track-level!) allows us to
    # recover mu+ and mu- from the SAME generated parent when both appear.
    # Store no row belonging to unselected events.
    samples = {src: defaultdict(lambda: defaultdict(list)) for src in SOURCES}
    counts = {}
    for src in SOURCES:
        filename = args.crossings / f"crossings_{src}_ms2_530.csv"
        if not filename.is_file():
            raise FileNotFoundError(filename)
        total,chosen = 0,0
        with filename.open(newline="",encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            needed = {
                "source","job","event","track","pdg","station","det_id",
                "x_in_cm","y_in_cm","z_in_cm",
                "x_out_cm","y_out_cm","z_out_cm",
                "chamber_center_y_cm","theta_u_deg","theta_w_deg",
                "gas_path_mm","p_GeV","layer",
            }
            if not needed.issubset(reader.fieldnames or []):
                raise ValueError(f"{filename} missing {needed-set(reader.fieldnames or [])}")
            for row in reader:
                total += 1
                if row["source"] != src:
                    raise ValueError(f"Unexpected source in {filename}")
                ev = row["job"],row["event"]
                if stable_u64(src,*ev) % args.sample_modulo != 0:
                    continue
                if abs(int(row["pdg"]))!=13:
                    continue
                key = (row["track"],row["pdg"])
                samples[src][ev][key].append(row)
                chosen += 1
        counts[src] = (total,chosen,len(samples[src]))
        print(f"{src}: read {total:,} projected gas crossings; "
              f"selected {len(samples[src])} events / {chosen} crossings",
              flush=True)
    return samples, counts


def build(args):
    if args.sample_modulo < 1 or args.max_events_per_source < 1:
        raise ValueError("Sampling controls must be positive")
    samples,stats = select_rows(args)
    if not all(samples[src] for src in SOURCES):
        raise RuntimeError(
            "No sampled events for one channel. Decrease --sample-modulo "
            "(e.g. 512) to increase the demonstrator size."
        )
    # Defer loading existing Garfield lookup until after selection; no new
    # microscopic transport is run.
    models = get_models(args)
    analog=bool(models is not None and not getattr(
        args, "no_anode_digitizer",False
    ))
    gain_path=getattr(args,"anode_gain_csv",None)
    gain_bank=empirical_gains(gain_path) if (analog and gain_path) else None
    sources = {}
    for src in SOURCES:
        ranked=[]
        for (job,event),group in samples[src].items():
            coverage = len({int(row["station"])
                            for crossings in group.values() for row in crossings})
            ncharges = len({int(pdg) for track,pdg in group})
            score = (coverage, ncharges, sum(map(len, group.values())))
            ranked.append((score,stable_u64("display",src,job,event),(job,event),group))
        # Curated interesting event selection for demonstration, NOT unbiased
        # statistics. Rank by number of hit stations, then both muon signs.
        ranked.sort(key=lambda entry:(-entry[0][0],-entry[0][1],
                                       -entry[0][2],entry[1]))
        events=[]
        for score,_,(job,event),group in ranked[:args.max_events_per_source]:
            tracks=[]
            for (track,pdg),rows in sorted(group.items()):
                points=[describe_row(row,models,analog,gain_bank)
                        for row in rows]
                points.sort(key=lambda p:(p["z"],p["det"]))
                if points:
                    tracks.append({
                        "id":str(track),"pdg":int(pdg),
                        "points":points,
                        "stations":sorted(set(p["st"] for p in points)),
                    })
            if tracks:
                events.append({
                    "job":str(job),"event":str(event),"tracks":tracks,
                    "hitStations":score[0],
                })
        if not events:
            raise RuntimeError(f"No displayable events: {src}")
        sources[src]=events
    payload={
        "sources":sources,
        "geometry":[
            [
                {"det":c.det_id,"x":round(c.x,4),"y":round(c.y,4),
                 "st":c.station,"layer":c.layer,
                 "row":c.row,"col":c.col}
                for c in station
            ] for station in GEOMETRY
        ],
        "meta":{
            "dataKind":getattr(args,"data_kind","historical"),
            "sample":(
                "CURATED true current-geometry Geant4 gas hits: parent chunk, "
                "MCKine-associated muon PDG, actual module gas entry/exit and Edep"
            ) if getattr(args,"data_kind","historical")=="native" else
                "SHA256 deterministic event sample of historical MS2=530cm gas projections; display is curated, NOT a random detector-physics sample",
            "sampleModulo":args.sample_modulo,
            "chargeModel":models is not None,
            "anodeSignalModel":analog,
            "anodeModel":(
                "STUDY ONLY: Stage-A unconditioned empirical gas gain"
                if gain_bank is not None else
                "STUDY ONLY: uncalibrated Polya gas gain plus zero-gain fraction"
            ) if analog else None,
            "anodeElectronicsParameters":asdict(AnodeConfig()) if analog else None,
            "anodeWaveformPhysics":(
                "PARAMETRIC anode fast electron + slow positive-ion current, "
                "nominal CR-RC amplifier and fixed seeded white-noise samples. "
                "NOT true anode Shockley-Ramo field or prototype calibration."
            ) if analog else None,
            "nSources":len(SOURCES),
            "stationZcm":[300.,340.,530.,590.,810.,850.],
            "geometryModel":"nominal checkerboard, finite MS0 77 wires / other stations 128, pitch 4mm",
            "trajectoryModel":(
                "Lines JOIN true native Geant4 gas-hit midpoints; connections "
                "between chambers are only visual guides, not step histories."
            ) if getattr(args,"data_kind","historical")=="native" else
                "line segments joining projected gas-crossing midpoints, NOT true Geant4 steps or exact bending between stations",
            "electronModel":"one deterministic event-seeded Heed Poisson/cluster and stored Garfield single-electron transport realization, B=0 reference, 6mm gap mapping",
            "electronicsModel":"NOT implemented: no avalanche gain, induced-current weighting field, shaping, threshold, cross talk or noise",
            "meaningMissingStation":"No projected GAS crossing for selected track; does not prove it did not physically traverse the station plane",
        },
    }
    marker="__PI_TRACK_DATA__"
    template=args.template.read_text(encoding="utf-8")
    if template.count(marker)!=1:
        raise ValueError("Track viewer template marker must occur once")
    data=json.dumps(payload,ensure_ascii=False,separators=(",",":")).replace("</","<\\/")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(template.replace(marker,data),encoding="utf-8")
    print(f"PASS: Track Explorer: {sum(map(len,sources.values()))} real "
          f"source/job/event groups, Garfield charge model={models is not None}, "
          f"anode surrogate waveforms={analog}")
    print(f"Open: {args.out.resolve()} ({args.out.stat().st_size/1048576:.1f} MiB)")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crossings",type=Path,default=Path(
        "test_runs/mwpc_stage_c/checkerboard_full"))
    ap.add_argument("--heed-tracks",type=Path,default=Path(
        "test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_tracks.csv"))
    ap.add_argument("--heed-clusters",type=Path,default=Path(
        "test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_clusters.csv"))
    ap.add_argument("--garfield",type=Path,default=None)
    ap.add_argument("--data-kind",choices=("historical","native"),
                    default="historical")
    ap.add_argument("--sample-modulo",type=int,default=2048)
    ap.add_argument("--max-events-per-source",type=int,default=44)
    ap.add_argument("--no-anode-digitizer",action="store_true",
                    help="Keep Garfield collected electron counts but no "
                         "parametric analog anode waveforms")
    ap.add_argument("--anode-gain-csv",type=Path,default=None,
                    help="Unconditioned Phase A single seed gain ensemble "
                         "with collected_on_wire column (including zeros). "
                         "Otherwise a named UNCALIBRATED Polya model is used.")
    ap.add_argument("--geometry-only",action="store_true",
                    help="Avoid loading Heed/Garfield and show geometric anode cells only")
    ap.add_argument("--template",type=Path,default=HERE/"pi_muon_track_explorer_template.html")
    ap.add_argument("--out",type=Path,default=Path(
        "test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_muon_track_explorer.html"))
    args=ap.parse_args()
    if args.garfield is None:
        a=Path("test_runs/mwpc_stage_c/"
               "c_transport_box_map_hybrid_boundary_2p4/"
               "c_transport_box_electrons.csv")
        b=Path("test_runs/mwpc_stage_c/"
               "c_transport_box_map_gm2p0_gp4p0/"
               "c_transport_box_electrons.csv")
        args.garfield=a if a.is_file() else b
    build(args)


if __name__=="__main__":
    main()

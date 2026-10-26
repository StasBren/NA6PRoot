#!/usr/bin/env python3
"""LUT500 J/psi parent-event up/down asymmetry for all six stations / DETs.

Independent, non-destructive sidecar for the full 2M-parent LUT20 dashboard.
Uses existing SHA-checked native 50k-parent crossing CSVs and p3-only
900-cell/500-trial Garfield maps. No Geant4/Heed/Garfield re-simulation:
the original stable-seeded Poisson-Heed sampler is used.

Parent (source, chunk, event) counts are summed per STATION and separately
per physical DET. Each parent/DET is a separate statistical unit; crossing
moments are never substituted for parent-event moments. A undefined if
N_up+N_down=0; such events are not recorded as A=0.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import time
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import build_chambers
from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord
from pi_six_geometry_native_timing import (
    HEED_TRACKS,HEED_CLUSTERS,birth_sample,load_map,
    response_for_births,stable_seed)
from pi_pitch_fieldwire_pilot import tag
from pi_lut_20_vs_500_convergence import ROOT500
from pi_event_asymmetry_500_ms0 import up_down
from run_pi_per_chamber_response import read_heed
from run_pi_geant4_campaign import REPO

HERE=Path(__file__).resolve().parent
TEMPLATE=HERE/"pi_lut500_event_asymmetry_root_template.html"
OUT=REPO/"test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1"
SOURCE="jpsi40-bore"
PITCH=3.
CONFIGS=(tag(PITCH,0.),tag(PITCH,75.))
N_PARENTS=50_000
A_EDGES=np.linspace(-1.,1.,42)
SCHEMA="PI_JPSI_SIX_STATIONS_PER_DET_PARENT_ASYM_P3_LUT500_ROOT_V1"
CSV_HEADER=("source","chunk","event","station","det_id","gas_crossings",
            "up_noFW","down_noFW","up_FW75","down_FW75")


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


class Distribution:
    """Equal-weight parent-event A distribution, NO electron weighting."""
    def __init__(self):
        self.parents_with_crossing=0
        self.valid=0
        self.undefined=0
        self.sum_a=0.
        self.sum_a2=0.
        self.sum_up=0
        self.sum_down=0
        self.hist=np.zeros(41,dtype=np.int64)

    def add(self,up,down):
        u=int(up);d=int(down)
        if min(u,d)<0:raise ValueError("Negative classified electron count")
        self.parents_with_crossing+=1
        self.sum_up+=u
        self.sum_down+=d
        if u+d==0:
            self.undefined+=1
            return
        self.valid+=1
        a=(u-d)/(u+d)
        j=min(40,max(0,int(np.searchsorted(A_EDGES,a,side="right")-1)))
        self.hist[j]+=1
        self.sum_a+=a
        self.sum_a2+=a*a

    def result(self):
        if self.parents_with_crossing!=self.valid+self.undefined or \
                int(self.hist.sum())!=self.valid:
            raise AssertionError("Event-level A conservation failure")
        mean=self.sum_a/self.valid if self.valid else None
        rms=math.sqrt(max(0.,self.sum_a2/self.valid-mean*mean)) \
            if self.valid else None
        return dict(parents_with_gas_crossing=self.parents_with_crossing,
                    valid_A=self.valid,undefined_A=self.undefined,
                    total_up=self.sum_up,total_down=self.sum_down,
                    mean_A=mean,rms_A=rms,hist=self.hist.tolist(),
                    electron_weighted_A=(
                        (self.sum_up-self.sum_down)/(self.sum_up+self.sum_down)
                        if self.sum_up+self.sum_down else None))


class Paired:
    def __init__(self):
        self.both_valid=0
        self.sum_delta=0.
        self.sum_delta2=0.
    def add(self,v):
        u0,d0,u1,d1=[int(x) for x in v]
        if u0+d0==0 or u1+d1==0:return
        a0=(u0-d0)/(u0+d0)
        a1=(u1-d1)/(u1+d1)
        delta=a1-a0
        self.both_valid+=1
        self.sum_delta+=delta
        self.sum_delta2+=delta*delta
    def result(self):
        n=self.both_valid
        mean=self.sum_delta/n if n else None
        return {
            "both_valid_parents":n,
            "mean_delta_FW_minus_noFW":mean,
            "rms_delta":(
                math.sqrt(max(0.,self.sum_delta2/n-mean*mean))
                if n else None)}


def update(accs,paired,counts):
    u0,d0,u1,d1=[int(x) for x in counts]
    accs[CONFIGS[0]].add(u0,d0)
    accs[CONFIGS[1]].add(u1,d1)
    paired.add(counts)


def fresh_acc():
    return {cfg:Distribution() for cfg in CONFIGS}


def summary_for_chunks(chunks,stations,max_crossings,out,make_html=True):
    if not 1<=chunks<=40 or any(st not in range(6) for st in stations) \
            or not stations or max_crossings<0:
        raise ValueError("Use chunks=1..40, station subset 0..5, max_crossings>=0")
    out=Path(out)
    targets=[out/"summary.json",out/"dashboard.html",
             out/"event_counts.csv.gz",out/"manifest.json"]
    if any(p.exists() for p in targets):
        raise FileExistsError("Dashboard output already exists: "+str(out))
    started=time.perf_counter()
    maps={};metadata={}
    for diameter,key in zip((0.,75.),CONFIGS):
        maps[key],metadata[key]=load_map(PITCH,diameter,ROOT500)
        if maps[key].shape!=(900,500,5):
            raise ValueError("Expected immutable 900 x 500 x 5 Garfield source")
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    layout=defaultdict(list)
    det_ids=defaultdict(set)
    for group in build_chambers():
        for c in group:
            if c.station not in stations:continue
            if c.det_id in det_ids[c.station]:
                raise ValueError("Duplicate physical DET in layout")
            det_ids[c.station].add(c.det_id)
            layout[c.station].append(dict(id=int(c.det_id),row=int(c.row),
                 col=int(c.col),layer=str(c.layer),
                 x=float(c.x),y=float(c.y),
                 hx=float(c.half_x),hy=float(c.half_y)))
    for st in stations:
        if not det_ids[st]:raise ValueError("Missing physical DET layout "+str(st))
    station_acc={st:fresh_acc() for st in stations}
    station_pair={st:Paired() for st in stations}
    chamber_acc={st:{d:fresh_acc() for d in det_ids[st]} for st in stations}
    chamber_pair={st:{d:Paired() for d in det_ids[st]} for st in stations}
    total_crossings=Counter()
    total_electrons=Counter()
    total_boundary=Counter()
    input_provenance=[]
    out.mkdir(parents=True,exist_ok=True)
    raw_tmp=out/"event_counts.csv.gz.incomplete"
    if raw_tmp.exists():raise FileExistsError("Incomplete old output: "+str(raw_tmp))
    if make_html and not TEMPLATE.is_file():raise FileNotFoundError(TEMPLATE)
    try:
        with gzip.open(raw_tmp,"wt",encoding="utf-8",newline="",
                       compresslevel=4) as gz:
            writer=csv.writer(gz)
            writer.writerow(CSV_HEADER)
            for chunk in range(1,chunks+1):
                path,info=crossing_provenance("Jpsi",chunk)
                input_provenance.append({"chunk":chunk,"crossing_sha256":
                                info["csv_sha256"],"csv":str(path)})
                station_counts={st:np.zeros((N_PARENTS,4),dtype=np.int64)
                                for st in stations}
                station_hits={st:np.zeros(N_PARENTS,dtype=np.int32)
                              for st in stations}
                det_counts={st:defaultdict(dict) for st in stations}
                det_hits={st:defaultdict(Counter) for st in stations}
                processed=Counter()
                with path.open(newline="",encoding="utf-8") as f:
                    for row in csv.DictReader(f):
                        st=int(row["station"])
                        if st not in stations or abs(int(row["pdg"]))!=13:
                            continue
                        if max_crossings and processed[st]>=max_crossings:
                            continue
                        det=int(row["det_id"])
                        if det not in det_ids[st]:
                            raise ValueError(f"Unknown physical DET {det} in MS{st}")
                        ev=int(row["event"])
                        if not 0<=ev<N_PARENTS:
                            raise ValueError("Invalid local Jpsi parent event")
                        station_hits[st][ev]+=1
                        det_hits[st][det][ev]+=1
                        length=float(row["path_length_mm"])
                        if length==0.:
                            if not coincident_zero_chord(row):
                                raise ValueError("Invalid zero-length native crossing")
                            positions=np.empty((0,3))
                            q=np.empty(0)
                        else:
                            seed=stable_seed(
                              "native-six-pitch-shared-Heed-v1",SOURCE,
                              f"chunk{chunk:04d}",row["event"],
                              row["hit_index"],row["track_id"])
                            positions,q=birth_sample(row,lam,marks,seed)
                        total_electrons[st]+=len(q)
                        values=[]
                        for cfg in CONFIGS:
                            finite,rawcol,wires,samples=response_for_births(
                                positions[:,0],positions[:,1],q,PITCH,st,maps[cfg])
                            up,down,b=up_down(samples,finite)
                            if up+down+b!=int(np.count_nonzero(finite)):
                                raise AssertionError("Non-conserved anode population")
                            total_boundary[st,cfg]+=b
                            values.extend((up,down))
                        station_counts[st][ev,:]+=values
                        old=det_counts[st][det].get(ev)
                        if old is None:
                            det_counts[st][det][ev]=list(values)
                        else:
                            for j,value in enumerate(values):old[j]+=value
                        processed[st]+=1
                        total_crossings[st]+=1
                        if sum(processed.values())%20000==0:
                            print("PROGRESS source Jpsi chunk",chunk,
                                  "gas crossings",sum(processed.values()),
                                  flush=True)
                for st in stations:
                    # Check event-wise conservation of module -> station counts.
                    check=np.zeros_like(station_counts[st])
                    for det,entries in det_counts[st].items():
                        for ev,counts in entries.items():
                            check[ev,:]+=counts
                    if not np.array_equal(check,station_counts[st]):
                        raise AssertionError(f"DET-to-station parent electron conservation MS{st}")
                    if sum(sum(h.values()) for h in det_hits[st].values()) != \
                            int(station_hits[st].sum()):
                        raise AssertionError("DET crossing sum mismatch")
                    for ev in np.flatnonzero(station_hits[st]):
                        row=station_counts[st][ev]
                        update(station_acc[st],station_pair[st],row)
                        writer.writerow((SOURCE,chunk,int(ev),st,-1,
                                         int(station_hits[st][ev]),*map(int,row)))
                    for det,hits in det_hits[st].items():
                        for ev,n_cross in hits.items():
                            vals=det_counts[st][det].get(ev,[0,0,0,0])
                            update(chamber_acc[st][det],chamber_pair[st][det],vals)
                            writer.writerow((SOURCE,chunk,int(ev),st,det,
                                             int(n_cross),*map(int,vals)))
                    print(f"PASS chunk {chunk}/{chunks} MS{st}: "
                          f"{processed[st]:,} muon gas crossings, "
                          f"{int(np.count_nonzero(station_hits[st])):,} "
                          "parents with gas crossing",flush=True)
        os.replace(raw_tmp,out/"event_counts.csv.gz")
        stations_result={}
        for st in stations:
            station_val={cfg:acc.result() for cfg,acc in station_acc[st].items()}
            chamber_val={}
            for det in sorted(det_ids[st]):
                dist={cfg:a.result() for cfg,a in chamber_acc[st][det].items()}
                for cfg in CONFIGS:
                    if station_val[cfg]["valid_A"]<dist[cfg]["valid_A"]:
                        raise AssertionError("Per-chamber valid parent > station valid")
                chamber_val[str(det)]={"configurations":dist,
                                       "paired":chamber_pair[st][det].result()}
            for cfg in CONFIGS:
                # Counts additive; parent A hist and mean/rms NOT additive.
                for field in ("total_up","total_down"):
                    if sum(row["configurations"][cfg][field]
                           for row in chamber_val.values())!=station_val[cfg][field]:
                        raise AssertionError("Station vs chamber electron totals differ")
            stations_result[str(st)]={
                "gas_crossings":int(total_crossings[st]),
                "heed_primaries":int(total_electrons[st]),
                "layout":layout[st],"station_parent":{
                    "configurations":station_val,"paired":station_pair[st].result()},
                "chambers":chamber_val,
                "n_physical_chambers":len(det_ids[st])}
        result={
          "schema":SCHEMA,"source":SOURCE,"pitch_mm":PITCH,"B_t":0,
          "generated_parent_slots_per_chunk":N_PARENTS,
          "number_of_source_chunks":chunks,
          "parents_generated_in_input_chunks":chunks*N_PARENTS,
          "partial_exposure":bool(max_crossings),
          "max_crossings_per_station_and_chunk":max_crossings,
          "n_stations_with_computed_results":len(stations),
          "computed_stations":list(stations),
          "configs":[dict(key=CONFIGS[0],label="No field wire",
                          field_wire_diameter_um=0),
                     dict(key=CONFIGS[1],label="Field wire 75 um",
                          field_wire_diameter_um=75)],
          "garfield_LUT500":metadata,
          "native_csv_inputs":input_provenance,
          "Heed_library":{"tracks":str(HEED_TRACKS),"clusters":str(HEED_CLUSTERS),
               "lambda_per_mm":float(lam),"mean_cluster_electrons":
                               float(np.mean(marks))},
          "asymmetry":{"definition":"A=(N_up-N_down)/(N_up+N_down)",
              "statistical_unit":"generated Jpsi parent per station OR per physical DET",
              "gas_crossings":"one physical chamber muon traversal",
              "up":"0<phi_hit_deg<180 at anode surface",
              "down":"-180<phi_hit_deg<0 at anode surface",
              "boundary":"phi_hit exactly 0 or +/-180 excluded",
              "undefined_policy":"N_up+N_down=0 => undefined, not A=0",
              "valid_hist_bins":41,"hist_edges":A_EDGES.tolist(),
              "mean_RMS":"arithmetic equal-weight parent events with finite A",
              "no_gain":True,"no_induced_signal":True,
              "no_electronics":True,"magnetic_field_t":0},
          "stations":stations_result,
          "boundary_electrons":{f"MS{st}/{cfg}":int(total_boundary[st,cfg])
                              for st in stations for cfg in CONFIGS},
          "elapsed_wall_seconds":time.perf_counter()-started,
          "limitations":[
              "LUT500 data exist for p=3mm ONLY. p=3.5/4mm plots "
              "on older dashboard use DIFFERENT LUT20 and are not mixed.",
              "The 2M Jpsi dashboard remains a separate LUT20 data product.",
              "Parent event A for a station sums classified electrons from "
              "all physical DETs traversed by that parent in that station.",
              "Parent event A for a DET is computed from only that DET, "
              "not averaged over gas crossings.",
              "A station histogram CANNOT be calculated by adding per-DET "
              "A histograms; counts must first be pooled per parent.",
              "Missing/undefined A never becomes numerical zero.",
              "One Jpsi parent can contribute to multiple stations/DETs; "
              "they are not independent across views.",
              "Heed responses draw finite 500-trial LUT electron outcomes; "
              "these are not newly independent Garfield drifts.",
              "No avalanche, weighting-field signal, or electronics."
          ]
        }
        # Require no partial source data when labeling output complete.
        if not max_crossings and not all(stations_result[str(st)]["gas_crossings"]
                                         for st in stations):
            raise ValueError("No station crossing data")
        (out/"summary.json").write_text(
            json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
            encoding="utf-8")
        if make_html:
            template=TEMPLATE.read_text(encoding="utf-8")
            if template.count("__PAYLOAD__")!=1:
                raise ValueError("ROOT HTML template must have one payload marker")
            payload=json.dumps(result,ensure_ascii=False,separators=(",",":"),
                               allow_nan=False).replace("</","<\\/")
            page=template.replace("__PAYLOAD__",payload)
            (out/"dashboard.html").write_text(page,encoding="utf-8")
        manifest={"schema":SCHEMA,
                  "complete":not bool(max_crossings) and stations==list(range(6)),
                  "computed_stations":list(stations),
                  "source_chunks":chunks,
                  "files_sha256":{p.name:sha(p) for p in targets[:-1] if p.exists()},
                  "source_csv_sha256":[x["crossing_sha256"]
                                       for x in input_provenance]}
        (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",
                                         encoding="utf-8")
        print(f"PASS ROOT-style eventwise LUT500 dashboard: {out/'dashboard.html'}")
        print(f"PASS parent+DET event counts: {out/'event_counts.csv.gz'}")
        return result
    finally:
        # A rejected or crashed run leaves only an explicit incomplete gzip.
        # We do not remove it silently; inspect the partial data before retry.
        pass


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunks",type=int,default=1)
    ap.add_argument("--stations",type=int,nargs="+",default=list(range(6)))
    ap.add_argument("--max-crossings",type=int,default=200,
                    help="Max muon gas crossings per station per 50k chunk; "
                         "0 = FULL and eligible for scientific report")
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--no-html",action="store_true")
    args=ap.parse_args()
    summary_for_chunks(args.chunks,sorted(set(args.stations)),
                       args.max_crossings,args.out,not args.no_html)


if __name__=="__main__":
    main()

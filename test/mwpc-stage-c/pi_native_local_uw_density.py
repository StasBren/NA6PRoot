#!/usr/bin/env python3
"""Real native Geant4 TGeo-local muon midpoint u-w histograms.

first: existing first chunk CSV per source (3000 events/channel).
all: each immutable original ROOT chunk re-exported as temporary CSV,
     with source/code-frozen resumable chunk summaries, cross-checked
     against independent full native angular quality analysis.
No Geant4 simulation, Garfield, Heed or cathode strip model.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
import shutil
from pathlib import Path
from audit_pi_geant4_campaign import audit, DEFAULT as CAMPAIGN
from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_native_full_campaign_qa import export_chunk, file_id, ANGLE, bin_of

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=REPO/"test_runs/mwpc_stage_c"
FIRST=BASE/"pi_native_anode_first_chunk/native_exports"
OUT=BASE/"pi_native_local_uw_density"
QA=BASE/"pi_native_full_quality"
FIELDS=("theta_u_in_deg","theta_w_in_deg","theta_u_out_deg","theta_w_out_deg")
SELECTORS=("positive","far10","near2")
NU,NW=16,24
GEOM={c.det_id:c for g in build_chambers() for c in g}
CHANNELS={"jpsi40-bore":"jpsi","phi40-bore":"phi","omega40-bore":"omega"}


def resolve_root_executable(build_dirs=None,env=None):
    """Locate ROOT binary without relying on interactive-shell PATH.

    Prefer PI_ROOT_EXECUTABLE, current PATH, then ROOTSYS and the
    ROOT_DIR/CMAKE_PREFIX_PATH from EXISTING NA6PRoot CMake build.
    """
    env=os.environ if env is None else env
    explicit=env.get("PI_ROOT_EXECUTABLE","").strip()
    if explicit:
        selected=shutil.which(explicit,path=env.get("PATH",""))
        if not selected and Path(explicit).is_file() and os.access(explicit,os.X_OK):
            selected=explicit
        if not selected:
            raise FileNotFoundError(
                "PI_ROOT_EXECUTABLE set but cannot be executed: "+explicit)
        return Path(selected).resolve()
    selected=shutil.which("root",path=env.get("PATH",""))
    if selected:
        return Path(selected).resolve()
    possible=[]
    def offer(value):
        if not value:
            return
        # Existing CMake caches may store semicolon-delimited prefixes.
        for token in value.split(";"):
            token=token.strip()
            if token and "$" not in token:
                possible.append(Path(os.path.expanduser(token)))
    for name in ("ROOTSYS","ROOT_DIR","ROOT_ROOT","CMAKE_PREFIX_PATH"):
        offer(env.get(name,""))
    if build_dirs is None:
        home=Path.home()
        build_dirs=(
            REPO/"build",home/"na6p/build",
            home/"na6p/build/NA6PRoot",
            home/"na6p/build/NA6PRoot-mwpc",
            home/"na6p/src/NA6PRoot-mwpc/build",
        )
    for directory in build_dirs:
        cache=Path(directory)/"CMakeCache.txt"
        if not cache.is_file():
            continue
        for line in cache.read_text(encoding="utf-8",errors="replace").splitlines():
            if "=" not in line or line.startswith(("//","#")):
                continue
            key,value=line.split("=",1)
            field=key.split(":",1)[0]
            if field in ("ROOT_DIR","ROOTSYS","ROOT_ROOT","CMAKE_PREFIX_PATH",
                         "ROOT_INCLUDE_DIR","ROOT_LIBRARY_DIR",
                         "ROOT_BIN_DIR","ROOT_BINARY_DIR","ROOT_DIR_HINT"):
                offer(value)
    for candidate in possible:
        for prefix in (candidate,*list(candidate.parents)[:7]):
            exe=prefix/"bin/root"
            if exe.is_file() and os.access(exe,os.X_OK):
                return exe.resolve()
    raise FileNotFoundError(
        "ROOT binary not found. NA6PRoot Stage-C is installed, but "
        "ROOT's bin directory is not on PATH and the CMake cache did "
        "not locate it. Find it with "
        "'find ~/na6p -path */bin/root -type f' "
        "then launch with PI_ROOT_EXECUTABLE=/absolute/path/bin/root "
        "bash test/mwpc-stage-c/run_pi_native_local_uw.sh all."
    )


def activate_root_for_export(build_dirs=None,env=None):
    """Expose native ROOT macro, rootcling and ROOT shared libraries."""
    env=os.environ if env is None else env
    exe=resolve_root_executable(build_dirs=build_dirs,env=env)
    bindir=exe.parent
    prefix=bindir.parent
    env["PI_ROOT_EXECUTABLE"]=str(exe)
    env["ROOTSYS"]=str(prefix)
    env["PATH"]=os.pathsep.join(
        [str(bindir)]+[p for p in env.get("PATH","").split(os.pathsep)
                       if p and p!=str(bindir)])
    additions=[str(p) for p in (prefix/"lib",prefix/"lib64") if p.is_dir()]
    env["LD_LIBRARY_PATH"]=os.pathsep.join(
        additions+[p for p in env.get("LD_LIBRARY_PATH","").split(os.pathsep)
                   if p and p not in additions])
    return str(exe)


def local_bounds(ch):
    return 10*ch.half_y,10*ch.half_x


def fresh():
    return {"n":0,"outside":0,"all":{},"selected":{},
            "angles":{f:[0]*(len(ANGLE)-1) for f in FIELDS}}


def increment(grid,index,value=1):
    key=str(index)
    grid[key]=grid.get(key,0)+value


def selections(theta):
    return (("positive",theta>=0),("far10",abs(theta)>=10),
            ("near2",-2<=theta<2))


def process_csv(path):
    result={}
    anomalies=defaultdict(int)
    with path.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f)
        needed={"pdg","det_id","station","path_length_mm","edep_eV",
                "u_in_mm","u_out_mm","w_in_mm","w_out_mm",*FIELDS}
        if not needed.issubset(rd.fieldnames or []):
            raise ValueError("Native ROOT export missing real local fields: "+
                             repr(needed-set(rd.fieldnames or [])))
        for row in rd:
            if abs(int(row["pdg"]))!=13:continue
            det=int(row["det_id"]);st=int(row["station"])
            if det not in GEOM or GEOM[det].station!=st:
                raise ValueError(f"Muon gas hit in invalid DET {det}, MS{st}")
            angular=[float(row[f]) for f in FIELDS]
            pathlen=float(row["path_length_mm"])
            edep=float(row["edep_eV"])
            if not all(math.isfinite(v) for v in angular+[pathlen,edep]) or (
                pathlen<=0 or edep<0):
                anomalies["invalid_native_hit"]+=1
                continue
            r=result.setdefault(str(det),fresh());r["n"]+=1
            for field,theta in zip(FIELDS,angular):
                r["angles"][field][bin_of(theta,ANGLE)]+=1
            u=.5*(float(row["u_in_mm"])+float(row["u_out_mm"]))
            w=.5*(float(row["w_in_mm"])+float(row["w_out_mm"]))
            hu,hw=local_bounds(GEOM[det])
            if not(math.isfinite(u) and math.isfinite(w)) or (
                abs(u)>hu+1.e-5 or abs(w)>hw+1.e-5
            ):
                r["outside"]+=1
                anomalies["outside_physical_module_or_nonfinite"]+=1
                continue
            iu=max(0,min(NU-1,int(math.floor((u+hu)/(2*hu)*NU))))
            iw=max(0,min(NW-1,int(math.floor((w+hw)/(2*hw)*NW))))
            cell=iw*NU+iu
            increment(r["all"],cell)
            for field,theta in zip(FIELDS,angular):
                target=r["selected"].setdefault(field,{})
                for category,accept in selections(theta):
                    if accept:increment(target.setdefault(category,{}),cell)
    for det,r in result.items():
        if any(sum(h)!=r["n"] for h in r["angles"].values()):
            raise AssertionError("Angular counts do not conserve native hits")
        if sum(r["all"].values())!=r["n"]-r["outside"]:
            raise AssertionError("Midpoint geometry 2D count conservation")
        for field in FIELDS:
            for category in SELECTORS:
                if sum(r["selected"].get(field,{}).get(category,{}).values())>(
                    r["n"]-r["outside"]):
                    raise AssertionError("Selected 2D density exceeds total")
    return result,dict(anomalies)


def merge(total,records):
    for det,r in records.items():
        t=total.setdefault(det,fresh())
        t["n"]+=r["n"];t["outside"]+=r["outside"]
        for field,arr in r["angles"].items():
            for k,val in enumerate(arr):t["angles"][field][k]+=val
        for key,count in r["all"].items():increment(t["all"],key,count)
        for field,selectors in r["selected"].items():
            q=t["selected"].setdefault(field,{})
            for category,grid in selectors.items():
                d=q.setdefault(category,{})
                for key,count in grid.items():increment(d,key,count)


def check_full_quality(total):
    qa_counts={}
    with (QA/"pi_native_all_chamber_quality.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            qa_counts[row["source"],int(row["det_id"])]=int(row["gas_hit_rows"])
    qa_angles={}
    with (QA/"pi_native_all_histograms.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            field=row["variable"]
            if field not in FIELDS:continue
            lo=float(row["bin_lower"]);hi=float(row["bin_upper"])
            ix=next((j for j in range(len(ANGLE)-1)
                     if abs(ANGLE[j]-lo)<1e-8 and
                        abs(ANGLE[j+1]-hi)<1e-8),None)
            if ix is None:raise ValueError("Frozen native QA angles mismatch")
            key=(row["source"],int(row["det_id"]),field)
            counts=qa_angles.setdefault(key,[None]*(len(ANGLE)-1))
            if counts[ix] is not None:raise ValueError("Repeated QA angle bin")
            counts[ix]=int(row["hit_count"])
    if len(qa_counts)!=3*298 or len(qa_angles)!=3*298*4:
        raise ValueError("Full quality analysis lacks all 894 chamber records")
    for src in SOURCES:
        for det in GEOM:
            item=total[src].get(str(det),fresh())
            if item["n"]!=qa_counts[src,det]:
                raise ValueError(f"Independent full QA count disagreement "
                                 f"{src} DET{det}")
            for field in FIELDS:
                if item["angles"][field]!=qa_angles[src,det,field]:
                    raise ValueError(f"Independent full QA angle histogram "
                                     f"mismatch {src} DET{det} {field}")


def build(args):
    args.out.mkdir(parents=True,exist_ok=True)
    exporter=HERE/"export_c_geant_muon_crossings.C"
    epoch=hashlib.sha256(Path(__file__).read_bytes()+
                         exporter.read_bytes()).hexdigest()
    total={src:{} for src in SOURCES}
    anomalies=defaultdict(int)
    if args.mode=="first":
        for src,channel in CHANNELS.items():
            path=args.first/f"{channel}_geant_gas_crossings.csv"
            if not path.is_file():raise FileNotFoundError(path)
            this,bad=process_csv(path)
            merge(total[src],this)
            for key,n in bad.items():anomalies[key]+=n
            print(f"PASS first {src}: {sum(v['n'] for v in this.values()):,} "
                  "valid native muon gas hit segments",flush=True)
    else:
        status=audit(args.campaign,False)
        if status["status"]!="COMPLETE":
            raise ValueError("Frozen Geant4 campaign incomplete")
        # Fail fast before 201 exports if Stage-C is active but ROOT is not.
        root_binary=activate_root_for_export()
        print("PI native ROOT executable:",root_binary,flush=True)
        print("ROOTSYS:",os.environ["ROOTSYS"],flush=True)
        checkpoint=args.out/"checkpoints";checkpoint.mkdir(exist_ok=True)
        current=args.out/"_export_current.csv"
        log=args.out/"_export_current.log"
        for src,channel in CHANNELS.items():
            for i in range(1,68):
                chunk=args.campaign/channel/f"chunk{i:04d}"
                inputs={name:file_id(chunk/name) for name in (
                    "HitsMuonSpec.root","MCKine.root","geometry.root")}
                target=checkpoint/f"{channel}_chunk{i:04d}.json"
                if target.is_file():
                    cp=json.loads(target.read_text(encoding="utf-8"))
                    if cp.get("epoch")!=epoch or cp.get("inputs")!=inputs:
                        raise ValueError(f"Frozen input/code mismatch: {target}; "
                                         "use a new --out directory")
                    this,bad=cp["records"],cp["anomalies"]
                    state="SKIP"
                else:
                    if current.exists():current.unlink()
                    print(f"RUN {src} chunk{i:04d}/67: native ROOT export...",
                          flush=True)
                    export_chunk(root_binary,chunk,current,log)
                    this,bad=process_csv(current)
                    tmp=target.with_suffix(".json.tmp")
                    tmp.write_text(json.dumps(
                        {"epoch":epoch,"inputs":inputs,"records":this,
                         "anomalies":bad},separators=(",",":")),
                        encoding="utf-8")
                    os.replace(tmp,target)
                    current.unlink()
                    state="PASS"
                merge(total[src],this)
                for key,n in bad.items():anomalies[key]+=n
                print(f"{state} {src} chunk{i:04d}/67",flush=True)
        check_full_quality(total)
    payload={
       "meta":{
          "scope":("FIRST_CHUNK_3000_EVENTS_PER_SOURCE" if args.mode=="first"
                   else "FULL_201_CHUNKS_200000_EVENTS_PER_SOURCE"),
          "generated_events_per_source":3000 if args.mode=="first" else 200000,
          "unit":"one valid true Geant4 muon gas hit segment, not unique track",
          "coordinates":"True TGeo-local u transverse and w along wires, mm",
          "method":"One midpoint (u_in+u_out)/2,(w_in+w_out)/2 per gas hit",
          "angleCuts":"positive: theta>=0; far10: |theta|>=10; "
                      "near2: -2<=theta<2; native in/out u/w angle fields",
          "nu":NU,"nw":NW,"u_fastest":True,
          "source_epoch":epoch,
          "anomalies":dict(anomalies),
          "full_independent_qa_passed":args.mode=="all",
          "not":"Not Heed/Garfield, not incoming station muons",
       },
       "modules":{},
    }
    for src in SOURCES:
        m={}
        for det in GEOM:
            item=total[src].get(str(det),fresh())
            m[str(det)]={
                "n":item["n"],"outside":item["outside"],
                "all":[[int(k),v] for k,v in item["all"].items()],
                "selected":{
                    field:{
                        category:[[int(k),v] for k,v in item["selected"]
                             .get(field,{}).get(category,{}).items()]
                        for category in SELECTORS}
                    for field in FIELDS
                }
            }
        payload["modules"][src]=m
    output=args.out/"pi_native_local_uw_density.json"
    temp=output.with_suffix(".json.tmp")
    temp.write_text(json.dumps(payload,ensure_ascii=False,separators=(",",":")),
                    encoding="utf-8")
    os.replace(temp,output)
    print(f"PASS: local u-w midpoint maps built from "
          f"{sum(v['n'] for m in total.values() for v in m.values()):,} "
          f"valid native muon hits; scope={args.mode}")
    print("Output:",output)
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("first","all"))
    p.add_argument("--first",type=Path,default=FIRST)
    p.add_argument("--campaign",type=Path,default=CAMPAIGN)
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    build(args)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Freeze 10,000 REAL gas crossings from original 2M Geant4 parents/channel.

For each of 40 verified original 50k-parent native export chunks, pick
250 independent generated parents using deterministic SHA ranking and one
SHA-ranked eligible crossing per parent. No new Geant4, Heed or Garfield.
Sample units are PARENT events. Not a full detector-efficiency denominator.
Requires original SHA-provenanced geant_muon_crossings.csv exports on WSL.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

CHANNEL_SOURCE = {"Jpsi": "jpsi40-bore", "Phi": "phi40-bore", "Omega": "omega40-bore"}
TAG = "NA60-PENNING057-ORIGINAL2M-PERCHUNK250-20261108-V1"
N_CHUNKS = 40
PER_CHUNK = 250
FIELDS = (
    "source","chunk","event","hit_index","track_id","station","det_id","pdg",
    "u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
    "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",
    "pu_in_GeV","pv_in_GeV","pw_in_GeV","momentum_GeV","heed_seed"
)
NEEDED = set(FIELDS) - {"source","chunk","momentum_GeV","heed_seed"}


def must(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def seed(*words):
    h = hashlib.sha256("\x1f".join(map(str,words)).encode()).digest()
    return 1 + int.from_bytes(h[:4], "little") % 2_000_000_000


def choose(row):
    try:
        if abs(int(row["pdg"]))!=13 or int(row["station"]) not in range(6):
            return None
        v={k:float(row[k]) for k in (
            "u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
            "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",
            "pu_in_GeV","pv_in_GeV","pw_in_GeV")}
        if not all(math.isfinite(x) for x in v.values()):
            return None
        L=math.sqrt(sum((v[f"{ax}_out_mm"]-v[f"{ax}_in_mm"])**2 for ax in "uvw"))
        if not .05<L<100. or abs(L-v["path_length_mm"])>max(.02,L*.005):
            return None
        if max(abs(v["v_in_mm"]),abs(v["v_out_mm"]))>3.05:
            return None
        half=153.6 if int(row["station"])==0 else 256.
        if max(abs(v["u_in_mm"]),abs(v["u_out_mm"])) >= half-16.:
            return None
        momentum=math.sqrt(sum(v[f"p{ax}_in_GeV"]**2 for ax in "uvw"))
        if not .11<momentum<100000.:
            return None
        dot=sum(((v[f"{ax}_out_mm"]-v[f"{ax}_in_mm"])/L)*
                (v[f"p{ax}_in_GeV"]/momentum) for ax in "uvw")
        if dot<.999 or max(abs(v["theta_u_in_deg"]),abs(v["theta_w_in_deg"]))>180.:
            return None
        # Nine-wire periodic Garfield patch: existing Stage-C source acceptance.
        if abs(v["u_out_mm"]-v["u_in_mm"])>8. or \
           abs(v["w_out_mm"]-v["w_in_mm"])>40.:
            return None
        return momentum
    except (ValueError,TypeError,OverflowError,KeyError):
        return None


def source_chunk(export_root, channel, i):
    folder=export_root/channel.lower()/f"chunk{i:04d}"
    path=folder/"geant_muon_crossings.csv"
    m=folder/"crossings_manifest.json"
    # A completed original Geant4 ROOT chunk may not have a native CSV
    # export yet. Use verified exports only, and report limited coverage.
    if not path.exists() and not m.exists():
        return None
    must(path.is_file() and m.is_file(),
         f"Partially exported original Geant4 chunk: {folder}. "
         "Do not use an unmanifested CSV.")
    doc=json.loads(m.read_text())
    info=doc["immutable_input"]
    must(info["channel"]==channel and info["chunk"]==i and
         info["source"]==CHANNEL_SOURCE[channel] and
         info["n_generated_parent_events"]==50000 and
         info["global_parent_index_interval"]==[(i-1)*50000,i*50000],
         "Original Geant4 source provenance mismatch: "+str(m))
    must(sha(path)==doc["csv_sha256"],
         "Original Geant4 native crossing CSV SHA mismatch: "+str(path))
    return path,doc


def select(export_root, channel, out):
    must(not out.exists(), "Refusing to overwrite existing selected input: "+str(out))
    allselected=[]
    scanned=[]
    available=[]
    missing=[]
    for i in range(1,N_CHUNKS+1):
        found=source_chunk(export_root,channel,i)
        if found is None:
            missing.append(i)
        else:
            path,doc=found
            available.append((i,path,doc))
    must(available, f"No original verified Geant4 native exports for {channel}")
    n_available=len(available)
    print(f"{channel}: {n_available}/40 verified original Geant4 CSV export chunks "
          f"available; absent exports: {missing}",flush=True)
    # Choose 10k unique muon-parent events with quotas balanced across the
    # AVAILABLE verified exports. Every selected parent is genuine original
    # Geant4 data, even if not all 40 source chunks have been CSV-exported.
    quota,remainder=divmod(10000,n_available)
    for ordinal,(i,path,doc) in enumerate(available):
        required=quota+(ordinal<remainder)
        unique={}
        eligible=0
        with path.open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f)
            must(reader.fieldnames is not None and NEEDED.issubset(reader.fieldnames),
                 "Original Geant4 native CSV lacks physical fields: "+str(path))
            for row in reader:
                p=choose(row)
                if p is None:
                    continue
                ev=int(row["event"])
                must(0<=ev<50000, "Original Geant4 local parent-event index outside 0..49999")
                eligible+=1
                ident=(CHANNEL_SOURCE[channel],i,ev,
                       int(row["hit_index"]),int(row["track_id"]),
                       int(row["station"]),int(row["det_id"]))
                parent=(CHANNEL_SOURCE[channel],i,ev)
                rank=hashlib.sha256((TAG+"|crossing|"+"|".join(map(str,ident))).encode()).hexdigest()
                prior=unique.get(parent)
                if prior is None or rank<prior[0]:
                    rec={k:row[k] for k in FIELDS
                         if k not in ("source","chunk","momentum_GeV","heed_seed")}
                    rec.update({"source":CHANNEL_SOURCE[channel],"chunk":i,
                                "momentum_GeV":f"{p:.12g}",
                                "heed_seed":seed(TAG,"heed",*ident)})
                    unique[parent]=(rank,rec)
        must(len(unique)>=required,
             f"{channel} chunk{i}: only {len(unique)} eligible unique "
             f"parent events; need {required}. Export more original Geant4 "
             "chunks with the existing ROOT exporter; no campaign created.")
        candidates=[]
        for parent,(_,rec) in unique.items():
            parenthash=hashlib.sha256((TAG+"|parent|"+"|".join(map(str,parent))).encode()).hexdigest()
            candidates.append((parenthash,parent,rec))
        candidates.sort()
        selected=[rec for _,_,rec in candidates[:required]]
        allselected.extend(selected)
        scanned.append({
            "chunk":i, "csv_sha256":doc["csv_sha256"],
            "root_sha256":doc["immutable_input"]["input_root_sha256"],
            "eligible_gas_crossings":eligible,
            "eligible_unique_generated_parents":len(unique),
            "n_selected":len(selected),
        })
        print(f"PASS {channel} source chunk {i:02d}/40: {len(unique)} eligible parents; "
              f"{required} independent parents selected", flush=True)
    must(len(allselected)==10000 and len({
        (r["source"],int(r["chunk"]),int(r["event"])) for r in allselected})==10000,
        "Not 10,000 globally independent generated parents")
    out.mkdir(parents=True,exist_ok=False)
    csvpath=out/"native_muon_crossings_input.csv"
    with csvpath.open("x",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
        w.writeheader()
        w.writerows(allselected)
    meta={
        "schema":"NA60_PENNING057_ORIGINAL_2M_SOURCE_SELECTION_V1",
        "channel":channel,"source":CHANNEL_SOURCE[channel],
        "generated_parent_events_in_source":2000000,
        "n_selected":10000,"n_source_chunks":n_available,
        "n_source_chunks_possible":N_CHUNKS,
        "missing_unexported_chunks":missing,
        "per_export_chunk_selection":"as balanced as possible among presently verified exports",
        "input_sha256":sha(csvpath),
        "selection_tag":TAG,
        "original_sources":scanned,
        "selection_note":"10,000 hash-ranked distinct parent events sampled across "
                         "the currently available manifest-verified original 50k-parent "
                         "Geant4 CSV exports. Equal per-export quotas. ONE gas crossing "
                         "per parent. Missing unexported chunks never fabricated. "
                         "Not an unbiased all-crossings sample, a full-2M sample, "
                         "or a detector-efficiency denominator.",
        "physics_note":"Native ORIGINAL Geant4 crossing only; later run full TrackHeed "
                       "and Penning r=0.57 Garfield on the ATLAS worker nodes.",
    }
    (out/"selection_manifest.json").write_text(json.dumps(meta,indent=2,sort_keys=True)+"\n")
    print(f"PASS 10000 unique original Geant4 parent crossings for {channel} "
          f"using {n_available} verified source chunks (of 40 total)")
    print("SELECTED CSV:",csvpath)
    print("SHA256:",meta["input_sha256"])
    print("SELECTION MANIFEST:",out/"selection_manifest.json")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--channel",required=True,choices=sorted(CHANNEL_SOURCE))
    ap.add_argument("--export-root",type=Path,
                    default=Path.home()/"na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1")
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    root=a.export_root.expanduser().resolve()
    dest=a.out or (Path.home()/"na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c"
                   /"pi_penning057_from_original2m_v1"/a.channel.lower())
    select(root,a.channel,dest)


if __name__=="__main__":
    main()

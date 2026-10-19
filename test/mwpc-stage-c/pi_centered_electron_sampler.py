#!/usr/bin/env python3
"""NEW centered-v=0 / 3+3 mm Garfield + Heed native gas-hit sampler.

No legacy +1 mm translation. Keeps the old 2+4mm code entirely untouched.
Electrons are sampled along the ACTUAL TGeo-local (u,v,w) true Geant4
gas-hit segment; independent unconditioned Heed primary ionisation. The
saved symmetric Garfield map gives 1-electron collection, wire-offset and
DRIFT time. This is NOT the true anode weighting signal or ASIC firing.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from pi_validate_centered_box_map import validate_map

GEOMETRY_ID = "MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1"
PITCH_MM = 4.0
NBOX = 600
NV = 30
N_PER_STATION = (77,128,128,128,128,128)
REQUIRED_HIT = {
    "event","hit_index","track_id","pdg","det_id","station",
    "path_length_mm","u_in_mm","u_out_mm","v_in_mm","v_out_mm",
    "w_in_mm","w_out_mm","theta_u_in_deg","theta_w_in_deg",
    "edep_eV",
}


def stable_seed(*parts):
    """An absolute per-hit seed, independent of pilot processing order."""
    token="\x1f".join(map(str,parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(token).digest()[:8],"little")


def n_anodes(station):
    if not 0<=int(station)<6:
        raise ValueError("Invalid station")
    return N_PER_STATION[int(station)]


def load_centered_map(path:Path, manifest_path:Path):
    """Fail closed for wrong geometry, wrong map bytes, or missing provenance."""
    path=Path(path)
    manifest_path=Path(manifest_path)
    if not manifest_path.is_file():
        raise FileNotFoundError("Missing centered Garfield provenance manifest")
    meta=json.loads(manifest_path.read_text(encoding="utf-8"))
    expected={
        "geometry_id":GEOMETRY_ID,"gap_minus_mm":3.,
        "gap_plus_mm":3.,"wire_pitch_mm":4.,
        "wire_diameter_um":30.,"anode_voltage_V":1800.,
        "B_local_u_T":0.,"u_bin_min_mm":-2.,"v_bin_min_mm":-3.,
        "du_mm":0.2,"dv_mm":0.2,"n_u":20,"n_v":30,
    }
    for key,value in expected.items():
        if key not in meta or meta[key]!=value:
            raise ValueError(f"Centered response geometry mismatch: {key}")
    if meta.get("provenance_log_checked") is not True:
        raise ValueError("Centred Garfield map not verified against run log")
    if not path.is_file():
        raise FileNotFoundError(path)
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    if sha!=meta.get("response_csv_sha256"):
        raise ValueError("Garfield CSV SHA256 does not match verified manifest")
    n=int(meta["electrons_per_box"])
    if n<=0 or n*600!=int(meta["total_response_rows"]):
        raise ValueError("Incompatible Garfield electron-per-box count")
    # This additionally checks 600 exact u/v bins, all sample identifiers,
    # nonnegative times, valid collection wire ids and total sample exposure.
    summary=validate_map(path,n,None)
    if summary["response_csv_sha256"]!=sha or (
        summary["collected_response_rows"]!=meta.get("collected_response_rows")
    ):
        raise ValueError("Garfield transport summary disagreement")
    pool=np.empty((NBOX,n,3),dtype=np.float64)
    seen=np.zeros((NBOX,n),dtype=bool)
    with path.open(newline="",encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            bid=int(row["box_id"])
            sample=int(row["sample"])
            if seen[bid,sample]:
                raise ValueError(f"Duplicate centered Garfield trial {bid}/{sample}")
            seen[bid,sample]=True
            yes=int(row["collected"])==1
            pool[bid,sample,0]=int(row["collected"])
            pool[bid,sample,1]=int(row["wire_index"]) if yes else 0
            try:
                pool[bid,sample,2]=float(row["drift_time_ns"]) if yes else math.nan
            except ValueError:
                pool[bid,sample,2]=math.nan
    if not seen.all():
        raise ValueError("Centered Garfield map has incomplete sample coverage")
    return pool,meta


def ionize_and_drift(row,rng,lam,marks,pool):
    """Generate one native hit incl n=0: NO rigid alignment offset in v."""
    station=int(row["station"])
    nw=n_anodes(station)
    if pool.ndim!=3 or pool.shape[0]!=600 or pool.shape[2]!=3:
        raise ValueError("Centered transport response has wrong shape")
    path=float(row["path_length_mm"])
    if not(math.isfinite(path) and 0<path<1000):
        raise ValueError(f"Invalid native gas segment path: {path}")
    a=np.array([float(row[f"{c}_in_mm"]) for c in "uvw"])
    b=np.array([float(row[f"{c}_out_mm"]) for c in "uvw"])
    if not(np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError("Nonfinite native TGeo-local coordinate")
    if max(abs(a[1]),abs(b[1]))>3.05:
        raise ValueError("Native gas hit outside centered physical 6-mm gap")
    if (lam<0 or not math.isfinite(lam) or len(marks)==0
        or np.any(np.asarray(marks)<1)):
        raise ValueError("Invalid independent Heed primary ionisation library")
    nc=int(rng.poisson(lam*path))
    mult=marks[rng.integers(len(marks),size=nc)] if nc else np.empty(0,dtype=int)
    frac=np.repeat(rng.random(nc),mult)
    xyz=a[None,:]+frac[:,None]*(b-a)[None,:] if len(frac) else np.empty((0,3))
    nprim=len(frac)
    if nprim==0:
        return dict(nclusters=nc,nprimary=0,ncollected=0,nrawcollected=0,
                    nedge=0,nuncollected=0,
                    u=np.empty(0),v=np.empty(0),w=np.empty(0),
                    u_cell=np.empty(0),anodes=np.empty(0,dtype=int),
                    times=np.empty(0),valid_times=np.empty(0,dtype=bool))
    u,v,w=xyz.T
    if np.any(np.abs(v)>3.05):
        raise ValueError("Heed clusters outside actual 3+3mm gas volume")
    # Fixed finite physical grid: 77 anodes MS0 / 128 all other stations.
    first=np.floor(u/PITCH_MM+(nw-1)/2+.5).astype(np.int64)
    ul=u-(first-(nw-1)/2)*PITCH_MM
    # Nominal wire index uses infinite grid BEFORE finite edge loss check.
    if np.max(np.abs(ul))>2.000001:
        raise AssertionError("Nearest-anode cell normalization failed")
    iu=np.clip(np.floor((ul+2)/.2).astype(np.int64),0,19)
    iv=np.clip(np.floor((v+3)/.2).astype(np.int64),0,29)
    ibox=iu*NV+iv
    sample=np.minimum((rng.random(nprim)*pool.shape[1]).astype(np.int64),
                      pool.shape[1]-1)
    chosen=pool[ibox,sample]
    collected=chosen[:,0]>.5
    target=first+chosen[:,1].astype(np.int64)
    inside=collected&(target>=0)&(target<nw)
    edge=int(np.count_nonzero(collected&~inside))
    good_time=np.isfinite(chosen[inside,2])&(chosen[inside,2]>=0)
    return dict(nclusters=nc,nprimary=nprim,
                ncollected=int(np.count_nonzero(inside)),
                nrawcollected=int(np.count_nonzero(collected)),
                nedge=edge,nuncollected=int(np.count_nonzero(~collected)),
                u=u[inside],v=v[inside],w=w[inside],u_cell=ul[inside],
                anodes=target[inside],times=chosen[inside,2],
                valid_times=good_time)


def multiplicities(hit,cutoffs=(50,100,200,400,800)):
    wires=hit["anodes"]
    times=hit["times"]
    valid=hit["valid_times"]
    counts=[int(np.unique(wires[valid & (times<=t)]).size) for t in cutoffs]
    counts.append(int(np.unique(wires).size))
    if not all(0<=counts[i]<=counts[i+1] for i in range(len(counts)-1)):
        raise AssertionError("Time-gated distinct anode multiplicity nonmonotonic")
    if counts[-1]>hit["ncollected"] or (
        hit["nprimary"] != hit["nrawcollected"] + hit["nuncollected"]
    ) or hit["ncollected"]+hit["nedge"]!=hit["nrawcollected"]:
        raise AssertionError("Primary/collected/edge electron conservation failed")
    return counts

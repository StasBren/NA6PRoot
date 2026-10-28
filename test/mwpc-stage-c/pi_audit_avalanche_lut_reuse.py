#!/usr/bin/env python3
"""READ-ONLY physical/schema audit of frozen LUT20 and LUT500 for ion-avalanche sidecars.

No Garfield execution; no modification to old drift maps. --export writes a
NEW, separate set of stratified *initial positions* to attach full Garfield
AvalancheElectron trials by (geometry, old_LUT_SHA256, box_id, sample).

Important: a fresh avalanche from the same start is NOT a literal continuation
of the old drift trajectory, and not necessarily paired to its old endpoint.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=REPO/"test_runs/mwpc_stage_c"
ROOTS={
    20:BASE/"pi_pitch_fieldwire_pilot_v1/pilot",
    500:None,
}
CONFIGS=(("p3_noFW",3.,0.),("p3_FW75",3.,75.),
         ("p4_noFW",4.,0.),("p4_FW75",4.,75.))
NEEDED=("iu","iv","box_id","sample","random_seed",
        "u_start_mm","v_start_mm","w_start_mm",
        "drift_ok","collected","field_wire_collected",
        "wire_index","u_end_mm","v_end_mm","w_end_mm",
        "drift_time_ns","phi_hit_deg",
        "last_step_u_mm","last_step_v_mm","last_step_energy_eV")
EXPORT_FIELDS=("lut_trials","geometry","original_lut_sha256",
               "box_id","iu","iv","sample","random_seed",
               "u_start_mm","v_start_mm","w_start_mm",
               "old_drift_ok","old_collected","old_field_wire_collected",
               "old_wire_index","old_phi_hit_deg","old_drift_time_ns",
               "old_last_step_energy_eV")

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def tag(pitch,fw):
    return f"p{pitch:.1f}".replace(".","p")+(
        f"_field{fw:.1f}um".replace(".","p") if fw else "_anode_only")

def locations(trials,pitch,fw):
    root=ROOTS[20] if trials==20 else BASE/(
        "pi_pitch_fieldwire_500_p3_v1/pilot" if pitch==3. else
        "pi_pitch_fieldwire_500_p4_v1/pilot")
    directory=root/tag(pitch,fw)
    return directory/"frozen_config.json", directory/"garfield_electrons.csv"

def inspect(trials,name,pitch,fw,export_rows):
    confp,raw=locations(trials,pitch,fw)
    if not confp.is_file() or not raw.is_file():
        raise FileNotFoundError(
            f"Frozen LUT {trials}/{name} missing; expected {confp} and {raw}")
    conf=json.loads(confp.read_text())
    expect={
        "mode":"pilot","pitch_mm":pitch,"field_wire_diam_um":fw,
        "gap_minus_mm":3.,"gap_plus_mm":3.,
        "anode_voltage_v":1800.,"anode_diameter_um":30.,
        "magnetic_field_t":0.,"gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
        "n_u":round(pitch/.1),"n_v":30,
        "du_mm":.1,"dv_mm":.2,"samples_per_box":trials,
    }
    for key,val in expect.items():
        if conf.get(key)!=val:
            raise ValueError(f"{name} LUT{trials} incompatible {key}: {conf.get(key)} vs {val}")
    expected=expect["n_u"]*30*trials
    covered=bytearray(expected)
    found=0;anode=0;fwcol=0
    chosenU={0,expect["n_u"]//2,expect["n_u"]-1}
    chosenV={0,7,14,15,22,29}
    chosenTrial={0,trials//2,trials-1}
    selected=[]
    with raw.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        if not set(NEEDED).issubset(r.fieldnames or ()):
            raise ValueError(f"Old LUT missing fields for avalanche reuse: {name}/{trials}")
        for row in r:
            iu=int(row["iu"]);iv=int(row["iv"]);j=int(row["sample"])
            bid=int(row["box_id"])
            if not (0<=iu<expect["n_u"] and 0<=iv<30 and
                    0<=j<trials and bid==iu*30+iv):
                raise ValueError(f"Invalid cell/trial mapping in {raw}")
            index=bid*trials+j
            if covered[index]:raise ValueError(f"Duplicate LUT key: {name}/{trials}/{bid}/{j}")
            covered[index]=1;found+=1
            col=int(row["collected"]);fwc=int(row["field_wire_collected"])
            if col not in (0,1) or fwc not in (0,1) or col+fwc>1:
                raise ValueError(f"Invalid mutually-exclusive endpoint labels {name}/{trials}")
            anode+=col;fwcol+=fwc
            u=float(row["u_start_mm"]);v=float(row["v_start_mm"])
            if not (math.isfinite(u) and math.isfinite(v) and
                    (-pitch/2+iu*.1-1.e-6)<=u<=( -pitch/2+(iu+1)*.1+1.e-6) and
                    (-3+iv*.2-1.e-6)<=v<=(-3+(iv+1)*.2+1.e-6)):
                raise ValueError(f"Saved exact start is not in expected box: {name}/{trials}/{bid}/{j}")
            if iu in chosenU and iv in chosenV and j in chosenTrial:
                selected.append({
                    "lut_trials":trials,"geometry":name,"box_id":bid,
                    "iu":iu,"iv":iv,"sample":j,
                    "random_seed":int(row["random_seed"]),
                    "u_start_mm":row["u_start_mm"],
                    "v_start_mm":row["v_start_mm"],
                    "w_start_mm":row["w_start_mm"],
                    "old_drift_ok":row["drift_ok"],
                    "old_collected":row["collected"],
                    "old_field_wire_collected":row["field_wire_collected"],
                    "old_wire_index":row["wire_index"],
                    "old_phi_hit_deg":row["phi_hit_deg"],
                    "old_drift_time_ns":row["drift_time_ns"],
                    "old_last_step_energy_eV":row["last_step_energy_eV"],
                })
    if found!=expected or 0 in covered:
        raise ValueError(f"Partial LUT {name}/{trials}: {found}/{expected}")
    if len(selected)!=len(chosenU)*len(chosenV)*len(chosenTrial):
        raise ValueError(f"Partial stratified exact-start selection {name}/{trials}: {len(selected)}")
    digest=sha(raw)
    for item in selected:item["original_lut_sha256"]=digest
    stat={"geometry":name,"trials":trials,"n_boxes":expect["n_u"]*30,
          "n_rows":found,"collected_on_anode":anode,
          "collected_on_field_wire":fwcol,
          "n_stratified_starts":len(selected),
          "lut_csv_sha256":digest,"frozen_config_sha256":sha(confp),
          "original_csv":str(raw),"mode":"READ_ONLY_AUDIT"}
    export_rows.extend(selected)
    return stat

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--trials",choices=("20","500","all"),default="all")
    a.add_argument("--export",action="store_true",
                   help="Write new immutable index of original LUT exact starts; NEVER physics results")
    args=a.parse_args()
    kinds=(20,500) if args.trials=="all" else (int(args.trials),)
    collected=[];stats=[]
    for k in kinds:
        for name,pitch,fw in CONFIGS:
            s=inspect(k,name,pitch,fw,collected)
            stats.append(s)
            print(f"PASS LUT{k:3d} {name:8s}: {s['n_rows']:7d} verified trials, "
                  f"{s['collected_on_anode']:7d} anode, "
                  f"{s['collected_on_field_wire']:7d} field-wire, "
                  f"{s['n_stratified_starts']} exact-start pilot keys",flush=True)
    print(f"PASS read-only LUT reuse audit; {len(collected)} indexed exact starts "
          "(NO Garfield avalanches run)")
    if not args.export:return
    out=BASE/f"pi_avalanche_lut_reuse_v1/frozen_exact_starts_{args.trials}"
    if out.exists():raise FileExistsError(f"Refusing to overwrite existing LUT sidecar index: {out}")
    out.mkdir(parents=True)
    csvout=out/"lut_exact_start_index.csv"
    with csvout.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=EXPORT_FIELDS);w.writeheader()
        w.writerows(collected)
    manifest={"schema":"NA60_PI_LUT_AVALANCHE_EXACT_START_AUDIT_V1",
       "physics_status":"NO_AVALANCHES_COMPUTED; source coordinates only",
       "do_not_infer":"New independent AvalancheElectron from same start IS NOT "
                      "a continuous drift-to-avalanche simulation",
       "sampling":"3 u × 6 v × 3 original trial IDs per geometry; NOT exhaustive or representative",
       "input_luts":stats,
       "index_csv_sha256":sha(csvout),"n_index_rows":len(collected)}
    (out/"index_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("PASS immutable index:",out)

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Select reproducible, strictly held-out Garfield microscopic LUT500 starts.

No Garfield/Geant4 invocation. Read both frozen full-resolution p=3mm CSVs;
validate all 450,000 rows/geometry, grid indices, frozen settings and pinned
source SHA256. Write ONE non-overwriting sample table with independent
train/validation examples from EVERY u/v cell. MC will reuse exactly these
birth coordinates; collision RNG streams are distinct and independent.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
SOURCE=REPO/"test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/pilot"
DEST=REPO/"test_runs/mwpc_stage_c/pi_fast_lut_reference_v1"
GEOMS={
 "p3p0_anode_only":("8b00c43680c8fb63abb0f83e1da6600d48c83cce211bd205bca8c9688740cdbf",0.),
 "p3p0_field75p0um":("0b3bb89883a798fe334195cc8948a198c6c5470d4e45ca1e13af48025360e797",75.),
}
FIELDS=("geometry","split","iu","iv","sample","random_seed","u_start_mm",
        "v_start_mm","reference_collected","reference_field_wire_collected",
        "reference_status","reference_wire_index","reference_phi_hit_deg",
        "reference_drift_time_ns")
NCELL=30*30
NSOURCE=500

def checksum(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda:stream.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()

def sample_plan(geom,nfit,nvalid,seed):
    chosen={}
    for cell in range(NCELL):
        material=f"{seed}:{geom}:{cell}".encode("utf-8")
        rng=random.Random(int.from_bytes(hashlib.sha256(material).digest()[:8],"big"))
        ss=rng.sample(range(NSOURCE),nfit+nvalid)
        for i,slot in enumerate(ss):
            split="fit" if i<nfit else "validation"
            # Distinct MC stream, never the frozen microscopic seed.
            tag=f"fastmc-start-v1:{geom}:{cell}:{slot}:{seed}"
            ms=int.from_bytes(hashlib.sha256(tag.encode()).digest()[:4],"big")&0x7fffffff
            chosen[cell*NSOURCE+slot]=(split,ms)
    return chosen

def validate_config(path,fw):
    manifest=json.loads(path.read_text(encoding="utf-8"))
    expected={
      "study":"PI_GARFIELD_PITCH_FIELDWIRE_V1","pitch_mm":3.,
      "field_wire_diam_um":fw,"samples_per_box":500,
      "mode":"pilot","magnetic_field_t":0.,"n_u":30,"n_v":30,
      "du_mm":.1,"dv_mm":.2,"gap_minus_mm":3.,"gap_plus_mm":3.,
      "anode_diameter_um":30.,"anode_voltage_v":1800.,
      "gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
    }
    for key,value in expected.items():
        if manifest.get(key)!=value:
            raise ValueError(f"Frozen manifest mismatch: {path} {key}: {manifest.get(key)!r}")
    return checksum(path)

def select(out,nfit,nvalid,seed):
    if nfit<1 or nvalid<1 or nfit+nvalid>NSOURCE:
        raise ValueError("Need >=1 fit and validation trial/cell; total <=500")
    if out.exists():raise FileExistsError(f"Refusing to overwrite frozen selection: {out}")
    manifest_path=out.with_suffix(".manifest.json")
    if manifest_path.exists():
        raise FileExistsError(f"Selection manifest already exists: {manifest_path}")
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_name(out.name+".tmp")
    if tmp.exists():raise FileExistsError(f"Stale partial selection, inspect first: {tmp}")
    provenance={}
    selected_all=0
    try:
        with tmp.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.DictWriter(stream,fieldnames=FIELDS)
            writer.writeheader()
            for geom,(pinned_sha,fw) in GEOMS.items():
                directory=SOURCE/geom
                original=directory/"garfield_electrons.csv"
                cfg_hash=validate_config(directory/"frozen_config.json",fw)
                actual=checksum(original)
                if actual!=pinned_sha:
                    raise ValueError(f"FROZEN SHA mismatch {geom}: expected {pinned_sha}, got {actual}")
                chosen=sample_plan(geom,nfit,nvalid,seed)
                seen=bytearray(NCELL*NSOURCE)
                rows=0
                exported=0
                with original.open(newline="",encoding="utf-8") as raw:
                    reader=csv.DictReader(raw)
                    required={"iu","iv","sample","box_id","u_start_mm",
                              "v_start_mm","collected","field_wire_collected",
                              "status","wire_index","phi_hit_deg","drift_time_ns"}
                    if not required.issubset(reader.fieldnames or []):
                        raise ValueError(f"Frozen reference missing columns: {geom}")
                    for row in reader:
                        iu,iv,trial=[int(row[c]) for c in ("iu","iv","sample")]
                        if not (0<=iu<30 and 0<=iv<30 and 0<=trial<NSOURCE):
                            raise ValueError(f"Bad LUT500 index in {geom}: {iu},{iv},{trial}")
                        cell=iu*30+iv
                        if int(row["box_id"])!=cell:
                            raise ValueError(f"Bad LUT cell ID: {geom} cell={cell}")
                        idx=cell*NSOURCE+trial
                        if seen[idx]:raise ValueError("Duplicate micro trial")
                        seen[idx]=1
                        u,v=float(row["u_start_mm"]),float(row["v_start_mm"])
                        if not (math.isfinite(u) and math.isfinite(v)
                           and -1.5+.1*iu-1e-8<=u<=-1.5+.1*(iu+1)+1e-8
                           and -3.+.2*iv-1e-8<=v<=-3.+.2*(iv+1)+1e-8):
                            raise ValueError("Micro birth outside declared cell")
                        c,fwc=int(row["collected"]),int(row["field_wire_collected"])
                        if c not in (0,1) or fwc not in (0,1) or c+fwc>1:
                            raise ValueError("Inconsistent frozen collection")
                        if idx in chosen:
                            split,mc_seed=chosen[idx]
                            writer.writerow({
                             "geometry":geom,"split":split,"iu":iu,"iv":iv,
                             "sample":trial,"random_seed":mc_seed,
                             "u_start_mm":row["u_start_mm"],
                             "v_start_mm":row["v_start_mm"],
                             "reference_collected":c,
                             "reference_field_wire_collected":fwc,
                             "reference_status":row["status"],
                             "reference_wire_index":row["wire_index"],
                             "reference_phi_hit_deg":row["phi_hit_deg"],
                             "reference_drift_time_ns":row["drift_time_ns"]})
                            exported+=1
                        rows+=1
                if rows!=NCELL*NSOURCE or sum(seen)!=rows:
                    raise ValueError(f"Corrupt incomplete micro source: {geom}, count={rows}")
                expected_count=NCELL*(nfit+nvalid)
                if exported!=expected_count:
                    raise ValueError(f"Wrong selected count {geom}: {exported}!={expected_count}")
                provenance[geom]={
                   "sha256":actual,"frozen_config_sha256":cfg_hash,
                   "source_rows":rows,"selected_rows":exported}
                selected_all+=exported
                print(f"VERIFIED frozen {geom}: {rows:,} rows, "
                      f"{exported:,} selected, sha256={actual[:16]}...",flush=True)
        import os
        os.replace(tmp,out)
        info={"schema":"PI_FAST_LUT_EXISTING_MICRO_STARTS_V1",
              "source_provenance":provenance,
              "fit_per_cell":nfit,"validation_per_cell":nvalid,
              "grid":"30x30 cells, 500 microscopic trials/cell",
              "seed":seed,"selection_csv_sha256":checksum(out),
              "selected_electrons":selected_all,
              "split_integrity":"Independent frozen trial indices in every cell; "
               "MC seeds deliberately independent of microscopic random seeds",
              "scope":"p3 B=0 only, no new microscopic trajectories"}
        manifest_path.write_text(json.dumps(info,indent=2)+"\n",encoding="utf-8")
        print(f"PASS selected {selected_all:,} reference starts -> {out}")
        print(f"PASS provenance -> {manifest_path}")
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=DEST/"selected_micro_reference.csv")
    p.add_argument("--fit-per-cell",type=int,default=4)
    p.add_argument("--validation-per-cell",type=int,default=4)
    p.add_argument("--seed",type=int,default=276014)
    a=p.parse_args()
    select(a.output,a.fit_per_cell,a.validation_per_cell,a.seed)

if __name__=="__main__":main()

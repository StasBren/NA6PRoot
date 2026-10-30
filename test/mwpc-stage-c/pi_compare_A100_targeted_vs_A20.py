#!/usr/bin/env python3
"""Re-evaluate the immutable SAME 120 real Geant4 muons with A100 targeted.

Original A20@D500 and full-Geant4/Heed/Garfield reference NEVER changed.
After 16k new A100 original-start microscopic avalanches, 50 cells per
geometry have 100 independent Garfield responses; other cells keep 20.

Uses only frozen Heed primaries; 200 NEW stochastic resamples/muon/field.
Compare coverage, gain, ion-sector asymmetry and unconditional fluctuations.
One physical Garfield response per muon remains a stochastic reference.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,random,statistics
from collections import Counter,defaultdict
from pathlib import Path

from pi_run_muon_lut_validation import (
    load_lut,PRIMARY_FIELDS,LUT_FIELDS,SEED_FIELDS,spatial_cell,
    read_csv,stable_seed,summarize_samples,asym,percent,id_of,
)
from pi_holdout120_A20_from_D500 import (
    REF,GEOMS,ROOT as HOLDOUT_ROOT,COMP as ORIGINAL_COMP,
    payload_rows,verify_full,write_csv,write_json,
)
from pi_A100_targeted_50cell_campaign import (
    ROOT as NEW,PLAN,CHUNK,TOTAL,N_CHUNKS,TOP,ADD,
    verify_plan,validate_chunk,chunk_dir,
)
from pi_prepare_muon_avalanche_b1 import FIELDS,sha
from pi_run_native_muon_b1 import CROSSING_COLUMNS
from pi_publish_run_to_drive import (
    prepare as drive_prepare,publish as drive_publish,EXPECTED_REMOTE,
)

OUTPUT=NEW/"compare_frozen_120_real_muons_A100_vs_original_A20"
N_REP=200
KEY=CROSSING_COLUMNS[0]

def validate_extra():
    m,rows,targets=verify_plan()
    man=NEW/"manifest.json"
    if not man.is_file():
        raise FileNotFoundError("Targeted full A100 Garfield not yet complete")
    final=json.loads(man.read_text())
    if (final.get("status")!="complete" or
        final.get("new_full_Garfield_avalanches")!=TOTAL or
        final.get("original_target_plan_sha256")!=sha(NEW/"plan_manifest.json")):
        raise ValueError("A100 original physics manifest invalid")
    for name,digest in final["files_sha256"].items():
        if sha(NEW/name)!=digest:
            raise ValueError("A100 physics aggregate SHA changed: "+name)
    for i in range(N_CHUNKS):
        chunk=rows[i*CHUNK:(i+1)*CHUNK]
        validate_chunk(chunk_dir(i),i,chunk,m["compiled_executable_sha256"])
        if sha(chunk_dir(i)/"complete.json")!=final["chunk_complete_manifest_sha256"][str(i)]:
            raise ValueError("A100 original Garfield chunk SHA modified")
    return m,final,targets,sha(man)

def load_extended(roots):
    base,hash20=load_lut() # validates entire immutable 84k original A20
    m,final,targets,sha100=roots
    target={(r["geometry"],int(r["box_id"])) for r in targets}
    for geom in GEOMS:
        f=NEW/f"trials_extra_{geom}.csv"
        rows=read_csv(f,LUT_FIELDS)
        if len(rows)!=TOP*ADD:
            raise ValueError("Original extra trials for geometry incomplete")
        seen=set()
        for r in rows:
            b=int(r["box_id"]);s=int(r["sample"])
            if (geom,b) not in target or r["configuration"]!=geom:
                raise ValueError("A100 physics outside selected original target cells")
            if (b,s) in seen:
                raise ValueError("Duplicate extra immutable physics response")
            seen.add((b,s))
            if (int(r["capped"]) or
                int(r["successful"]) not in (0,1)):
                raise ValueError("Unexpected cap/physics status")
            base[geom][b].append((
                int(r["ion_births"]),int(r["ion_up"]),
                int(r["ion_down"]),int(r["ion_top"]),int(r["ion_bottom"])))
        if len(seen)!=TOP*ADD:
            raise ValueError("Not all 80 extra starts per selected cell")
    for geom in GEOMS:
        for b,arr in base[geom].items():
            want=100 if (geom,b) in target else 20
            if len(arr)!=want:
                raise RuntimeError(f"Wrong A100 source empirical pool size {geom} box={b}")
    return base,hash20,sha100,target

def evaluate(reps=N_REP):
    if reps!=N_REP:
        raise ValueError("Fix 200 reps for provenance-comparable original 120 muon benchmark")
    originals=validate_extra()
    rows,stat,ref_sha=verify_full()
    if OUTPUT.exists():
        raise FileExistsError("Existing A100 comparison preserved: "+str(OUTPUT))
    table,oldsha,newsha,target=load_extended(originals)
    orig=read_csv(ORIGINAL_COMP/"muon_comparison.csv",
       ("crossing_id","configuration","station","theta_u_deg","theta_w_deg",
        "path_length_mm","n_Heed_clusters","n_Heed_primaries",
        "full_ion_births","lut_mean_ion_births","lut_sd_ion_births",
        "lut_p05","lut_p50","lut_p95","full_in_lut_90pct",
        "full_zero_birth_Heed_seeds","lut_mean_zero_birth_Heed_seeds",
        "full_A_hemi","lut_mean_A_hemi","lut_A_hemi_p05",
        "lut_A_hemi_p95","full_A_TB","lut_mean_A_TB",
        "lut_A_TB_p05","lut_A_TB_p95",
        "A_hemi_in_lut_90pct","A_TB_in_lut_90pct"))
    orig_map={(r["crossing_id"],r["configuration"]):r for r in orig}
    physical=defaultdict(list)
    for r in read_csv(REF/"muon_reference_heed_primaries.csv",PRIMARY_FIELDS):
        physical[r["crossing_id"]].append(r)
    crossing={(r[KEY],r["configuration"]):r for r in
             read_csv(REF/"muon_reference_crossing.csv",CROSSING_COLUMNS)}
    more=[];replays=[]
    for j,mu in enumerate(rows):
        mid=id_of(mu)
        prim=physical[mid]
        for geom in GEOMS:
            counts=Counter()
            for e in prim:
                b,_,_=spatial_cell(
                    geom,int(mu["station"]),float(mu["u_in_mm"]),
                    float(e["u_rel_cm"]),float(e["v_abs_cm"]))
                counts[b]+=1
            touched=sum(n for b,n in counts.items() if (geom,b) in target)
            rng=random.Random(stable_seed("A100_TARGETED_120MUONS_NOT_OLD_A20_V1",mid,geom))
            ys=[];hemi=[];tb=[]
            for r in range(reps):
                ion=up=down=top=bottom=0
                for b,n in counts.items():
                    arr=table[geom][b]
                    for t in range(n):
                        x=arr[rng.randrange(len(arr))]
                        ion+=x[0];up+=x[1];down+=x[2];top+=x[3];bottom+=x[4]
                ys.append(ion);hemi.append(asym(up,down));tb.append(asym(top,bottom))
            q=summarize_samples(ys)
            ah=summarize_samples([x for x in hemi if x is not None])
            at=summarize_samples([x for x in tb if x is not None])
            old=orig_map[mid,geom]
            full=crossing[mid,geom]
            birth=int(full["ion_births"])
            fhemi=asym(int(full["ion_up"]),int(full["ion_down"]))
            ftb=asym(int(full["ion_top"]),int(full["ion_bottom"]))
            if birth!=int(old["full_ion_births"]):
                raise ValueError("Historic original full Garfield muon reference changed")
            more.append({
                "crossing_id":mid,"configuration":geom,
                "n_real_Heed_electrons":len(prim),
                "n_electrons_in_selected_A100_cells":touched,
                "full_reference_ion_births":birth,
                "A20_original_mean":float(old["lut_mean_ion_births"]),
                "A100_targeted_mean":q["mean"],
                "A20_original_90pct_covered":int(old["full_in_lut_90pct"]),
                "A100_targeted_90pct_covered":int(q["p05"]<=birth<=q["p95"]),
                "A20_old_p05":float(old["lut_p05"]),
                "A20_old_p95":float(old["lut_p95"]),
                "A100_p05":q["p05"],"A100_p95":q["p95"],
                "A20_original_A_hemi_covered":int(old["A_hemi_in_lut_90pct"]),
                "A100_A_hemi_covered":int(ah["p05"]<=fhemi<=ah["p95"]),
                "A20_original_A_TB_covered":int(old["A_TB_in_lut_90pct"]),
                "A100_A_TB_covered":int(at["p05"]<=ftb<=at["p95"]),
                "full_A_hemi":fhemi,"A100_A_hemi_mean":ah["mean"],
                "full_A_TB":ftb,"A100_A_TB_mean":at["mean"],
            })
        if (j+1)%10==0:
            print(f"PASS physical Heed A100 conditional LUT fold {j+1}/120",flush=True)
    if len(more)!=480:raise ValueError("Lost a held-out original muon field")
    bygeom={}
    bootstrap=random.Random(20261030)
    for geom in GEOMS:
        g=[r for r in more if r["configuration"]==geom]
        n=len(g)
        delta=[r["A100_targeted_90pct_covered"]-r["A20_original_90pct_covered"]
               for r in g]
        gains=[r["A100_targeted_mean"] for r in g]
        full=[r["full_reference_ion_births"] for r in g]
        boot=[]
        for t in range(1500):
            inds=[bootstrap.randrange(n) for x in range(n)]
            boot.append(sum(delta[i] for i in inds)/n)
        boot.sort()
        bygeom[geom]={
            "n_independent_Geant4_muons":n,
            "n_real_primary_electrons":sum(r["n_real_Heed_electrons"] for r in g),
            "fraction_primary_electrons_in_A100_selected_cells":
                sum(r["n_electrons_in_selected_A100_cells"] for r in g)/
                sum(r["n_real_Heed_electrons"] for r in g),
            "original_A20_gain_90pct_coverage":
                sum(r["A20_original_90pct_covered"] for r in g)/n,
            "targeted_A100_gain_90pct_coverage":
                sum(r["A100_targeted_90pct_covered"] for r in g)/n,
            "coverage_change_paired_bootstrap95":
                [percent(boot,.025),percent(boot,.975)],
            "original_A20_A_hemi_coverage":
                sum(r["A20_original_A_hemi_covered"] for r in g)/n,
            "targeted_A100_A_hemi_coverage":
                sum(r["A100_A_hemi_covered"] for r in g)/n,
            "original_A20_A_TB_coverage":
                sum(r["A20_original_A_TB_covered"] for r in g)/n,
            "targeted_A100_A_TB_coverage":
                sum(r["A100_A_TB_covered"] for r in g)/n,
            "targeted_A100_over_full_Garfield_mean_gain":
                sum(gains)/sum(full),
            "targeted_A100_ion_hemi_absolute_mean_error":
                statistics.mean(abs(r["A100_A_hemi_mean"]-r["full_A_hemi"]) for r in g),
            "targeted_A100_ion_TB_absolute_mean_error":
                statistics.mean(abs(r["A100_A_TB_mean"]-r["full_A_TB"]) for r in g),
        }
    OUTPUT.mkdir(parents=True)
    out=OUTPUT/"muon_A20_vs_targeted_A100.csv"
    write_csv(out,tuple(more[0]),more)
    summary={
        "schema":"PI_A100_TARGETED_ON_SAME_120_HEED_MUONS_V1",
        "status":"complete","n_true_heldout_Geant4_muons":120,
        "targeted_original_lut_cells":TOP*4,
        "original_A20_manifest_sha256":oldsha,
        "targeted_A100_full_Garfield_manifest_sha256":newsha,
        "original_Geant4_120_full_Garfield_manifest_sha256":ref_sha,
        "n_new_original_Garfield_avalanche_trials":TOTAL,
        "reps_per_muon_and_geometry":reps,
        "by_geometry":bygeom,
        "important_caveats":"This is a position/exposure-driven A100 tuned to the 120 original Heed source positions, NOT a new blinded external 120-muon sample. The original full-Garfield ion outcomes were not used to select target cells. 200 conditional replay replicas share 20 or 100 original Garfield training trials/cell. The present pure empirical LUT does not add an exact metallic-start veto; 11 physically embedded starts in original p4_FW75 remain important sensitivity. No B field, ion drift or strip signal.",
    }
    write_json(OUTPUT/"summary.json",summary)
    summary["outputs_sha256"]={
        p.name:sha(p) for p in OUTPUT.iterdir() if p.is_file()}
    write_json(OUTPUT/"manifest.json",summary)
    for geom,entry in bygeom.items():
        print("OUTCOME",geom,entry,flush=True)
    print("PASS A100 target-cell sensitivity analysis, ORIGINAL 120 full Garfield "
          "physics and 84k baseline A20 untouched, SHA-sealed",flush=True)

def publish_results():
    manifest=OUTPUT/"manifest.json"
    if not manifest.is_file():
        raise FileNotFoundError("No SHA-sealed A20 vs A100 comparison yet")
    data=json.loads(manifest.read_text())
    if data.get("status")!="complete":
        raise RuntimeError("Incomplete A100 comparison")
    for name,d in data.get("outputs_sha256",{}).items():
        if sha(OUTPUT/name)!=d:
            raise ValueError("Original A100 comparison data changed: "+name)
    bundle=drive_prepare(
        OUTPUT,"mwpc_A20_vs_A100_targeted50x4_same120muons_v1")
    drive_publish(bundle,EXPECTED_REMOTE,upload=True)
    print("PASS published unchanged, SHA-checked A20 versus "
          "targeted A100 comparison on original 120 Heed muons",flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("compare","publish"))
    a=p.parse_args()
    if a.mode=="compare":evaluate()
    elif a.mode=="publish":publish_results()

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Immutable pilot of a separate Garfield++ macroscopic drift + Polya proxy.

Mode prepare: SHA-validate the original 16k targeted microscopic campaign and
select 25 independent original exact-start positions per geometry, WITHOUT
looking at their ion-birth outcomes. Mode buildcheck: validate new executable
source/binary and C++ CSV parser without running gas physics. Mode gas: create
new Ar70/CO2 30 Magboltz macroscopic transport table B=0. Mode run: execute
new separate RKF drift + *illustrative* zero-inflated Polya proxy (not an
ion birth census). Mode compare: compare sampled proxy and full microscopic
Garfield physical results on exact same starts, no overwrite of any source.

No Geant4 or TrackHeed; original A20, A100, muon120 immutable.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics,subprocess,random
from pathlib import Path
from collections import defaultdict
from pi_A100_targeted_50cell_campaign import ROOT as A100,INPUT,GEOMS,verify_plan,validate_chunk,chunk_dir,N_CHUNKS,CHUNK
from pi_compare_A100_targeted_vs_A20 import validate_extra
from pi_prepare_muon_avalanche_b1 import REPO,sha

ROOT=A100.parent/'pi_RKF_Polya_separate_B0_pilot_v1'
SOURCE=REPO/'test/mwpc-avalanche/phaseC_rkf_polya_pilot.cxx'
EXEC=REPO/'test/mwpc-avalanche/build/mwpc_phase_c_rkf_polya_pilot'
INPUT_SUBSET=ROOT/'exact_original_D500_pilot_100_starts.csv'
PLAN=ROOT/'pilot_manifest.json'
GAS=ROOT/'Ar70_CO2_30_293K_760T_B0_Magboltz35E.gas'
OUTPUT=ROOT/'rkf_polya_proxy_100_starts.csv'
SUMMARY=ROOT/'summary.json'
REPORT=ROOT/'report_manifest.json'
N_PER_GEOM=25
COLS=('lut_trials','geometry','original_lut_sha256','box_id','iu','iv','sample','random_seed','u_start_mm','v_start_mm','w_start_mm','old_drift_ok','old_collected','old_field_wire_collected','old_wire_index','old_phi_hit_deg','old_drift_time_ns','old_last_step_energy_eV')
RKF_COLUMNS=('geometry','box_id','sample','source_u_mm','source_v_mm','source_w_mm','physical_metal','rkf_transport_ok','rkf_mean_gain','rkf_loss','rkf_final_u_mm','rkf_final_v_mm','rkf_final_w_mm','rkf_end_status','polya_proxy_ion_births','polya_theta','p0_proxy','seconds_drift')

def read_csv(path,fields):
 with path.open(newline='') as f:
  rd=csv.DictReader(f)
  if tuple(rd.fieldnames or ())!=fields:raise ValueError('CSV columns do not match: '+str(path))
  return list(rd)

def json_out(p,obj):p.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def prepare():
 if ROOT.exists():raise FileExistsError('Existing independent RKF pilot preserved: '+str(ROOT))
 original,complete,targets,a100_sha=validate_extra()  # all 40 chunk hashes
 entries=read_csv(INPUT,COLS)
 by_geom=defaultdict(lambda: defaultdict(list))
 for r in entries:by_geom[r['geometry']][int(r['box_id'])].append(r)
 selected=[]
 for geom in GEOMS:
  cells=sorted(by_geom[geom]);rnd=random.Random(int(hashlib.sha256(('PHASEC_RKF_FIXED_20261031_'+geom).encode()).hexdigest()[:16],16))
  rnd.shuffle(cells)
  if len(cells)!=50:raise ValueError('Not exactly 50 A100 cells per geometry')
  for box in cells[:N_PER_GEOM]:
   candidates=sorted(by_geom[geom][box],key=lambda r:int(r['sample']))
   if len(candidates)!=80:raise ValueError('Not 80 truly distinct new samples per cell')
   selected.append(candidates[rnd.randrange(len(candidates))])
 if len(selected)!=100 or len({(r['geometry'],r['box_id'],r['sample']) for r in selected})!=100:raise ValueError('Pilot cardinality/uniqueness wrong')
 ROOT.mkdir(parents=True)
 with INPUT_SUBSET.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=COLS);w.writeheader();w.writerows(selected)
 json_out(PLAN,{'schema':'MWPC_PHASEC_RKF_POLYA_PILOT_PLAN_V1','status':'prepared','original_A100_manifest_sha256':a100_sha,
               'original_D500_selected_starts_sha256':sha(INPUT),'original_source_cpp_sha256':sha(SOURCE),
               'pilot_starts_sha256':sha(INPUT_SUBSET),'n_pilot':100,'n_per_geometry':N_PER_GEOM,
               'selection':'25 randomly chosen distinct cells per geometry; 1 uniformly selected unused original D500 start per chosen cell, deterministic seed; full micro Garfield ion outcomes NOT inspected for selection',
               'original_A20_and_120_full_microscopic_immutable':True})
 print('PASS RKF PILOT PREPARED: 25 original unused D500 full microscopic starts x 4 geometries; NO physics writes',flush=True)

def verify_plan():
 if not PLAN.is_file():raise FileNotFoundError('Run prepare first')
 m=json.loads(PLAN.read_text())
 if m.get('schema')!='MWPC_PHASEC_RKF_POLYA_PILOT_PLAN_V1' or m.get('n_pilot')!=100 or sha(SOURCE)!=m.get('original_source_cpp_sha256') or sha(INPUT_SUBSET)!=m.get('pilot_starts_sha256') or sha(INPUT)!=m.get('original_D500_selected_starts_sha256') or sha(A100/'manifest.json')!=m.get('original_A100_manifest_sha256'):
  raise ValueError('Original pilot/physics source SHA altered; preserve all sources')
 rows=read_csv(INPUT_SUBSET,COLS)
 if len(rows)!=100 or len({(r['geometry'],int(r['box_id']),int(r['sample'])) for r in rows})!=100:raise ValueError('Pilot input corrupt')
 if any(sum(r['geometry']==g for r in rows)!=25 for g in GEOMS):raise ValueError('Unequal geometry balance')
 return rows,m

def check_executable():
 rows,m=verify_plan()
 if not EXEC.is_file():raise FileNotFoundError('Build target mwpc_phase_c_rkf_polya_pilot first: '+str(EXEC))
 proc=subprocess.run([str(EXEC),'--check',str(INPUT_SUBSET)],capture_output=True,text=True)
 if proc.returncode or 'PASS C++ ORIGINAL D500 exact-start parser (100 rows)' not in proc.stdout:
  raise RuntimeError('Actual compiled C++ preflight failed: '+proc.stdout+'\n'+proc.stderr)
 print(proc.stdout.strip())
 return rows,m

def make_gas(ncoll=5):
 check_executable()
 if GAS.exists():raise FileExistsError('Refuse to overwrite gas table: '+str(GAS))
 subprocess.run([str(EXEC),'--make-gas',str(GAS),str(ncoll)],check=True)
 print('PASS RKF new Magboltz gas table SHA256:',sha(GAS),flush=True)

def run(theta,pzero):
 rows,m=check_executable()
 if OUTPUT.exists():raise FileExistsError('Existing RKF pilot output preserved: '+str(OUTPUT))
 if not GAS.is_file():raise FileNotFoundError('Generate new B=0 gas table via gas mode first')
 if not (0.<=theta<=10. and 0.<=pzero<1.):raise ValueError('Invalid proxy parameters')
 subprocess.run([str(EXEC),'--run',str(INPUT_SUBSET),str(GAS),str(OUTPUT),str(theta),str(pzero)],check=True)
 parsed=read_csv(OUTPUT,RKF_COLUMNS)
 if len(parsed)!=len(rows):raise RuntimeError('Incomplete RKF pilot output')
 got={(r['geometry'],int(r['box_id']),int(r['sample'])) for r in parsed}
 want={(r['geometry'],int(r['box_id']),int(r['sample'])) for r in rows}
 if got!=want:raise ValueError('RKF row identities changed')
 json_out(ROOT/'rkf_run_receipt.json',{'schema':'MWPC_PHASEC_RKF_PILOT_RECEIPT_V1','status':'complete','pilot_input_sha256':sha(INPUT_SUBSET),'rkf_file_sha256':sha(OUTPUT),'Magboltz_table_sha256':sha(GAS),'source_cpp_sha256':sha(SOURCE),'n_processed':len(parsed),'theta_proxy':theta,'zero_gate_proxy':pzero,'notes':'RKF GetGain is NOT microscopic ion birth census. Positive Gamma/Polya proxy and Bernoulli zero gate are uncalibrated, explicitly separate from Garfield SetGainFluctuationsPolya which is signal-scaling API.'})
 print('PASS RKF PILOT 100 exact same original D500 starts, no original Garfield rewrite',flush=True)

def compare():
 rows,m=verify_plan()
 receipt=ROOT/'rkf_run_receipt.json'
 if not receipt.is_file():raise FileNotFoundError('No completed RKF pilot run')
 rec=json.loads(receipt.read_text());
 if rec.get('status')!='complete' or rec.get('pilot_input_sha256')!=sha(INPUT_SUBSET) or rec.get('rkf_file_sha256')!=sha(OUTPUT) or rec.get('Magboltz_table_sha256')!=sha(GAS) or rec.get('source_cpp_sha256')!=sha(SOURCE):
  raise ValueError('RKF source/output/gas table SHA mismatch')
 if SUMMARY.exists() or REPORT.exists():raise FileExistsError('RKF pilot comparison already exists; preserve original outputs')
 rkf=read_csv(OUTPUT,RKF_COLUMNS)
 truth={}
 a100=json.loads((A100/'manifest.json').read_text());assert a100['status']=='complete' and a100['new_full_Garfield_avalanches']==16000
 for geom in GEOMS:
  path=A100/f'trials_extra_{geom}.csv'
  if sha(path)!=a100['files_sha256'][path.name]:raise ValueError('A100 trial SHA mismatch')
  with path.open(newline='') as f:
   rd=csv.DictReader(f)
   for r in rd:
    key=(r['configuration'],int(r['box_id']),int(r['sample']))
    if key in truth:raise ValueError('Duplicate A100 physics')
    truth[key]=int(r['ion_births'])
 out={}
 for geom in GEOMS:
  rr=[r for r in rkf if r['geometry']==geom]
  if len(rr)!=25:raise ValueError('Partial geometry rows')
  good=[r for r in rr if r['rkf_transport_ok']=='1']
  metals=[r for r in rr if r['physical_metal']=='1']
  fails=[r for r in rr if r['rkf_transport_ok']!='1' and r['physical_metal']!='1']
  meanrkf=statistics.mean(float(r['rkf_mean_gain']) for r in good) if good else None
  meanmicro=statistics.mean(truth[r['geometry'],int(r['box_id']),int(r['sample'])] for r in rr)
  meanproxy=statistics.mean(float(r['polya_proxy_ion_births']) for r in rr if r['polya_proxy_ion_births']) if not fails else None
  # An event-by-event randomized macro gain cannot be expected to match a microscopic single draw.
  out[geom]={'n_exact_original_D500_starts':len(rr),'n_rkf_transport_ok':len(good),'n_metal_zero_geometry':len(metals),'n_nonmetal_transport_failure':len(fails),
             'mean_microscopic_ion_births_reference':meanmicro,'mean_RKF_GetGain_not_ions':meanrkf,
             'mean_separately_sampled_Polya_proxy':meanproxy,
             'median_seconds_drift':statistics.median(float(r['seconds_drift']) for r in rr),
             'runtime_total_seconds':sum(float(r['seconds_drift']) for r in rr)}
 summ={'schema':'MWPC_PHASEC_RKF_POLYA_EXACT_D500_PILOT_COMPARISON_V1','status':'complete','n_exact_pairs':100,'geometries':out,
       'original_A100_manifest_sha256':m['original_A100_manifest_sha256'],'gas_table_sha256':sha(GAS),'rkf_output_sha256':sha(OUTPUT),
       'critical_interpretation':'Mean RKF gain is not the same observable as Garfield Microscopic ion_births. Zero-gate and theta are provisional; distinct physical ion-birth census, ion sectors and signal are NOT produced. Numerical comparison only exploratory, not validation nor a claim of accuracy. Separate transport failures explicitly counted.'}
 json_out(SUMMARY,summ)
 json_out(REPORT,{'status':'complete','files_sha256':{p.name:sha(p) for p in (PLAN,INPUT_SUBSET,GAS,OUTPUT,SUMMARY)}})
 print(json.dumps(summ,indent=2,ensure_ascii=False),flush=True)
 print('PASS RKF pilot diagnostic, SHA-sealed, original 100k microscopic Garfield immutable',flush=True)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('mode',choices=['prepare','buildcheck','gas','run','compare','status'])
 ap.add_argument('--ncoll',type=int,default=5)
 ap.add_argument('--theta',type=float,default=0.4,help='ILLUSTRATIVE Pólya shape; not calibrated to RKF gain')
 ap.add_argument('--pzero',type=float,default=0.17,help='ILLUSTRATIVE independent zero gate; not calibrated for RKF path')
 a=ap.parse_args()
 if a.mode=='prepare':prepare()
 elif a.mode=='buildcheck':check_executable()
 elif a.mode=='gas':make_gas(a.ncoll)
 elif a.mode=='run':run(a.theta,a.pzero)
 elif a.mode=='compare':compare()
 else:
  print('RKF_ROOT',ROOT,'exists',ROOT.exists())
  print('PILOT_PLAN',PLAN.exists(),'GAS',GAS.exists(),'CSV',OUTPUT.exists(),'COMPARISON',SUMMARY.exists())
if __name__=='__main__':main()
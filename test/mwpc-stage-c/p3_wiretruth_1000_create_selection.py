#!/usr/bin/env python3
"""Freeze 1000 real Geant4 parent-muon crossings from the immutable 30k export.
No physics simulation or reweighting; one randomly ranked parent per channel.
"""
import argparse,csv,gzip,hashlib,io,json,zipfile
from collections import Counter
from pathlib import Path

NATIVE_HEADER=['source','chunk','event','hit_index','track_id','station','det_id','pdg','u_in_mm','v_in_mm','w_in_mm','u_out_mm','v_out_mm','w_out_mm','theta_u_in_deg','theta_w_in_deg','path_length_mm','edep_eV','pu_in_GeV','pv_in_GeV','pw_in_GeV','momentum_GeV','heed_seed']
COUNTS={'Jpsi':334,'Phi':333,'Omega':333}
SALT='NA60-P3-FWTRUTH-20261114-v1;strict-1000-from-30k'

def main():
 p=argparse.ArgumentParser();p.add_argument('--source',default='/mnt/data/NA60_30k_farm_five_raw_exports.zip');p.add_argument('--out',default='/mnt/data/p3_wiretruth_1000/selection_1000.csv');a=p.parse_args()
 source=Path(a.source);out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
 by={k:[] for k in COUNTS}
 with zipfile.ZipFile(source) as z, z.open('crossings_3d.csv.gz') as f:
  for row in csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=f))):
   if row['channel'] in by:by[row['channel']].append(row)
 selected=[]
 for channel,n in COUNTS.items():
  rows=by[channel]
  assert len(rows)==10000, (channel,len(rows))
  assert len({r['parent_index'] for r in rows})==len(rows)
  ranked=sorted(rows,key=lambda r:hashlib.sha256((SALT+'|'+channel+'|'+r['parent_index']+'|'+r['crossing_id']).encode()).digest())[:n]
  for r in ranked:
   seed=1+int.from_bytes(hashlib.sha256((SALT+'|heed|'+r['crossing_id']).encode()).digest()[:8],'big')%1_999_999_998
   selected.append((channel,int(r['parent_index']),r,seed))
 # Deterministic interleaving of 334/333/333 so every 20-event native chunk samples channels.
 ordered=[]
 for index in range(334):
  for channel in ('Jpsi','Phi','Omega'):
   xs=sorted((x for x in selected if x[0]==channel),key=lambda x:x[1])
   if index<len(xs):ordered.append(xs[index])
 assert len(ordered)==1000
 seen=set()
 with out.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=NATIVE_HEADER,lineterminator='\n');writer.writeheader()
  for c,index,r,seed in ordered:
   assert r['crossing_id'] not in seen;seen.add(r['crossing_id'])
   d={k:r[k] for k in NATIVE_HEADER[:-1]};d['heed_seed']=str(seed);writer.writerow(d)
 meta={'protocol':'p3_noFW_penning057_wiretruth_1000_v1','salt':SALT,'source_archive_name':source.name,
       'source_archive_sha256':hashlib.file_digest(source.open('rb'),'sha256').hexdigest(),
       'selected_csv_sha256':hashlib.file_digest(out.open('rb'),'sha256').hexdigest(),
       'count':len(ordered),'by_channel':dict(Counter(x[0] for x in ordered)),
       'source_input':'30k original Geant4 crossing metadata; new TrackHeed samples and Garfield avalanche runs; do not claim exact 30k Heed primary replay',
       'sampling':'SHA256 pseudo-random, equal channel allocation, no angle-based enrichment; 20 muons per chunk',
       'seed_derivation':'SHA256 crossing_id+salt; independent new Heed streams',
       'native_header':NATIVE_HEADER}
 out.with_suffix('.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
 print(json.dumps(meta,indent=2))

if __name__=='__main__':main()

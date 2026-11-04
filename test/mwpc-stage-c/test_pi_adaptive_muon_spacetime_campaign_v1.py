#!/usr/bin/env python3
"""Synthetic isolated source/parser/weighted parent-bootstrap regression. No physics."""
import ast
import contextlib
import csv
import io
import json
import math
import pathlib
import tempfile
from collections import Counter, defaultdict
import numpy as np

SOURCE=pathlib.Path(__file__).with_name('pi_adaptive_muon_spacetime_campaign_v1.py')
tree=ast.parse(SOURCE.read_text())
functions={x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}

def attach(names,ns):
    for name in names:
        exec(compile(ast.Module(body=[functions[name]],type_ignores=[]),
                     str(SOURCE), 'exec'),ns)


def main():
    fields=('source|chunk|event|hit_index|track_id|station|det_id',
            'configuration','G_eff_per_Heed_primary','A_ion_hemi','A_ion_TB',
            'collected_avalanche_electrons','heed_primaries',
            'ion_births','zero_seeds','ion_top','ion_bottom')
    geoms=('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
    parents=[]
    for st in range(6):
        for side in (-1,1):
            for _ in range(10):
                j=len(parents)
                parents.append({'station':st,'theta_u_in_deg':side*3,
                                'id':f'jpsi|2|{j}|0|1|{st}|1'})
    def id_of(row):return row['id']
    ns={'Counter':Counter,'defaultdict':defaultdict,'np':np,'math':math,
        'GEOMS':geoms,'GROUP':10,'id_of':id_of,'csv':csv,'json':json}
    attach(('parent_key','select','source_split','selftest'),ns)
    with contextlib.redirect_stdout(io.StringIO()):ns['selftest']()
    with tempfile.TemporaryDirectory() as tmp:
        base=pathlib.Path(tmp)
        def chunk_dir(i):return base/f'chunk_{i:04d}'
        def read_completed():
            return parents,{'dev_parents':120,'n_parents':120,
              'eligible_source_counts':{f'{st}_{side}':100+st*12+side
               for st in range(6) for side in (-1,1)}},list(range(12))
        def load_json(_):return {'n_full_Garfield_avalanches':1500}
        def write_json(p,o):p.write_text(json.dumps(o))
        for i in range(12):
            d=chunk_dir(i);d.mkdir()
            with (d/'muon_reference_crossing.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
                for parent in parents[i*10:(i+1)*10]:
                    ix=int(parent['id'].split('|')[2])
                    for gnum,g in enumerate(geoms):
                        w.writerow({fields[0]:parent['id'],
                          'configuration':g,
                          'G_eff_per_Heed_primary':1000*(gnum+1)+ix*5,
                          'A_ion_hemi':.15*math.sin(ix*1.7)+.04*gnum,
                          'A_ion_TB':.2*math.sin(ix*1.7)/math.sqrt(gnum+1),
                          'collected_avalanche_electrons':gnum*1000+100,
                          'heed_primaries':50,'ion_births':gnum*1000+120,
                          'zero_seeds':gnum+2,
                          'ion_top':40+gnum*300,'ion_bottom':35+gnum*200})
        ns.update({'REPORTS':base/'reports','read_completed':read_completed,
                   'chunk_dir':chunk_dir,'load_json':load_json,
                   'write_json':write_json})
        attach(('quantile_interval','report'),ns)
        with contextlib.redirect_stdout(io.StringIO()):
            ns['report'](60,.03,.1,.03)
        p=base/'reports/development_0120_parents.json'
        result=json.loads(p.read_text())
        assert result['n_parent_events']==120
        assert result['n_true_Garfield_avalanche_realizations']==18000
        assert result['source_proxy_effective_sample_size']<120
        assert set(result['per_geometry'])==set(geoms)
        assert all(result['per_geometry'][g]['gain']['defined_count']==120
                   for g in geoms)
        assert all('source_proxy_moments' in result['per_geometry'][g]
                   for g in geoms)
        assert all(result['FW_minus_noFW'][key]['source_proxy']['CI95']
                   for key in result['FW_minus_noFW'])
        # No overwrite of an existing report.
        old=p.read_bytes()
        with contextlib.redirect_stdout(io.StringIO()):
            ns['report'](60,.03,.1,.03)
        assert p.read_bytes()==old
    print('PASS independent parent selection + 12-stratum paired bootstrap + source-proxy weights + immutable report; SYNTHETIC ONLY')

if __name__=='__main__':main()

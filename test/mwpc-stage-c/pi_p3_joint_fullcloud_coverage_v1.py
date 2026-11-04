#!/usr/bin/env python3
"""P3 joint birth-cloud bank: physically guarded coverage in original D500 grid.

SHA-verifies full-cloud native physics and original 500-start-per-cell drift
source. Audits deterministic reserved source starts (passes 140...455), never
trains or launches Garfield. Does not change source files or existing bank.
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict

from pi_audit_avalanche_lut_reuse import CONFIGS, locations, sha
from pi_run_lut500_overnight import offset_for, STRIDE, source_maps
import pi_p3_joint_fullcloud_bank_v1 as bankmod

CORE = bankmod.CORE
GEOMS = bankmod.GEOMS
# Ten independent original D500 start sites per spatial box; disjoint from
# old original A20/A100, first pilot and reserved 100..139 extended pilot.
PASSES = (140, 175, 210, 245, 280, 315, 350, 385, 420, 455)


def require(ok, why):
    if not ok:
        raise ValueError(why)


def indices(g, box):
    require(g in GEOMS and 0 <= box < 900, 'Bad original D500 box ID')
    result = [(offset_for(g,box)+STRIDE*j)%500 for j in PASSES]
    require(len(set(result)) == len(PASSES), 'D500 permutation repeats')
    return result


def selftest():
    require(PASSES[0] >= 140 and max(PASSES) < 500, 'Overlapping D500 passes')
    for g in GEOMS:
        for box in (0,23,456,899):
            allseen = {(offset_for(g,box)+STRIDE*p)%500 for p in range(140)}
            require(len(set(indices(g,box))) == len(PASSES) and
                    not allseen.intersection(indices(g,box)),
                    'Validation source overlaps previous full-cloud and A20/A100')
    print('PASS deterministic D500 coverage source selftest')


def coverage():
    bankmod.audit()
    b=bankmod.Bank()
    source,index_sha=source_maps()  # SHA-protected original D500 provenance.
    print('P3 FULL-CLOUD SOURCE COVERAGE MAP: 900 boxes per geometry;')
    print('10 archived independent reserved D500 start sites per box; NOT a muon-flux distribution')
    print('ABSTAIN = no eligible physically guarded bank cloud within 5x5 cells.')
    for geom,pitch,fw in CONFIGS:
        if geom not in GEOMS:
            continue
        _, raw=locations(500,pitch,fw)
        require(sha(raw)==source[geom]['lut_csv_sha256'],
                'Original D500 source changed: '+geom)
        targets={box:set(indices(geom,box)) for box in range(900)}
        checked=set()
        stats=defaultdict(lambda:{'n':0,'gas':0,'metal':0,'supported':0,
                                  'local_supported':0,'eligible_distinct_cells':[],
                                  'eligible_trials':[]})
        represented={r['box_id'] for r in b.items[geom]}
        bybox=defaultdict(lambda: {'n':0,'gas':0,'metal':0,'supported':0})
        with raw.open(newline='') as fh:
            reader=csv.DictReader(fh)
            require({'box_id','sample','u_start_mm','v_start_mm','w_start_mm'}.issubset(
                reader.fieldnames or []), 'Original D500 columns missing')
            for row in reader:
                box=int(row['box_id'])
                sample=int(row['sample'])
                if not (0<=box<900) or sample not in targets[box]:
                    continue
                k=(box,sample)
                require(k not in checked, 'Duplicate original D500 site')
                checked.add(k)
                u=float(row['u_start_mm']);v=float(row['v_start_mm'])
                w=float(row['w_start_mm'])
                require(all(map(math.isfinite,(u,v,w))) and abs(w)<1.e-8 and
                        -1.5 <= u < 1.5 and -3 <= v < 3,
                        'Out-of-domain original D500 start')
                label=CORE.stratum(u,v,geom)
                reg=label[0] if label is not None else 'metal'
                status=stats[reg];status['n']+=1
                bybox[box]['n']+=1
                if label is None:
                    status['metal']+=1
                    bybox[box]['metal']+=1
                    continue
                status['gas']+=1
                bybox[box]['gas']+=1
                choices=b.candidates(u,v,geom,mode='guarded')
                if choices:
                    status['supported']+=1
                    bybox[box]['supported']+=1
                    status['eligible_trials'].append(len(choices))
                    status['eligible_distinct_cells'].append(
                        len({r['box_id'] for r,p in choices}))
                if b.candidates(u,v,geom,mode='local'):
                    status['local_supported']+=1
        require(len(checked)==900*len(PASSES),
                'Original D500 source has missing reserved sites')
        require(all(x['n']==len(PASSES) for x in bybox.values()) and
                len(bybox)==900,'Not all 900 D500 boxes covered')
        gas=sum(s['gas'] for s in stats.values())
        supported=sum(s['supported'] for s in stats.values())
        local=sum(s['local_supported'] for s in stats.values())
        unavailable_cells=[c for c, val in bybox.items()
                           if val['gas']>0 and val['supported']==0]
        print('\nGEOMETRY',geom,
              '\n bank training avalanches',len(b.items[geom]),
              'represented cells',len(represented),
              '\n original D500 test starts',len(checked),'gas',gas,
              'in metal',len(checked)-gas,
              '\n guarded source coverage',f'{supported}/{gas}',
              f'({100*supported/max(gas,1):.2f}%)',
              '\n local-only source coverage',f'{local}/{gas}',
              f'({100*local/max(gas,1):.2f}%)',
              '\n gas cells with NO supported starts',len(unavailable_cells),
              '\n cells WITHOUT any training cloud',900-len(represented))
        for reg in ('near_anode','near_field','transition','bulk','metal'):
            if reg not in stats:continue
            s=stats[reg]
            mean_neigh=(sum(s['eligible_distinct_cells'])/
                        len(s['eligible_distinct_cells'])
                        if s['eligible_distinct_cells'] else None)
            print(' ',reg,'gas',s['gas'],'supported',s['supported'],
                  'local_supported',s['local_supported'],
                  'avg eligible distinct source cells',mean_neigh)
        # Only counts, not new prescriptions based on currently selected data.
        print('NO support box IDs (up to first 35):',sorted(unavailable_cells)[:35])
        print('Reference future uniform campaign: 5 independent new Garfield'
              ' starts per 900 cells = 4500 trials in this geometry;'
              ' initial pilot budget, NOT commissioned')
    print('\nPASS read-only, no new Garfield, no changes to any frozen archive')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('selftest','coverage'))
    a=p.parse_args()
    (selftest if a.action=='selftest' else coverage)()


if __name__=='__main__':
    main()

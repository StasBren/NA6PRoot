#!/usr/bin/env python3
"""Read-only SHA-verified exact Heed wire exposure and ion-birth tail audit.

Replaces misleading earlier median-zero 60/60 near-wire fraction split.
Does not re-run Garfield, rewrite artifacts or optimize the frozen LUT.
"""
import argparse,csv,json,math
from collections import defaultdict
import numpy as np
import pi_external_new120_readonly_diagnostics_v1 as prior

def distances(u,v,g):
    if g not in prior.GEOMS:raise ValueError('Unknown geometry')
    p=3 if g.startswith('p3_') else 4
    da=min(math.hypot(u-k*p,v) for k in range(-4,5))
    df=(min(math.hypot(u-(k+.5)*p,v) for k in range(-4,4))
        if 'FW75' in g else float('inf'))
    return da,df

def source():
    out=defaultdict(lambda:[0,0,0,0,0,0,0])
    seen=set()
    with prior.SEEDS.open(newline='') as f:
        for r in csv.DictReader(f):
            key=(r['crossing_id'],r['configuration'],int(r['primary_index']))
            if key in seen:raise ValueError('Duplicate full Garfield electron')
            seen.add(key)
            a,b=distances(float(r['seed_u_mm']),float(r['seed_v_mm']),key[1])
            z=out[key[:2]]
            n=int(r['ion_births'])
            if n<0:raise ValueError('Negative microscopic ions')
            z[0]+=1
            z[1]+=int(min(a,b)<.35)
            z[2]+=int(a<.015)
            z[3]+=int(b<.0375)
            z[4]+=int(n==0)
            z[5]+=n
            z[6]+=int(int(r['successful'])==0)
    if len(out)!=480:raise ValueError('Expected 480 reference parent/geometry rows')
    return out

def summary(rr):
    if not rr:return {'n':0}
    src=prior.measure(rr)
    return {'n':len(rr),
        'coverage_A20_pct':src['A20']['coverage_pct'],
        'coverage_guarded_A20_pct':src['guarded_A20']['coverage_pct'],
        'MAE_A20':src['A20']['mean_absolute_error'],
        'MAE_guarded_A20':src['guarded_A20']['mean_absolute_error'],
        'PI_width_A20':src['A20']['mean_90PI_width'],
        'PI_width_guarded_A20':src['guarded_A20']['mean_90PI_width']}

def subgroup_ci(rr):
    if len(rr)<2:return None
    a=np.array([100*(int(r['guarded_A20_covered'])-int(r['A20_covered']))
                for r in rr],float)
    b=np.array([float(r['A20_MAE'])-float(r['guarded_A20_MAE'])
                for r in rr],float)
    rng=np.random.default_rng(20261104)
    ii=rng.integers(len(rr),size=(3000,len(rr)))
    return {'paired_coverage_gain_pp':[round(float(x),3) for x in
            (np.mean(a),*np.quantile(a[ii].mean(axis=1),[.025,.975]))],
            'paired_MAE_reduction_ions':[round(float(x),3) for x in
            (np.mean(b),*np.quantile(b[ii].mean(axis=1),[.025,.975]))]}

def analyze(rows,orig):
    bykey={(r['parent'],r['geometry']):r for r in rows}
    if len(bykey)!=480 or set(bykey)!=set(orig):
        raise ValueError('New 120 parent/geometry identity mismatch')
    for k,r in bykey.items():
        n,near,anode,field,zero,ions,fail=orig[k]
        if n!=int(r['n_primaries']) or ions!=int(r['full_Garfield_ions']):
            raise ValueError('Full true ion births/Heed source mismatch: '+str(k))
    result={}
    for g in prior.GEOMS:
        rr=[r for r in rows if r['geometry']==g]
        exposed=[r for r in rr if orig[r['parent'],g][1]>0]
        unexposed=[r for r in rr if orig[r['parent'],g][1]==0]
        tail=sorted(rr,key=lambda r:(-int(r['full_Garfield_ions']),r['parent']))
        alln=sum(orig[r['parent'],g][0] for r in rr)
        result[g]={'total_exact_Heed_starts':alln,
            'within_0_35mm_of_wire':sum(orig[r['parent'],g][1] for r in rr),
            'inside_anode_metal':sum(orig[r['parent'],g][2] for r in rr),
            'inside_field_wire_metal':sum(orig[r['parent'],g][3] for r in rr),
            'microscopic_zero_ion_primary_fraction':round(
                sum(orig[r['parent'],g][4] for r in rr)/alln,5),
            'native_unsuccessful_primaries':sum(orig[r['parent'],g][6] for r in rr),
            'at_least_one_near_wire_primary':summary(exposed),
            'near_wire_bootstrap_exploratory':subgroup_ci(exposed),
            'no_near_wire_primary':summary(unexposed),
            'far_wire_bootstrap_exploratory':subgroup_ci(unexposed),
            'top_12_OBSERVED_Garfield_ion_birth_muons':summary(tail[:12]),
            'other_108_observed_muons':summary(tail[12:])}
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument('mode',choices=['selftest','audit'])
    a=p.parse_args()
    if a.mode=='selftest':
        for g in prior.GEOMS:
            x,y=distances(.00237593484574722,.0107657602328344,g)
            assert x<.015 and y>0
            assert distances(.02,.03,g)[0]>.015
            if 'FW75' in g:
                pitch=3 if g.startswith('p3_') else 4
                assert distances(pitch/2,0,g)[1]==0
        print('PASS strict physical 30um anode and 75um field wire distance selftest; no physics')
        return
    rr,_=prior.data() # verifies complete original native SHA seals
    z=source()
    print(json.dumps({'schema':'NA60_EXTERNAL_NEW120_WIRE_AND_TAIL_POSTHOC_V2',
        'status':'read_only','per_geometry':analyze(rr,z),
        'notes':'TRUE near-wire subgroup means >=1 full Garfield Heed seed within 0.35mm. Top12 tail is selected on OBSERVED Garfield truth, not predictive calibration. Fixed original LUT and full physics unchanged.'},
        indent=2))
    print('PASS SHA-verified real-geometry wire exposure / tail posthoc; no physics or original edits')

if __name__=='__main__':main()

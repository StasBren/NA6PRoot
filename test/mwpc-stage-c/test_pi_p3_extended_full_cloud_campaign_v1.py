#!/usr/bin/env python3
"""Synthetic integration test of extended 3mm full-cloud campaign selection.

No Garfield++ calls, source archive access, modifications of canonical pilot,
or network requests. Uses temporary 2x900x40 D500 records and asserts that
480 and 960 per geometry are balanced by regime/sign and disjoint from the
original A20/A100 pass interval and excluded old pilot IDs.
"""
from __future__ import annotations
import csv
import tempfile
from pathlib import Path
from unittest.mock import patch

import pi_p3_extended_full_cloud_campaign_v1 as campaign
from pi_run_lut500_overnight import offset_for, STRIDE
from pi_audit_avalanche_lut_reuse import sha


def run():
    campaign.selftest()
    with tempfile.TemporaryDirectory(prefix='na60_p3_extended_synthetic_') as d:
        tmp=Path(d)
        physical_index={}
        excluded=set()
        for geom in campaign.GEOMS:
            out=tmp/(geom+'.csv')
            cols=('box_id','sample','iu','iv','u_start_mm','v_start_mm',
                  'w_start_mm','random_seed','drift_ok','collected',
                  'field_wire_collected','wire_index','phi_hit_deg',
                  'drift_time_ns','last_step_energy_eV')
            with out.open('w',newline='') as fh:
                writer=csv.DictWriter(fh,fieldnames=cols)
                writer.writeheader()
                for box in range(900):
                    iu,iv=divmod(box,30)
                    for passage in range(campaign.PASS_START,campaign.PASS_END):
                        sample=(offset_for(geom,box)+STRIDE*passage)%500
                        if box==450 and passage==100:
                            excluded.add((geom,box,sample))
                        writer.writerow({
                            'box_id':box,'sample':sample,'iu':iu,'iv':iv,
                            'u_start_mm':-1.5+(iu+.5)*.1,
                            'v_start_mm':-3.+(iv+.5)*.2,
                            'w_start_mm':0,
                            'random_seed':314159+box*500+sample,
                            'drift_ok':1,'collected':1,
                            'field_wire_collected':0,'wire_index':0,
                            'phi_hit_deg':0,'drift_time_ns':1,
                            'last_step_energy_eV':1})
            physical_index[geom]={'lut_csv_sha256':sha(out)}
        def get_source(_n,pitch,wire_diameter):
            geom='p3_FW75' if wire_diameter else 'p3_noFW'
            return None,tmp/(geom+'.csv')
        campaign.ROOT=tmp/'extended'
        campaign.INPUT=campaign.ROOT/'input'
        campaign.CHUNKS=campaign.ROOT/'full_cloud_chunks'
        campaign.configure_core()
        with patch.object(campaign,'locations',side_effect=get_source),\
             patch.object(campaign,'source_maps',return_value=(physical_index,'SYNTHETIC-INDEX-SHA')),\
             patch.object(campaign,'old_pilot_keys',return_value=(excluded,'SYNTHETIC-OLD48-SHA')):
            for n in (480,960):
                rows,tags,shas,indexsha,counts,oldsha=campaign.selection_from_sources(n)
                assert len(rows)==len(tags)==2*n
                assert not set((r['geometry'],r['box_id'],r['sample'])
                               for r in rows)&excluded
                assert len({(r['geometry'],r['box_id'],r['sample'])
                            for r in rows})==2*n
                for geom in campaign.GEOMS:
                    assert sum(x['chosen'] for x in counts[geom]['strata'].values())==n
                    near=counts[geom]['strata']["('near_anode', 'positive')"]
                    assert near['chosen_distinct_boxes']>1
                print(f'PASS synthetic source {n} distinct starts PER geometry; stratum/box diversity, old-48 exclusions')
            campaign.prepare(480)
            rows,tags,manifest=campaign.check_plan()
            assert len(rows)==len(tags)==960
            assert tuple(manifest['caps'])==campaign.CAPS
            assert manifest['pass_range_reserved']==[100,140]
            assert all(tag['pass_index']>=100 for tag in tags)
            print('PASS synthetic prepare->SHA-check-plan for 480+480 starts, no physics')

if __name__=='__main__':run()

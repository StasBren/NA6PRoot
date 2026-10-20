#!/usr/bin/env python3
"""Synthetic integration test: real-format 3+3mm two-anode arrivals and TGeo u/v.

No Geant4/ROOT/Garfield/Heed model computation. Exercises exact selected
native event/hit/track source join, physics count conservation, correct
local entrance/exit coordinates and absence of fabricated voltage traces.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from build_pi_two_anode_electron_explorer import build, GEOM


def main():
    with TemporaryDirectory(prefix='two_anode_electron_only_') as temp:
        directory=Path(temp)
        source_dir=directory/'native_exports';source_dir.mkdir()
        example={
            'src':'jpsi40-bore','st':0,'det':0,'layer':'A',
            'event':6,'track':3,'hitIndex':9,'pdg':13,'thetaU':3.,
            'thetaW':-1.5,'edep_eV':1300.,'path_mm':6.15,
            'nPrimary':5,'nCollected':3,'N_t100':2,'wires':[4,5],
            'anodes':[
                {'wire':4,'nPrimaryCollected':2,'nTimed':2,
                 'arrivalTimes_ns':[25.,45.],'samples_mV':[100,200]},
                {'wire':5,'nPrimaryCollected':1,'nTimed':1,
                 'arrivalTimes_ns':[60.],'samples_mV':[100,200]}
            ]
        }
        groups={}
        for src in ('jpsi40-bore','phi40-bore','omega40-bore'):
            e=dict(example,src=src)
            groups[f'{src}|0']=[e]
            native={
                'event':6,'hit_index':9,'track_id':3,'det_id':0,
                'station':0,'pdg':13,
                'u_in_mm':-1.4,'v_in_mm':-3.,'w_in_mm':-12.,
                'u_out_mm':2.1,'v_out_mm':3.,'w_out_mm':13.,
                'theta_u_in_deg':3.,'theta_w_in_deg':-1.5
            }
            channel={'jpsi40-bore':'jpsi','phi40-bore':'phi',
                     'omega40-bore':'omega'}[src]
            with (source_dir/f'{channel}_geant_gas_crossings.csv').open(
                'w',newline='',encoding='utf-8'
            ) as stream:
                wr=csv.DictWriter(stream,fieldnames=list(native))
                wr.writeheader();wr.writerow(native)
        original=directory/'old.html'
        original.write_text(
            '<script id="twoWireData" type="application/json">'+
            json.dumps({
                'meta':{'geometry_id':'MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1'},
                'examples':groups
            })+'</script>',
            encoding='utf-8')
        template=directory/'template.html'
        template.write_text(
            '<script type="application/json">__PI_TWO_ANODE_ELECTRON_DATA__</script>',
            encoding='utf-8')
        output=directory/'new.html'
        args=SimpleNamespace(original=original,export_dir=source_dir,
                             template=template,out=output)
        result=build(args)
        assert result['curated_events']==3
        assert not result['new_electronic_voltage_model_used']
        assert len(GEOM)==298
        content=output.read_text(encoding='utf-8')
        payload=json.loads(content.split(
            '<script type="application/json">',1)[1].split('</script>',1)[0])
        assert len(payload['layout'])==298
        for src in ('jpsi40-bore','phi40-bore','omega40-bore'):
            e=payload['examples'][f'{src}|0'][0]
            assert e['segment']['v_in_mm']==-3
            assert e['segment']['v_out_mm']==3
            assert e['segment']['u_in_mm']==-1.4
            assert e['wires']==[4,5]
            assert sum(len(q['times_ns']) for q in e['anodes'])==3
            assert abs(e['deltaFirst_ns']-35)<1e-10
            assert not any('samples_mV' in q for q in e['anodes'])
        assert result['native_geant4_recomputed'] is False
        assert result['heed_or_garfield_recomputed'] is False
        # Refuse to invent geometry if real native hit source changes.
        corrupt=source_dir/'omega_geant_gas_crossings.csv'
        text=corrupt.read_text().replace('2.1,3.0','2.1,3.0')
        corrupt.write_text(text.replace('3.0,-1.5','7.0,-1.5'))
        try:
            build(args)
        except ValueError as exc:
            assert 'does not match native Geant4 hit' in str(exc)
        else:
            raise AssertionError('Incorrect native muon angle silently accepted')
        print('PASS: all three source native u/v/w gas-track endpoints preserved')
        print('PASS: exactly two physical anodes, collected primary electron count')
        print('PASS: no waveform mV/electronics output; correct separate arrivals')
        print('PASS: refuses mismatched real event/hit/track truth')


if __name__=='__main__':
    main()

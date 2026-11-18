"""Synthetic end-to-end smoke: no claim of new Geant4/Garfield physics."""
from __future__ import annotations
import csv
import math
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import run_comparison as r

class TestPipeline(unittest.TestCase):
    def test_selection_and_comparison_read_write_contract(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t);root=t/'repo';(root/'test/mwpc-stage-c').mkdir(parents=True)
            raw=t/'raw';raw.mkdir()
            paths={}
            for ch in r.CHANNELS:
                file=raw/f'{ch}.csv';paths[ch]=file
                # One deterministic muon/station; additionally one overlapping
                # chamber hit of same physical muon, which must be de-duplicated.
                with file.open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=list(r.NATIVE_REQUIRED));w.writeheader()
                    for st in r.STATIONS:
                        for track in (1,2,3):
                            u0,u1=-5.,5.;v0,v1=-2.9,2.9
                            length=math.hypot(u1-u0,v1-v0)
                            rec={k:'0' for k in r.NATIVE_REQUIRED}
                            rec.update(event=str(track),track_id=str(track),
                                station=str(st),pdg='13',det_id=str(st*10+2),
                                hit_index=str(st*10+track),
                                path_length_mm=str(length),
                                u_in_mm=str(u0),u_out_mm=str(u1),
                                v_in_mm=str(v0),v_out_mm=str(v1),
                                w_in_mm='0',w_out_mm='0',
                                pu_in_GeV=str(10*(u1-u0)/length),
                                pv_in_GeV=str(10*(v1-v0)/length),pw_in_GeV='0',
                                theta_u_in_deg='15',theta_w_in_deg='0')
                            w.writerow(rec)
                            if track==1:
                                dup=dict(rec);dup['det_id']=str(st*10+3)
                                dup['hit_index']=str(st*10+track+20)
                                w.writerow(dup)
            def fake_provenance(channel,chunk):
                self.assertEqual(chunk,1)
                return paths[channel],{'csv_sha256':r.sha(paths[channel])}
            def birth(row,lam,marks,seed):
                pos=np.array([[-6.,0.,0.],[-3.,0.,0.],[0.,0.,0.],
                              [3.,0.,0.],[6.,0.,0.]])
                return pos,np.full(len(pos),0.5)
            ht=t/'tracks.csv';hc=t/'clusters.csv';ht.write_text('dummy');hc.write_text('dummy')
            def native(repo):
                return fake_provenance,{x:x for x in r.CHANNELS},birth,ht,hc,lambda a,b:(5.,np.array([1]))
            target=t/'comparison'
            with patch.object(r,'native_interface',native):
                args=SimpleNamespace(repo=root,out=target,chunks=1,n=2,overwrite=False)
                r.select(args)
                with zipfile.ZipFile(t/'lut500.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
                    payload=['iu,iv,collected,wire_index\n']
                    for iu in range(30):
                        for iv in range(30):
                            payload.extend(f'{iu},{iv},1,0\n' for _ in range(500))
                    z.writestr('NA60_LUT500/garfield_electrons.csv',''.join(payload))
                cargs=SimpleNamespace(repo=root,out=target,selected_n=2,limit=2,
                       lut=t/'lut500.zip',calibration=None,weights='unit',overwrite=False)
                r.compare(cargs)
                folder=target/'comparison_2_unit'
                self.assertEqual(len(list(folder.glob('summary_*.json'))),9)
                self.assertTrue((folder/'comparison_table.csv').is_file())
                with (folder/'comparison_table.csv').open(newline='') as f:
                    table=list(csv.DictReader(f))
                self.assertEqual(len(table),9)
                self.assertTrue(all(int(x['n_crossings'])==2 for x in table))
                self.assertEqual(len(list(target.glob('*selected_2.csv.gz'))),9)
                import json
                sel=json.loads((target/'selection_manifest.json').read_text())
                self.assertTrue(all(x['redundant_module_crossings']>=1 for x in sel['samples'].values()))

if __name__=='__main__':unittest.main()
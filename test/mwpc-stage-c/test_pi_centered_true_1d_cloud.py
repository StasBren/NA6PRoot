#!/usr/bin/env python3
"""No-Geant4 no-Garfield microscopic trajectory histogram regression.

Proves that 1D cloud histograms use intermediate Garfield coordinates,
drop collected/stopped electrons after final t, and conserve electron
counts at each requested snapshot (not birth-vs-arrival time plots).
"""
from __future__ import annotations
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

from pi_centered_true_1d_cloud import (
    analyze,load_micro,position_at,TIMES,EDGES,
)


def main():
    with TemporaryDirectory(prefix="pi_true_micro_drift_fixture_") as tmp:
        p=Path(tmp)/"micro.csv"
        columns=("electron_id","point_index","t_rel_birth_ns",
                 "u_cell_mm","v_mm","w_mm","end_status","number_points")
        tracks={
            0:[(0,0.,-1.,-2.,0.),(1,5.,1.,-1.,.15),
               (2,10.,2.,0.,.25),(3,20.,3.,1.,.5)],
            1:[(0,0.,0.,2.,0.),(1,6.,-.5,0.,.1)],
            2:[(0,0.,1.,1.,0.),(1,2.,1.2,.5,.05)],
        }
        with p.open('w',newline='',encoding='utf-8') as stream:
            wr=csv.DictWriter(stream,fieldnames=columns);wr.writeheader()
            for e,points in tracks.items():
                for ip,t,u,v,w in points:
                    wr.writerow(dict(electron_id=e,point_index=ip,
                        t_rel_birth_ns=t,u_cell_mm=u,v_mm=v,w_mm=w,
                        end_status=0,number_points=len(points)))
        raw,status=load_micro(p)
        assert len(raw)==3 and len(status)==3
        assert position_at(raw[0],5.)==(1.,-1.,.15)
        assert position_at(raw[0],7.5)==(1.5,-.5,.2)
        assert position_at(raw[1],10.) is None
        assert position_at(raw[2],5.) is None
        meta=[dict(electron_id=i,source="jpsi40-bore",station=0,
                   angle_group=0,wire_center_u_mm=4.,
                   w_physical_birth_mm=10.) for i in range(3)]
        report=analyze(meta,raw,{})
        b=report["jpsi40-bore|0|0"]
        assert b["born"]==3
        assert b["alive"][0]==3
        assert b["alive"][TIMES.index(5.)]==2
        assert b["alive"][TIMES.index(10.)]==1
        assert b["alive"][TIMES.index(20.)]==1
        assert b["alive"][TIMES.index(40.)]==0
        for hist in b["hist"].values():
            assert [sum(row) for row in hist]==b["alive"]
        assert len(b["hist"])==len(EDGES)==4
        print("PASS: true MICROSCOPIC intermediate Garfield drift points")
        print("PASS: t=5/10ns includes only surviving electron trajectories")
        print("PASS: 1D u(t),v(t),w(t) normalizations conserve survivors")
        print("PASS: no fake extrapolation from birth to endpoint")


if __name__=="__main__":
    main()

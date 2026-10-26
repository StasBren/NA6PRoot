#!/usr/bin/env python3
"""Deterministic in-memory regression for joint phi/time diagnostic."""
from __future__ import annotations

import math
import numpy as np
from pi_fast_lut_transport_joint_audit import (
    PHI_EDGES,angle_hemi,common,grouped_stats,hemi_cluster_bootstrap,
    joint_hist,total_variation
)

def main():
    pairs=[]
    for iu in range(30):
        for iv in range(30):
            for sample in range(4):
                angle=90. if iv>=15 else -90.
                t=2.0+.21*iu+1.3*iv+.04*sample
                ref={"iu":str(iu),"iv":str(iv),
                     "reference_collected":"1",
                     "reference_phi_hit_deg":str(angle),
                     "reference_drift_time_ns":str(t)}
                mc={"anode_collected":"1","phi_hit_deg":str(angle),
                    "drift_time_ns":str(t)}
                pairs.append((ref,mc))
    assert len(pairs)==3600
    x=common(pairs)
    assert len(x)==3600
    time_edges=np.asarray([0.,10.,18.,26.,34.,np.inf])
    assert len(PHI_EDGES)==13
    micro=joint_hist(x,"micro",time_edges)
    fast=joint_hist(x,"mc",time_edges)
    assert micro.shape==(12,5) and micro.sum()==3600
    assert total_variation(micro,fast)==0.
    hemi=angle_hemi(x,"micro")
    assert hemi["P_phi_ge_0"]==.5
    ci=hemi_cluster_bootstrap(x,100,261027)
    assert abs(ci[0])<1e-12 and abs(ci[1])<1e-12
    zones=grouped_stats(x)
    assert len(zones["u_block"])==6
    assert len(zones["v_block"])==6
    assert len(zones["u_edge_vs_inner"])==2
    assert all(abs(z["delta_P_upper"])<1e-12 for z in
              zones["v_block"].values())
    # A changed MC angle in one half of the sample must be detectable.
    for r,m in x[:500]:
        m["phi_hit_deg"]="-90."
    fast2=joint_hist(x,"mc",time_edges)
    assert total_variation(micro,fast2)>0.01
    print("PASS joint phi/time transport: true-zero control, 12x5 bins,"
          " 900-cell clustered hemisphere intervals, spatial zones,"
          " nonzero changed-angle detection.")

if __name__=="__main__":
    main()

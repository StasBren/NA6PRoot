#!/usr/bin/env python3
"""Regression of transport-only statistics on a synthetic 30x30x4 cohort.

No Garfield input or frozen LUT access; validates time bias, wire PDFs,
36-bin phi, 6x6 spatial coverage, and attachment-off fail-closed behavior.
"""
from __future__ import annotations
import math
from pi_fast_lut_transport_only_audit import group_report,ecdf_ks,tv

def synthetic_pairs():
    records=[]
    for iu in range(30):
        for iv in range(30):
            for k in range(4):
                survivor=int(k!=0) # 3/4 microscopic hits; MC all
                t=17. + .2*iv + .1*iu + .01*k
                phi=90. if iv>=15 else -90.
                ref={
                    "iu":iu,"iv":iv,
                    "reference_collected":survivor,
                    "reference_wire_index":0,
                    "reference_phi_hit_deg":phi if survivor else "nan",
                    "reference_drift_time_ns":t if survivor else "nan"}
                mc={
                    "anode_collected":1,"mc_attachment":"off",
                    "wire_index":0,"phi_hit_deg":phi,
                    "drift_time_ns":t+1.5,"drift_wall_us":700.}
                records.append((ref,mc))
    return records

def main():
    assert ecdf_ks([1.,1.,2.],[1.,2.,2.])>0.
    assert tv([-90.,90.],[-90.,90.])==0.
    pairs=synthetic_pairs()
    report=group_report(pairs,123,100)
    assert report["n_birth_sites"]==3600
    assert report["measured_mc_drift_wall_us_per_start"]["mean"]==700.
    assert report["n_micro_anode_survivors"]==2700
    assert report["n_mc_anode_survivors"]==3600
    assert report["n_common_survivor_birth_sites"]==2700
    assert report["wire_index"]["wire_pdf_TV"]==0.
    assert report["wire_index"]["same_wire_on_same_birth_and_both_survived_fraction"]==1.
    assert report["common_birth_survivor_angular_pdf_36bins"]["TV"]==0.
    dt=report["common_birth_survivor_time_ns"]["mean_difference_mc_minus_micro_ns"]
    assert math.isclose(dt,1.5,abs_tol=1.e-12)
    ci=report["common_birth_time_cell_bootstrap"]["cell_bootstrap_95pct_interval_ns"]
    assert all(math.isclose(t,1.5,abs_tol=1e-10) for t in ci)
    blocks=report["coarse_spatial_blocks"]
    assert len(blocks)==36
    assert all(b["n_starts"]==100 and b["n_joint_survivors"]==75
               for b in blocks)
    assert all(math.isclose(b["common_survivors_mean_time_mc_minus_micro_ns"],
                            1.5,abs_tol=1.e-12) for b in blocks)
    pairs[0][1]["mc_attachment"]="on"
    try:
        group_report(pairs,123,100)
    except ValueError as e:
        assert "attachment OFF" in str(e)
    else:
        raise AssertionError("Transport-only report accepted native attachment-on")
    print("PASS transport-only audit regression: common birth, drift bias, "
          "bootstrap, 36-bin angles, 6x6 blocks and off-only enforcement")

if __name__=="__main__":
    main()

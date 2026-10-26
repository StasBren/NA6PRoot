#!/usr/bin/env python3
"""Synthetic CSV regression for the isolated Garfield three-method report."""
import csv,tempfile
from pathlib import Path
from pi_garfield_drift_methods_report import read,stats,GEOMS,METHODS

def main():
    columns=("geometry,method,mc_step_um,iu,iv,trial,random_seed,"
             "u_start_mm,v_start_mm,drift_ok,status,u_end_mm,v_end_mm,"
             "drift_time_ns,anode_collected,field_wire_collected,other_end,"
             "wire_index,phi_hit_deg,drift_wall_us").split(",")
    with tempfile.TemporaryDirectory(prefix="pi_garfield_methods_test_") as temp:
        path=Path(temp)/"fake.csv"
        with path.open("w",newline="",encoding="utf-8") as handle:
            writer=csv.DictWriter(handle,fieldnames=columns)
            writer.writeheader()
            for geometry in GEOMS:
                for method in METHODS:
                    for trial in range(8):
                        anode=trial%3!=0
                        writer.writerow(dict(
                            geometry=geometry,method=method,mc_step_um=5,
                            iu=14,iv=11,trial=trial,random_seed=410000+trial,
                            u_start_mm=0.01*trial,v_start_mm=-1.2,
                            drift_ok=1,status=0,u_end_mm=0.,
                            v_end_mm=0.015,drift_time_ns=17.,
                            anode_collected=int(anode),field_wire_collected=0,
                            other_end=int(not anode),wire_index=0,
                            phi_hit_deg=(90 if anode else "nan"),
                            drift_wall_us=250. if method=="microscopic" else 5.)))
        groups=read(path)
        assert len(groups)==len(GEOMS)*len(METHODS)
        for rows in groups.values():
            x=stats(rows)
            assert x["starts"]==8 and x["anode"]==5
            assert sum(x["phi_probability_per_30deg"])==1.
            assert x["mean_drift_ns"]==17.
    print("PASS: matched birth sites, exclusive collection, finite "
          "timings and normalized coarse angular PDFs")

if __name__=="__main__":
    main()

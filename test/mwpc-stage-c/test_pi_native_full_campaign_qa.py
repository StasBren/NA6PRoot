#!/usr/bin/env python3
"""No-Geant4 synthetic conservation regression for full native PI QA."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv

from pi_native_full_campaign_qa import (
    process_csv, bin_of, EDEP, ANGLE, DELTA, EXPECTED,
)


def main():
    assert len(EXPECTED)==298
    assert EXPECTED[0]==0 and EXPECTED[32]==1
    assert bin_of(1200.,EDEP)>=0
    assert bin_of(-0.001,ANGLE)>=0
    assert bin_of(180.,DELTA)==len(DELTA)-2

    fields=[
        "event","hit_index","track_id","det_id","station","pdg",
        "edep_eV","path_length_mm","theta_u_in_deg",
        "theta_u_out_deg","theta_w_in_deg","theta_w_out_deg",
    ]
    base=dict(event=0,track_id=18,pdg=13,edep_eV=1200.,
              path_length_mm=6.,"theta_u_in_deg":1.,
              "theta_u_out_deg":1.05,
              "theta_w_in_deg":0.5,"theta_w_out_deg":.7)
    with TemporaryDirectory(prefix="native_full_pi_qa_") as temporary:
        p=Path(temporary)/"mock.csv"
        records=[
            dict(base,hit_index=0,station=0,det_id=0),
            dict(base,hit_index=1,station=0,det_id=1),
            dict(base,hit_index=2,station=0,det_id=0), # re-entry SAME det
            dict(base,hit_index=3,station=1,det_id=32),
            dict(base,hit_index=4,event=1,track_id=77,
                 station=5,det_id=278,edep_eV=0),
            dict(base,hit_index=5,event=1,track_id=79,
                 station=5,det_id=278,pdg=11),       # excluded electron
            dict(base,hit_index=6,event=1,track_id=78,
                 station=0,det_id=0,path_length_mm=0.), # flagged zero path
        ]
        with p.open("w",newline="",encoding="utf-8") as stream:
            wr=csv.DictWriter(stream,fieldnames=fields)
            wr.writeheader()
            wr.writerows(records)
        r=process_csv(p)
        assert r["total_detector_hits"]==7
        assert r["muon_hits"]==6 and r["nonmuon_hits"]==1
        assert r["invalid_muon_hits"]==1
        assert r["unique_tracks_with_gas"]==3
        assert r["tracks_crossing_gas_in_all_six"]==0
        assert r["station_unique_gas_tracks"]["0"]==2
        assert r["station_overlap_tracks"]["0"]==1
        assert r["station_extra_chamber_hits"]["0"]==1
        assert r["station_repeated_same_module_segments"]["0"]==1
        assert r["station"]["0"]["n"]==3
        assert r["module"]["0"]["n"]==2
        assert r["module"]["1"]["n"]==1
        assert r["module"]["278"]["positive_edep"]==0
        assert sum(r["module"]["0"]["h"]["edep_eV"])==2
        assert r["module"]["0"]["anomalies"]["nonphysical_edep_or_path"]==1
        # Duplicate physical hit index must never be accepted silently.
        records[1]["hit_index"]=0
        with p.open("w",newline="",encoding="utf-8") as stream:
            wr=csv.DictWriter(stream,fieldnames=fields)
            wr.writeheader()
            wr.writerows(records)
        try:
            process_csv(p)
        except ValueError as e:
            assert "Non-unique" in str(e)
        else:
            raise AssertionError("Duplicate hit index incorrectly allowed")
    print("PASS: native Geant4 gas-hit/source/station/chamber conservation")
    print("PASS: overlap vs same-module repeated segments counted separately")
    print("PASS: nonmuon exclusion, invalid geometry flag, angle/Edep spectra")
    print("PASS: unsafe duplicate hit-index collision rejected")


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Synthetic no-Garfield tests of the centered-field provenance and map schema."""
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

from pi_validate_centered_box_map import validate_map


def main():
    with TemporaryDirectory(prefix="centered_anode_map_test_") as directory:
        root=Path(directory)/"c_transport_box_map_gm3p0_gp3p0"
        root.mkdir()
        table=root/"c_transport_box_electrons.csv"
        log=root/"scan_progress.log"
        log.write_text(
            "gap (-,+) [mm] : 3.0, 3.0\n"
            "B [T] : 0\n"
            "Stage-C transport box map COMPLETE\n"
        )
        columns=("iu","iv","box_id","u_center_mm","v_center_mm",
                 "sample","collected","drift_time_ns","wire_index")
        with table.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=columns)
            writer.writeheader()
            for iu in range(20):
                for iv in range(30):
                    for sample in range(2):
                        writer.writerow({
                            "iu":iu,"iv":iv,"box_id":iu*30+iv,
                            "u_center_mm":-1.9+.2*iu,
                            "v_center_mm":-2.9+.2*iv,
                            "sample":sample,
                            "collected":sample,
                            "drift_time_ns":75. if sample else "nan",
                            "wire_index":0 if sample else 999999,
                        })
        info=validate_map(table,2,log)
        assert info["total_response_rows"]==1200
        assert info["collected_response_rows"]==600
        assert info["gap_minus_mm"]==info["gap_plus_mm"]==3.
        assert info["provenance_log_checked"] is True
        log.write_text("gap (-,+) [mm] : 2.0, 4.0\nB [T] : 0\nCOMPLETE\n")
        try:
            validate_map(table,2,log)
        except ValueError as e:
            assert "non-centered" in str(e)
        else:
            raise AssertionError("Legacy 2+4 mm Garfield field accepted")
        try:
            validate_map(table,3,None)
        except ValueError as e:
            assert "Map incomplete" in str(e)
        else:
            raise AssertionError("Incomplete response library accepted")

    print("PASS: all 600 centered Garfield boxes, exact u/v cell centres")
    print("PASS: provenance 3+3 gap, B0 and fixed response sample numbers")
    print("PASS: reject legacy 2+4 electric fields and incomplete response maps")


if __name__=="__main__":
    main()

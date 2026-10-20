#!/usr/bin/env python3
"""No-Garfield synthetic tests for six-geometry transport pilot validation."""
import csv
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from pi_pitch_fieldwire_pilot import (
    DU, NV, PITCHES, expected_rows, mean_rms, summarise, tag,
    validate_column,
)


def main():
    assert PITCHES == (3., 3.5, 4.)
    assert [round(p/DU) for p in PITCHES] == [30, 35, 40]
    assert expected_rows(3., 2, "smoke") == 60
    assert expected_rows(3.5, 20, "pilot") == 21000
    assert tag(3.5, 0.) == "p3p5_anode_only"
    assert tag(4., 100.) == "p4p0_field100p0um"
    moments = mean_rms([10., 20., 10., 20.], 123, reps=60)
    assert moments["n"] == 4
    assert moments["mean"] == 15. and moments["rms"] == 5.
    assert moments["mean_error"] > 0
    assert mean_rms([], 1)["mean"] is None

    with TemporaryDirectory(prefix="pi_pitch_fieldwire_") as t:
        path = Path(t)/"column.csv"
        pitch, sample_count = 3., 2
        iu = 15
        fields = (
            "iu","iv","box_id","sample","u_center_mm","v_center_mm",
            "collected","field_wire_collected","field_wire_index",
            "drift_time_ns","phi_hit_deg","incidence_normal_deg",
        )
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for iv in range(NV):
                for sample in range(sample_count):
                    c = (sample == 0)
                    field = (sample == 1 and iv % 2 == 0)
                    writer.writerow({
                        "iu":iu,"iv":iv,"box_id":iu*NV+iv,
                        "sample":sample,
                        "u_center_mm":-pitch/2+(iu+.5)*DU,
                        "v_center_mm":-3+(iv+.5)*.2,
                        "collected":int(c),
                        "field_wire_collected":int(field),
                        "field_wire_index":0 if field else 999999,
                        "drift_time_ns":10+iv if c else "nan",
                        "phi_hit_deg":30 if c else "nan",
                        "incidence_normal_deg":30+iv if c else "nan",
                    })
        assert validate_column(path, iu, pitch, 2) == 60
        summary = summarise(path, pitch, 100., "smoke", 2)
        assert summary["n_electrons"] == 60
        assert summary["n_anode_collected"] == 30
        assert summary["n_field_wire_collected"] == 15
        assert summary["n_not_on_wires"] == 15
        assert summary["drift_time_ns"]["n"] == 30
        assert math.isclose(summary["drift_time_ns"]["mean"], 24.5)
        assert summary["drift_time_ns"]["mean_error"] > 0
        assert summary["electron_incidence_normal_deg"]["rms"] > 0
    print("PASS: six pitch tags, exact grids, anode/field-wire separation, "
          "mean/RMS bootstrap errors and conservative smoke-scope labeling")


if __name__ == "__main__":
    main()

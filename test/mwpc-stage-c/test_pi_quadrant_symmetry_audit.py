#!/usr/bin/env python3
"""Synthetic unit tests for the read-only full Geant4 quadrant preflight."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_native_full_campaign_qa import ANGLE
from pi_quadrant_symmetry_audit import (
    FIELDS, geometry_summary, quadrant, total_variation, run,
)


def main():
    geom = geometry_summary()
    assert len(geom) == 6
    assert [g["modules"] for g in geom] == [32, 20, 42, 42, 72, 90]
    assert all(g["four_module_orbits"] > 0 for g in geom)
    assert all(not g["exact_quadrant_chamber_geometry"] for g in geom)
    assert total_variation([1, 3], [1, 3]) == 0
    assert total_variation([1, 0], [0, 1]) == 1
    assert total_variation([0, 0], [0, 0]) is None

    with TemporaryDirectory(prefix="pi_quadrant_preflight_") as tmp:
        root = Path(tmp)
        data = root / "qa"
        data.mkdir()
        out = root / "out"
        (data / "pi_native_all_QA_report.json").write_text(
            json.dumps({"full_campaign_analyzed": True,
                        "processed_chunks_per_source": 67}),
            encoding="utf-8")
        chambers = [c for station in build_chambers() for c in station]
        with (data / "pi_native_all_chamber_quality.csv").open(
                "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=(
                "source", "station", "det_id", "gas_hit_rows"))
            w.writeheader()
            for source in SOURCES:
                for c in chambers:
                    w.writerow(dict(source=source, station=c.station,
                                    det_id=c.det_id, gas_hit_rows=4))
        with (data / "pi_native_all_histograms.csv").open(
                "w", newline="", encoding="utf-8") as f:
            cols = ("source", "station", "det_id", "variable",
                    "bin_lower", "bin_upper", "hit_count",
                    "gas_hit_exposure")
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for source in SOURCES:
                for c in chambers:
                    for field in FIELDS:
                        for i in range(len(ANGLE)-1):
                            w.writerow(dict(
                                source=source, station=c.station,
                                det_id=c.det_id, variable=field,
                                bin_lower=ANGLE[i], bin_upper=ANGLE[i+1],
                                hit_count=4 if i == 0 else 0,
                                gas_hit_exposure=4))
        summary, rows = run(data, out)
        assert not summary["exact_four_quadrant_geometry"]
        assert len(rows) == len(SOURCES) * 6 * 5
        # An empty AXIS category is valid for even-by-even station MS0.
        # The per-chamber rate is undefined (None), not zero or four.
        assert all(row["segments_per_chamber"] == 4
                   for row in rows if row["chambers"] > 0)
        assert all(row["segments_per_chamber"] is None and
                   row["gas_hit_segments"] == 0 and
                   row["ratio_to_Q1_per_chamber"] is None
                   for row in rows if row["chambers"] == 0)
        assert sum(row["chambers"] == 0 for row in rows) == len(SOURCES)
        assert all(row["station"] == 0 and row["quadrant"] == "AXIS"
                   for row in rows if row["chambers"] == 0)
        assert all(row["ratio_to_Q1_per_chamber"] == 1
                   for row in rows if row["quadrant"] in ("Q1","Q2","Q3","Q4"))
        assert sum(row["gas_hit_segments"] for row in rows) == \
               len(SOURCES) * len(chambers) * 4
        assert all(quadrant(c) in ("Q1","Q2","Q3","Q4","AXIS")
                   for c in chambers)
        assert (out / "pi_quadrant_native_counts.csv").is_file()
        print("PASS: geometry parity, all 894 x 2 angle distributions, "
              "quadrant/axis count conservation, no false fourfold symmetry")


if __name__ == "__main__":
    main()

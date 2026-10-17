#!/usr/bin/env python3
"""Offline synthetic regression of per-module Garfield/Heed simulation and resume.

Uses deliberately tiny synthetic files, not the user's ROOT production.
"""
import argparse
import csv
import json
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from run_pi_per_chamber_response import (
    SOURCES, read_garfield, read_heed, run_source, simulate, summaries,
)


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def make_fixture(base):
    heed_tracks = base / "heed_tracks.csv"
    heed_clusters = base / "heed_clusters.csv"
    garfield = base / "map.csv"
    crossings = base / "check_full"
    crossings.mkdir()

    write_csv(
        heed_tracks, ("clusters",), [{"clusters": 20.}, {"clusters": 22.}]
    )
    write_csv(
        heed_clusters, ("electrons",),
        [{"electrons": 1}, {"electrons": 2}, {"electrons": 5}],
    )
    fields = (
        "iu", "iv", "u_center_mm", "v_center_mm",
        "collected", "wire_index", "drift_time_ns",
    )
    # Complete 20x30 response map; each box contains one collected
    # electron and one electron lost; boundary wire shift is deterministic.
    with garfield.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for iu in range(20):
            for iv in range(30):
                for collected in (0, 1):
                    w.writerow({
                        "iu": iu, "iv": iv,
                        "u_center_mm": -1.9 + .2 * iu,
                        "v_center_mm": -1.9 + .2 * iv,
                        "collected": collected,
                        "wire_index": 0 if iu < 10 else 1,
                        "drift_time_ns": 80 if collected else "",
                    })
    source = SOURCES[0]
    records = []
    for index, det in enumerate((52, 53, 52)):
        records.append({
            "source": source, "job": "chunk001",
            "event": index // 2, "track": index % 2,
            "station": 2, "det_id": det, "gas_path_mm": 6.3,
            "gas_center_z_cm": 530,
            "theta_u_deg": 7, "theta_w_deg": -3,
            "y_in_cm": -.2, "y_out_cm": .2,
            "z_in_cm": 529.7, "z_out_cm": 530.3,
        })
    fields = tuple(records[0].keys())
    write_csv(
        crossings / f"crossings_{source}_ms2_530.csv", fields, records
    )
    modules = []
    for det in range(298):
        modules.append({
            "source": source, "scenario_ms2_z_cm": 530,
            "station": 2 if 52 <= det < 94 else 0,
            "det_id": det, "muon_crossings": 2 if det == 52 else
                              1 if det == 53 else 0,
        })
    write_csv(
        crossings / "modules.csv", tuple(modules[0].keys()), modules
    )
    return source, heed_tracks, heed_clusters, garfield, crossings


def main():
    with TemporaryDirectory(prefix="pi_chamber_test_") as temp:
        base = Path(temp)
        source, ht, hc, garfield, crossings = make_fixture(base)
        lam, marks = read_heed(ht, hc)
        boxes = read_garfield(garfield)

        rng = np.random.default_rng(3)
        with (crossings / f"crossings_{source}_ms2_530.csv").open() as f:
            row = next(csv.DictReader(f))
        sim = simulate(row, rng, lam, marks, boxes)
        assert sim["geom_wire_cells"] >= 1
        assert 0 <= sim["collected_electrons"] <= sim["primary_electrons"]
        assert sim["wires_ge1e"] >= sim["wires_ge2e"]
        assert sim["wires_ge2e"] >= sim["wires_ge5e"]
        assert sim["wires_ge5e"] >= sim["wires_ge10e"]
        assert sim["wires_ge1e_100ns"] <= sim["wires_ge1e_200ns"]
        print("PASS: Poisson/Heed/Garfield single segment sanity checks")

        settings = {
            "crossings": crossings,
            "out": base / "resumed",
            "seed": 4123, "max_per_source": 1,
            "checkpoint_every": 1, "started": time.monotonic(),
        }
        args = argparse.Namespace(**settings)
        args.out.mkdir()
        deadline = time.monotonic() + 60
        assert run_source(source, args, deadline, lam, marks, boxes) is False
        assert json.loads(
            (args.out / f"checkpoint_{source}.json").read_text()
        )["processed"] == 1
        args.max_per_source = 0
        assert run_source(source, args, deadline, lam, marks, boxes) is True
        assert (args.out / f"complete_{source}.ok").is_file()
        args.heed_tracks, args.heed_clusters, args.garfield = ht, hc, garfield
        summaries(args)

        args2 = argparse.Namespace(**{
            **settings, "out": base / "uninterrupted", "max_per_source": 0,
        })
        args2.out.mkdir()
        assert run_source(
            source, args2, time.monotonic() + 60, lam, marks, boxes
        ) is True
        file1 = args.out / f"wire_crossings_{source}.csv"
        file2 = args2.out / f"wire_crossings_{source}.csv"
        assert file1.read_bytes() == file2.read_bytes(), (
            "Restart changed RNG sequence or lost CSV rows"
        )
        print("PASS: exact checkpoint/resume reproduces uninterrupted CSV")
        with (args.out / "pi_wire_multiplicity_per_chamber.csv").open() as f:
            rows = list(csv.DictReader(f))
        row = next(r for r in rows if int(r["det_id"]) == 52)
        assert int(row["processed_chamber_crossings"]) == 2
        assert int(row["expected_all_crossings"]) == 2
        assert abs(
            sum(float(row[f"P_nwire_{k}"])
                for k in ("0", "1", "2", "3", "4", "5plus")) - 1.
        ) < 1.e-12
        assert len(rows) == 298
        print("PASS: separate probabilities and crossing counts for 298 modules")


if __name__ == "__main__":
    main()

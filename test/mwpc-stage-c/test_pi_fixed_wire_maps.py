#!/usr/bin/env python3
"""Small synthetic regression for fixed-wire PI atlas, zero external data."""
import csv
import math
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_fixed_wire_occupancy import (
    CHAMBERS, LOOKUP, calc, electronic_response, geometric_cells,
    n_wires, nominal_wire_index, positions,
)
from plot_pi_fixed_wire_atlas import load_wire_csv
from prepare_pi_wire_geometry_controls import produce
from plot_pi_wire_physical_validation_v2 import (
    load_geometric, station_figure, contact_figure, detailed_figure,
)


def save_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_fake_transport(path):
    fields = ("iu", "iv", "collected", "wire_index", "drift_time_ns",
              "u_center_mm", "v_center_mm")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for iu in range(20):
            for iv in range(30):
                for collected in (0, 1):
                    writer.writerow({
                        "iu": iu, "iv": iv, "collected": collected,
                        "wire_index": 0, "drift_time_ns": 95 if collected else "",
                        "u_center_mm": -1.9 + .2*iu,
                        "v_center_mm": -1.9 + .2*iv,
                    })


def main():
    assert len(LOOKUP) == 298
    assert sum(n_wires(ch.station) for ch in LOOKUP.values()) == 36512
    assert positions(0).size == 77
    assert positions(1).size == 128
    assert abs(positions(0)[0] + 152.) < 1.e-10
    assert abs(positions(1)[0] + 254.) < 1.e-10
    assert int(nominal_wire_index(0., 0, 0)) == 38
    assert int(nominal_wire_index(0., 1, 0)) == 64

    # Direct micro-response invariant: collected electrons never assigned
    # to nonexistent wires; phase-specific wire indices are stable.
    offsets = np.arange(600, dtype=np.int64)*2
    sizes = np.full(600, 2, dtype=np.int64)
    response = np.tile([[1, 0, 75], [1, 1, 150]], (600, 1)).astype(float)
    u = np.asarray([0., 0., -255.9, 255.9])
    v = np.zeros(4)
    q = np.asarray([.1, .9, .2, .8])
    counts, early, missing = electronic_response(
        u, v, 1, 0., (offsets, sizes, response), q
    )
    assert int(counts.sum()) + missing == 4
    assert np.all(early <= counts)
    assert int(early.sum()) == 2
    print("PASS: finite wire edges and timing, 36,512-wire geometry")

    with TemporaryDirectory(prefix="pi_fixed_wire_") as folder:
        base = Path(folder)
        inp = base / "crossings"
        inp.mkdir()
        heed_tracks = base / "tracks.csv"
        heed_clusters = base / "clusters.csv"
        garfield = base / "map.csv"
        output = base / "results"
        save_csv(heed_tracks, ("clusters",),
                 [{"clusters": 21}, {"clusters": 23}])
        save_csv(heed_clusters, ("electrons",),
                 [{"electrons": 1}, {"electrons": 2}, {"electrons": 3}])
        make_fake_transport(garfield)
        module_records = []
        for src in SOURCES:
            for det, ch in LOOKUP.items():
                module_records.append({
                    "source": src, "scenario_ms2_z_cm": 530.,
                    "station": ch.station, "det_id": det,
                    "muon_crossings": int(det in (0, 52)),
                })
            records = []
            for det in (0, 52):
                ch = LOOKUP[det]
                z = 530. if ch.station == 2 else 300.
                z += ch.z_offset
                records.append({
                    "source": src, "station": ch.station,
                    "det_id": det, "chamber_center_y_cm": ch.y,
                    "y_in_cm": ch.y-.02, "y_out_cm": ch.y+.02,
                    "z_in_cm": z-.3, "z_out_cm": z+.3,
                    "gas_center_z_cm": z,
                    "gas_path_mm": 6.01,
                    "theta_u_deg": 3.8,
                    "x_in_cm": ch.x-.05, "x_out_cm": ch.x+.05,
                })
            save_csv(
                inp / f"crossings_{src}_ms2_530.csv",
                tuple(records[0]), records,
            )
        save_csv(inp / "modules.csv", tuple(module_records[0]), module_records)
        args = Namespace(
            crossings=inp, heed_tracks=heed_tracks,
            heed_clusters=heed_clusters, garfield=garfield,
            out=output, max_per_source=0, seed=9,
        )
        calc(args)
        vals, expo, counts = load_wire_csv(
            output / "pi_fixed_wire_counts.csv"
        )
        assert vals.shape == (3, 298, 128, 3)
        assert np.all(expo[:, 0] == 1)
        assert np.all(expo[:, 52] == 1)
        assert np.all(expo[:, 1] == 0)
        assert np.isnan(vals[0, 1, 0, 0])
        # A chamber crossing may activate two or more DIFFERENT wires.
        # Consequently the SUM of per-wire activations is NOT bounded by
        # the number of chamber crossings. Each individual wire's count is.
        for det in (0, 52):
            n = n_wires(LOOKUP[det].station)
            active = vals[:, det, :n, :]
            assert np.all(np.isfinite(active))
            assert np.all((active >= 0.) & (active <= 1.))
            # Both stricter conditions are subsets of >=1 collected e-.
            assert np.all(active[:, :, 1] <= active[:, :, 0])
            assert np.all(active[:, :, 2] <= active[:, :, 0])
            # A single crossing can contribute to as many as n wires.
            assert np.all(counts[:, det, :] <= expo[:, det, None] * n)
        # The core physical invariant is already independently enforced
        # for EVERY wire by load_wire_csv: 0 <= n_wire <= chamber crossings.
        with (output / "pi_fixed_wire_station_multiplicity.csv").open() as f:
            rows = list(csv.DictReader(f))
        for source in SOURCES:
            total = [
                int(r["crossings"]) for r in rows
                if r["source"] == source and r["condition"] == "ge1e"
                and float(r["phase_mm"]) == 0
            ]
            assert sum(total) == 2
        with (output / "pi_fixed_wire_module_multiplicity.csv").open() as f:
            per_module = list(csv.DictReader(f))
        assert len(per_module) == 3 * 298 * 3 * 4
        for row in per_module:
            n = int(row["crossings"])
            assert n == int(row["reference_expected_crossings"])
            bins = [int(row["N" + category]) for category in
                    ("0", "1", "2", "3", "4", "5plus")]
            assert sum(bins) == n
            probabilities = [float(row["P" + category]) for category in
                             ("0", "1", "2", "3", "4", "5plus")]
            if n:
                assert abs(sum(probabilities)-1.) < 1e-12
            else:
                assert all(math.isnan(p) for p in probabilities)
        print("PASS: PI per-module P0/P1/P2/P3/P4/P5+ for all 298 chambers")
        print("PASS: full 3-channel/298-detector accounting, wire CSV integrity")
        # Regression for the improved physics presentation: XY trajectory
        # control and independently computed geometric wire profiles.
        control = output / "physical_validation"
        produce(inp, control)
        geometric = load_geometric(control / "pi_geom_wire_counts.csv", expo)
        assert geometric.shape == (3, 298, 128)
        assert geometric[0, 0].sum() >= 1
        assert geometric[0, 52].sum() >= 1
        assert geometric[0, 1].sum() == 0
        xy = np.load(control / "pi_crossing_density_xy.npz")
        assert sum(xy[f"{SOURCES[0]}_MS{st}"].sum() for st in range(6)) == 2
        # One station and chamber rendering in the same plotting functions
        # used by the real atlas; avoids 894 untested PNGs after long run.
        reconstructed = np.zeros((3, 298, 128, 3), dtype=np.int64)
        for src_i in range(3):
            for det, chamber in LOOKUP.items():
                n = expo[src_i, det]
                if n:
                    nw = n_wires(chamber.station)
                    reconstructed[src_i, det, :nw, :] = np.rint(
                        vals[src_i, det, :nw, :] * n
                    ).astype(np.int64)
        testing = control / "atlas_v2"
        for dirname in ("station_validation", "contact_sheets", "chamber_profiles"):
            (testing / dirname).mkdir(parents=True)
        station_figure(testing / "station_validation", 0, 0, expo,
                       geometric, reconstructed, xy, 200., 200.)
        contact_figure(testing / "contact_sheets", 0, 0, expo,
                       geometric, reconstructed)
        detailed_figure(testing / "chamber_profiles", 0, LOOKUP[0],
                        expo, geometric, reconstructed)
        assert len(list(testing.rglob("*.png"))) == 3
        print("PASS: XY density + geometric/collected overlay and real PNG rendering")
        print("PASS: three phase runs, all 3 response conditions, smoke output")


if __name__ == "__main__":
    main()

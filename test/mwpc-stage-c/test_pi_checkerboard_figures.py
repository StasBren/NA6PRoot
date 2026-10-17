#!/usr/bin/env python3
"""Regression checks for the PI gas-crossing hitmaps and angle selection."""
from __future__ import annotations

import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize, PowerNorm

from make_pi_checkerboard_figures import (
    GEOMETRY, allocate, consume_crossings, validate,
    build_color_norm, plot_modules, plot_spatial,
)
from mwpc_checkerboard_projection import NX


def check_color_norm():
    values = np.asarray([0., 1., 4., 9., 1000.])
    norm, clipped = build_color_norm(values, "sqrt", .75)
    assert isinstance(norm, PowerNorm) and clipped
    assert norm(1.) < norm(4.) < norm(9.)
    linear, clipped = build_color_norm(values, "linear", 1.)
    assert isinstance(linear, Normalize) and not clipped
    assert abs(float(linear(1000.)) - 1.) < 1.e-12
    log_norm, clipped = build_color_norm(values, "log", .75)
    assert isinstance(log_norm, LogNorm) and clipped
    assert log_norm(1.) < log_norm(4.) < log_norm(9.)
    empty, clipped = build_color_norm(np.zeros(10), "sqrt", .995)
    assert isinstance(empty, Normalize) and not clipped

    sample = allocate(2.)[2]
    sample["spatial"][10, 10] = 200
    sample["spatial"][10, 11] = 1
    sample["chamber_counts"][2, 2] = 300
    sample["chamber_counts"][2, 3] = 1
    fig, axes = plt.subplots(1, 2)
    for mode in ("linear", "sqrt", "log"):
        spatial, clipped_spatial = plot_spatial(
            axes[0], sample, 2, mode, .5
        )
        modules, clipped_modules = plot_modules(
            axes[1], sample, 2, mode, .5
        )
        assert spatial is not None and clipped_spatial
        assert modules is not None and clipped_modules
        fig.canvas.draw()
        axes[0].cla()
        axes[1].cla()
    plt.close(fig)
    print("PASS: sqrt, log and linear normalization + clipping render correctly")


def main():
    check_color_norm()
    # Two neighboring modules in MS2; same muon crosses both.
    station = 2
    a = GEOMETRY[station][3 * NX[station] + 2]
    b = GEOMETRY[station][3 * NX[station] + 3]
    assert a.layer != b.layer
    # Their sensitive X intervals overlap by 3 cm.
    x_overlap_midpoint = (a.x + b.x) / 2.
    y = a.y
    assert abs(x_overlap_midpoint - a.x) < 66.56 / 2.
    assert abs(x_overlap_midpoint - b.x) < 66.56 / 2.

    fields = (
        "job", "event", "track", "station", "det_id", "row", "col",
        "layer", "x_in_cm", "y_in_cm", "x_out_cm", "y_out_cm",
        "theta_u_deg", "theta_w_deg",
    )
    records = [
        {
            "job": "chunk001", "event": 4, "track": 1, "station": 2,
            "det_id": chamber.det_id, "row": chamber.row,
            "col": chamber.col, "layer": chamber.layer,
            "x_in_cm": x_overlap_midpoint,
            "x_out_cm": x_overlap_midpoint, "y_in_cm": y,
            "y_out_cm": y, "theta_u_deg": 3., "theta_w_deg": 11.,
        }
        for chamber in (a, b)
    ]
    records.append({
        **records[0],
        "event": 5, "track": 2, "theta_u_deg": 5.,
    })
    with TemporaryDirectory() as directory:
        path = Path(directory) / "synthetic.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            w = csv.DictWriter(stream, fieldnames=fields)
            w.writeheader()
            w.writerows(records)
        entries = allocate(2.)
        n = consume_crossings(path, entries)
        validate("jpsi40-bore", entries, Path(directory), 530., False)
        assert n == 2, n
        data = entries[2]
        assert data["n_muons"] == 2
        assert data["n_crossings"] == 3
        assert data["n_multi"] == 1
        assert data["multiplicity"][1] == 1
        assert data["multiplicity"][2] == 1
        assert data["chamber_counts"][a.row, a.col] == 2
        assert data["chamber_counts"][b.row, b.col] == 1
        assert data["spatial"].sum() == 3
        assert data["angles"]["theta_u_deg"] == [3., 5.]
        assert data["angles"]["theta_w_deg"] == [11., 11.]

    print("PASS: overlapping checkerboard modules each receive a spatial hit")
    print("PASS: two muons, three chamber crossings, one double crossing")
    print("PASS: one angle entry per muon/station, not per chamber crossing")


if __name__ == "__main__":
    main()

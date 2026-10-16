#!/usr/bin/env python3
"""Quick deterministic checks for checkerboard geometry / 3-D gas intersection."""
from __future__ import annotations

import math
import random

import mwpc_checkerboard_projection as geom


def test_layout():
    groups = geom.build_chambers()
    assert [len(g) for g in groups] == [32, 20, 42, 42, 72, 90]
    assert [g[0].det_id for g in groups] == [0, 32, 52, 94, 136, 208]
    assert groups[0][0].layer == "A"
    assert groups[0][1].layer == "B"
    assert groups[0][4].layer == "C"
    assert groups[0][5].layer == "D"
    assert all((g[i].det_id + 1 == g[i+1].det_id)
               for g in groups for i in range(len(g)-1))
    return groups


def test_segment_length(groups):
    ch = groups[1][0]
    zs = 340.
    centre_z = zs + ch.z_offset
    ray = geom.intersect_gas(ch.x, ch.y, zs, 0., 0., ch, zs)
    assert ray is not None
    assert abs(ray[0] - (centre_z - 0.3)) < 1.e-10
    assert abs(ray[1] - (centre_z + 0.3)) < 1.e-10
    assert abs(10. * (ray[1] - ray[0]) - 6.) < 1.e-9
    # Forward line inside gas crosses either opposite z faces or a side.
    for sy in (0.0, 0.1, 0.5):
        s = geom.intersect_gas(ch.x, ch.y, zs, 0., sy, ch, zs)
        assert s is not None
        assert s[0] >= centre_z - 0.3 - 1.e-8
        assert s[1] <= centre_z + 0.3 + 1.e-8


def test_candidate_pruning(groups, trials_per_station=1500):
    rng = random.Random(298040)
    checked = 0
    for station, group in enumerate(groups):
        station_z = 530. if station == 2 else geom.Z_STATION[station]
        for _ in range(trials_per_station):
            x0 = rng.uniform(-310., 310.)
            y0 = rng.uniform(-260., 260.)
            sx = rng.uniform(-1.0, 1.0)
            sy = rng.uniform(-0.7, 0.7)
            z0 = 550. if station == 2 else station_z
            fast = {c.det_id for c in geom.candidate_chambers(
                group, station, station_z, x0, y0, z0, sx, sy
            ) if geom.intersect_gas(x0, y0, z0, sx, sy, c, station_z)}
            brute = {c.det_id for c in group
                     if geom.intersect_gas(x0, y0, z0, sx, sy, c, station_z)}
            assert fast == brute, (station, x0, y0, sx, sy, fast, brute)
            checked += 1
    print(f"PASS: candidate pruning equals brute force for {checked} rays")


if __name__ == "__main__":
    groups = test_layout()
    test_segment_length(groups)
    test_candidate_pruning(groups)
    print("PASS: 298 modules, checkerboard IDs/layers, gas thickness")

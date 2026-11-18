"""Synthetic no-physics tests for fine-angle postrun audit primitives."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from inspect_angle_tail import EDGES, REQUIRED, angle_bin, add, zero_bucket


class TestPostRunAngleAudit(unittest.TestCase):
    def test_boundaries_and_finite_angles(self):
        self.assertEqual(angle_bin(0), (0, 5))
        self.assertEqual(angle_bin(4.999), (0, 5))
        self.assertEqual(angle_bin(30), (30, 35))
        self.assertEqual(angle_bin(35), (35, 40))
        self.assertEqual(angle_bin(179.99), (90, 180))
        self.assertEqual(angle_bin(180), (90, 180))
        for value in (-0.001, 180.001, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                angle_bin(value)

    def test_multiplicities(self):
        r = {"n_source_cells": "4", "geometry_n_gt5pct": "3", "lut_n_gt5pct": "2"}
        a = zero_bucket()
        add(a, r)
        self.assertEqual(a["n"], 1)
        self.assertEqual(a["geometry_ge2"], 1)
        self.assertEqual(a["geometry_ge3"], 1)
        self.assertEqual(a["lut_ge2"], 1)
        self.assertEqual(a["lut_ge3"], 0)
        self.assertEqual(a["source_cells_ge3"], 1)
        self.assertEqual(a["geometry_ge4"], 0)


if __name__ == "__main__":
    unittest.main()

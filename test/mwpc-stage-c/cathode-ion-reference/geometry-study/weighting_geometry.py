"""Readout-strip (weighting) geometry for the frozen P3 spectral weighting model, without editing frozen files.

The frozen modules (wire branch, test/mwpc-stage-c/numerical-weighting) hard-code the readout geometry as module
constants that their functions read at CALL time:
  p3_fourier_mfs_mode_pilot.A                    stereo tan in the Fourier-mode slab solution exp(-i k A u), |k|sqrt(1+A^2)
  p3_fourier_strip_reconstruction_pilot.A, .P    stereo tan and strip width (= pitch) in the exact bare-strip term
  p3_spectral_weighting_lut.A, .P                phase eta = w - A u, strip-indicator Fourier coefficients, meta, checks
  p3_replay_full_spectral_muon.A, .P             strip offsets (w - P k) and the landing-strip assignment
`apply()` imports those modules under their own names and sets these constants ONCE per process, before any use, to
the values of a geometry file; every function then runs unchanged. A second, different geometry in the same process
is refused. With the baseline geometry the constants are re-set to their frozen values (no-op), which the tests use to
show that the overlay reproduces the frozen LUT bit for bit.

Parameters that may change (weighting only, transport fixed):  stereo_tan (both planes, U = +A on +v, V = -A on -v,
as in the frozen model), strip_pitch_mm (contiguous strips: width = pitch).
Parameters that may NOT change here (they enter the 2100 V drift field / ion transport or the 9-wire model):
anode pitch 3 mm, anode radius 0.015 mm, anode-cathode gap 3 mm, nine wires, voltage, gas, B = 0.
Strip width != pitch (interstrip gaps) is NOT supported: the frozen model has ideal contiguous strips, and a gap needs
a boundary-condition model (ground, linear, or dielectric) that does not exist in the frozen code.
"""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

FROZEN = {"stereo_tan": 0.1, "strip_pitch_mm": 1.7, "strip_width_mm": 1.7,
          "anode_pitch_mm": 3.0, "anode_radius_mm": 0.015, "gap_mm": 3.0, "n_wires": 9,
          "plane_stereo_signs": "U=+A on +v cathode, V=-A on -v cathode",
          "model": "P3_9wire_infinite_plates_ideal_contiguous_strips"}
TRANSPORT_FIXED = ("anode_pitch_mm", "anode_radius_mm", "gap_mm", "n_wires", "model", "plane_stereo_signs")
MODULES = ("p3_fourier_mfs_mode_pilot", "p3_fourier_strip_reconstruction_pilot", "p3_spectral_weighting_lut",
           "p3_replay_full_spectral_muon")
_APPLIED = None


@dataclass(frozen=True)
class Geometry:
    name: str
    stereo_tan: float
    strip_pitch_mm: float
    strip_width_mm: float
    anode_pitch_mm: float
    anode_radius_mm: float
    gap_mm: float
    n_wires: int
    plane_stereo_signs: str
    model: str

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text())
        g = cls(**{**FROZEN, **d})
        g.check()
        return g

    def physics(self):
        d = asdict(self); d.pop("name"); return d

    @property
    def geometry_id(self):
        return hashlib.sha256(json.dumps(self.physics(), sort_keys=True).encode()).hexdigest()

    @property
    def is_baseline(self):
        return all(getattr(self, k) == FROZEN[k] for k in FROZEN)

    def check(self):
        for k in TRANSPORT_FIXED:
            if getattr(self, k) != FROZEN[k]:
                raise ValueError(f"{k} = {getattr(self, k)} differs from the frozen value {FROZEN[k]}: it changes the "
                                 "drift field / ion transport, which this overlay must not touch")
        if self.strip_width_mm != self.strip_pitch_mm:
            raise ValueError("strip width != pitch (interstrip gaps) is not representable in the frozen ideal "
                             "contiguous-strip model")
        if not (0 < self.stereo_tan < 0.5 and 0.5 <= self.strip_pitch_mm <= 6.0):
            raise ValueError("geometry outside the range this overlay has been reasoned about")


def frozen_sha256(numw_dir):
    out = {}
    for m in MODULES:
        p = Path(numw_dir) / f"{m}.py"
        out[m] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def apply(geometry: Geometry, numw_dir):
    """Import the frozen modules from numw_dir and set their readout-geometry constants. Once per process."""
    global _APPLIED
    if _APPLIED is not None:
        if _APPLIED.geometry_id != geometry.geometry_id:
            raise RuntimeError("a different weighting geometry is already applied in this process")
        return _APPLIED_MODULES
    numw_dir = str(Path(numw_dir).resolve())
    for m in MODULES:
        if m in sys.modules and not str(getattr(sys.modules[m], "__file__", "")).startswith(numw_dir):
            raise RuntimeError(f"{m} already imported from elsewhere")
    sys.path.insert(0, numw_dir)
    mods = {m: importlib.import_module(m) for m in MODULES}
    A, P = geometry.stereo_tan, geometry.strip_pitch_mm
    # values as the frozen modules define them; asserted so a frozen-code change is detected
    for m in MODULES:
        if hasattr(mods[m], "A"):
            assert mods[m].A == FROZEN["stereo_tan"], f"{m}.A is not the frozen value: frozen code changed?"
        if hasattr(mods[m], "P"):
            assert mods[m].P == FROZEN["strip_pitch_mm"], f"{m}.P is not the frozen value: frozen code changed?"
    assert mods["p3_fourier_mfs_mode_pilot"].PITCH == FROZEN["anode_pitch_mm"]
    assert mods["p3_fourier_mfs_mode_pilot"].R == FROZEN["anode_radius_mm"]
    assert mods["p3_fourier_mfs_mode_pilot"].GAP == FROZEN["gap_mm"]
    mods["p3_fourier_mfs_mode_pilot"].A = A
    mods["p3_fourier_strip_reconstruction_pilot"].A = A
    mods["p3_fourier_strip_reconstruction_pilot"].P = P
    mods["p3_spectral_weighting_lut"].A = A
    mods["p3_spectral_weighting_lut"].P = P
    mods["p3_replay_full_spectral_muon"].A = A
    mods["p3_replay_full_spectral_muon"].P = P
    _APPLIED = geometry
    globals()["_APPLIED_MODULES"] = mods
    return mods


def check_table(path, geometry: Geometry):
    """The frozen WeightingLUT checks strip_mm but NOT stereo_tan: a table built for another stereo angle would load
    silently. Every geometry-study tool therefore calls this before using a table."""
    import numpy as np
    with np.load(path, allow_pickle=False) as f:
        meta = json.loads(f["meta"].item())
    for key, val in (("stereo_tan", geometry.stereo_tan), ("strip_mm", geometry.strip_pitch_mm),
                     ("pitch_mm", geometry.anode_pitch_mm), ("radius_mm", geometry.anode_radius_mm),
                     ("gap_mm", geometry.gap_mm), ("model", geometry.model)):
        if meta.get(key) != val:
            raise ValueError(f"LUT {path}: meta {key} = {meta.get(key)!r}, geometry {geometry.name} needs {val!r}")
    return meta


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def describe(geometry: Geometry, numw_dir):
    return {"geometry": asdict(geometry), "geometry_id": geometry.geometry_id, "baseline": geometry.is_baseline,
            "frozen_module_sha256": frozen_sha256(numw_dir)}

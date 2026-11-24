#!/usr/bin/env python3
"""Run one LUT-dependent frozen tool under a readout-strip geometry (weighting_geometry.apply), files untouched.

  geom_run.py GEOMETRY.json NUMW_DIR lut-build OUT.npz [frozen build options]
  geom_run.py GEOMETRY.json NUMW_DIR replay    LUT.npz EXPORT OUT [--half-strips 5 ...]   (frozen replay main)
  geom_run.py GEOMETRY.json NUMW_DIR prompt    LUT.npz EXPORT OUT --replay-dir NUMW_DIR [--half-strips 5]
  geom_run.py GEOMETRY.json NUMW_DIR green     LUT.npz GRID OUT --replay-dir NUMW_DIR --d-max-mm X [--d-sub 10]
  geom_run.py GEOMETRY.json NUMW_DIR describe

Each sub-command calls the unchanged main() of the frozen module (or of the S/T-ladder tool) with the remaining
arguments, after apply() set the strip constants. Before any LUT is used, its meta is checked against the geometry
(check_table: the frozen loader does not check stereo_tan). Writes <OUT>.geometry.json (or OUT/geometry.json for a
directory output) with geometry, geometry_id and the sha256 of the frozen modules. The baseline geometry is a no-op.
"""
import json
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import weighting_geometry as wg  # noqa: E402

TOOLS = {"replay": None, "prompt": HERE.parent / "timing-ablation" / "t_prompt_electron_term.py",
         "green": HERE.parent / "spatial-ablation" / "s_green_table.py"}


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    geom = wg.Geometry.load(sys.argv[1])
    numw = Path(sys.argv[2])
    cmd, rest = sys.argv[3], sys.argv[4:]
    mods = wg.apply(geom, numw)
    info = wg.describe(geom, numw)
    if cmd == "describe":
        print(json.dumps(info, indent=1)); return 0
    if cmd == "lut-build":
        out = Path(rest[0])
        sys.argv = ["p3_spectral_weighting_lut.py", "build"] + rest
        rc = mods["p3_spectral_weighting_lut"].main() or 0
        meta = wg.check_table(out, geom)
        info.update(lut=str(out), lut_sha256=wg.sha256(out), lut_meta=meta, build_args=rest)
        Path(str(out) + ".geometry.json").write_text(json.dumps(info, indent=1) + "\n")
        print(f"weighting_id (LUT sha256) {info['lut_sha256']}  geometry_id {geom.geometry_id}")
        return rc
    if cmd not in TOOLS:
        raise SystemExit(f"unknown command {cmd}")
    table, out = Path(rest[0]), Path(rest[2])
    info.update(lut=str(table), lut_sha256=wg.sha256(table), lut_meta=wg.check_table(table, geom), args=rest)
    if cmd == "replay":
        sys.argv = ["p3_replay_full_spectral_muon.py"] + rest
        rc = mods["p3_replay_full_spectral_muon"].main() or 0
    else:
        sys.argv = [str(TOOLS[cmd])] + rest
        try:
            runpy.run_path(str(TOOLS[cmd]), run_name="__main__")
            rc = 0
        except SystemExit as e:
            rc = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    # the frozen modules must still carry this geometry (nothing re-imported them)
    for m, mod in mods.items():
        if hasattr(mod, "A") and mod.A != geom.stereo_tan or hasattr(mod, "P") and mod.P != geom.strip_pitch_mm:
            raise SystemExit(f"{m}: strip constants changed during the run")
    if rc == 0 and out.is_dir():
        (out / "geometry.json").write_text(json.dumps(info, indent=1) + "\n")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

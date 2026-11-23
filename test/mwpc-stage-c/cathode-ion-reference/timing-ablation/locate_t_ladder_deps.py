#!/usr/bin/env python3
"""Find and verify the EXISTING dependencies of the real 2100 V T-ladder on this machine. Read-only; builds nothing.

  python3 locate_t_ladder_deps.py [ROOT ...]  (default roots: ~/na6p ~/NA6PRoot-signal40-2100 ~/NA6PRoot*)
Writes ./t_ladder_deps_report.txt and ./t_ladder_deps.env (candidate exports; review before sourcing).

Looks for
 1. exporter executable  mwpc_stage_c_p3_export_full_spectral_paths
 2. ion mobility table   IonMobility_Ar+_Ar.txt (the provisional Ar+ in Ar table named in the 2100 V reference docs)
 3. spectral LUT         *.npz with arrays near/bulk/meta and meta.model == P3_9wire_infinite_plates_ideal_contiguous_strips
 4. replay code          p3_replay_full_spectral_muon.py (+ its imports) in a numerical-weighting directory
 5. previous frozen path exports (ions.csv + vertices.f64 + complete.txt) and frozen replay outputs
    (waveforms.csv + channels.csv + summary.json), with rows / exact_ions so they can be paired with a muon
For each LUT the full meta and the SHA256 are printed; a LUT rebuilt with the default `build` settings has the meta
printed at the end (identical meta = identical physics model and grid; bytes may differ across numpy/BLAS builds).
"""
import hashlib
import json
import os
import sys
from pathlib import Path

DEFAULT_META = {"model": "P3_9wire_infinite_plates_ideal_contiguous_strips", "L_mm": 48.0, "cutoff_per_mm": 6.0,
                "sources": 24, "image_extra": 4, "phase_max_mm": 12.0, "radius_mm": 0.015, "pitch_mm": 3.0,
                "strip_mm": 1.7, "gap_mm": 3.0, "stereo_tan": 0.1, "near_max_mm": 0.75, "blend_mm": [0.5, 0.65],
                "u_max_mm": 13.5, "near_shape": [5, 32, 48, 241], "bulk_shape": [91, 41, 241],
                "field_accuracy_detector_validated": False}
SKIP = {".git", "node_modules", "__pycache__", ".cache", "proc", "sys"}


def sha(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def walk(roots):
    seen = set()
    for r in roots:
        r = Path(r).expanduser()
        if not r.exists():
            continue
        for d, dirs, files in os.walk(r, followlinks=False):
            dirs[:] = [x for x in dirs if x not in SKIP and not d.startswith("/mnt/")]
            rd = os.path.realpath(d)
            if rd in seen:
                dirs[:] = []
                continue
            seen.add(rd)
            yield Path(d), set(files)


def seal(d):
    try:
        return dict(l.split("=", 1) for l in (d / "complete.txt").read_text().split() if "=" in l)
    except Exception:
        return {}


def main():
    roots = sys.argv[1:] or [str(p) for p in [Path("~/na6p").expanduser(), Path("~/NA6PRoot-signal40-2100").expanduser()]
                             + sorted(Path.home().glob("NA6PRoot*"))]
    exe, mob, luts, numw, exports, replays = [], [], [], [], [], []
    for d, files in walk(roots):
        if "mwpc_stage_c_p3_export_full_spectral_paths" in files and os.access(d / "mwpc_stage_c_p3_export_full_spectral_paths", os.X_OK):
            exe.append(d / "mwpc_stage_c_p3_export_full_spectral_paths")
        if "IonMobility_Ar+_Ar.txt" in files:
            mob.append(d / "IonMobility_Ar+_Ar.txt")
        if "p3_replay_full_spectral_muon.py" in files and "p3_spectral_weighting_lut.py" in files:
            numw.append(d)
        for f in files:
            if f.endswith(".npz"):
                luts.append(d / f)
        if {"ions.csv", "vertices.f64", "complete.txt"} <= files:
            exports.append(d)
        if {"waveforms.csv", "channels.csv", "summary.json"} <= files:
            replays.append(d)
    out = [f"# T-ladder dependency scan, roots: {roots}"]
    out.append("\n## 1 exporter executables")
    out += [f"  {p}  sha256={sha(p)}" for p in exe] or ["  NONE FOUND -> build target mwpc_stage_c_p3_export_full_spectral_paths "
                                                       "from feature/p3-2100v-wire-aware-stereo-weighting-20261120 (test/mwpc-stage-c);"
                                                       " compiling unchanged source does not change physics"]
    out.append("\n## 2 mobility tables")
    out += [f"  {p}  sha256={sha(p)}" for p in mob] or ["  NONE FOUND -> must be the file used by the frozen replays; do not substitute"]
    out.append("\n## 3 spectral LUTs (npz with near/bulk/meta)")
    good = []
    try:
        import numpy as np
        for p in luts:
            try:
                with np.load(p, allow_pickle=False) as z:
                    if not {"near", "bulk", "meta"} <= set(z.files):
                        continue
                    meta = json.loads(z["meta"].item())
            except Exception:
                continue
            same = meta == DEFAULT_META
            good.append((p, same))
            out.append(f"  {p}  sha256={sha(p)}  meta==default_build: {same}")
            if not same:
                out.append(f"    meta={json.dumps(meta)}")
    except ImportError:
        out.append("  numpy missing: cannot inspect npz files")
    if not good:
        out.append("  NONE FOUND -> rebuild: (cd NUMW && python3 p3_spectral_weighting_lut.py build NEW.npz); default "
                   "settings reproduce the documented model/grid (meta below); ~4 min, 21 MiB")
    out.append("\n## 4 replay code directories")
    out += [f"  {d}  replay sha256={sha(d / 'p3_replay_full_spectral_muon.py')}" for d in numw] or ["  NONE FOUND"]
    out.append("\n## 5a path exports (rows, exact_ions, vertices)")
    for d in exports:
        s = seal(d)
        out.append(f"  {d}  rows={s.get('rows')} exact_ions={s.get('exact_ions')} vertices={s.get('vertices')}"
                   + ("  [T-ladder variant]" if (d / "ablation.json").exists() else ""))
    out.append("\n## 5b replay outputs (rows, exact_ions, wall_s)")
    for d in replays:
        try:
            s = json.loads((d / "summary.json").read_text())
            out.append(f"  {d}  rows={s.get('rows')} exact_ions={s.get('exact_ions')} wall_s={s.get('wall_s', 0):.1f} "
                       f"electron_current_included={s.get('electron_current_included')}")
        except Exception:
            out.append(f"  {d}  (summary.json unreadable)")
    out.append("\n## default LUT build meta\n  " + json.dumps(DEFAULT_META))
    rep = "\n".join(out) + "\n"
    Path("t_ladder_deps_report.txt").write_text(rep)
    env = ["# candidates found by locate_t_ladder_deps.py -- REVIEW before sourcing"]
    if exe: env.append(f"export MWPC_EXPORT_BIN='{exe[0]}'")
    if mob: env.append(f"export MWPC_ION_MOBILITY='{mob[0]}'")
    if good: env.append(f"export MWPC_SPECTRAL_LUT='{sorted(good, key=lambda x: not x[1])[0][0]}'")
    if numw: env.append(f"export MWPC_NUMW_DIR='{numw[0]}'")
    Path("t_ladder_deps.env").write_text("\n".join(env) + "\n")
    print(rep)
    print("wrote t_ladder_deps_report.txt and t_ladder_deps.env")


if __name__ == "__main__":
    main()

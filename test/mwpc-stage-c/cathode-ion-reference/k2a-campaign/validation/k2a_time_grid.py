#!/usr/bin/env python3
"""Require ONE time grid for every compared waveform (S0 and K2A variants), explicitly; FAIL on any mismatch.

  k2a_time_grid.py OUT.json --window-ns W NAME=DIR [NAME=DIR ...]

Each DIR has replay/waveforms.csv and prompt/prompt_waveforms.csv (first column time_center_ns). Required, for every
variant and both files, against the FIRST variant's replay grid:
  * at least 2 bins, strictly uniform spacing (every interval equal to the first one, exactly)
  * the same sampling interval dt and the same first bin centre (exact floating-point equality)
  * at least nb = W / dt bins (W must be an integer number of bins), i.e. the comparison window is fully covered
  * the first nb time values identical (exact equality) to the reference's
  * replay and prompt of the same variant have identical full time columns and the same number of data columns
The comparison window W is that of the fast model (k2a_strip_eval window_ns); variants with more bins (the frozen S0
replay clock) are compared on the first nb bins only, and the fraction of their |ion current| beyond the window is
recorded (information). Nothing is resampled or interpolated. Exit 0 PASS, 1 FAIL (reasons listed in OUT.json).
check() is also used by k2a_label_union.py.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

FILES = ("replay/waveforms.csv", "prompt/prompt_waveforms.csv")


def read(p):
    with open(p) as f:
        h = f.readline().rstrip("\n").split(",")
    a = np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)
    return h, a


def check(variants, window_ns):
    """variants: list of (name, Path). Returns (ok, report dict, nb)."""
    bad, rep, ref = [], {"window_ns": window_ns, "variants": {}}, None
    for n, d in variants:
        vr = {"dir": str(d)}
        rep["variants"][n] = vr
        cols = {}
        for f in FILES:
            p = Path(d) / f
            if not p.is_file():
                bad.append(f"{n}: missing {f}"); continue
            try:
                h, a = read(p)
            except ValueError as e:
                bad.append(f"{n}/{f}: unreadable ({e})"); continue
            if a.shape[1] != len(h):
                bad.append(f"{n}/{f}: {a.shape[1]} columns but {len(h)} header fields"); continue
            t = a[:, 0]
            cols[f] = (t, a)
            vr[f] = {"bins": int(len(t)), "t0_ns": float(t[0]) if len(t) else None,
                     "dt_ns": float(t[1] - t[0]) if len(t) > 1 else None, "t_last_ns": float(t[-1]) if len(t) else None}
            if len(t) < 2:
                bad.append(f"{n}/{f}: fewer than 2 bins"); continue
            dts = np.diff(t)
            if not np.all(np.isfinite(t)) or not np.all(dts == dts[0]) or not dts[0] > 0:
                bad.append(f"{n}/{f}: time axis not strictly uniform (intervals {np.unique(dts)[:5].tolist()})")
        if len(cols) == 2:
            (tr, ar), (tp, ap_) = cols[FILES[0]], cols[FILES[1]]
            if len(tr) != len(tp) or not np.array_equal(tr, tp):
                bad.append(f"{n}: replay and prompt time axes differ ({len(tr)} vs {len(tp)} bins)")
            if ar.shape[1] != ap_.shape[1]:
                bad.append(f"{n}: replay and prompt have different strip counts")
        if FILES[0] in cols and ref is None:
            ref = (n, cols[FILES[0]][0])
        if FILES[0] in cols:
            vr["_t"] = cols[FILES[0]][0]; vr["_I"] = cols[FILES[0]][1][:, 1:]
    nb = None
    if ref is not None and len(ref[1]) > 1:
        rn, rt = ref
        dt, t0 = rt[1] - rt[0], rt[0]
        x = window_ns / dt
        if not (x > 0 and abs(x - round(x)) < 1e-9):
            bad.append(f"window {window_ns} ns is not an integer number of {dt} ns bins")
        else:
            nb = int(round(x))
            rep["dt_ns"], rep["t0_ns"], rep["bins_compared"], rep["reference"] = float(dt), float(t0), nb, rn
            for n, d in variants:
                vr = rep["variants"][n]
                t = vr.pop("_t", None); I = vr.pop("_I", None)
                if t is None or len(t) < 2:
                    continue
                if t[1] - t[0] != dt:
                    bad.append(f"{n}: sampling interval {t[1] - t[0]} ns != {dt} ns ({rn})")
                if t[0] != t0:
                    bad.append(f"{n}: first bin centre {t[0]} ns != {t0} ns ({rn})")
                if len(t) < nb:
                    bad.append(f"{n}: {len(t)} bins do not cover the {window_ns} ns comparison window ({nb} bins)")
                elif not np.array_equal(t[:nb], rt[:nb]):
                    k = int(np.flatnonzero(t[:nb] != rt[:nb])[0])
                    bad.append(f"{n}: time value of bin {k} is {t[k]} ns, {rn} has {rt[k]} ns")
                tot = float(np.abs(I).sum())
                vr["ion_abs_current_beyond_window_frac"] = float(np.abs(I[nb:]).sum()) / tot if tot > 0 and len(t) > nb else 0.0
    for v in rep["variants"].values():
        v.pop("_t", None); v.pop("_I", None)
    rep["status"] = "FAIL" if bad else "PASS"
    rep["failures"] = bad
    return not bad, rep, nb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--window-ns", type=float, required=True); ap.add_argument("variants", nargs="+")
    a = ap.parse_args()
    vs = [(x.split("=", 1)[0], Path(x.split("=", 1)[1])) for x in a.variants]
    ok, rep, nb = check(vs, a.window_ns)
    rep["tool"] = "k2a_time_grid.py"
    a.out.write_text(json.dumps(rep, indent=1) + "\n")
    if ok:
        print(f"PASS time grid: {len(vs)} variants share dt {rep['dt_ns']:g} ns, first bin {rep['t0_ns']:g} ns, "
              f"{nb} bins in the {a.window_ns:g} ns comparison window (" + ", ".join(f"{n} {v['replay/waveforms.csv']['bins']}" for n, v in rep["variants"].items()) + " bins)")
        return 0
    print("FAIL time grid: " + "; ".join(rep["failures"][:4]) + (" ..." if len(rep["failures"]) > 4 else ""))
    return 1


if __name__ == "__main__":
    sys.exit(main())

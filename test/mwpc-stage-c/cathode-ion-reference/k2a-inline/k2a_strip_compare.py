#!/usr/bin/env python3
"""Compare two absolute-strip K2A responses in replay format (exit 0 PASS, 1 FAIL, 3 INCONCLUSIVE).

  k2a_strip_compare.py REF_DIR TEST_DIR OUT.json [--tol 1e-12] [--label-shift M_U M_V] [--name TEXT]

REF_DIR / TEST_DIR   replay/waveforms.csv (ion current, fC/ns) and prompt/prompt_waveforms.csv (electron step), as
                     written by k2a_recentered_k2a.py, the frozen replay format and k2a_strip_eval
--label-shift        TEST labels are compared with REF labels + M (U_i <-> U_{i+M_U}): displacement symmetry test
Metrics, computed here from the waveforms only (independent of both evaluators), per plane and for ion / electron:
  labels      the absolute strip sets are identical (after the label shift)                         required
  time grid   identical time_center_ns column                                                        required
  L2          ||test - ref|| / ||ref|| over all strips and bins                                      <= tol
  strip L2    max over strips carrying > 1e-9 of the plane's L2 norm of the per-strip relative L2     <= tol
  charge      per-strip window charge q_s = sum_t i_s(t) dt: max |dq_s| / sum_s |q_s|, and the plane total
              relative difference                                                                    <= tol
  U/V         ratio of the plane window charges, relative difference                                  <= tol
Default tol 1e-12: the compiled adapter and the Python engine evaluate the same kernel_eval.cc arithmetic and differ
only by last-ulp rounding of the per-cloud coefficients (numpy/BLAS vs C++ summation order of exactly-zero dot products
of the hemisphere constants, ~1e-15 relative); 1e-12 leaves a factor > 100 margin and is 9 orders of magnitude below the
0.5 % physics warning level of the E2 comparisons.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def load(p):
    with open(p) as f:
        h = f.readline().strip().split(",")
    return h[1:], np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref", type=Path); ap.add_argument("test", type=Path); ap.add_argument("out", type=Path)
    ap.add_argument("--tol", type=float, default=1e-12); ap.add_argument("--label-shift", type=int, nargs=2, default=(0, 0))
    ap.add_argument("--name", default="")
    a = ap.parse_args()
    rep = {"tool": "k2a_strip_compare.py", "name": a.name, "ref": str(a.ref), "test": str(a.test), "tol": a.tol,
           "label_shift": list(a.label_shift), "checks": {}, "metrics": {}}
    files = [("ion", "replay/waveforms.csv"), ("electron", "prompt/prompt_waveforms.csv")]
    for _, f in files:
        for d in (a.ref, a.test):
            if not (d / f).is_file():
                rep["status"] = "INCONCLUSIVE"; rep["missing"] = str(d / f)
                a.out.write_text(json.dumps(rep, indent=1) + "\n"); print(f"INCONCLUSIVE {a.name}: missing {d / f}"); return 3
    ok = True
    q = {}
    for term, f in files:
        hr, R = load(a.ref / f); ht, T = load(a.test / f)
        shift = {"U": a.label_shift[0], "V": a.label_shift[1]}
        mapped = [f"{l[0]}_{int(l.split('_')[1]) - shift[l[0]]}" for l in ht]   # test label -> reference label
        same_labels = sorted(mapped) == sorted(hr) and len(set(mapped)) == len(mapped)
        same_t = R.shape[0] == T.shape[0] and np.array_equal(R[:, 0], T[:, 0])
        rep["checks"][f"{term}_labels"] = bool(same_labels); rep["checks"][f"{term}_time_grid"] = bool(same_t)
        ok &= same_labels and same_t
        if not (same_labels and same_t):
            rep["metrics"][term] = {"ref_labels": hr, "test_labels_mapped": mapped}
            continue
        col = {l: i for i, l in enumerate(hr)}
        Tm = np.zeros_like(R[:, 1:])
        for j, l in enumerate(mapped):
            Tm[:, col[l]] = T[:, 1 + j]
        Rv = R[:, 1:]
        dt = float(R[1, 0] - R[0, 0]) if len(R) > 1 else 1.0
        m = {}
        for s in ("U", "V"):
            idx = [i for i, l in enumerate(hr) if l[0] == s]
            r, t = Rv[:, idx], Tm[:, idx]
            nr = np.linalg.norm(r)
            m[f"{s}_L2"] = float(np.linalg.norm(t - r) / nr) if nr > 0 else float(np.linalg.norm(t - r))
            ns = np.linalg.norm(r, axis=0)
            big = ns > 1e-9 * nr
            m[f"{s}_strip_L2_max"] = float(np.max(np.linalg.norm(t - r, axis=0)[big] / ns[big])) if big.any() else 0.0
            qr, qt = r.sum(0) * dt, t.sum(0) * dt
            Qa = np.abs(qr).sum()
            m[f"{s}_strip_charge_max_rel"] = float(np.max(np.abs(qt - qr)) / Qa) if Qa > 0 else float(np.max(np.abs(qt - qr), initial=0))
            m[f"{s}_total_charge_ref_fC"] = float(qr.sum()); m[f"{s}_total_charge_test_fC"] = float(qt.sum())
            m[f"{s}_total_charge_rel"] = float(abs(qt.sum() - qr.sum()) / abs(qr.sum())) if qr.sum() != 0 else float(abs(qt.sum()))
            m[f"{s}_strips"] = len(idx)
            m[f"{s}_bitwise_equal"] = bool(np.array_equal(r, t))
            q[(term, s)] = (qr.sum(), qt.sum())
        if q[(term, "V")][0] != 0 and q[(term, "V")][1] != 0:
            ur, ut = q[(term, "U")][0] / q[(term, "V")][0], q[(term, "U")][1] / q[(term, "V")][1]
            m["UV_ratio_rel"] = float(abs(ut - ur) / abs(ur)) if ur != 0 else float(abs(ut))
        for k_, v in m.items():
            if k_.endswith(("_L2", "_L2_max", "_charge_max_rel", "_total_charge_rel", "UV_ratio_rel")):
                good = v <= a.tol
                rep["checks"][f"{term}_{k_}"] = bool(good); ok &= good
        rep["metrics"][term] = m
    rep["status"] = "PASS" if ok else "FAIL"
    a.out.write_text(json.dumps(rep, indent=1) + "\n")
    mi = rep["metrics"].get("ion", {})
    print(f"{rep['status']} {a.name}: ion L2 U {mi.get('U_L2', float('nan')):.2e} V {mi.get('V_L2', float('nan')):.2e}, "
          f"strip charge U {mi.get('U_strip_charge_max_rel', float('nan')):.2e} V {mi.get('V_strip_charge_max_rel', float('nan')):.2e}, "
          f"electron L2 U {rep['metrics'].get('electron', {}).get('U_L2', float('nan')):.2e} (tol {a.tol:g})"
          + ("" if ok else "; failed: " + ", ".join(k for k, v in rep["checks"].items() if not v)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Diagnose stereo-coordinate outliers in a t_shaped_compare run (e.g. the 200 ns TW/TC u outlier).

  t_centroid_diagnostic.py COMPARE_DIR --tp 200 --variant TW [--abs-thr-fC 0.5,1,2] [--muon muon_NNNN]

Reads COMPARE_DIR/per_muon_Tp<tp>.jsonl (written by t_shaped_compare.py at or after the commit that adds U/V_peaks to
every row). Without --muon it takes the muon with the largest |u(variant) - u(reference)| under the DEFAULT centroid
(peak-weighted over strips with peak >= 2% of the reference leading peak, the original definition; nothing is
changed or hidden). Writes centroid_diagnostic_Tp<tp>_<variant>.md next to the jsonl:
  1. the muon and its original u, w shifts
  2. per plane: every strip's reference and variant peak, distance from the 2% inclusion threshold, membership change
  3. u, w with membership fixed from the reference vs threshold-dependent membership
  4. u, w for alternative centroids: lead +-1, lead +-2, all positive strips, and absolute thresholds (fC-equivalent,
     unit-peak stand-in shaper) given with --abs-thr-fC
Stereo: U = w - A u, V = w + A u, u = pitch (cV - cU) / (2A) with pitch 1.7 mm, A = 0.1, so 1 strip of centroid
difference is 8.5 mm of u.
"""
import argparse
import json
from pathlib import Path

import numpy as np

P, A = 1.7, 0.1


def cen(pk, idx, sel):
    s = np.sum(pk[sel])
    return float(np.sum(pk[sel] * idx[sel]) / s) if s > 0 else float("nan")


def uw(cu, cv):
    return P * (cv - cu) / (2 * A), 0.5 * P * (cu + cv)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("compare_dir", type=Path)
    ap.add_argument("--tp", default="200")
    ap.add_argument("--variant", default="TW")
    ap.add_argument("--muon", default=None)
    ap.add_argument("--incl-frac", type=float, default=0.02, help="original inclusion fraction of the reference lead")
    ap.add_argument("--abs-thr-fC", default="0.5,1,2,5")
    a = ap.parse_args()
    rows = [json.loads(l) for l in open(a.compare_dir / f"per_muon_Tp{a.tp}.jsonl")]
    refname = next(r["variant"] for r in rows if "thresholds_fCeq" in r)
    ref = {r["muon"]: r for r in rows if r["variant"] == refname}
    var = {r["muon"]: r for r in rows if r["variant"] == a.variant}
    if not var:
        raise SystemExit(f"variant {a.variant} not in the jsonl")
    if any("U_peaks" not in r for r in list(ref.values()) + list(var.values())):
        raise SystemExit("jsonl lacks per-strip peaks: re-run t_shaped_compare.py from this commit (no new simulation needed)")
    ranked = sorted(var, key=lambda m: -abs(var[m]["u_mm"]))
    m = a.muon or ranked[0]
    r0, r1 = ref[m], var[m]
    L = [f"# Centroid diagnostic: {a.variant} vs {refname}, Tp = {a.tp} ns (provisional CR-RC^3 stand-in)", "",
         f"Source: `{a.compare_dir}/per_muon_Tp{a.tp}.jsonl`; {len(var)} muons.", "",
         "Largest |du| under the original centroid: " + ", ".join(f"{k} {var[k]['u_mm']:+.4f} mm" for k in ranked[:5]), "",
         f"## 1. Muon {m}", "",
         f"Original definition (strips >= {a.incl_frac:g} x reference lead peak, per plane): du = {r1['u_mm']:+.5f} mm, "
         f"dw = {r1['w_mm']:+.5f} mm; reference u = {r0['u_mm']:.4f} mm, w = {r0['w_mm']:.4f} mm.", ""]
    cents = {}
    for s in ("U", "V"):
        idx = np.array(r0[f"{s}_strip_index"])
        p0, p1 = np.array(r0[f"{s}_peaks"]), np.array(r1[f"{s}_peaks"])
        thr = a.incl_frac * p0.max()
        m0, m1 = p0 >= thr, p1 >= thr
        L += [f"## 2. Plane {s}: inclusion threshold {thr:.5g} fC-eq ({a.incl_frac:g} x reference lead {p0.max():.5g})", "",
              "| strip | ref peak | var peak | ref - thr | var - thr | ref in | var in |", "|---|---|---|---|---|---|---|"]
        for i in range(len(idx)):
            flag = " **changed**" if m0[i] != m1[i] else ""
            L.append(f"| {idx[i]} | {p0[i]:.5g} | {p1[i]:.5g} | {p0[i] - thr:+.4g} | {p1[i] - thr:+.4g} | {int(m0[i])} | {int(m1[i])}{flag} |")
        l0 = int(np.argmax(p0))
        near = lambda h: np.abs(np.arange(len(idx)) - l0) <= h
        cents[s] = {"original (threshold-dependent)": (cen(p0, idx, m0), cen(p1, idx, m1)),
                    "membership fixed from reference": (cen(p0, idx, m0), cen(p1, idx, m0)),
                    "lead +-1": (cen(p0, idx, near(1) & (p0 > 0)), cen(p1, idx, near(1) & (p1 > 0))),
                    "lead +-2": (cen(p0, idx, near(2) & (p0 > 0)), cen(p1, idx, near(2) & (p1 > 0))),
                    "all positive strips": (cen(p0, idx, p0 > 0), cen(p1, idx, p1 > 0))}
        for q in [float(x) for x in a.abs_thr_fC.split(",") if x]:
            cents[s][f"absolute {q:g} fC-eq"] = (cen(p0, idx, p0 >= q), cen(p1, idx, p1 >= q))
        L.append("")
    L += ["## 3-4. Reconstructed coordinates", "",
          "| centroid definition | u ref [mm] | du [mm] | dw [mm] |", "|---|---|---|---|"]
    for k in cents["U"]:
        (cu0, cu1), (cv0, cv1) = cents["U"][k], cents["V"][k]
        u0, w0 = uw(cu0, cv0); u1, w1 = uw(cu1, cv1)
        L.append(f"| {k} | {u0:.4f} | {u1 - u0:+.5f} | {w1 - w0:+.5f} |")
    L += ["", "Interpretation aid: a strip of peak p at distance d (strips) from the centroid entering a plane sum S "
          "moves that plane's centroid by about p d / S strips, i.e. 8.5 p d / S mm of u."]
    out = a.compare_dir / f"centroid_diagnostic_Tp{a.tp}_{a.variant}.md"
    out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

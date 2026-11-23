#!/usr/bin/env python3
"""APPROXIMATE avalanche-electron term: an instantaneous step at the ion birth time.

What is exact and what is not
* The frozen replay (p3_replay_full_spectral_muon.py) is ION-ONLY and, for the frozen RKF paths and spectral LUT, it is
  the time-dependent Shockley-Ramo ion signal: Q_ion,k(t) = +q [phi_k(r_ion(t)) - phi_k(r_birth)] (chargeSign = +1).
* The true electron signal is also time dependent. In the same convention (carrier charge -q), each avalanche
  electron gives Q_e,k(t) = -q [phi_k(r_e(t)) - phi_k(r_e(t_create))], following its own path to the anode; the
  primary electron's drift adds one more such term. Summed over the avalanche and taken to completion this is
  +q phi_k(r_birth) per ion/electron pair (phi_k = 0 on the anode; electron created where its partner ion is born). The archives do not contain avalanche-electron
  trajectories, so this integral is all that can be formed from them.
* THIS helper replaces i_e(t) by its integral, deposited as a step in the 10 ns bin of the ion birth time. Neglected:
  the electron transit/avalanche development time (sub-ns to ns scale, not measured here), the primary-electron drift
  current (one electron, not multiplied), and any shape of the electron current inside the bin.
  Use it only AFTER the ion-only comparison, and report both. Sum ion + step = +q phi_k(r_ion(t)) at bin resolution.

Uses the replay's own LUT, query() and channel frame (imported from --replay-dir), so labels and strip bases match
waveforms.csv.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", type=Path)
    ap.add_argument("export", type=Path)
    ap.add_argument("out", type=Path, help="NEW directory; writes prompt_waveforms.csv and prompt_summary.json")
    ap.add_argument("--replay-dir", type=Path, required=True, help="directory with p3_replay_full_spectral_muon.py")
    ap.add_argument("--half-strips", type=int, default=5)
    a = ap.parse_args()
    sys.path.insert(0, str(a.replay_dir.resolve()))
    import p3_replay_full_spectral_muon as rp  # noqa: E402  (frozen replay, unchanged)

    rows, vertices, weights, exact = rp.read_export(a.export)
    offsets = np.arange(-a.half_strips, a.half_strips + 1)
    nc = len(offsets)
    bases = np.array([int(rows[0]["upper_index"]), int(rows[0]["lower_index"])])
    lut = rp.WeightingLUT(a.table)
    a.out.mkdir(parents=True, exist_ok=False)
    total = np.zeros((rp.NBINS, 2, nc))
    sum_birth = np.zeros((2, nc))
    for r, wt in zip(rows, weights):
        if [int(r["upper_index"]), int(r["lower_index"])] != bases.tolist():
            raise ValueError("inconsistent physical strip frame")
        p = rp.path_for(r, vertices)
        phi0 = rp.query(lut, p, np.array([p[0, 0]]), offsets)[0]
        q = wt * rp.EFC
        b = int(np.floor(p[0, 0] / rp.DT))
        if b < rp.NBINS:
            total[b] += q * phi0
        sum_birth += q * phi0
    labels = [f"{s}_{base + int(k)}" for s, base in zip(("U", "V"), bases) for k in offsets]
    np.savetxt(a.out / "prompt_waveforms.csv",
               np.column_stack((rp.EDGES[:-1] + rp.DT / 2, total.reshape(rp.NBINS, -1) / rp.DT)),
               delimiter=",", header="time_center_ns," + ",".join(labels), comments="")
    summ = {"rows": len(rows), "exact_ions": exact, "prompt_charge_fC": dict(zip(labels, sum_birth.ravel().tolist())),
            "approximation": "instantaneous electron step: integral +q phi_k(r_birth) placed in the 10 ns bin of birth_ns; true time-dependent electron current not available from the archives"}
    (a.out / "prompt_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(f"instantaneous electron-step approximation: {len(rows)} rows, plane sums U={sum_birth[0].sum():.6g} fC V={sum_birth[1].sum():.6g} fC")


if __name__ == "__main__":
    main()

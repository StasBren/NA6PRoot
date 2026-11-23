#!/usr/bin/env python3
"""Avalanche-electron (prompt) term for a path export, on the replay's 10 ns grid and channel set.

The frozen replay (p3_replay_full_spectral_muon.py) is ION-ONLY: Q_ion,k(t) = +q [phi_k(r_ion(t)) - phi_k(r_birth)].
The avalanche electron paired with each ion moves from r_birth to the anode (phi_k(anode) = 0 for a cathode strip)
in well under 1 ns, adding Q_e,k = +q phi_k(r_birth) as a step at the birth time (frozen sign convention,
chargeSign = +1). Total = +q phi_k(r_ion(t)). Moving a birth time therefore moves BOTH terms; the T-ladder compare
uses ion + prompt.

Uses the replay's own LUT, query() and channel frame (imported from --replay-dir), so the channel labels and strip
bases are identical to waveforms.csv. Approximation: the electron term is a step at birth_ns (electron transit time
neglected; it is below the 10 ns bin).
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
            "convention": "Q_e,k = +q phi_k(r_birth) at birth_ns; replay ion term + this = +q phi_k(r_ion(t))"}
    (a.out / "prompt_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(f"prompt term: {len(rows)} rows, plane sums U={sum_birth[0].sum():.6g} fC V={sum_birth[1].sum():.6g} fC")


if __name__ == "__main__":
    main()

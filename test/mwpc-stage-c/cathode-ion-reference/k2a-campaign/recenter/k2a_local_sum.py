#!/usr/bin/env python3
"""Sum per-group frozen-replay outputs onto ABSOLUTE strip labels (linear superposition; nothing rescaled).

  k2a_local_sum.py GROUPS_DIR S0GROUPS_DIR OUT_VARIANT_DIR

GROUPS_DIR      k2a_group_export.py output (groups.json)
S0GROUPS_DIR    contains <group>/replay/{waveforms.csv,channels.csv,summary.json} and <group>/prompt/prompt_waveforms.csv
                (frozen p3_replay_full_spectral_muon.py and t_prompt_electron_term.py run on each sub-export)
OUT_VARIANT_DIR replay/waveforms.csv + prompt/prompt_waveforms.csv on the union of labels (U then V, ascending index),
                replay/channels.csv summed per absolute strip, local_sum.json with the census check:
                sum over groups of collected (landing) ion charge == exact ions of the source x e (all within windows).
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

EFC = 1.602176634e-4


def load(p):
    with open(p) as f:
        h = f.readline().strip().split(",")
    return h[1:], np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)


def main():
    G, S, O = (Path(x) for x in sys.argv[1:4])
    if O.exists():
        raise SystemExit(f"refusing: {O} exists")
    meta = json.loads((G / "groups.json").read_text())
    acc = {"replay": defaultdict(lambda: 0.0), "prompt": defaultdict(lambda: 0.0)}
    tgrid = None
    chan = defaultdict(lambda: np.zeros(6))
    exact = 0; landed = 0.0
    for name, g in meta["groups"].items():
        for sub, fn in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv")):
            lab, d = load(S / name / sub / fn)
            if tgrid is None:
                tgrid = d[:, 0]
            elif not np.array_equal(tgrid, d[:, 0]):
                raise SystemExit("time grids differ between groups")
            for j, l in enumerate(lab):
                acc[sub][l] = acc[sub][l] + d[:, 1 + j]
        for r in csv.DictReader(open(S / name / "replay" / "channels.csv", newline="")):
            k = f"{r['side']}_{r['strip_index']}"
            chan[k] += np.array([float(r[c]) for c in ("window_charge_fC", "endpoint_charge_fC", "unrecorded_tail_fC",
                                                        "collected_ion_charge_fC", "first100ns_charge_fC", "waveform_L2_fC_per_ns")])
            landed += float(r["collected_ion_charge_fC"])
        exact += g["exact_ions"]
    order = sorted(acc["replay"], key=lambda l: (l[0] != "U", int(l.split("_")[1])))
    for sub, fn in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv")):
        (O / sub).mkdir(parents=True)
        M = np.column_stack([tgrid] + [acc[sub][l] if not np.isscalar(acc[sub][l]) else np.zeros_like(tgrid) for l in order])
        np.savetxt(O / sub / fn, M, delimiter=",", header="time_center_ns," + ",".join(order), comments="")
    with open(O / "replay" / "channels.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["side", "strip_index", "window_charge_fC", "endpoint_charge_fC", "unrecorded_tail_fC",
                                       "collected_ion_charge_fC", "first100ns_charge_fC", "sum_of_group_L2"])
        for l in order:
            s, i = l.split("_"); w.writerow([s, i] + [repr(float(x)) for x in chan[l]])
    census = {"exact_ions": exact, "expected_collected_fC": exact * EFC, "collected_fC": landed,
              "rel_diff": abs(landed - exact * EFC) / (exact * EFC), "labels": order, "groups": list(meta["groups"])}
    (O / "local_sum.json").write_text(json.dumps(census, indent=1) + "\n")
    ok = census["rel_diff"] < 1e-9
    print(("PASS" if ok else "FAIL") + f" local S0 sum: {len(meta['groups'])} groups, {len(order)} absolute strips, "
          f"collected ion charge {landed:.9g} fC vs census {exact * EFC:.9g} fC (rel {census['rel_diff']:.2g})")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

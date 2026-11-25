#!/usr/bin/env python3
"""K2A with per-primary strip recentering, accumulated on ABSOLUTE strip indices, using the UNCHANGED engine.

  k2a_recentered_k2a.py KERNEL_DIR MOMENTS.csv OUT_VARIANT_DIR [--groups GROUPS.json] [--model K2A] [--impl vector]
                        [--half-strips 5] [--window-ns 7000] [--kernel-id ID]

For each group of clouds with the same (m_U, m_V) the engine is called with D'_U = D_U - m_U P, D'_V = D_V - m_V P and
base_U + m_U, base_V + m_V (a Moments copy; the kernel, interpolation and evaluation code are the frozen engine's).
Responses are added per absolute label. m per primary comes from GROUPS.json (k2a_group_export.py, so the local S0
reference uses exactly the same windows) or, without it, from the moments: floor(Dbar/P + 1/2) per primary,
charge-weighted (the same rule). No clamping: a cloud outside the kernel after recentering raises as in the engine.
Writes replay/waveforms.csv, prompt/prompt_waveforms.csv (union labels, U then V ascending), recentered.json.
"""
import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "spatial-ablation" / "engine"))
import cathode_engine as ce  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kernel", type=Path); ap.add_argument("moments", type=Path); ap.add_argument("out", type=Path)
    ap.add_argument("--groups", type=Path); ap.add_argument("--model", default="K2A"); ap.add_argument("--impl", default="vector")
    ap.add_argument("--half-strips", type=int, default=5); ap.add_argument("--window-ns", type=float, default=7000.0)
    ap.add_argument("--kernel-id")
    ap.add_argument("--per-group-out", type=Path, help="also write each group's response to DIR/<g_mU_mV>/ (same format)")
    a = ap.parse_args()
    if a.out.exists() or (a.per_group_out and a.per_group_out.exists()):
        raise SystemExit(f"refusing: {a.out} or the per-group directory exists")
    k = ce.KernelTable(a.kernel, expected_id=a.kernel_id)
    eng = ce.CathodeResponseEngine(k, model=a.model, half_strips=a.half_strips, window_ns=a.window_ns)
    m = ce.read_moments_csv(a.moments)
    X = m.cols; P = k.pitch_mm
    prim = X["primary"].astype(int)
    if a.groups:
        ps = {int(p): tuple(v) for p, v in json.loads(a.groups.read_text())["primary_shift"].items()}
    else:
        ps = {}
        for p in np.unique(prim):
            s = prim == p; w = X["N_fC"][s] / X["N_fC"][s].sum()
            ps[int(p)] = (int(np.floor((w * X["D_U_mm"][s]).sum() / P + .5)), int(np.floor((w * X["D_V_mm"][s]).sum() / P + .5)))
    key = np.array([ps[int(p)] for p in prim])
    t0 = time.perf_counter()
    ion, ele, tc = defaultdict(lambda: 0.0), defaultdict(lambda: 0.0), None
    groups = {}
    for mu, mv in sorted({tuple(x) for x in key.tolist()}):
        s = (key[:, 0] == mu) & (key[:, 1] == mv)
        cols = {c: v[s].copy() for c, v in X.items()}
        cols["D_U_mm"] = cols["D_U_mm"] - P * mu; cols["D_V_mm"] = cols["D_V_mm"] - P * mv
        r = eng.evaluate(ce.Moments(cols, m.base_U + mu, m.base_V + mv, f"group {mu},{mv}"), impl=a.impl)
        tc = r.t_center_ns
        for j, l in enumerate(r.labels):
            si, o = (0 if l[0] == "U" else 1), j % len(eng.offsets)
            ion[l] = ion[l] + r.ion[:, si, o]; ele[l] = ele[l] + r.electron[:, si, o]
        groups[f"g_{mu}_{mv}"] = int(s.sum())
        if a.per_group_out:
            gd = a.per_group_out / f"g_{mu}_{mv}"
            for sub, arr, fn in (("replay", r.ion, "waveforms.csv"), ("prompt", r.electron, "prompt_waveforms.csv")):
                (gd / sub).mkdir(parents=True)
                np.savetxt(gd / sub / fn, np.column_stack([tc] + [arr[:, 0 if l[0] == "U" else 1, j % len(eng.offsets)]
                                                                  for j, l in enumerate(r.labels)]),
                           delimiter=",", header="time_center_ns," + ",".join(r.labels), comments="")
    wall = time.perf_counter() - t0
    order = sorted(ion, key=lambda l: (l[0] != "U", int(l.split("_")[1])))
    for sub, acc, fn in (("replay", ion, "waveforms.csv"), ("prompt", ele, "prompt_waveforms.csv")):
        (a.out / sub).mkdir(parents=True)
        np.savetxt(a.out / sub / fn, np.column_stack([tc] + [acc[l] for l in order]), delimiter=",",
                   header="time_center_ns," + ",".join(order), comments="")
    (a.out / "recentered.json").write_text(json.dumps({"kernel_id": k.kernel_id, "model": a.model, "impl": a.impl,
        "half_strips": a.half_strips, "window_ns": a.window_ns, "clouds_per_group": groups, "labels": order,
        "evaluate_wall_s": wall, "shifts_from": str(a.groups) if a.groups else "moments"}, indent=1) + "\n")
    print(f"PASS recentered {a.model}: {m.n} clouds in {len(groups)} groups -> {len(order)} absolute strips, {wall * 1e3:.2f} ms")


if __name__ == "__main__":
    main()

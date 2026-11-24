#!/usr/bin/env python3
"""Plan the S-ladder expansion: pick representative non-sealed development muons and estimate cost before running.

  s_plan_expansion.py ARCHIVE DEV.csv SEALED.csv RUN --tag TAG [--n 5] [--include 0000] [--nr 14 --ntheta 72]

READ-ONLY. Uses only the exact per-(primary, wire) census `muon_NNNN/microscopic_wire_seed.csv` of each candidate
(no ROOT, no export) and the measured costs of the existing RUN (green_TAG/axes.json, timing.txt, S0 replay
summary.json, s_synthesize_*.json). Sealed indices are removed before any file is opened.

Selection (deterministic): always the --include muons (default 0000, already prepared), then greedily the muon that
maximizes the distance to the already chosen ones in the standardized space of
  log(total ion births), number of anodes with >= 5% of the births, share of births off the dominant anode,
  number of anodes with ion births,
so that events with births around different / several anodes are represented. The table prints every candidate.

Cost model (scaled from the measured RUN, labelled as an estimate):
  K256 rows per muon            = sum_p min(256, ion births of primary p)
  S0 replay time                ~ rows x (measured s per row)
  Green table bytes             = 2 sides x wires x nr x ntheta x nD x ntau x 4 B (nD, ntau as in the existing table)
  Green build time              ~ grid paths x (measured s per grid path), RKF paths ~ grid paths x measured
  synthesis memory per variant  ~ rows x contributions per ion x 40 B (S1-S3 type: 2 x ntheta x 2; G/K type: 8)
"""
import argparse
import csv
import glob
import json
from pathlib import Path

import numpy as np


def idx(path):
    with open(path, newline="") as f:
        return [f"{int(r['index']):04d}" for r in csv.DictReader(f) if r.get("index", "").strip()]


def census(path):
    births, coll = {}, {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            p, g = int(r["primary_index"]), int(r["wire_index_global"])
            births.setdefault(p, {}).setdefault(g, 0)
            births[p][g] += int(r["ion_births"])
    perwire, rows = {}, 0
    for p, d in births.items():
        rows += min(256, sum(d.values())) if sum(d.values()) > 0 else 0
        for g, n in d.items():
            perwire[g] = perwire.get(g, 0) + n
    tot = sum(perwire.values())
    shares = sorted((n / tot for n in perwire.values() if n > 0), reverse=True) if tot else [0]
    return dict(births=tot, rows=rows, wires_any=sum(1 for n in perwire.values() if n > 0),
                wires_5pct=sum(1 for s in shares if s >= 0.05), off_dominant=1 - shares[0],
                wire_list=sorted(g for g, n in perwire.items() if n > 0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", type=Path); ap.add_argument("dev", type=Path); ap.add_argument("sealed", type=Path)
    ap.add_argument("run", type=Path); ap.add_argument("--tag", required=True)
    ap.add_argument("--n", type=int, default=5); ap.add_argument("--include", default="0000")
    ap.add_argument("--nr", type=int, default=14); ap.add_argument("--ntheta", type=int, default=72)
    ap.add_argument("--wires", type=int, default=3, help="anodes expected in the shared Green table")
    a = ap.parse_args()
    sealed = set(idx(a.sealed))
    cands = [m for m in idx(a.dev) if m not in sealed]
    info = {}
    for m in cands:
        c = a.archive / f"muon_{m}" / "microscopic_wire_seed.csv"
        if c.exists():
            info[m] = census(c)
    keys = sorted(info)
    F = np.array([[np.log(info[m]["births"] + 1), info[m]["wires_5pct"], info[m]["off_dominant"], info[m]["wires_any"]]
                  for m in keys])
    Z = (F - F.mean(0)) / np.where(F.std(0) > 0, F.std(0), 1)
    chosen = [m for m in a.include.split(",") if m in info]
    while len(chosen) < min(a.n, len(keys)):
        cz = Z[[keys.index(m) for m in chosen]] if chosen else np.zeros((1, Z.shape[1]))
        d = np.min(np.linalg.norm(Z[:, None, :] - cz[None, :, :], axis=2), axis=1)
        d[[keys.index(m) for m in chosen]] = -1
        chosen.append(keys[int(np.argmax(d))])
    print(f"# candidates: {len(keys)} non-sealed development muons with a census ({len(sealed)} sealed excluded)")
    print("muon  births     K256rows  anodes(any/>=5%)  off-dominant  chosen")
    for m in keys:
        i = info[m]
        print(f"{m}  {i['births']:9d}  {i['rows']:8d}  {i['wires_any']:3d}/{i['wires_5pct']:<3d}          "
              f"{i['off_dominant']:.4f}        {'*' if m in chosen else ''}")
    print(f"\nchosen: {','.join(chosen)}")
    # measured costs of the existing run
    ax = json.loads((a.run / f"green_{a.tag}" / "axes.json").read_text())
    shp = ax["shape"]                     # (2, k, r, theta, D, tau)
    paths_old = shp[1] * shp[2] * shp[3]
    timing = {}
    tf = a.run / "timing.txt"
    if tf.exists():
        for l in tf.read_text().split():
            if "=" in l:
                k_, v_ = l.split("=", 1); timing[k_] = float(v_)
    s0 = sorted(glob.glob(str(a.run / "muon_*" / "S0" / "replay" / "summary.json")))
    s0rate = []
    for f in s0:
        js = json.loads(Path(f).read_text())
        if js.get("rows") and js.get("wall_s"):
            s0rate.append(js["wall_s"] / js["rows"])
    syn = sorted(glob.glob(str(a.run / "muon_*" / "s_synthesize_*.json")))
    synrate = {}
    for f in syn:
        js = json.loads(Path(f).read_text())
        for v, d in js["variants"].items():
            synrate.setdefault(v, []).append(d["wall_s"] / js["rows"])
    rows_new = sum(info[m]["rows"] for m in chosen)
    paths_new = a.wires * a.nr * a.ntheta
    bytes_new = 2 * a.wires * a.nr * a.ntheta * shp[4] * shp[5] * 4
    print("\n# cost estimate (scaled from the measured run; ESTIMATE, not a measurement)")
    print(f"existing table green_{a.tag}: shape {shp}, {paths_old} grid paths, build {ax.get('build_s', float('nan')):.0f} s")
    print(f"shared table for the chosen muons: {a.wires} anodes x {a.nr} r x {a.ntheta} theta = {paths_new} grid paths, "
          f"{bytes_new / 2**20:.0f} MiB (+ PHI0 {bytes_new / shp[5] * 2 / 2**20:.1f} MiB); harmonic kernels "
          f"{7 / a.ntheta * bytes_new / 2**20:.0f} MiB")
    tb = ax.get("build_s")
    if tb:
        print(f"  table build ~ {tb * paths_new / paths_old / 60:.1f} min (LUT evaluation, scales with grid paths)")
    gp = timing.get(f"green_paths_{a.tag}_s")
    if gp:
        print(f"  RKF grid paths ~ {gp * paths_new / paths_old / 60:.1f} min")
    if s0rate:
        r_ = float(np.median(s0rate))
        print(f"S0 replays: {rows_new} K256 rows in total, ~ {rows_new * r_ / 60:.1f} min at {r_ * 1e3:.2f} ms/row measured"
              f" (+ export and T1 collapse per muon)")
    for v, lst in synrate.items():
        r_ = float(np.median(lst))
        print(f"synthesis {v:5s}: ~ {rows_new * r_:.0f} s for all chosen muons")
    big = max(info[m]["rows"] for m in chosen)
    print(f"synthesis peak memory (largest chosen muon, {big} rows): S1/S2/S3-type ~ {big * 2 * a.ntheta * 2 * 40 / 2**20:.0f} MiB "
          f"per plane; G/K-type ~ {big * 8 * 40 / 2**20:.1f} MiB; the Green table is memory-mapped")


if __name__ == "__main__":
    main()

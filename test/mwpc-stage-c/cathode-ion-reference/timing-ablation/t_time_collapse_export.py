#!/usr/bin/env python3
"""T-ladder timing ablation on a FROZEN RKF path export (p3_export_full_spectral_paths output).

Writes a NEW export directory in the identical format (ions.csv, vertices.f64, complete.txt) in which every ion path is
shifted rigidly in time so that all ions of a group start at the group's effective birth time:

  identity : no change (gate: output must be byte-identical to the source)
  primary  : group = Heed primary                -> T1 (one effective time per primary, spatial cloud kept)
  wire     : group = nearest-birth anode         -> timing part of T4 (one time per wire, cloud kept)
  muon     : group = whole crossing               -> one time per muon (extreme control)

Effective time = ion-weight mean of birth_ns over the group (weights = inverse inclusion weights).
Why the shifted paths are exact: the exporter drifts ions with Garfield RKF in a static 2D field with a field-dependent
but time-independent mobility, so a path depends only on its start point; shifting its start time shifts every vertex
time by the same amount. Positions and weights are not touched. The frozen replay then evaluates the real spectral
U/V weighting potential on the shifted paths, unchanged.
Read-only on the source; refuses to overwrite.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

PITCH_MM = 3.0


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def read_seal(d):
    return dict(line.split("=", 1) for line in (d / "complete.txt").read_text().splitlines() if "=" in line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path, help="frozen export directory (ions.csv, vertices.f64, complete.txt)")
    ap.add_argument("out", type=Path, help="NEW output directory")
    ap.add_argument("--mode", choices=["identity", "primary", "wire", "muon"], required=True)
    a = ap.parse_args()
    src, out = a.source, a.out
    if out.exists():
        raise SystemExit(f"refusing: {out} exists")
    seal = read_seal(src)
    nv = int(seal["vertices"])
    if (src / "vertices.f64").stat().st_size != 40 * nv:
        raise SystemExit("vertices.f64 size does not match complete.txt")
    with (src / "ions.csv").open(newline="") as f:
        rd = csv.DictReader(f)
        header = rd.fieldnames
        rows = list(rd)
    need = ["row", "primary", "weight", "birth_ns", "end_ns", "first_vertex", "vertices", "upper_index", "lower_index"]
    if header != need:
        raise SystemExit(f"unexpected ions.csv header {header}")
    V = np.memmap(src / "vertices.f64", dtype="<f8", mode="r", shape=(nv, 5))
    n = len(rows)
    w = np.array([float(r["weight"]) for r in rows])
    tb = np.array([float(r["birth_ns"]) for r in rows])
    first = np.array([int(r["first_vertex"]) for r in rows])
    cnt = np.array([int(r["vertices"]) for r in rows])
    prim = np.array([int(r["primary"]) for r in rows])
    if np.any(np.abs(V[first, 0] - tb) > 1e-7):
        raise SystemExit("first vertex time != birth_ns")
    wire = np.rint(V[first, 1] / PITCH_MM).astype(int)     # nearest anode to the ion birth (u in mm, patch frame)
    key = {"identity": np.arange(n), "primary": prim, "wire": wire, "muon": np.zeros(n, int)}[a.mode]
    shift = np.zeros(n)
    if a.mode != "identity":
        for g in np.unique(key):
            m = key == g
            shift[m] = np.sum(w[m] * tb[m]) / np.sum(w[m]) - tb[m]
    if np.any(tb + shift < 0):
        raise SystemExit("negative shifted birth time")
    out.mkdir(parents=True)
    with (out / "vertices.f64").open("wb") as f:
        for j in range(n):
            blk = np.array(V[first[j]:first[j] + cnt[j]], dtype="<f8")
            if shift[j] != 0.0:
                blk[:, 0] += shift[j]
            f.write(blk.tobytes())
    with (out / "ions.csv").open("w", newline="") as f:
        if a.mode == "identity":
            f.write((src / "ions.csv").read_text())
        else:
            f.write(",".join(need) + "\n")
            for r, s in zip(rows, shift):
                r = dict(r)
                r["birth_ns"] = repr(float(float(r["birth_ns"]) + s))
                r["end_ns"] = repr(float(float(r["end_ns"]) + s))
                f.write(",".join(r[k] for k in need) + "\n")
    (out / "complete.txt").write_text((src / "complete.txt").read_text())
    wsum = np.sum(w)
    meta = {"tool": "t_time_collapse_export.py", "mode": a.mode, "source": str(src.resolve()),
            "source_vertices_sha256": sha256(src / "vertices.f64"), "source_ions_sha256": sha256(src / "ions.csv"),
            "rows": n, "groups": int(len(np.unique(key))),
            "weighted_rms_shift_ns": float(np.sqrt(np.sum(w * shift ** 2) / wsum)),
            "max_abs_shift_ns": float(np.max(np.abs(shift))),
            "nearest_birth_wires": sorted(set(int(x) for x in wire)),
            "assumption": "static field, time-independent mobility: path depends only on start point"}
    (out / "ablation.json").write_text(json.dumps(meta, indent=2) + "\n")
    if a.mode == "identity":
        same = sha256(out / "vertices.f64") == meta["source_vertices_sha256"] and sha256(out / "ions.csv") == meta["source_ions_sha256"]
        print(("PASS" if same else "FAIL") + " identity export is byte-identical to the frozen source")
        if not same:
            raise SystemExit(3)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

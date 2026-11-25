#!/usr/bin/env python3
"""SYNTHETIC: the streaming k2a_group_export.py writes byte-identical exports to the previous in-memory version
(commit 9af160a), in both modes, and its peak memory stays bounded on a large export. Not NA60+ data.

  test_group_export_equivalence.py WORKDIR [--old-ref 9af160a] [--vertices 4000000]

1. a synthetic export (integer weight per primary, strip columns spread over several strips, two strip-index pairs)
2. old (git show OLD_REF:...) and new tool, modes `primary` and `shift --m-U 1 --m-V -2`: every output file
   (vertices.f64, ions.csv, complete.txt, groups.json) compared byte for byte
3. peak RSS of each run (separate processes); the new tool must stay below the export size, the old one is reported
"""
import argparse
import csv
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
NEW = HERE / "recenter" / "k2a_group_export.py"
REL = "test/mwpc-stage-c/cathode-ion-reference/k2a-campaign/recenter/k2a_group_export.py"


def make_export(d, nvert, rng):
    d.mkdir(parents=True)
    nrow = max(200, nvert // 200)
    counts = rng.integers(50, 350, nrow); counts = (counts * nvert / counts.sum()).astype(int) + 1
    nv = int(counts.sum())
    prim = np.sort(rng.integers(0, max(20, nrow // 50), nrow))
    V = np.empty((nv, 5), "<f8")
    off = 0; rows = []
    shift = {p: rng.uniform(-14, 21) for p in np.unique(prim)}       # strip columns spread like the far-D events
    for j in range(nrow):
        n = int(counts[j]); t0 = rng.uniform(0, 300)
        V[off:off + n, 0] = t0 + np.cumsum(rng.uniform(0.1, 50, n))
        V[off:off + n, 1] = rng.normal(3 * rng.integers(-2, 3), 0.05, n)
        V[off:off + n, 2] = np.linspace(rng.uniform(-0.1, 0.1), rng.choice([-3.0, 3.0]), n)
        V[off:off + n, 3] = shift[prim[j]] + rng.normal(0, 0.05)
        V[off:off + n, 4] = -shift[prim[j]] + rng.normal(0, 0.05)
        rows.append({"row": j, "primary": int(prim[j]), "weight": 0.0, "birth_ns": t0, "end_ns": float(V[off + n - 1, 0]),
                     "first_vertex": off, "vertices": n, "upper_index": 7, "lower_index": -4})
        off += n
    for p in np.unique(prim):                                        # integer census per primary: weight = N / saved
        sel = [r for r in rows if r["primary"] == p]; N = int(rng.integers(len(sel), 50 * len(sel)))
        for r in sel:
            r["weight"] = N / len(sel)
    V.tofile(d / "vertices.f64")
    with open(d / "ions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    ws = sum(r["weight"] for r in rows)
    (d / "complete.txt").write_text(f"rows={nrow}\nexact_ions={int(round(ws))}\nweight_sum={ws!r}\nvertices={nv}\nstep_mm=0.005\n")
    return nv


def run(tool, args):
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    subprocess.run([sys.executable, str(tool)] + args, check=True, stdout=subprocess.DEVNULL)
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, before


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("work", type=Path)
    ap.add_argument("--old-ref", default="9af160a"); ap.add_argument("--vertices", type=int, default=4_000_000)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=False)
    old = a.work / "k2a_group_export_old.py"
    repo = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--show-toplevel"], capture_output=True, text=True, check=True).stdout.strip()
    old.write_text(subprocess.run(["git", "-C", repo, "show", f"{a.old_ref}:{REL}"], capture_output=True, text=True, check=True).stdout)
    rng = np.random.default_rng(17)
    src = a.work / "export"; nv = make_export(src, a.vertices, rng)
    size = (src / "vertices.f64").stat().st_size
    print(f"synthetic export: {nv} vertices, {size / 2**20:.1f} MiB")
    rss = {}
    for mode, extra in (("primary", []), ("shift", ["--m-U", "1", "--m-V", "-2"])):
        n_new = a.work / f"new_{mode}"
        r_new, _ = run(NEW, [str(src), str(n_new), "--mode", mode] + extra)
        r_old, _ = run(old, [str(src), str(a.work / f"old_{mode}"), "--mode", mode] + extra)
        rss[mode] = r_new
        files = sorted(str(p.relative_to(n_new)) for p in n_new.rglob("*") if p.is_file())
        files_old = sorted(str(p.relative_to(a.work / f"old_{mode}")) for p in (a.work / f"old_{mode}").rglob("*") if p.is_file())
        assert files == files_old, (files, files_old)
        for f in files:
            assert (n_new / f).read_bytes() == (a.work / f"old_{mode}" / f).read_bytes(), f"{mode}: {f} differs"
        print(f"PASS {mode}: {len(files)} files byte-identical to {a.old_ref} "
              f"({', '.join(sorted({p.split('/')[0] for p in files}))})")
    # memory: measure the new tool alone in a fresh process (RUSAGE_CHILDREN is a running maximum)
    for mode, extra in (("primary", []), ("shift", ["--m-U", "1", "--m-V", "-2"])):
        out = subprocess.run([sys.executable, "-c", (
            "import resource, subprocess, sys; subprocess.run(sys.argv[1:], check=True, stdout=subprocess.DEVNULL); "
            "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)"),
            sys.executable, str(NEW), str(src), str(a.work / f"mem_new_{mode}"), "--mode", mode] + extra,
            capture_output=True, text=True, check=True).stdout.split()[-1]
        outo = subprocess.run([sys.executable, "-c", (
            "import resource, subprocess, sys; subprocess.run(sys.argv[1:], check=True, stdout=subprocess.DEVNULL); "
            "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)"),
            sys.executable, str(old), str(src), str(a.work / f"mem_old_{mode}"), "--mode", mode] + extra,
            capture_output=True, text=True, check=True).stdout.split()[-1]
        new_mib, old_mib = int(out) / 1024, int(outo) / 1024
        print(f"peak RSS {mode}: new {new_mib:.0f} MiB, old {old_mib:.0f} MiB, export {size / 2**20:.0f} MiB")
        assert new_mib < 0.5 * size / 2**20 + 100, "streaming export is not memory-bounded"
    print("PASS streaming group export: byte-identical outputs, bounded memory")


if __name__ == "__main__":
    main()

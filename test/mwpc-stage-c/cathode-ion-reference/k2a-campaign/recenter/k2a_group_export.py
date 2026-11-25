#!/usr/bin/env python3
"""Split / relabel an archived ion-path export into strip-recentered sub-exports (same geometry; READ-ONLY on source).

  k2a_group_export.py SRC_EXPORT OUT_DIR --mode primary [--pitch-mm 1.7]
  k2a_group_export.py SRC_EXPORT OUT_DIR --mode shift --m-U M --m-V M

primary: every Heed primary p gets m_U(p) = floor(Dbar_U(p)/P + 1/2), m_V(p) likewise, with Dbar the ion-weight mean of
         the first-vertex strip columns (= the charge-weighted cloud D of the moment pipeline, combined per primary).
         Primaries with the same (m_U, m_V) form one sub-export OUT_DIR/g_<mU>_<mV>/ in which upper/lower_index are
         base + m and vertex columns 3/4 are D - m P. Groups are whole primaries, so each sub-export's weight sum is an
         integer exact census (the frozen replay requires it). OUT_DIR/groups.json lists primary -> (m_U, m_V).
shift:   the WHOLE export relabelled by one (m_U, m_V): physically the same event, window moved by m strips (symmetry
         test of the frozen replay / LUT strip mapping).
Transport columns (t, u, v), path lengths, weights and birth/end times are byte-identical; their sha256 is recorded.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def seal(d):
    return dict(l.split("=", 1) for l in (d / "complete.txt").read_text().splitlines() if "=" in l)


def write_export(dst, rows, V, idx_rows, mU, mV, P, src_seal):
    dst.mkdir(parents=True)
    out_rows, off, blocks = [], 0, []
    for new_j, (j, mu, mv) in enumerate(zip(idx_rows, mU, mV)):
        r = rows[j]
        s, n = int(r["first_vertex"]), int(r["vertices"])
        b = np.array(V[s:s + n], dtype="<f8")
        if mu:
            b[:, 3] -= P * mu
        if mv:
            b[:, 4] -= P * mv
        blocks.append(b)
        rr = dict(r); rr["row"] = str(new_j); rr["first_vertex"] = str(off)
        rr["upper_index"] = str(int(r["upper_index"]) + mu); rr["lower_index"] = str(int(r["lower_index"]) + mv)
        out_rows.append(rr); off += n
    allv = np.vstack(blocks)
    allv.astype("<f8").tofile(dst / "vertices.f64")
    with open(dst / "ions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(out_rows)
    wsum = sum(float(r["weight"]) for r in out_rows)
    exact = int(round(wsum))
    if not np.isclose(wsum, exact, rtol=1e-9, atol=1e-6):        # the frozen replay's own census test
        raise SystemExit(f"{dst}: group weight sum {wsum} is not an integer census (groups must be whole primaries)")
    lines = [f"rows={len(out_rows)}", f"exact_ions={exact}", f"weight_sum={wsum!r}", f"vertices={len(allv)}",
             f"step_mm={src_seal['step_mm']}"]
    extra = {k: v for k, v in src_seal.items() if k not in ("rows", "exact_ions", "weight_sum", "vertices", "step_mm")}
    lines += [f"{k}={v}" for k, v in extra.items()]
    (dst / "complete.txt").write_text("\n".join(lines) + "\n")
    return hashlib.sha256(np.ascontiguousarray(allv[:, :3]).tobytes()).hexdigest(), len(out_rows), exact


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", type=Path); ap.add_argument("out", type=Path)
    ap.add_argument("--mode", choices=("primary", "shift"), required=True)
    ap.add_argument("--m-U", type=int, default=0); ap.add_argument("--m-V", type=int, default=0)
    ap.add_argument("--pitch-mm", type=float, default=1.7)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    S = seal(a.src)
    rows = list(csv.DictReader(open(a.src / "ions.csv", newline="")))
    V = np.memmap(a.src / "vertices.f64", dtype="<f8", mode="r", shape=(int(S["vertices"]), 5))
    first = np.array([int(r["first_vertex"]) for r in rows])
    p = np.array([int(r["primary"]) for r in rows]); w = np.array([float(r["weight"]) for r in rows])
    P = a.pitch_mm
    src_tuv = hashlib.sha256(np.ascontiguousarray(np.asarray(V)[:, :3]).tobytes()).hexdigest()
    a.out.mkdir(parents=True)
    meta = {"tool": "k2a_group_export.py", "source": str(a.src), "mode": a.mode, "pitch_mm": P,
            "source_tuv_sha256": src_tuv, "groups": {}}
    if a.mode == "shift":
        idx = list(range(len(rows)))
        tuv, n, ex = write_export(a.out / f"g_{a.m_U}_{a.m_V}", rows, V, idx, [a.m_U] * len(rows), [a.m_V] * len(rows), P, S)
        if tuv != src_tuv:
            raise SystemExit("transport columns changed")
        meta["groups"][f"g_{a.m_U}_{a.m_V}"] = {"rows": n, "exact_ions": ex, "m_U": a.m_U, "m_V": a.m_V}
    else:
        D = np.asarray(V[first, 3:5])
        mprim = {}
        for pp in np.unique(p):
            s = p == pp
            dU, dV = np.sum(w[s] * D[s, 0]) / w[s].sum(), np.sum(w[s] * D[s, 1]) / w[s].sum()
            mprim[int(pp)] = (int(np.floor(dU / P + 0.5)), int(np.floor(dV / P + 0.5)))
        grp = defaultdict(list)
        for j, pp in enumerate(p):
            grp[mprim[int(pp)]].append(j)
        for (mu, mv), idx in sorted(grp.items()):
            name = f"g_{mu}_{mv}"
            tuv, n, ex = write_export(a.out / name, rows, V, idx, [mu] * len(idx), [mv] * len(idx), P, S)
            meta["groups"][name] = {"rows": n, "exact_ions": ex, "m_U": mu, "m_V": mv,
                                    "primaries": sorted({int(p[j]) for j in idx})}
        meta["primary_shift"] = {str(k): list(v) for k, v in sorted(mprim.items())}
        if sum(g["exact_ions"] for g in meta["groups"].values()) != int(S["exact_ions"]):
            raise SystemExit("group censuses do not add up to the source census")
    (a.out / "groups.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"GROUPED {a.src} -> {len(meta['groups'])} sub-exports: "
          + ", ".join(f"{k} ({v['rows']} rows)" for k, v in meta["groups"].items()))


if __name__ == "__main__":
    main()

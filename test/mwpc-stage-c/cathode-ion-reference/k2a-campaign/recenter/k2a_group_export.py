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
Memory: vertex blocks are streamed to disk one ion path at a time (a 1.33 GiB T1 export is never held in memory); the
source is read with explicit bounded reads (no mmap). Disk: shift mode writes one full copy, primary mode one copy in total.
"""
import argparse
import csv
import hashlib
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np


def seal(d):
    return dict(l.split("=", 1) for l in (d / "complete.txt").read_text().splitlines() if "=" in l)


def write_export(dst, rows, V, idx_rows, mU, mV, P, src_seal):
    """Stream the selected rows' vertex blocks to dst/vertices.f64, one ion path at a time (memory = one path). The
    arithmetic per block is the same as the original in-memory implementation (copy as <f8, column 3 -= P m_U,
    column 4 -= P m_V), so the bytes are identical; the (t, u, v) hash is accumulated block by block over the same
    row-major bytes, so it equals the hash of the whole array."""
    dst.mkdir(parents=True)
    out_rows, off = [], 0
    h_tuv = hashlib.sha256()
    with open(dst / "vertices.f64", "wb", buffering=1 << 22) as fv:
        for new_j, (j, mu, mv) in enumerate(zip(idx_rows, mU, mV)):
            r = rows[j]
            s, n = int(r["first_vertex"]), int(r["vertices"])
            b = np.array(V[s:s + n], dtype="<f8")
            if mu:
                b[:, 3] -= P * mu
            if mv:
                b[:, 4] -= P * mv
            fv.write(b.tobytes())
            h_tuv.update(np.ascontiguousarray(b[:, :3]).tobytes())
            rr = dict(r); rr["row"] = str(new_j); rr["first_vertex"] = str(off)
            rr["upper_index"] = str(int(r["upper_index"]) + mu); rr["lower_index"] = str(int(r["lower_index"]) + mv)
            out_rows.append(rr); off += n
    if (dst / "vertices.f64").stat().st_size != 40 * off:
        raise SystemExit(f"{dst}: vertices.f64 size mismatch")
    with open(dst / "ions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(out_rows)
    wsum = sum(float(r["weight"]) for r in out_rows)
    exact = int(round(wsum))
    if not np.isclose(wsum, exact, rtol=1e-9, atol=1e-6):        # the frozen replay's own census test
        raise SystemExit(f"{dst}: group weight sum {wsum} is not an integer census (groups must be whole primaries)")
    lines = [f"rows={len(out_rows)}", f"exact_ions={exact}", f"weight_sum={wsum!r}", f"vertices={off}",
             f"step_mm={src_seal['step_mm']}"]
    extra = {k: v for k, v in src_seal.items() if k not in ("rows", "exact_ions", "weight_sum", "vertices", "step_mm")}
    lines += [f"{k}={v}" for k, v in extra.items()]
    (dst / "complete.txt").write_text("\n".join(lines) + "\n")
    return h_tuv.hexdigest(), len(out_rows), exact


class VertexFile:
    """Read-only access to an export's vertices.f64 (n x 5 little-endian float64) by explicit reads, never mmap,
    so the process holds only the block it asked for (mmap'd pages would count toward RSS)."""

    def __init__(self, path, n):
        self.f = open(path, "rb"); self.n = n
        if os.fstat(self.f.fileno()).st_size != 40 * n:
            raise SystemExit(f"{path}: size does not match the seal")

    def __len__(self):
        return self.n

    def __getitem__(self, sl):
        if isinstance(sl, tuple):                     # V[a:b, cols]
            rs, cs = sl
            return self[rs][:, cs]
        if not isinstance(sl, slice) or sl.step not in (None, 1):
            raise TypeError("VertexFile supports contiguous row slices only")
        a, b = max(0, sl.start or 0), min(self.n, self.n if sl.stop is None else sl.stop)
        if b <= a:
            return np.empty((0, 5), "<f8")
        buf = os.pread(self.f.fileno(), 40 * (b - a), 40 * a)
        return np.frombuffer(buf, dtype="<f8").reshape(b - a, 5)

    def rows(self, idx):
        """vertex rows idx (any order), one 40-byte read each"""
        return np.array([np.frombuffer(os.pread(self.f.fileno(), 40, 40 * int(i)), dtype="<f8") for i in idx]).reshape(-1, 5)


def tuv_sha_streaming(V, chunk=1 << 20):
    h = hashlib.sha256()
    for s in range(0, len(V), chunk):
        h.update(np.ascontiguousarray(V[s:s + chunk, :3]).tobytes())
    return h.hexdigest()


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
    V = VertexFile(a.src / "vertices.f64", int(S["vertices"]))
    first = np.array([int(r["first_vertex"]) for r in rows])
    p = np.array([int(r["primary"]) for r in rows]); w = np.array([float(r["weight"]) for r in rows])
    P = a.pitch_mm
    src_tuv = tuv_sha_streaming(V)                    # bounded memory, same bytes as the whole-array hash
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
        D = V.rows(first)[:, 3:5]
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

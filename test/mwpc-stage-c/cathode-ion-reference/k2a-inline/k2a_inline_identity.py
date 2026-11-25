#!/usr/bin/env python3
"""T1 identity checks between runs of the inline generator (exit 0 PASS, 1 FAIL / MISMATCH, 3 INCONCLUSIVE).

  k2a_inline_identity.py moments A_PREFIX B_PREFIX OUT.json
      PREFIX.moments.csv, .clouds.csv and .primaries.csv byte-identical; .inline.json identical except the
      run-dependent keys (prefix paths, storage, process peak RSS). Used for: rerun reproducibility (R1 vs R2), the
      production moments-only path vs the full-archive path (R3 vs R1), weighted-archive run vs full (R4 vs R1).
  k2a_inline_identity.py dumps A_DUMP B_DUMP OUT.json
      births.csv and primaries.csv of two archive dumps byte-identical (rerun reproducibility of the archive content).
  k2a_inline_identity.py subset FULL_DUMP WEIGHTED_DUMP OUT.json [--primaries N] [--label TEXT]
      the weighted (reservoir) archive is a sample of the SAME avalanches: for its primaries < N (default all of the
      full dump) garfield_seed and ion_births_exact equal the full archive's, every weighted row (primary,
      callback_ordinal) exists in the full archive with bit-identical u, v, w, t, ion_births_saved = min(N, K) and
      each weight = exact / saved. With --label frozen-production the same test against the frozen production K256
      archive of this muon is reported as MATCH / MISMATCH (environment-dependent provenance diagnostic).
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

RUN_DEPENDENT = {"moments_file", "storage", "process_peak_rss_kib_at_finish"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def finish(out, rep, ok, label="identity"):
    rep["status"] = rep.get("status") or ("PASS" if ok else "FAIL")
    Path(out).write_text(json.dumps(rep, indent=1) + "\n")
    print(f"{rep['status']} {label}: " + "; ".join(rep.get("notes", [])))
    return 0 if rep["status"] in ("PASS", "MATCH") else (3 if rep["status"] == "INCONCLUSIVE" else 1)


def moments(a, b, out):
    rep = {"mode": "moments", "a": a, "b": b, "files": {}, "notes": []}
    ok = True
    for s in (".moments.csv", ".clouds.csv", ".primaries.csv", ".inline.json"):
        if not (Path(a + s).is_file() and Path(b + s).is_file()):
            rep["status"] = "INCONCLUSIVE"; rep["notes"].append(f"missing {s}")
            return finish(out, rep, False, "moments identity")
    for s in (".moments.csv", ".clouds.csv", ".primaries.csv"):
        ha, hb = sha(a + s), sha(b + s)
        rep["files"][s] = {"a": ha, "b": hb, "identical": ha == hb}
        ok &= ha == hb
        rep["notes"].append(f"{s} {'identical' if ha == hb else 'DIFFERENT'}")
    ja, jb = json.loads(Path(a + ".inline.json").read_text()), json.loads(Path(b + ".inline.json").read_text())
    diff = sorted(k for k in set(ja) | set(jb) if k not in RUN_DEPENDENT and ja.get(k) != jb.get(k))
    rep["inline_json_differences"] = {k: [ja.get(k), jb.get(k)] for k in diff}
    rep["rss_kib"] = [ja.get("process_peak_rss_kib_at_finish"), jb.get("process_peak_rss_kib_at_finish")]
    ok &= not diff
    rep["notes"].append("inline.json identical outside run-dependent keys" if not diff else f"inline.json differs: {diff}")
    return finish(out, rep, ok, "moments identity")


def dumps(a, b, out):
    rep = {"mode": "dumps", "a": a, "b": b, "files": {}, "notes": []}
    ok = True
    for s in ("births.csv", "primaries.csv"):
        pa, pb = Path(a) / s, Path(b) / s
        if not (pa.is_file() and pb.is_file()):
            rep["status"] = "INCONCLUSIVE"; rep["notes"].append(f"missing {s}")
            return finish(out, rep, False, "archive identity")
        ha, hb = sha(pa), sha(pb)
        rep["files"][s] = {"a": ha, "b": hb, "identical": ha == hb}
        ok &= ha == hb
        rep["notes"].append(f"{s} {'identical' if ha == hb else 'DIFFERENT'}")
    return finish(out, rep, ok, "archive identity")


def subset(full, weighted, out, nprim=None, label=None):
    rep = {"mode": "subset", "full": full, "weighted": weighted, "label": label, "notes": [], "problems": []}
    for p in (Path(full) / "births.csv", Path(full) / "primaries.csv", Path(weighted) / "births.csv", Path(weighted) / "primaries.csv",
              Path(weighted) / "muon.json"):
        if not p.is_file():
            rep["status"] = "INCONCLUSIVE"; rep["notes"].append(f"missing {p}")
            return finish(out, rep, False, "subset")
    fp, wp = rows(Path(full) / "primaries.csv"), rows(Path(weighted) / "primaries.csv")
    K = json.loads((Path(weighted) / "muon.json").read_text())["reservoir"]
    N = len(fp) if nprim is None else nprim
    fmap = {int(r["primary_index"]): r for r in fp}
    wmap = {int(r["primary_index"]): r for r in wp}
    P = sorted(i for i in fmap if i < N)
    prob = rep["problems"]
    for i in P:
        f, w = fmap[i], wmap.get(i)
        if w is None:
            prob.append(f"primary {i} absent from the weighted archive"); continue
        if f["garfield_seed"] != w["garfield_seed"]:
            prob.append(f"primary {i}: garfield_seed {f['garfield_seed']} vs {w['garfield_seed']}")
        ne = int(f["ion_births_exact"])
        if ne != int(w["ion_births_exact"]):
            prob.append(f"primary {i}: census full {ne} vs weighted {w['ion_births_exact']}")
        if int(w["ion_births_saved"]) != min(ne, K):
            prob.append(f"primary {i}: saved {w['ion_births_saved']} != min({ne}, {K})")
    full_rows = {}
    with open(Path(full) / "births.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if int(r["primary_index"]) < N:
                full_rows[(int(r["primary_index"]), int(r["callback_ordinal"]))] = (r["u_patch_cm"], r["v_absolute_cm"], r["w_relative_cm"], r["t_birth_ns"])
    nw = 0
    with open(Path(weighted) / "births.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            i = int(r["primary_index"])
            if i >= N:
                continue
            nw += 1
            key = (i, int(r["callback_ordinal"]))
            vals = (r["u_patch_cm"], r["v_absolute_cm"], r["w_relative_cm"], r["t_birth_ns"])
            if full_rows.get(key) != vals:   # text of %.17g doubles: equal text <=> equal bits
                if len(prob) < 50:
                    prob.append(f"weighted row {key} {'absent from' if key not in full_rows else 'differs from'} the full archive")
            exp = int(wmap[i]["ion_births_exact"]) / int(wmap[i]["ion_births_saved"])
            if abs(float(r["inverse_inclusion_weight"]) - exp) > 1e-12 * exp:
                prob.append(f"weighted row {key}: weight {r['inverse_inclusion_weight']} != {exp}")
    rep.update({"primaries_compared": len(P), "weighted_rows_compared": nw, "reservoir": K})
    rep["notes"].append(f"{len(P)} primaries, {nw} reservoir rows compared, {len(prob)} problems")
    ok = not prob
    if label == "frozen-production":
        rep["status"] = "MATCH" if ok else "MISMATCH"
        rep["notes"].append("environment-dependent provenance diagnostic (not part of the T1 verdict)")
    return finish(out, rep, ok, f"subset ({label or 'same environment'})")


def main():
    a = sys.argv[1:]
    if len(a) >= 4 and a[0] == "moments":
        return moments(a[1], a[2], a[3])
    if len(a) >= 4 and a[0] == "dumps":
        return dumps(a[1], a[2], a[3])
    if len(a) >= 4 and a[0] == "subset":
        n = int(a[a.index("--primaries") + 1]) if "--primaries" in a else None
        lab = a[a.index("--label") + 1] if "--label" in a else None
        return subset(a[1], a[2], a[3], n, lab)
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())

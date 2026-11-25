#!/usr/bin/env python3
"""T1: in-memory full-birth K2A moments vs two independent full-population references of the SAME avalanche.

  k2a_inline_compare.py --inline PREFIX --dump DUMPDIR --batch BATCH.moments.csv --input input.csv
                        --max-heed-primaries N --out OUTDIR

PREFIX   inline generator outputs (PREFIX.inline.json, .moments.csv, .clouds.csv, .primaries.csv) of a run that ALSO
         wrote the full-mode archive (every callback birth, weight 1) from the same birth stream
DUMPDIR  k2a_archive_dump of that full archive (raw TTree read; no shared code with the batch path)
BATCH    reference (a): spatial-ablation p3_birth_moments on the full archive (frozen reader + BirthMoments.hh)
reference (b), computed here: numpy per-ion values + math.fsum (correctly rounded sums) per cloud, two-pass variances.

Checks (each PASS/FAIL; any FAIL -> FAIL, exit 1; missing/unreadable input -> INCONCLUSIVE, exit 3; else PASS, exit 0):
  CENSUS    per primary: archive rows == ion_births_exact == ion_births_saved == inline callback census == births
            received by the accumulator; callback ordinals 0..N-1 in order; all weights exactly 1; zero-ion primaries
            agree; total births and per-cloud births exact and integer; sum of cloud births == primary births
  SEEDS     primary indices 0..P-1 with P = min(N, generated Heed primaries); garfield_seed == 1 + (heed_seed +
            (i+1) 104729) mod 2e9 on every primary and ion row (original indices and seeds retained)
  IDENTITY  the (primary, wire) sets and order of inline, batch and reference are identical; wire = rint(u_mm/3)
  FRAME     patch_shift_cm, u_in_mm, w_in_mm bitwise equal to the archive's; base_U/base_V equal everywhere
  STATUS    inline event usable, status OK (ZERO_ION_PRIMARY informational), no orphan births, moments file present
  MOMENTS   every cloud field within its error bound of BOTH references (below); bitwise agreement with the batch
            path is counted separately per field (expected for all spatial moments when the two binaries are built
            without floating-point contraction)
  CHARGE    sum of cloud N_fC == total births x e within (n+1) u
Error bounds, u = 2^-53, n = births of the cloud (n_p of the primary), M1 = mean |x| over the cloud's births:
  mean of a per-ion quantity x (D_U, D_V, C1, S1, C2, S2, A, log r):  (4 n + 16) u M1
      recursive summation and the division: 2 gamma_n M1 (Higham 2002, eq. 4.4, |mean| <= M1); per-ion libm value
      (glibc vs numpy, assumed <= 4 ulp): 8 u M1; reference rounding u; factor 2 margin on the summation term.
      A (signs): plus 2 n_amb / n with n_amb = births whose |sin theta| < 1e-14 (sign not reproducible across libm).
  r_eff: relative, (4 n + 16) u mean|log r| + 4 u;  min_clearance: 4 u r_min (hypot, <= 2 ulp)
  N_fC: (n + 1) u N_fC (n equal positive terms)
  T1 time t_ns (= primary weight-mean birth time): vs reference (4 n_p + 16) u mean_p|t|; vs batch (which forms
      t + (mean - t) per ion and averages again): 4 (n + 2) u T, T = max |t| of the cloud and of the T1 time
  t_rms_ns: inline exactly 0 (collapsed, by definition); batch rounding noise <= sqrt(8 (n + 2) u) T
  raw timing (t_raw_mean_ns, t_raw_rms_ns per cloud and per primary): mean (4 n + 16) u mean|t|; RMS
      min(sqrt(dV), dV / rms_ref), dV = 2 (4 n + 16) u <t^2>
Diagnostic (reported, never a failure): births whose generator census wire llround(u/0.3) differs from the K2A wire
rint(u_mm/3) - possible only on a half-cell boundary; each such birth's distance to the boundary is listed.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

U = 2.0 ** -53
EFC = 1.602176634e-4
SPATIAL = ["D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A"]


def llround(x):
    a = abs(x); q = math.floor(a); r = q + 1 if a - q >= 0.5 else q
    return int(math.copysign(r, x)) if r else 0


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def bits_equal(a, b):
    return np.float64(a).tobytes() == np.float64(b).tobytes()


class Report:
    def __init__(self):
        self.groups = {}

    def check(self, group, ok, what):
        g = self.groups.setdefault(group, {"pass": 0, "fail": 0, "failures": []})
        if ok:
            g["pass"] += 1
        else:
            g["fail"] += 1
            if len(g["failures"]) < 50:
                g["failures"].append(what)
        return ok

    def status(self):
        return "FAIL" if any(g["fail"] for g in self.groups.values()) else "PASS"


def reference(dump, frame, base_U, base_V):
    """Independent reference (b) from the dumped births."""
    d = np.genfromtxt(dump / "births.csv", delimiter=",", names=True, dtype=None, encoding=None)
    if d.size == 0:
        return d, {}, {}
    d = np.atleast_1d(d)
    p = d["primary_index"].astype(np.int64)
    u = 10 * d["u_patch_cm"]; v = 10 * d["v_absolute_cm"]
    k = np.rint(u / 3).astype(np.int64)
    du = u - 3 * k
    r = 1000 * np.hypot(du, v)
    th = np.arctan2(v, du)
    wIn = .1 * frame["w_in_mm"]
    wu = 10. * (d["w_relative_cm"] + wIn - .1 * frame["patch_shift_cm"] - .17 * base_U)
    wv = 10. * (d["w_relative_cm"] + wIn + .1 * frame["patch_shift_cm"] - .17 * base_V)
    w = d["inverse_inclusion_weight"].astype(float); t = d["t_birth_ns"].astype(float)
    s = np.sin(th)
    per = {"D_U_mm": wu, "D_V_mm": wv, "C1": np.cos(th), "S1": s, "C2": np.cos(2 * th), "S2": np.sin(2 * th),
           "A": np.sign(s), "logr": np.log(r)}
    fs = math.fsum
    prim = {}
    for pi in np.unique(p):
        m = p == pi
        W = fs(w[m]); tm = fs(w[m] * t[m]) / W
        prim[int(pi)] = {"n": int(m.sum()), "t_mean": tm, "mean_abs_t": fs(w[m] * np.abs(t[m])) / W,
                         "t_raw_rms": math.sqrt(fs(w[m] * (t[m] - tm) ** 2) / W), "mean_t2": fs(w[m] * t[m] ** 2) / W,
                         "max_abs_t": float(np.max(np.abs(t[m])))}
    clouds = {}
    order = np.lexsort((k, p))
    keys = np.stack([p[order], k[order]], axis=1)
    bounds = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
    for g in np.split(order, bounds):
        pi, ki = int(p[g[0]]), int(k[g[0]])
        W = fs(w[g]); q = w[g] * EFC
        c = {"n": len(g), "N_fC": fs(q), "weight_sum": W}
        for f, x in per.items():
            c[f] = fs(w[g] * x[g]) / W
            c["M1_" + f] = fs(w[g] * np.abs(x[g])) / W
        c["r_eff_um"] = math.exp(c["logr"])
        c["r_min"] = float(np.min(r[g])); c["min_clearance_um"] = c["r_min"] - 15.0
        tm = fs(w[g] * t[g]) / W
        c["t_raw_mean_ns"] = tm; c["t_raw_rms_ns"] = math.sqrt(fs(w[g] * (t[g] - tm) ** 2) / W)
        c["mean_abs_t"] = fs(w[g] * np.abs(t[g])) / W; c["mean_t2"] = fs(w[g] * t[g] ** 2) / W
        c["max_abs_t"] = float(np.max(np.abs(t[g])))
        c["n_amb"] = int(np.sum(np.abs(s[g]) < 1e-14))
        clouds[(pi, ki)] = c
    return d, prim, clouds


def rms_tol(n, mean_t2, rms_ref):
    dV = 2 * (4 * n + 16) * U * mean_t2
    return min(math.sqrt(dV), dV / rms_ref if rms_ref > 0 else math.inf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inline", dest="inl", required=True)
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--batch", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--max-heed-primaries", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    a.out.mkdir(parents=True)
    R = Report()
    need = [Path(a.inl + s) for s in (".inline.json", ".clouds.csv", ".primaries.csv")] + \
           [a.dump / "births.csv", a.dump / "primaries.csv", a.dump / "muon.json", a.batch, a.input]
    missing = [str(p) for p in need if not p.is_file()]
    if missing:
        (a.out / "compare.json").write_text(json.dumps({"status": "INCONCLUSIVE", "missing": missing}, indent=1) + "\n")
        print("INCONCLUSIVE T1 compare: missing", missing)
        return 3
    J = json.loads(Path(a.inl + ".inline.json").read_text())
    M = json.loads((a.dump / "muon.json").read_text())
    inp = read_csv(a.input)
    heed = int(inp[0]["heed_seed"])
    # ---------- STATUS ----------
    R.check("STATUS", J["usable_for_cathode_signal"] is True, f"inline event not usable: {J['status_names']}")
    R.check("STATUS", (J["status"] & ~32) == 0, f"inline status {J['status_names']}")
    R.check("STATUS", J["orphan_births"] == 0, "orphan births")
    R.check("STATUS", Path(a.inl + ".moments.csv").is_file(), "engine-facing moments file missing")
    R.check("STATUS", M["storage_mode"] == "full" and M["reservoir"] == 0, f"archive is not full mode: {M['storage_mode']}")
    if not Path(a.inl + ".moments.csv").is_file():
        (a.out / "compare.json").write_text(json.dumps({"status": "FAIL", "groups": R.groups}, indent=1) + "\n")
        print("FAIL T1 compare: no moments (blocked event)")
        return 1
    # ---------- FRAME ----------
    frame = {k: M[k] for k in ("patch_shift_cm", "u_in_mm", "w_in_mm")}
    for k in frame:
        R.check("FRAME", bits_equal(J[k], M[k]), f"{k}: inline {J[k]!r} archive {M[k]!r}")
    base_U = llround((.1 * frame["w_in_mm"] - .1 * (.1 * frame["u_in_mm"])) / .17)
    base_V = llround((.1 * frame["w_in_mm"] + .1 * (.1 * frame["u_in_mm"])) / .17)
    batch = read_csv(a.batch)
    R.check("FRAME", J["base_U"] == base_U and J["base_V"] == base_V, f"inline base {J['base_U']},{J['base_V']} vs ref {base_U},{base_V}")
    if batch:
        R.check("FRAME", int(batch[0]["base_U"]) == base_U and int(batch[0]["base_V"]) == base_V, "batch base_U/V")
    # ---------- references ----------
    d, rprim, rcl = reference(a.dump, frame, base_U, base_V)
    dp = read_csv(a.dump / "primaries.csv")
    ip = read_csv(a.inl + ".primaries.csv")
    ic = read_csv(a.inl + ".clouds.csv")
    im = read_csv(a.inl + ".moments.csv")
    # ---------- SEEDS ----------
    P = min(a.max_heed_primaries, J["n_heed_primaries_generated"]) if a.max_heed_primaries > 0 else J["n_heed_primaries_generated"]
    R.check("SEEDS", [int(r["primary_index"]) for r in dp] == list(range(P)), f"archive primary indices != 0..{P - 1}")
    R.check("SEEDS", [int(r["primary"]) for r in ip] == list(range(P)), "inline primary indices != 0..P-1")
    R.check("SEEDS", int(J["heed_seed"]) == heed and (M["heed_seed"] & 0xffffffff) == heed, "heed seed (input / inline / archive)")
    seed_of = {}
    for r in dp:
        i = int(r["primary_index"]); exp = 1 + (heed + (i + 1) * 104729) % 2_000_000_000
        seed_of[i] = int(r["garfield_seed"])
        R.check("SEEDS", int(r["garfield_seed"]) == exp, f"primary {i}: garfield_seed {r['garfield_seed']} != {exp}")
    gs = d["garfield_seed"].astype(np.int64) if d.size else np.array([], dtype=np.int64)
    pp = d["primary_index"].astype(np.int64) if d.size else np.array([], dtype=np.int64)
    R.check("SEEDS", all(int(x) in seed_of and bool(np.all(gs[pp == x] == seed_of[int(x)])) for x in np.unique(pp)),
            "ion rows carry their primary's garfield_seed")
    # ---------- CENSUS ----------
    ipm = {int(r["primary"]): r for r in ip}
    census_rows = []
    tot = 0
    for r in dp:
        i = int(r["primary_index"]); N = int(r["ion_births_exact"]); n_rows = int(np.sum(pp == i))
        ordv = d["callback_ordinal"][pp == i] if d.size else []
        q = ipm.get(i, {})
        cb = int(q.get("callback_births", -1)); seen = int(q.get("births_seen", -1))
        ok = N == int(r["ion_births_saved"]) == n_rows == cb == seen and int(r["successful"]) == 1
        R.check("CENSUS", ok, f"primary {i}: exact {N} saved {r['ion_births_saved']} rows {n_rows} inline callback {cb} seen {seen}")
        R.check("CENSUS", np.array_equal(np.asarray(ordv, dtype=np.int64), np.arange(n_rows)), f"primary {i}: callback ordinals not 0..N-1 in order")
        R.check("CENSUS", (N == 0) == ("ZERO_ION_PRIMARY" in q.get("status_names", "")), f"primary {i}: zero-ion flag")
        tot += N
        rp = rprim.get(i)
        census_rows.append({"primary": i, "garfield_seed": r["garfield_seed"], "archive_births": N, "inline_callback_births": cb,
                            "inline_births_seen": seen, "inline_clouds": q.get("n_clouds"), "inline_status": q.get("status_names"),
                            "T1_inline_ns": q.get("t_mean_ns"), "T1_reference_ns": repr(rp["t_mean"]) if rp else "",
                            "raw_rms_inline_ns": q.get("t_raw_rms_ns"), "raw_rms_reference_ns": repr(rp["t_raw_rms"]) if rp else ""})
        if rp:
            T1 = float(q["t_mean_ns"])
            R.check("TIMING", abs(T1 - rp["t_mean"]) <= (4 * rp["n"] + 16) * U * rp["mean_abs_t"], f"primary {i}: T1 {T1!r} vs ref {rp['t_mean']!r}")
            tr = float(q["t_raw_rms_ns"])
            R.check("TIMING", abs(tr - rp["t_raw_rms"]) <= rms_tol(rp["n"], rp["mean_t2"], rp["t_raw_rms"]), f"primary {i}: raw rms {tr!r} vs {rp['t_raw_rms']!r}")
    R.check("CENSUS", J["births"] == tot == len(d), f"total births inline {J['births']} archive {tot} rows {len(d)}")
    R.check("CENSUS", d.size == 0 or bool(np.all(d["inverse_inclusion_weight"] == 1.0)), "full-mode weights not all exactly 1")
    R.check("CENSUS", J["n_primaries_zero_ions"] == sum(int(r["ion_births_exact"]) == 0 for r in dp), "zero-ion primary count")
    # ---------- IDENTITY ----------
    kin = [(int(r["primary"]), int(r["wire"])) for r in ic]
    kmo = [(int(r["primary"]), int(r["wire"])) for r in im]
    kba = [(int(r["primary"]), int(r["wire"])) for r in batch]
    kre = sorted(rcl)
    R.check("IDENTITY", kin == kmo == kba == kre, f"cloud keys differ: inline {len(kin)} moments {len(kmo)} batch {len(kba)} ref {len(kre)}")
    per_p = {}
    for r in ic:
        per_p[int(r["primary"])] = per_p.get(int(r["primary"]), 0) + int(r["births"])
    R.check("CENSUS", all(per_p.get(i, 0) == int(r["ion_births_exact"]) for i, r in ((int(x["primary_index"]), x) for x in dp)), "sum of cloud births != primary births")
    # ---------- MOMENTS ----------
    bmap = {(int(r["primary"]), int(r["wire"])): r for r in batch}
    mmap = {(int(r["primary"]), int(r["wire"])): r for r in im}
    fields = ["N_fC", "r_eff_um", "min_clearance_um"] + SPATIAL + ["t_ns", "t_rms_ns", "t_raw_rms_ns"]
    bit = {f: [0, 0] for f in fields}
    worst = {}
    rows = []
    totN = 0.0
    for r in ic:
        key = (int(r["primary"]), int(r["wire"]))
        c = rcl.get(key); b = bmap.get(key); m = mmap.get(key)
        if c is None or b is None or m is None:
            continue
        n = int(r["births"]); totN += float(r["N_fC"])
        R.check("CENSUS", n == c["n"] == int(b["rows"]) == int(m["rows"]), f"{key}: cloud births inline {n} ref {c['n']} batch {b['rows']}")
        rp = rprim[key[0]]
        T = max(c["max_abs_t"], abs(rp["t_mean"]))
        tol = {"N_fC": (n + 1) * U * c["N_fC"],
               "r_eff_um": ((4 * n + 16) * U * c["M1_logr"] + 4 * U) * c["r_eff_um"],
               "min_clearance_um": 4 * U * c["r_min"],
               "t_raw_rms_ns": rms_tol(n, c["mean_t2"], c["t_raw_rms_ns"])}
        for f in SPATIAL:
            tol[f] = (4 * n + 16) * U * c["M1_" + f] + (2 * c["n_amb"] / n if f == "A" else 0)
        ref = {f: c[f] for f in ["N_fC", "r_eff_um", "min_clearance_um", "t_raw_rms_ns"] + SPATIAL}
        row = {"primary": key[0], "wire": key[1], "births": n}
        for f in fields:
            vi = float(m[f]); vb = float(b[f])
            bit[f][1] += 1
            be = bits_equal(vi, vb); bit[f][0] += be
            if f == "t_ns":
                tr_ = (4 * rp["n"] + 16) * U * rp["mean_abs_t"]; tb = 4 * (n + 2) * U * T
                R.check("TIMING", abs(vi - rp["t_mean"]) <= tr_, f"{key}: t_ns {vi!r} vs reference T1 {rp['t_mean']!r}")
                R.check("TIMING", vi == float(ipm[key[0]]["t_mean_ns"]), f"{key}: cloud t_ns != primary T1 time")
                R.check("TIMING", abs(vi - vb) <= tb, f"{key}: t_ns vs batch {vb!r}")
                dr, db, tl = abs(vi - rp["t_mean"]), abs(vi - vb), max(tr_, tb)
                row.update({f + "_inline": repr(vi), f + "_batch": repr(vb), f + "_ref": repr(rp["t_mean"])})
            elif f == "t_rms_ns":
                R.check("TIMING", vi == 0.0, f"{key}: collapsed t_rms_ns not exactly 0")
                tl = math.sqrt(8 * (n + 2) * U) * T
                R.check("TIMING", abs(vb) <= tl, f"{key}: batch t_rms {vb!r} exceeds rounding bound {tl:.3g}")
                dr, db = 0.0, abs(vb)
                row.update({f + "_inline": repr(vi), f + "_batch": repr(vb)})
            else:
                tl = tol[f]
                dr, db = abs(vi - ref[f]), abs(vi - vb)
                grp = "TIMING" if f == "t_raw_rms_ns" else "MOMENTS"
                R.check(grp, dr <= tl, f"{key} {f}: inline {vi!r} vs reference {ref[f]!r} |d| {dr:.3g} > tol {tl:.3g}")
                R.check(grp, db <= tl, f"{key} {f}: inline {vi!r} vs batch {vb!r} |d| {db:.3g} > tol {tl:.3g}")
                row.update({f + "_inline": repr(vi), f + "_batch": repr(vb), f + "_ref": repr(ref[f])})
            row[f + "_tol"] = f"{tl:.3e}"; row[f + "_bitwise_batch"] = int(be)
            worst[f] = max(worst.get(f, 0.0), (max(dr, db) / tl) if tl > 0 else (0.0 if max(dr, db) == 0 else math.inf))
        # raw cloud timing mean (inline diagnostics file only)
        tmi = float(r["t_raw_mean_ns"])
        R.check("TIMING", abs(tmi - c["t_raw_mean_ns"]) <= (4 * n + 16) * U * c["mean_abs_t"], f"{key}: t_raw_mean {tmi!r} vs {c['t_raw_mean_ns']!r}")
        rows.append(row)
    R.check("CHARGE", abs(totN - tot * EFC) <= (tot + 1) * U * tot * EFC * 2, f"sum N_fC {totN!r} vs births x e {tot * EFC!r}")
    # ---------- wire-rule diagnostic ----------
    diag = []
    if d.size:
        u = 10 * d["u_patch_cm"]; k = np.rint(u / 3).astype(np.int64)
        lw = d["local_wire_index"].astype(np.int64)
        for j in np.flatnonzero(lw != k):
            x = u[j] / 3
            diag.append({"row": int(j), "primary": int(pp[j]), "u_mm": repr(float(u[j])), "k2a_wire": int(k[j]), "census_wire": int(lw[j]),
                         "distance_to_half_cell_in_cells": abs(abs(x - math.floor(x)) - 0.5)})
    # ---------- write ----------
    with open(a.out / "clouds_comparison.csv", "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(a.out / "primary_census.csv", "w", newline="") as f:
        if census_rows:
            w = csv.DictWriter(f, fieldnames=list(census_rows[0])); w.writeheader(); w.writerows(census_rows)
    st = R.status()
    rep = {"tool": "k2a_inline_compare.py (T1)", "status": st, "groups": R.groups, "n_clouds": len(rows), "n_primaries": len(dp),
           "births": tot, "bitwise_inline_vs_batch": {f: {"equal": v[0], "compared": v[1]} for f, v in bit.items()},
           "worst_diff_over_tolerance": worst, "wire_rule_differences": diag, "frame": frame, "base_U": base_U, "base_V": base_V,
           "inline_state_bytes_peak": J.get("accumulator_state_bytes_peak"), "heed_seed": heed, "max_heed_primaries": a.max_heed_primaries}
    (a.out / "compare.json").write_text(json.dumps(rep, indent=1) + "\n")
    L = [f"# T1 full-birth moment comparison: **{st}**", "",
         f"{len(dp)} Heed primaries (original indices 0..{len(dp) - 1}), {tot} births, {len(rows)} clouds.", "",
         "| check group | pass | fail |", "|---|---|---|"]
    L += [f"| {g} | {v['pass']} | {v['fail']} |" for g, v in R.groups.items()]
    L += ["", "| field | bitwise inline = batch | worst max(|d_ref|,|d_batch|)/tolerance |", "|---|---|---|"]
    L += [f"| {f} | {bit[f][0]}/{bit[f][1]} | {worst.get(f, 0):.3g} |" for f in fields]
    L += ["", f"Wire-rule differences (generator llround vs K2A rint, diagnostic): {len(diag)}"]
    for g, v in R.groups.items():
        for x in v["failures"][:10]:
            L.append(f"- FAIL {g}: {x}")
    (a.out / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0 if st == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())

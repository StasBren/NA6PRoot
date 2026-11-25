#!/usr/bin/env python3
"""Phase 1b, READ-ONLY: every event's birth-record cloud moments against the support of the frozen K2A kernel.
No waveform is computed; nothing is clamped or extrapolated; every violation is listed with its exact cause.

  k2a_domain_check.py --kernel KERNEL_DIR [--kernel-id ID] --inventory INVENTORY_OUT --selection SELECTION_DIR
                      --moments DIR [DIR ...] [--sealed-moments DIR ...] --out NEW_DIR
                      [--half-strips 5 --window-ns 7000 --dt-ns 10] [--reference-run RUN --tag TAG]

KERNEL_DIR      s_kernel.py build output (kernel.json with green_axes; K.npy only hashed if --kernel-id is given)
INVENTORY_OUT   k2a_inventory.py output (slots.csv: roles, status)
--moments       p3_birth_moments --batch output dirs (<id>.moments.csv + .frame.json + .event.json, or <id>.error.txt)

Domain of the engine (CathodeResponseEngine._prepare / evaluate, reproduced exactly):
  WIRE      cloud wire k = rint(u_mm / 3) must be one of the kernel's wires
  R_EFF     rn[0] <= r_eff <= rn[-1]                      (log-r interpolation, no extrapolation)
  D_U, D_V  floor((D - D0)/dD) - d_sub*H >= 0 and floor((D - D0)/dD) + 1 + d_sub*H <= nD - 1   (all 2H+1 strips)
  T_WINDOW  0 <= t < window: a cloud born at t >= window is dropped SILENTLY by the engine (its electron step is not
            binned and its ion current starts after the window); listed here as out of domain
  FINITE    all moment columns finite, N > 0
Ion-level, relevant only where G/S0 are synthesized on the same grid (Phase 2), reported separately (not a K2A cause):
  G_ION_R_INNER   innermost ion of the cloud (min clearance + 15 um) below rn[0]
Warning (not a K2A cause): PATCH_EDGE_WIRE  cloud on |k| = 4, the outermost wire of the frozen 9-wire field patch.
Event-level: READ_FAILED (converter error, verbatim), NO_IONS (no cloud: zero-ion event), META_MISMATCH (event.json vs
the audited archive_inventory.csv: channel, station, det_id, pdg, Heed primaries, exact ion births).

Sealed validation: moments given with --sealed-moments are checked with the same code, but per-event results go ONLY
to NEW_DIR/SEALED_DO_NOT_OPEN/ and the main report shows counts. Extension requirements and costs are computed from
NON-SEALED events only, so no grid/kernel design decision uses the sealed set.

Outputs: per_event.csv, cloud_violations.csv, coverage.md, coverage.json (by event / ion weight (= charge) / cloud,
overall, per channel, station, role), extension_estimate.json/.md (what a kernel covering all non-sealed events needs
and what it costs, scaled from the measured kernel/Green build). Exit 0 always unless inputs are inconsistent (1).
"""
import argparse
import csv
import hashlib
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

R_WIRE_UM = 15.0
COLS = ("t_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A", "min_clearance_um")


def channel(s):
    s = s.lower()
    for c in ("jpsi", "omega", "phi"):
        if c in s:
            return c
    return "unknown:" + s


def sha256(p, chunk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


class Domain:
    def __init__(self, kd, H, window):
        self.meta = json.loads((kd / "kernel.json").read_text())
        ax = self.meta["green_axes"]
        self.ax = ax
        self.wires = sorted(int(x) for x in ax["wires"])
        self.rn = np.array(ax["r_um"], float)
        self.D = np.array(ax["D_mm"], float)
        self.tau = np.array(ax["tau_ns"], float)
        self.dD, self.dsub = float(ax["dD_mm"]), int(ax["d_sub"])
        self.H, self.window = H, window
        self.nD = len(self.D)
        # exact engine test, expressed as the supported D interval
        self.D_lo = self.D[0] + self.dD * self.dsub * H                     # jd - dsub H >= 0
        self.D_hi = self.D[0] + self.dD * (self.nD - 1 - self.dsub * H)     # jd + 1 + dsub H <= nD - 1  <=>  x < this
        if self.tau[-1] < window:
            raise SystemExit(f"kernel tau range {self.tau[-1]} ns shorter than the window {window} ns")

    def d_ok(self, Dv):
        x = (Dv - self.D[0]) / self.dD
        jd = np.floor(x).astype(np.int64)
        return (jd - self.dsub * self.H >= 0) & (jd + 1 + self.dsub * self.H <= self.nD - 1)


def read_moments(p):
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {c: np.zeros(0) for c in COLS + ("wire", "primary")}, None, None
    X = {c: np.array([float(r[c]) for r in rows]) for c in COLS}
    X["wire"] = np.array([int(r["wire"]) for r in rows]); X["primary"] = np.array([int(r["primary"]) for r in rows])
    return X, int(rows[0]["base_U"]), int(rows[0]["base_V"])


def check_event(eid, mdir, dom, inv_row):
    rec = {"event": eid, "status": "IN_DOMAIN", "causes": "", "clouds": 0, "clouds_out": 0, "charge_fC": 0.0,
           "charge_out_fC": 0.0, "rows": 0, "wires": "", "n_wires": 0, "n_wires_5pct": 0, "dominant_wire_share": 0.0,
           "r_eff_min_um": "", "r_eff_max_um": "", "ion_r_min_um": "", "D_U_min": "", "D_U_max": "", "D_V_min": "",
           "D_V_max": "", "t_min_ns": "", "t_max_ns": "", "g_ion_r_inner_clouds": 0, "meta_mismatch": ""}
    viol = []
    err = mdir / f"{eid}.error.txt"
    mp = mdir / f"{eid}.moments.csv"
    if err.exists():
        rec.update(status="READ_FAILED", causes="READ_FAILED: " + err.read_text().strip().replace("\n", " | "))
        return rec, viol, None
    if not mp.exists():
        rec.update(status="READ_FAILED", causes="READ_FAILED: no moments output (converter crashed or not run)")
        return rec, viol, None
    ev = json.loads((mdir / f"{eid}.moments.csv.event.json").read_text()) if (mdir / f"{eid}.moments.csv.event.json").exists() else {}
    X, bu, bv = read_moments(mp)
    n = len(X["t_ns"])
    rec["clouds"] = n
    rec["rows"] = int(ev.get("rows_read", 0))
    for k in ("source_channel", "station", "det_id", "pdg", "theta_u_deg", "theta_w_deg", "momentum_GeV",
              "n_primaries", "ion_births_exact", "n_primaries_zero_ions"):
        rec[k] = ev.get(k, "")
    rec["channel"] = channel(ev["source_channel"]) if "source_channel" in ev else (inv_row or {}).get("channel", "")
    mism = []
    if inv_row and ev:
        for key, evk, conv in (("channel", None, None), ("station", "station", int), ("det_id", "det_id", int),
                               ("pdg", "pdg", int), ("heed_primaries", "n_primaries", int),
                               ("ion_births_exact", "ion_births_exact", lambda x: int(float(x)))):
            if key == "channel":
                if rec["channel"] != inv_row["channel"]:
                    mism.append(f"channel {rec['channel']}!={inv_row['channel']}")
            elif conv(inv_row[key]) != conv(ev[evk]):
                mism.append(f"{key} {ev[evk]}!={inv_row[key]}")
    elif not ev:
        mism.append("no event.json")
    rec["meta_mismatch"] = ";".join(mism)
    if n == 0:
        rec.update(status="NO_IONS", causes="NO_IONS: no ion birth in the archive (no cloud)")
        return rec, viol, X
    q = X["N_fC"]; rec["charge_fC"] = float(q.sum())
    wires, wq = np.unique(X["wire"], return_inverse=True)
    share = np.bincount(wq, q) / q.sum()
    rec.update(wires=" ".join(map(str, wires)), n_wires=len(wires), n_wires_5pct=int(np.sum(share >= .05)),
               dominant_wire_share=float(share.max()),
               r_eff_min_um=float(X["r_eff_um"].min()), r_eff_max_um=float(X["r_eff_um"].max()),
               ion_r_min_um=float((X["min_clearance_um"] + R_WIRE_UM).min()),
               D_U_min=float(X["D_U_mm"].min()), D_U_max=float(X["D_U_mm"].max()),
               D_V_min=float(X["D_V_mm"].min()), D_V_max=float(X["D_V_mm"].max()),
               t_min_ns=float(X["t_ns"].min()), t_max_ns=float(X["t_ns"].max()))
    bad = np.zeros(n, bool)
    tests = [
        ("NONFINITE", ~np.all(np.isfinite(np.stack([X[c] for c in COLS])), axis=0) | (q <= 0), None, "finite, N>0"),
        ("WIRE", ~np.isin(X["wire"], dom.wires), X["wire"], f"wires {dom.wires}"),
        ("R_EFF_LOW", X["r_eff_um"] < dom.rn[0], X["r_eff_um"], f">= {dom.rn[0]!r} um"),
        ("R_EFF_HIGH", X["r_eff_um"] > dom.rn[-1], X["r_eff_um"], f"<= {dom.rn[-1]!r} um"),
        ("D_U", ~dom.d_ok(X["D_U_mm"]), X["D_U_mm"], f"[{dom.D_lo:.4f}, {dom.D_hi:.4f}) mm with {dom.H} half-strips"),
        ("D_V", ~dom.d_ok(X["D_V_mm"]), X["D_V_mm"], f"[{dom.D_lo:.4f}, {dom.D_hi:.4f}) mm with {dom.H} half-strips"),
        ("T_NEGATIVE", X["t_ns"] < 0, X["t_ns"], ">= 0 ns"),
        ("T_WINDOW", X["t_ns"] >= dom.window, X["t_ns"], f"< {dom.window} ns (engine drops later clouds silently)"),
    ]
    counts = {}
    for name, m, val, lim in tests:
        if m.any():
            counts[name] = (int(m.sum()), float(q[m].sum()))
            bad |= m
            for i in np.flatnonzero(m):
                viol.append({"event": eid, "primary": int(X["primary"][i]), "wire": int(X["wire"][i]), "cause": name,
                             "value": "" if val is None else repr(float(val[i])), "supported": lim,
                             "N_fC": float(q[i])})
    edge = np.abs(X["wire"]) >= 4                      # outermost wires of the 9-wire analytic patch (model edge)
    rec["patch_edge_wire_clouds"] = int(edge.sum()); rec["patch_edge_wire_charge_frac"] = float(q[edge].sum() / q.sum())
    gin = (X["min_clearance_um"] + R_WIRE_UM) < dom.rn[0]
    rec["g_ion_r_inner_clouds"] = int(gin.sum())
    rec["clouds_out"] = int(bad.sum()); rec["charge_out_fC"] = float(q[bad].sum())
    if bad.any():
        rec["status"] = "OUT_OF_DOMAIN"
        rec["causes"] = "; ".join(f"{k}: {c} clouds, {Q / q.sum():.4%} of charge" for k, (c, Q) in counts.items())
    if rec["meta_mismatch"]:
        rec["status"] = "META_MISMATCH" if rec["status"] == "IN_DOMAIN" else rec["status"]
        rec["causes"] = (rec["causes"] + "; " if rec["causes"] else "") + "META_MISMATCH: " + rec["meta_mismatch"]
    return rec, viol, X


def coverage(recs, key=None):
    groups = defaultdict(list)
    for r in recs:
        groups[r.get(key, "") if key else "all"].append(r)
    out = {}
    for g, rs in sorted(groups.items(), key=lambda kv: str(kv[0])):
        ok = [r for r in rs if r["status"] == "IN_DOMAIN"]
        ch = sum(r["charge_fC"] for r in rs); cl = sum(r["clouds"] for r in rs)
        out[str(g)] = {"events": len(rs), "events_in_domain": len(ok),
                       "event_fraction": len(ok) / len(rs) if rs else float("nan"),
                       "clouds": cl, "clouds_in_domain": cl - sum(r["clouds_out"] for r in rs),
                       "cloud_fraction": (cl - sum(r["clouds_out"] for r in rs)) / cl if cl else float("nan"),
                       "charge_fC": ch, "charge_fraction_in_domain": (ch - sum(r["charge_out_fC"] for r in rs)) / ch if ch else float("nan"),
                       "status_counts": {s: sum(1 for r in rs if r["status"] == s) for s in sorted({r["status"] for r in rs})}}
    return out


def extension(dom, recs, X_all, kd, reference_run, tag):
    """What a kernel covering every non-sealed in-file cloud needs, and its cost scaled from the measured build."""
    ws = sorted({int(w) for X in X_all for w in X["wire"]})
    allX = {c: np.concatenate([X[c] for X in X_all]) for c in ("r_eff_um", "D_U_mm", "D_V_mm", "t_ns", "min_clearance_um")} if X_all else {}
    need = {"wires": ws, "missing_wires": [w for w in ws if w not in dom.wires]}
    shape = dom.meta.get("shape")                                       # (basis, side, wire, r, D, tau)
    kbytes = int(dom.meta.get("bytes", 0))
    ax = dom.ax
    nW, nr, nD, nt = len(dom.wires), len(dom.rn), dom.nD, len(dom.tau)
    step = dom.rn[1] / dom.rn[0]
    if allX:
        rlo = float(allX["r_eff_um"].min()); rhi = float(allX["r_eff_um"].max())
        rion = float((allX["min_clearance_um"] + R_WIRE_UM).min())
        dlo = float(min(allX["D_U_mm"].min(), allX["D_V_mm"].min())); dhi = float(max(allX["D_U_mm"].max(), allX["D_V_mm"].max()))
        tmax = float(allX["t_ns"].max())
    else:
        rlo = rhi = rion = dlo = dhi = tmax = float("nan")
    add_lo = max(0, math.ceil(math.log(dom.rn[0] / rlo) / math.log(step))) if rlo < dom.rn[0] else 0
    add_lo_g = max(0, math.ceil(math.log(dom.rn[0] / rion) / math.log(step))) if rion < dom.rn[0] else 0
    add_hi = max(0, math.ceil(math.log(rhi / dom.rn[-1]) / math.log(step))) if rhi > dom.rn[-1] else 0
    nlo = math.ceil((dom.D_lo - dlo) / dom.dD) + 1 if dlo < dom.D_lo else 0      # +1: floor() at the node
    nhi = math.floor((dhi - dom.D_hi) / dom.dD) + 1 if dhi >= dom.D_hi else 0
    nt_new = nt if tmax < dom.window else None
    need.update(r_eff_min_um=rlo, r_eff_max_um=rhi, ion_r_min_um=rion, D_min_mm=dlo, D_max_mm=dhi, t_max_ns=tmax,
                supported={"wires": dom.wires, "r_um": [float(dom.rn[0]), float(dom.rn[-1])],
                           "D_mm_cloud": [float(dom.D_lo), float(dom.D_hi)], "window_ns": dom.window,
                           "tau_max_ns": float(dom.tau[-1])},
                radial_nodes_to_add={"inner_K2A": add_lo, "inner_ion_level_G": add_lo_g, "outer": add_hi},
                D_nodes_to_add={"low": nlo, "high": nhi},
                window_extension_needed=nt_new is None,
                tau_axis_covers_t_max=bool(tmax < dom.tau[-1]))
    W2 = nW + len(need["missing_wires"]); r2 = nr + max(add_lo, add_lo_g) + add_hi; D2 = nD + nlo + nhi
    f_table = (W2 / nW) * (r2 / nr) * (D2 / nD)
    f_paths_new = ((W2 * r2) - (nW * nr)) / (nW * nr)          # NEW grid paths relative to the existing grid paths
    gdir = Path(dom.meta.get("green", ""))
    gbytes = sum(p.stat().st_size for p in gdir.glob("*.npy")) if gdir.is_dir() else None
    timing = {}
    if reference_run and (reference_run / "timing.txt").exists():
        for t in (reference_run / "timing.txt").read_text().split():
            if "=" in t:
                k, v = t.split("=", 1)
                try:
                    timing[k] = float(v)
                except ValueError:
                    pass
    paths_s = timing.get(f"green_paths_{tag}_s") if tag else None
    gpdir = reference_run / f"green_paths_{tag}" if (reference_run and tag) else None
    gpbytes = sum(p.stat().st_size for p in gpdir.rglob("*") if p.is_file()) if gpdir and gpdir.is_dir() else None
    cost = {"factor_table_and_kernel_size": f_table, "factor_new_grid_paths": f_paths_new,
            "current": {"wires": nW, "r_nodes": nr, "D_nodes": nD, "tau_nodes": nt, "kernel_bytes": kbytes,
                        "kernel_shape": shape, "green_table_bytes": gbytes, "green_build_s": ax.get("build_s"),
                        "kernel_build_s": dom.meta.get("build_s"), "grid_paths_s": paths_s, "grid_paths_bytes": gpbytes},
            "extended": {"wires": W2, "r_nodes": r2, "D_nodes": D2, "tau_nodes": nt,
                         "kernel_bytes": kbytes * f_table if kbytes else None,
                         "green_table_bytes": gbytes * f_table if gbytes else None,
                         "green_build_s": ax.get("build_s") * f_table if ax.get("build_s") else None,
                         "kernel_build_s": dom.meta.get("build_s") * f_table if dom.meta.get("build_s") else None,
                         "new_grid_paths_s": paths_s * f_paths_new if paths_s else None,
                         "new_grid_paths_bytes": gpbytes * f_paths_new if gpbytes else None,
                         "engine_resident_bytes": kbytes * f_table if kbytes else None},
            "basis": "EXPECTED (linear scaling of the measured build; not a measurement). Table/kernel size and build "
                     "time scale with wires x r nodes x D nodes; grid paths (RKF) scale with new (wire, r) start points; "
                     "a D or tau extension needs no new paths; a window beyond tau_max needs longer tau axes."}
    return need, cost


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", type=Path, required=True); ap.add_argument("--kernel-id")
    ap.add_argument("--inventory", type=Path, required=True); ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--moments", type=Path, nargs="+", required=True); ap.add_argument("--sealed-moments", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--half-strips", type=int, default=5); ap.add_argument("--window-ns", type=float, default=7000.0)
    ap.add_argument("--reference-run", type=Path); ap.add_argument("--tag")
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    t0 = time.time()
    kid = None
    if a.kernel_id:
        h = hashlib.sha256()
        for f in ("K.npy", "KPHI0.npy", "kernel.json"):
            h.update(sha256(a.kernel / f).encode())
        kid = h.hexdigest()
        if kid != a.kernel_id:
            raise SystemExit(f"kernel id mismatch: {kid} != {a.kernel_id}")
    dom = Domain(a.kernel, a.half_strips, a.window_ns)
    slots = {r["slot"]: r for r in csv.DictReader(open(a.inventory / "slots.csv", newline=""))}
    inv = {f"{int(r['index']):04d}": r for r in csv.DictReader(open(a.selection / "archive_inventory.csv", newline=""))}
    sealed_ids = {f"muon_{s}" for s, r in slots.items() if r["role"] == "sealed_validation"}

    def scan(dirs, expect_sealed):
        recs, viols, Xs = [], [], []
        seen = set()
        for d in dirs:
            ids = sorted({re.sub(r"\.(moments\.csv|error\.txt)$", "", p.name) for p in list(d.rglob("*.moments.csv")) + list(d.rglob("*.error.txt"))})
            for eid in ids:
                if eid in seen:
                    raise SystemExit(f"event {eid} appears twice in the moments dirs")
                seen.add(eid)
                if (eid in sealed_ids) != expect_sealed:
                    raise SystemExit(f"{eid}: sealed/non-sealed separation violated ({'sealed' if eid in sealed_ids else 'non-sealed'} event in the "
                                     f"{'sealed' if expect_sealed else 'non-sealed'} moments)")
                mdir = next(p.parent for p in list(d.rglob(f"{eid}.moments.csv")) + list(d.rglob(f"{eid}.error.txt")))
                slot = eid.replace("muon_", "")
                rec, v, X = check_event(eid, mdir, dom, inv.get(slot))
                rec["role"] = slots.get(slot, {}).get("role", "unknown"); rec["slot_status"] = slots.get(slot, {}).get("status", "unknown")
                if not rec.get("channel"):
                    rec["channel"] = inv.get(slot, {}).get("channel", "")
                if rec.get("station", "") == "":
                    rec["station"] = inv.get(slot, {}).get("station", "")
                recs.append(rec); viols += v
                if X is not None and len(X["t_ns"]):
                    Xs.append(X)
        return recs, viols, Xs, seen

    recs, viols, Xs, seen = scan(a.moments, False)
    expected = {f"muon_{s}" for s, r in slots.items() if r["status"] == "USABLE" and r["role"] != "sealed_validation"}
    missing = sorted(expected - seen)
    for eid in missing:   # usable but never converted: visible, not dropped
        slot = eid[5:]
        recs.append({"event": eid, "status": "NOT_CONVERTED", "causes": "NOT_CONVERTED: no moments/error output",
                     "role": slots[slot]["role"], "channel": inv.get(slot, {}).get("channel", ""),
                     "station": inv.get(slot, {}).get("station", ""), "clouds": 0, "clouds_out": 0, "charge_fC": 0.0,
                     "charge_out_fC": 0.0})
    a.out.mkdir(parents=True)
    fields = []
    for r in recs:
        fields += [k for k in r if k not in fields]
    with open(a.out / "per_event.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n", restval=""); w.writeheader(); w.writerows(recs)
    vf = ["event", "primary", "wire", "cause", "value", "supported", "N_fC"]
    with open(a.out / "cloud_violations.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=vf, lineterminator="\n"); w.writeheader(); w.writerows(viols)
    cov = {"overall": coverage(recs)["all"], "by_channel": coverage(recs, "channel"), "by_station": coverage(recs, "station"),
           "by_role": coverage(recs, "role")}
    # sealed: identical checks, blinded output
    sealed_summary = None
    if a.sealed_moments:
        srecs, sviols, _, sseen = scan(a.sealed_moments, True)
        sd = a.out / "SEALED_DO_NOT_OPEN"; sd.mkdir()
        sf = []
        for r in srecs:
            sf += [k for k in r if k not in sf]
        with open(sd / "per_event.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=sf, lineterminator="\n", restval=""); w.writeheader(); w.writerows(srecs)
        with open(sd / "cloud_violations.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=vf, lineterminator="\n"); w.writeheader(); w.writerows(sviols)
        (sd / "README.txt").write_text("Sealed-validation domain check. Do not open before the model, grid coverage and acceptance "
                                       "criteria are frozen (Phase 3).\n")
        sexp = {f"muon_{s}" for s, r in slots.items() if r["status"] == "USABLE" and r["role"] == "sealed_validation"}
        sealed_summary = {"events": len(srecs), "expected": len(sexp), "not_converted": len(sexp - sseen),
                          "status_counts": {s: sum(1 for r in srecs if r["status"] == s) for s in sorted({r["status"] for r in srecs})}}
    need, cost = extension(dom, recs, Xs, a.kernel, a.reference_run, a.tag)
    summary = {"tool": "k2a_domain_check.py", "kernel": str(a.kernel), "kernel_id_verified": kid,
               "transport_id": dom.ax.get("grid_vertices_sha256"), "weighting_id": dom.ax.get("lut_sha256"),
               "half_strips": a.half_strips, "window_ns": a.window_ns, "non_sealed": cov, "not_converted": missing,
               "sealed_blinded": sealed_summary, "cloud_violations_by_cause": {c: sum(1 for v in viols if v["cause"] == c) for c in sorted({v["cause"] for v in viols})},
               "g_ion_level_inner_radius_events": sum(1 for r in recs if r.get("g_ion_r_inner_clouds", 0)),
               "meta_mismatch_events": [r["event"] for r in recs if r.get("meta_mismatch")],
               "wall_s": time.time() - t0}
    (a.out / "extension_estimate.json").write_text(json.dumps({"requirement_non_sealed": need, "cost": cost}, indent=1) + "\n")
    pct = lambda x: "n/a" if x != x else f"{100 * x:.2f} %"
    L = ["# K2A kernel domain preflight (read-only; no waveform computed)", "",
         f"kernel `{a.kernel}`; kernel_id {'verified ' + kid[:16] if kid else 'not verified (no --kernel-id)'}; "
         f"wires {dom.wires}; r_eff [{dom.rn[0]:.6g}, {dom.rn[-1]:.6g}] um; cloud D [{dom.D_lo:.3f}, {dom.D_hi:.3f}) mm "
         f"(±{a.half_strips} strips); window {a.window_ns:g} ns (tau max {dom.tau[-1]:g} ns)", "",
         "## Coverage, non-sealed events", "",
         "| group | events | in domain (events) | clouds in domain | ion weight (charge) in domain | statuses |", "|---|---|---|---|---|---|"]
    for title, block in (("all", {"all": cov["overall"]}), ("channel", cov["by_channel"]), ("station", cov["by_station"]),
                         ("role", cov["by_role"])):
        for g, c in block.items():
            L.append(f"| {title if title == 'all' else title + ' ' + g} | {c['events']} | {c['events_in_domain']} "
                     f"({pct(c['event_fraction'])}) | {pct(c['cloud_fraction'])} | {pct(c['charge_fraction_in_domain'])} | "
                     + ", ".join(f"{k} {v}" for k, v in c["status_counts"].items()) + " |")
    L += ["", "## Every non-sealed event not in domain", "", "| event | role | channel | station | status | causes |", "|---|---|---|---|---|---|"]
    for r in recs:
        if r["status"] != "IN_DOMAIN":
            L.append(f"| {r['event']} | {r.get('role')} | {r.get('channel')} | {r.get('station')} | {r['status']} | {r.get('causes')} |")
    L += ["", f"Cloud-level violations by cause: {summary['cloud_violations_by_cause'] or 'none'} (full list: cloud_violations.csv)",
          f"Ion-level inner radius below the Green grid (matters only for G/S0 synthesis on this grid): "
          f"{summary['g_ion_level_inner_radius_events']} events",
          f"Clouds on the outermost wires (|k| = 4) of the 9-wire analytic patch (a FIELD-MODEL edge, shared by S0, G and "
          f"K2A; warning, not a kernel cause): {sum(1 for r in recs if r.get('patch_edge_wire_clouds', 0))} events", ""]
    # input-feature representativeness of the development set (no waveform result is involved)
    def q(rs, k, p):
        v = [float(r[k]) for r in rs if r.get(k, "") not in ("", None)]
        return f"{np.quantile(v, p):.3g}" if v else "n/a"
    L += ["## Development set vs population (input features only)", "",
          "| group | events | clouds p50 / p95 / max | charge fC p50 / p95 | wires>=5% p50 / max | dominant-wire share p5 | "
          "max abs D_U p95 | t_max ns p95 |", "|---|---|---|---|---|---|---|---|"]
    ok = [r for r in recs if r.get("clouds", 0)]
    for r in ok:
        r["_absD"] = max(abs(float(r["D_U_min"])), abs(float(r["D_U_max"])))
    groups = [("development (incl. references)", [r for r in ok if r.get("role") in ("development", "dev_reference")])]
    groups += [(f"non-selected {c}", [r for r in ok if r.get("role") in ("pool", "previously_inspected") and r.get("channel") == c])
               for c in ("jpsi", "phi", "omega")]
    for name, rs in groups:
        L.append(f"| {name} | {len(rs)} | {q(rs, 'clouds', .5)} / {q(rs, 'clouds', .95)} / {q(rs, 'clouds', 1)} | "
                 f"{q(rs, 'charge_fC', .5)} / {q(rs, 'charge_fC', .95)} | {q(rs, 'n_wires_5pct', .5)} / {q(rs, 'n_wires_5pct', 1)} | "
                 f"{q(rs, 'dominant_wire_share', .05)} | {q(rs, '_absD', .95)} | {q(rs, 't_max_ns', .95)} |")
    for r in ok:
        r.pop("_absD", None)
    L.append("")
    nonusable = [(k, v) for k, v in sorted(slots.items()) if v["status"] != "USABLE"]
    L += ["## Slots not processed (from the inventory check)", "",
          "| slot | status | sha | listed excluded |", "|---|---|---|---|"]
    L += [f"| {k} | {v['status']} | {v['sha_check']} | {v['in_excluded_list']} |" for k, v in nonusable] + [""]
    summary["slots_not_processed"] = {k: v["status"] for k, v in nonusable}
    if sealed_summary:
        L += ["## Sealed validation (blinded: counts only; details in SEALED_DO_NOT_OPEN/)", "",
              f"{sealed_summary['events']} of {sealed_summary['expected']} converted; statuses {sealed_summary['status_counts']}", ""]
    ext = cost["extended"]; cur = cost["current"]
    mib = lambda b: "n/a" if not b else f"{b / 2**20:.0f} MiB"
    sec = lambda x: "n/a" if x is None else f"{x:.0f} s"
    L += ["## Kernel extension needed to cover every non-sealed cloud (EXPECTED cost, scaled from the measured build)", "",
          f"missing wires {need['missing_wires'] or 'none'}; radial nodes to add {need['radial_nodes_to_add']}; "
          f"D nodes to add {need['D_nodes_to_add']}; latest cloud {need['t_max_ns']:.1f} ns: window extension needed "
          f"{need['window_extension_needed']} (tau axis covers it: {need['tau_axis_covers_t_max']})", "",
          "| quantity | current | extended |", "|---|---|---|",
          f"| wires x r x D x tau nodes | {cur['wires']} x {cur['r_nodes']} x {cur['D_nodes']} x {cur['tau_nodes']} | "
          f"{ext['wires']} x {ext['r_nodes']} x {ext['D_nodes']} x {ext['tau_nodes']} |",
          f"| kernel K.npy (engine resident RAM) | {mib(cur['kernel_bytes'])} | {mib(ext['kernel_bytes'])} |",
          f"| Green table | {mib(cur['green_table_bytes'])} | {mib(ext['green_table_bytes'])} |",
          f"| Green table build | {sec(cur['green_build_s'])} | {sec(ext['green_build_s'])} |",
          f"| kernel build | {sec(cur['kernel_build_s'])} | {sec(ext['kernel_build_s'])} |",
          f"| grid paths (RKF) | {sec(cur['grid_paths_s'])}, {mib(cur['grid_paths_bytes'])} | NEW: {sec(ext['new_grid_paths_s'])}, "
          f"{mib(ext['new_grid_paths_bytes'])} (factor {cost['factor_new_grid_paths']:.3g} of the existing grid) |",
          "", cost["basis"], ""]
    (a.out / "coverage.md").write_text("\n".join(L) + "\n")
    (a.out / "coverage.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

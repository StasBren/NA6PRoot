#!/usr/bin/env python3
"""Deterministic selection of INDEPENDENT non-sealed muons for the K2A-vs-S0 validation (input metadata only).

  k2a_validation_select.py --pf PREFLIGHT_DIR --recenter RECOVERY_DIR --out NEW_DIR
                           [--exclude-run RUN] [--exclude-slots 0000,0858,0560,0997,0133] [--pilot 5] [--full 16]
                           [--include-previously-inspected] [--max-rows 0]

Reads only non-sealed preflight products: PF/inventory/slots.csv (roles), PF/domain/per_event.csv,
RECOVERY_DIR/per_event_recenter.csv (k2a_recenter_recovery.py: domain AFTER per-primary recentering, shifts) and
PF/moments_nonsealed/**/muon_*.moments.csv (for the hemisphere asymmetry and the per-primary shifts). The sealed
moments directory is never opened; the tool refuses if a sealed id appears in any input it reads.

Independence. Candidates are roles 'pool' and 'development' (and 'previously_inspected' only with the flag), slot
status USABLE, a converted event with clouds. Excluded BEFORE any feature is used:
  * sealed_validation (never read), excluded / failed / unusable slots;
  * every muon prepared in --exclude-run RUN (RUN/provenance_prepare_<M>.txt): these built the Green table / kernel;
  * --exclude-slots: muons already used for development decisions (E2a/E2b references, deferred cases).

Features (no model output is used): station, channel, |theta_u|, |theta_w|, log ion yield per ionising primary (gain),
charge-weighted mean |A| (local hemisphere asymmetry), wires with >= 5 % of the charge, dominant-wire share, maximum
|m| of the per-primary strip shifts and the charge share outside the base window (displacement), clouds, domain status
after recentering (causes), patch-edge charge share. Cost proxy: K256 rows read (export + replay scale with it) and the
number of distinct (m_U, m_V) shift groups (one frozen replay each).

Picks (deterministic; each the best-scoring candidate of its target, the cheapest of the five best on ties of intent):
  pilot  typical          nearest to the candidate median (standardized), in domain, max |m| <= 1
         multiwire_asym   >= 2 wires with >= 5 %, then largest mean |A|, in domain
         displaced        max |m| >= 2, in domain after recentering, largest charge outside the base window
         gain_tail        ion yield per primary >= candidate p95, in domain
         out_of_domain    out of domain after recentering (explicit-reporting test), cheapest
  full   + every station and channel not yet covered (in domain, nearest median), |theta_u| p90, |theta_w| p90,
         gain p5, clouds p95, patch-edge charge (if any), then greedy max-min fill in the standardized space
Writes NEW_DIR/selection.csv (order, set, slot, reason, features, cost), selection.md (coverage of each stratum vs the
candidate population), selection.json (inputs sha256, exclusions, parameters), pilot_slots.txt, full_slots.txt.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

P = 1.7      # strip pitch (mm) of the frozen readout; shift rule m = floor(Dbar/P + 1/2)
H = 5        # half strips of the frozen output window


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tie(slot):
    return int(hashlib.sha256(f"k2a-validation-{slot}".encode()).hexdigest()[:12], 16)


def f(x, d=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def moments_features(path):
    """Charge-weighted mean |A|, per-primary shifts (production rule) -> distinct groups, max |m|."""
    rows = list(csv.DictReader(open(path, newline="")))
    if not rows:
        return {}
    by = defaultdict(list)
    for r in rows:
        by[int(float(r["primary"]))].append(r)
    q = np.array([f(r["N_fC"]) for r in rows]); A = np.array([f(r["A"]) for r in rows])
    groups, mmax = set(), 0
    for p, cl in by.items():
        N = np.array([f(r["N_fC"]) for r in cl])
        w = N / N.sum()
        mu = math.floor(float((w * np.array([f(r["D_U_mm"]) for r in cl])).sum()) / P + .5)
        mv = math.floor(float((w * np.array([f(r["D_V_mm"]) for r in cl])).sum()) / P + .5)
        groups.add((mu, mv)); mmax = max(mmax, abs(mu), abs(mv))
    return {"mean_abs_A": float((q * np.abs(A)).sum() / q.sum()), "event_A": float((q * A).sum() / q.sum()),
            "shift_groups": len(groups), "max_abs_m": mmax, "ionising_primaries": len(by)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pf", type=Path, required=True); ap.add_argument("--recenter", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--exclude-run", type=Path)
    ap.add_argument("--exclude-slots", default="0000,0858,0560,0997,0133")
    ap.add_argument("--pilot", type=int, default=5); ap.add_argument("--full", type=int, default=16)
    ap.add_argument("--include-previously-inspected", action="store_true")
    ap.add_argument("--max-rows", type=float, default=0, help="skip candidates with more K256 rows (0 = no limit)")
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    if not 3 <= a.pilot <= 5 or a.full < a.pilot:
        raise SystemExit("--pilot must be 3..5 and --full >= --pilot")
    slots_p = a.pf / "inventory" / "slots.csv"; dom_p = a.pf / "domain" / "per_event.csv"
    rec_p = a.recenter / "per_event_recenter.csv"
    slots = {r["slot"]: r for r in csv.DictReader(open(slots_p, newline=""))}
    sealed = {s for s, r in slots.items() if r["role"] == "sealed_validation"}
    dom = {r["event"]: r for r in csv.DictReader(open(dom_p, newline=""))}
    rec = {r["event"]: r for r in csv.DictReader(open(rec_p, newline=""))}
    for name, d in (("domain/per_event.csv", dom), ("per_event_recenter.csv", rec)):
        bad = [e for e in d if e.replace("muon_", "") in sealed]
        if bad:
            raise SystemExit(f"refusing: sealed event(s) {bad[:3]} in {name}")
    roles = {"pool", "development"} | ({"previously_inspected"} if a.include_previously_inspected else set())
    excluded = {}
    for s in (x.strip() for x in a.exclude_slots.split(",")):
        if s:
            excluded[f"{int(s):04d}"] = "--exclude-slots (used for development decisions)"
    if a.exclude_run:
        for p in sorted(a.exclude_run.glob("provenance_prepare_*.txt")):
            excluded[p.stem.split("_")[-1]] = f"prepared in {a.exclude_run.name} (kernel / Green-table muon)"
    mfiles = {p.name.replace(".moments.csv", ""): p for p in (a.pf / "moments_nonsealed").rglob("muon_*.moments.csv")}
    cands, skipped = [], defaultdict(int)
    for s, r in sorted(slots.items()):
        e = f"muon_{s}"
        if r["role"] == "sealed_validation":
            continue                                    # never read
        if s in excluded:
            skipped["excluded: " + excluded[s].split(" (")[0]] += 1; continue
        if r["role"] not in roles or r.get("status", "USABLE") != "USABLE":
            skipped[f"role/status {r['role']}/{r.get('status')}"] += 1; continue
        d, rr = dom.get(e), rec.get(e)
        if not d or not rr or d["status"] in ("READ_FAILED", "NO_IONS", "NOT_CONVERTED", "META_MISMATCH") or e not in mfiles:
            skipped["not converted / no ions / metadata"] += 1; continue
        mf = moments_features(mfiles[e])
        if not mf:
            skipped["empty moments"] += 1; continue
        n_ion_p = mf["ionising_primaries"]
        c = {"slot": s, "event": e, "role": r["role"], "channel": d.get("channel", ""), "station": int(f(d.get("station"), -1)),
             "abs_theta_u": abs(f(d.get("theta_u_deg"))), "abs_theta_w": abs(f(d.get("theta_w_deg"))),
             "log_gain": math.log(f(d.get("ion_births_exact")) / max(n_ion_p, 1)), "mean_abs_A": mf["mean_abs_A"],
             "n_wires_5pct": int(f(d.get("n_wires_5pct"), 0)), "dominant_wire_share": f(d.get("dominant_wire_share")),
             "max_abs_m": mf["max_abs_m"], "outside_base_window": f(rr.get("charge_local_strip_outside_base_window"), 0.0),
             "clouds": int(f(d.get("clouds"), 0)), "patch_edge_charge": f(d.get("patch_edge_wire_charge_frac"), 0.0),
             "domain_after": rr.get("status_after", "?"), "causes_after": rr.get("causes_after", ""),
             "rows": int(f(d.get("rows"), 0)), "shift_groups": mf["shift_groups"], "heed_primaries": int(f(d.get("n_primaries"), 0))}
        if a.max_rows and c["rows"] > a.max_rows:
            skipped["--max-rows"] += 1; continue
        cands.append(c)
    if len(cands) < a.full:
        raise SystemExit(f"only {len(cands)} candidates for --full {a.full}")
    FE = ["abs_theta_u", "abs_theta_w", "log_gain", "mean_abs_A", "n_wires_5pct", "max_abs_m", "clouds"]
    X = np.array([[c[k] for k in FE] for c in cands], float)
    mu, sd = np.nanmedian(X, 0), np.nanstd(X, 0); sd[sd == 0] = 1
    Z = np.nan_to_num((X - mu) / sd)
    st = sorted({c["station"] for c in cands})
    Z = np.hstack([Z, np.array([[1.0 if c["station"] == s else 0.0 for s in st] for c in cands])])
    ind = lambda c: c["domain_after"] == "IN_DOMAIN"
    picked, order = [], []

    def take(reason, ok, score):
        """best score among ok candidates not yet picked; among the 5 best (deterministic), the cheapest by rows."""
        pool = [i for i, c in enumerate(cands) if i not in picked and ok(c)]
        if not pool:
            order.append({"reason": reason, "slot": None}); return
        pool.sort(key=lambda i: (score(i), tie(cands[i]["slot"])))
        best = min(pool[:5], key=lambda i: (cands[i]["rows"], tie(cands[i]["slot"]))) if reason != "typical" else pool[0]
        picked.append(best); order.append({"reason": reason, "slot": cands[best]["slot"], "i": best})

    gq = np.nanpercentile([c["log_gain"] for c in cands], [5, 95])
    med = lambda i: float(np.sum(Z[i, :len(FE)] ** 2))
    pilot_targets = [
        ("typical", lambda c: ind(c) and c["max_abs_m"] <= 1, med),
        ("multiwire_asym", lambda c: ind(c) and c["n_wires_5pct"] >= 2, lambda i: -cands[i]["mean_abs_A"]),
        ("displaced", lambda c: ind(c) and c["max_abs_m"] >= 2, lambda i: -cands[i]["outside_base_window"] - 1e-3 * cands[i]["max_abs_m"]),
        ("gain_tail", lambda c: ind(c) and c["log_gain"] >= gq[1], lambda i: -cands[i]["log_gain"]),
        ("out_of_domain", lambda c: c["domain_after"] == "OUT_OF_DOMAIN", lambda i: cands[i]["rows"]),
    ][:a.pilot]
    for r_, ok, sc in pilot_targets:
        take(r_, ok, sc)
    n_pilot = len(picked)
    covered = lambda key, val: any(cands[i][key] == val for i in picked)
    for s in st:
        if len(picked) < a.full and not covered("station", s):
            take(f"station_{s}", lambda c, s=s: ind(c) and c["station"] == s, med)
    for ch in sorted({c["channel"] for c in cands}):
        if len(picked) < a.full and not covered("channel", ch):
            take(f"channel_{ch}", lambda c, ch=ch: ind(c) and c["channel"] == ch, med)
    extra = [("theta_u_p90", lambda c: ind(c), lambda i: -cands[i]["abs_theta_u"]),
             ("theta_w_p90", lambda c: ind(c), lambda i: -cands[i]["abs_theta_w"]),
             ("gain_p5", lambda c: ind(c) and c["log_gain"] <= gq[0], lambda i: cands[i]["log_gain"]),
             ("clouds_p95", lambda c: ind(c), lambda i: -cands[i]["clouds"]),
             ("patch_edge", lambda c: ind(c) and c["patch_edge_charge"] > 0, lambda i: -cands[i]["patch_edge_charge"])]
    for r_, ok, sc in extra:
        if len(picked) < a.full:
            take(r_, ok, sc)
    while len(picked) < a.full:                                   # greedy max-min fill
        rest = [i for i in range(len(cands)) if i not in picked]
        dmin = [min(float(np.sum((Z[i] - Z[j]) ** 2)) for j in picked) for i in rest]
        best = rest[int(np.lexsort(([tie(cands[i]["slot"]) for i in rest], -np.array(dmin)))[0])]
        picked.append(best); order.append({"reason": "maxmin_fill", "slot": cands[best]["slot"], "i": best})
    a.out.mkdir(parents=True)
    rows = []
    k = 0
    for o in order:
        if o["slot"] is None:
            rows.append({"order": "", "set": "", "slot": "", "reason": o["reason"] + " (no candidate)"}); continue
        k += 1
        c = cands[o["i"]]
        rows.append({"order": k, "set": "pilot" if k <= n_pilot else "full", "reason": o["reason"], **c})
    cols = ["order", "set", "slot", "reason", "event", "role", "channel", "station", "abs_theta_u", "abs_theta_w", "log_gain",
            "mean_abs_A", "n_wires_5pct", "dominant_wire_share", "max_abs_m", "outside_base_window", "clouds",
            "patch_edge_charge", "domain_after", "causes_after", "rows", "shift_groups", "heed_primaries"]
    with open(a.out / "selection.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", restval=""); w.writeheader(); w.writerows(rows)
    sel = [r for r in rows if r.get("slot")]
    (a.out / "pilot_slots.txt").write_text(" ".join(r["slot"] for r in sel if r["set"] == "pilot") + "\n")
    (a.out / "full_slots.txt").write_text(" ".join(r["slot"] for r in sel) + "\n")
    meta = {"tool": "k2a_validation_select.py", "inputs": {str(p): sha(p) for p in (slots_p, dom_p, rec_p)},
            "candidates": len(cands), "skipped": dict(skipped), "excluded": excluded, "roles": sorted(roles),
            "pilot": [r["slot"] for r in sel if r["set"] == "pilot"], "full": [r["slot"] for r in sel],
            "missing_targets": [r["reason"] for r in rows if not r.get("slot")], "max_rows": a.max_rows}
    (a.out / "selection.json").write_text(json.dumps(meta, indent=1) + "\n")
    L = ["# K2A-vs-S0 validation selection (independent, non-sealed)", "",
         f"{len(cands)} candidates (roles {', '.join(sorted(roles))}); excluded before reading features: "
         + ", ".join(f"{s} ({v.split(' (')[0]})" for s, v in sorted(excluded.items())), "",
         "| # | set | slot | reason | channel | station | abs th_u | abs th_w | ln gain | mean abs A | wires>=5% | max abs m | clouds | domain after | rows | groups |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if not r.get("slot"):
            L.append(f"| | | | {r['reason']} | | | | | | | | | | | | |"); continue
        L.append(f"| {r['order']} | {r['set']} | {r['slot']} | {r['reason']} | {r['channel']} | {r['station']} | {r['abs_theta_u']:.1f} | "
                 f"{r['abs_theta_w']:.1f} | {r['log_gain']:.2f} | {r['mean_abs_A']:.3f} | {r['n_wires_5pct']} | {r['max_abs_m']} | "
                 f"{r['clouds']} | {r['domain_after']} | {r['rows']} | {r['shift_groups']} |")
    L += ["", "## Coverage: selection vs candidate population (min / median / max)", "", "| feature | population | selection |", "|---|---|---|"]
    for key in FE + ["station"]:
        pv = np.array([c[key] for c in cands], float); sv = np.array([float(r[key]) for r in sel])
        L.append(f"| {key} | {np.nanmin(pv):.3g} / {np.nanmedian(pv):.3g} / {np.nanmax(pv):.3g} | "
                 f"{np.nanmin(sv):.3g} / {np.nanmedian(sv):.3g} / {np.nanmax(sv):.3g} |")
    L.append(f"| channels | {sorted({c['channel'] for c in cands})} | {sorted({r['channel'] for r in sel})} |")
    L.append(f"| domain after recentering | {dict((s, sum(c['domain_after'] == s for c in cands)) for s in {c['domain_after'] for c in cands})} | "
             f"{dict((s, sum(r['domain_after'] == s for r in sel)) for s in {r['domain_after'] for r in sel})} |")
    if meta["missing_targets"]:
        L.append(f"\nTargets without a candidate: {', '.join(meta['missing_targets'])}")
    (a.out / "selection.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:4 + len(rows) + 2]))
    print(f"pilot: {meta['pilot']}\nfull:  {meta['full']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Deterministic supplementary DEVELOPMENT selection of non-sealed phi and omega muons (input metadata only).

  k2a_supplementary_selection.py --pf PREFLIGHT_DIR --recenter RECENTER_DIR --out NEW_DIR [--per-channel 8]
                                 [--include-previously-inspected]

Candidates: role 'pool' (optionally 'previously_inspected') of k2a_inventory.py, channel phi or omega, a converted
event with clouds (not READ_FAILED / NO_IONS / META_MISMATCH / NOT_CONVERTED). The sealed 10 phi + 10 omega and the
20 J/psi development muons are excluded by role BEFORE any feature is read; the tool refuses if a sealed id appears.
Features come from domain/per_event.csv and k2a_recenter_recovery.py's per_event_recenter.csv (cloud moments and
muon metadata only; no waveform, no result of any model).

Per channel, in this order (each pick is the candidate nearest to the target, ties broken by a fixed hash of the slot):
  typical       nearest to the channel median in the standardized feature space
  clouds_p95    cloud count nearest the channel p95            (large multiplicity, not the single extreme)
  clouds_max    largest cloud count                            (stress)
  multiwire     most anodes with >= 5 % of the charge, then smallest dominant-wire share
  far_D         largest share of charge more than 5 strips from the entry strip (needs recentering), then max |D|
  ood_<cause>   one event per cause still out of domain AFTER recentering (WIRE, R_EFF_LOW, R_EFF_HIGH, T_WINDOW, ...)
  edge_wire     largest charge share on |k| = 4 (field-patch edge), if any
then a greedy max-min fill in the standardized space (with station one-hot) up to --per-channel.
Features: log clouds, log charge, wires >= 5 %, dominant-wire share, |theta_u|, |theta_w|, D_U span, union strips,
K256 stored fraction. Writes supplementary_dev.csv (index, channel, reasons, features), .md, .json, SHA256SUMS.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

FEATS = ("log_clouds", "log_charge", "n_wires_5pct", "dominant_wire_share", "abs_theta_u", "abs_theta_w",
         "D_U_span", "union_strips_U", "saved_fraction")


def tie(slot):
    return int(hashlib.sha256(f"k2a-supp-{slot}".encode()).hexdigest()[:12], 16)


def fnum(x, d=np.nan):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pf", type=Path, required=True); ap.add_argument("--recenter", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--per-channel", type=int, default=8)
    ap.add_argument("--include-previously-inspected", action="store_true")
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    slots = {f"muon_{r['slot']}": r for r in csv.DictReader(open(a.pf / "inventory" / "slots.csv", newline=""))}
    roles = {"pool"} | ({"previously_inspected"} if a.include_previously_inspected else set())
    allowed = {k for k, v in slots.items() if v["role"] in roles and v["status"] == "USABLE"}
    sealed = {k for k, v in slots.items() if v["role"] == "sealed_validation"}
    pe = {r["event"]: r for r in csv.DictReader(open(a.pf / "domain" / "per_event.csv", newline="")) if r["event"] in allowed}
    rc = {r["event"]: r for r in csv.DictReader(open(a.recenter / "per_event_recenter.csv", newline="")) if r["event"] in allowed}
    if (set(pe) | set(rc)) & sealed:
        raise SystemExit("sealed id among the candidates: refusing")
    cand = []
    for e, r in pe.items():
        if r["status"] in ("READ_FAILED", "NO_IONS", "META_MISMATCH", "NOT_CONVERTED") or r.get("channel") not in ("phi", "omega"):
            continue
        x = rc.get(e, {})
        if not x or not x.get("clouds"):
            continue
        cand.append({"event": e, "slot": e[5:], "channel": r["channel"], "station": r.get("station", ""),
                     "status_before": r["status"], "status_after": x.get("status_after", ""),
                     "causes_after": x.get("causes_after", ""),
                     "log_clouds": np.log(fnum(r["clouds"])), "log_charge": np.log(fnum(r["charge_fC"])),
                     "n_wires_5pct": fnum(r["n_wires_5pct"]), "dominant_wire_share": fnum(r["dominant_wire_share"]),
                     "abs_theta_u": abs(fnum(r.get("theta_u_deg"))), "abs_theta_w": abs(fnum(r.get("theta_w_deg"))),
                     "D_U_span": fnum(x.get("D_U_span")), "union_strips_U": fnum(x.get("union_strips_U")),
                     "saved_fraction": fnum(x.get("saved_fraction"), 1.0),
                     "far_charge": fnum(x.get("charge_local_strip_outside_base_window"), 0.0),
                     "absD": max(abs(fnum(x.get(k), 0)) for k in ("D_U_min", "D_U_max", "D_V_min", "D_V_max")),
                     "edge": fnum(r.get("patch_edge_wire_charge_frac"), 0.0), "clouds": int(fnum(r["clouds"]))})
    a.out.mkdir(parents=True)
    chosen = []
    for ch in ("phi", "omega"):
        C = sorted([c for c in cand if c["channel"] == ch], key=lambda c: tie(c["slot"]))
        if not C:
            continue
        M = np.array([[c[f] for f in FEATS] for c in C], float)
        M = np.where(np.isfinite(M), M, np.nanmedian(M, 0))
        sd = M.std(0); sd[sd == 0] = 1
        Z = np.nan_to_num((M - M.mean(0)) / sd)          # a feature missing for every candidate contributes 0
        st = sorted({c["station"] for c in C})
        Zs = np.hstack([Z, np.array([[1.0 if c["station"] == s else 0.0 for s in st] for c in C])])
        picked, reasons = [], {}

        def take(i, why):
            if i is None:
                return
            if i in picked:
                reasons[i].append(why); return
            if len(picked) < a.per_channel:
                picked.append(i); reasons[i] = [why]

        def argbest(key, pool=None):
            pool = range(len(C)) if pool is None else pool
            pool = list(pool)
            return min(pool, key=lambda i: (key(i), tie(C[i]["slot"]))) if pool else None

        med = np.median(Z, 0)
        take(argbest(lambda i: float(np.sum((Z[i] - med) ** 2))), "typical")
        p95 = np.quantile([c["clouds"] for c in C], .95)
        take(argbest(lambda i: abs(C[i]["clouds"] - p95)), "clouds_p95")
        take(argbest(lambda i: -C[i]["clouds"]), "clouds_max")
        take(argbest(lambda i: (-C[i]["n_wires_5pct"], C[i]["dominant_wire_share"])), "multiwire")
        far = [i for i in range(len(C)) if C[i]["far_charge"] > 0]
        take(argbest(lambda i: (-C[i]["far_charge"], -C[i]["absD"]), far), "far_D")
        causes = sorted({k.split(":")[0].strip() for c in C for k in c["causes_after"].split(";") if k.strip()})
        for k in causes:
            take(argbest(lambda i: -C[i]["clouds"], [i for i in range(len(C)) if k in C[i]["causes_after"]]), f"ood_{k}")
        edge = [i for i in range(len(C)) if C[i]["edge"] > 0]
        take(argbest(lambda i: -C[i]["edge"], edge), "edge_wire")
        while len(picked) < min(a.per_channel, len(C)):
            d = np.array([min(float(np.sum((Zs[i] - Zs[j]) ** 2)) for j in picked) if picked else 0.0 for i in range(len(C))])
            rest = [i for i in range(len(C)) if i not in picked]
            best = max(rest, key=lambda i: (d[i], -tie(C[i]["slot"])))
            take(best, "maxmin_fill")
        for i in picked:
            chosen.append({**C[i], "reasons": "+".join(reasons[i])})
    fields = ["slot", "event", "channel", "station", "reasons", "status_before", "status_after", "causes_after", "clouds",
              "far_charge", "absD", "edge"] + list(FEATS)
    with open(a.out / "supplementary_dev.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore", lineterminator="\n"); w.writeheader()
        for c in chosen:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in c.items()})
    L = ["# Supplementary development selection (non-sealed phi / omega; input metadata only)", "",
         f"candidates: {sum(1 for c in cand if c['channel'] == 'phi')} phi, {sum(1 for c in cand if c['channel'] == 'omega')} omega "
         f"(roles {sorted(roles)}); per channel {a.per_channel}; sealed set excluded by role before reading features", "",
         "| slot | channel | station | reasons | clouds | wires>=5% | far-D charge | max abs D | status before -> after recentering |",
         "|---|---|---|---|---|---|---|---|---|"]
    for c in chosen:
        L.append(f"| {c['slot']} | {c['channel']} | {c['station']} | {c['reasons']} | {c['clouds']} | {c['n_wires_5pct']:.0f} | "
                 f"{c['far_charge']:.3g} | {c['absD']:.2f} | {c['status_before']} -> {c['status_after']} |")
    (a.out / "supplementary_dev.md").write_text("\n".join(L) + "\n")
    (a.out / "supplementary_dev.json").write_text(json.dumps({"tool": "k2a_supplementary_selection.py", "per_channel": a.per_channel,
        "roles": sorted(roles), "indices": [c["slot"] for c in chosen],
        "inputs_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                          (a.pf / "inventory" / "slots.csv", a.pf / "domain" / "per_event.csv", a.recenter / "per_event_recenter.csv")}},
        indent=1) + "\n")
    (a.out / "SHA256SUMS").write_text("".join(f"{hashlib.sha256((a.out / n).read_bytes()).hexdigest()}  {n}\n"
                                              for n in ("supplementary_dev.csv", "supplementary_dev.md", "supplementary_dev.json")))
    print("\n".join(L))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""READ-ONLY, input metadata only: what per-primary strip recentering would recover, and why D is large.

  k2a_recenter_recovery.py --pf PREFLIGHT_DIR --kernel KERNEL_DIR --out NEW_DIR [--half-strips 5 --window-ns 7000]

Uses only the NON-SEALED preflight products (moments_nonsealed/**: *.moments.csv, *.event.json, *.frame.json;
inventory/slots.csv; domain/per_event.csv). Refuses if a sealed event id appears. No waveform, no archive access, no
SHA, no conversion.

Strip frame (exporter / BirthMoments.hh): base_U = llround((w_in - A u_in)/P) from the muon ENTRY point,
  D_U = w_phys - A patch_shift - P base_U     (mm; the cloud's longitudinal coordinate relative to the entry strip,
                                               WITHOUT the -A u_patch stereo term, which the kernel/LUT adds per wire)
so the strip coordinate of a cloud is eta_U = D_U - A u_patch, eta_V = D_V + A u_patch (u_patch ~ 3 k).
Recentering (per Heed primary, so that the frozen S0 replay can be run per group with integer weight sums):
  m_U = round(Dbar_U / P), m_V = round(Dbar_V / P)       (charge-weighted over the primary's clouds)
  D'_U = D_U - m_U P on strips base_U + m_U + o, o = -H..H (likewise V)  -> accumulate on ABSOLUTE strip indices.
In the frozen model (2D field, infinite contiguous identical strips, kernel/LUT a function of eta - jP only) this is an
exact relabelling; it is checked against the frozen S0 by k2a_recenter_experiment.sh, not assumed here.

Reports: (1) D anatomy: D span vs the track's w extent across the gap (6 mm |tan theta_w|), entry-offset term,
top events; charge whose local strip lies outside the current base +-H output (signal the current 11-strip output
cannot contain); (2) recovery: per out-of-domain event, the causes left after recentering (WIRE, R_EFF_LOW/HIGH,
T_WINDOW, D after recentering, NONFINITE); (3) radial grid: inner node clearance, log step, why one inward log step is
invalid, and the proposed clearance-geometric inner nodes; (4) K256 reservoir fraction from event.json;
(5) provenance of the preflight (SOURCE_GIT_SHA.txt, which domain-check version wrote coverage.md).
"""
import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from k2a_domain_check import Domain, read_moments, R_WIRE_UM  # noqa: E402

A, P = 0.1, 1.7


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    rx, ry = np.argsort(np.argsort(x[ok])), np.argsort(np.argsort(y[ok]))
    return float(np.corrcoef(rx, ry)[0, 1])


def primary_shifts(X):
    """m_U, m_V per Heed primary: nearest integer of the charge-weighted D / P (so |D'| <= P/2 + the primary's own
    spread; the stereo term A u_patch, <= 1.35 mm, stays inside the kernel, as now). Returned per cloud."""
    q = X["N_fC"]; u = 3.0 * X["wire"]
    mU = np.zeros(len(q), int); mV = np.zeros(len(q), int)
    for p in np.unique(X["primary"]):
        s = X["primary"] == p
        w = q[s] / q[s].sum()
        dU, dV, ub = (w * X["D_U_mm"][s]).sum(), (w * X["D_V_mm"][s]).sum(), (w * u[s]).sum()
        mU[s] = int(np.floor(dU / P + 0.5)); mV[s] = int(np.floor(dV / P + 0.5))
    return mU, mV


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pf", type=Path, required=True); ap.add_argument("--kernel", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--half-strips", type=int, default=5); ap.add_argument("--window-ns", type=float, default=7000.0)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    H = a.half_strips
    dom = Domain(a.kernel, H, a.window_ns)
    slots = {f"muon_{r['slot']}": r for r in csv.DictReader(open(a.pf / "inventory" / "slots.csv", newline=""))}
    sealed = {k for k, v in slots.items() if v["role"] == "sealed_validation"}
    pe = {r["event"]: r for r in csv.DictReader(open(a.pf / "domain" / "per_event.csv", newline=""))}
    if set(pe) & sealed:
        raise SystemExit("sealed event in domain/per_event.csv: refusing")
    files = sorted((a.pf / "moments_nonsealed").rglob("muon_*.moments.csv"))
    rows, causes_after, cloud_r = [], Counter(), {"low": [], "high": []}
    sf = defaultdict(list)
    for f in files:
        eid = re.sub(r"\.moments\.csv$", "", f.name)
        if eid in sealed:
            raise SystemExit(f"sealed event {eid} found among the non-sealed moments: refusing")
        X, bu, bv = read_moments(f)
        evp = Path(str(f) + ".event.json"); frp = Path(str(f) + ".frame.json")
        ev = json.loads(evp.read_text()) if evp.exists() else {}
        fr = json.loads(frp.read_text()) if frp.exists() else {}
        role = slots.get(eid, {}).get("role", "?"); ch = pe.get(eid, {}).get("channel", slots.get(eid, {}).get("channel", ""))
        rec = {"event": eid, "role": role, "channel": ch, "station": pe.get(eid, {}).get("station", ""),
               "status_before": pe.get(eid, {}).get("status", "?"), "causes_before": pe.get(eid, {}).get("causes", "")}
        if ev.get("ion_rows_saved") not in (None, "") and float(ev.get("ion_births_exact") or 0) > 0:
            rec["saved_fraction"] = float(ev["ion_rows_saved"]) / float(ev["ion_births_exact"])
            sf[ch].append(rec["saved_fraction"])
        n = len(X["t_ns"])
        if n == 0:
            rec["status_after"] = rec["status_before"]; rows.append(rec); continue
        q = X["N_fC"]
        mU, mV = primary_shifts(X)
        DpU, DpV = X["D_U_mm"] - P * mU, X["D_V_mm"] - P * mV
        thw = float(ev.get("theta_w_deg", "nan")); thu = float(ev.get("theta_u_deg", "nan"))
        entry_term = A * (float(fr.get("u_in_mm", "nan")) - 10 * float(fr.get("patch_shift_cm", "nan")))
        rec.update(clouds=n, charge_fC=float(q.sum()), D_U_min=float(X["D_U_mm"].min()), D_U_max=float(X["D_U_mm"].max()),
                   D_V_min=float(X["D_V_mm"].min()), D_V_max=float(X["D_V_mm"].max()),
                   D_U_span=float(np.ptp(X["D_U_mm"])), gap_w_extent_mm=6.0 * abs(math.tan(math.radians(thw))) if thw == thw else "",
                   theta_w_deg=thw, theta_u_deg=thu, path_length_mm=ev.get("path_length_mm", ""),
                   momentum_GeV=ev.get("momentum_GeV", ""), entry_offset_term_mm=entry_term,
                   m_U_min=int(mU.min()), m_U_max=int(mU.max()), m_V_min=int(mV.min()), m_V_max=int(mV.max()),
                   charge_local_strip_outside_base_window=float(q[(np.abs(mU) > H) | (np.abs(mV) > H)].sum() / q.sum()),
                   charge_local_strip_beyond_3=float(q[(np.abs(mU) > 3) | (np.abs(mV) > 3)].sum() / q.sum()),
                   union_strips_U=int(mU.max() - mU.min() + 2 * H + 1), union_strips_V=int(mV.max() - mV.min() + 2 * H + 1),
                   Dp_abs_max=float(max(np.abs(DpU).max(), np.abs(DpV).max())))
        # domain after recentering (same tests as the engine, D on the recentered values)
        left = {}
        tests = {"WIRE": ~np.isin(X["wire"], dom.wires), "R_EFF_LOW": X["r_eff_um"] < dom.rn[0],
                 "R_EFF_HIGH": X["r_eff_um"] > dom.rn[-1], "D_U_RECENTERED": ~dom.d_ok(DpU), "D_V_RECENTERED": ~dom.d_ok(DpV),
                 "T_WINDOW": (X["t_ns"] >= a.window_ns) | (X["t_ns"] < 0),
                 "NONFINITE": ~np.all(np.isfinite(np.stack([X[c] for c in ("t_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm")])), 0)}
        for k, m in tests.items():
            if m.any():
                left[k] = (int(m.sum()), float(q[m].sum() / q.sum()))
        cloud_r["low"] += [float(x) for x in X["r_eff_um"][tests["R_EFF_LOW"]]]
        cloud_r["high"] += [float(x) for x in X["r_eff_um"][tests["R_EFF_HIGH"]]]
        rec["ion_r_min_um"] = float((X["min_clearance_um"] + R_WIRE_UM).min())
        rec["wires"] = " ".join(map(str, sorted(set(X["wire"].tolist()))))
        rec["causes_after"] = "; ".join(f"{k}: {c} clouds, {fq:.4%} charge" for k, (c, fq) in left.items())
        if rec["status_before"] in ("READ_FAILED", "NOT_CONVERTED", "META_MISMATCH"):
            rec["status_after"] = rec["status_before"]
        else:
            rec["status_after"] = "IN_DOMAIN" if not left else "OUT_OF_DOMAIN"
        for k in left:
            causes_after[k] += 1
        rows.append(rec)
    a.out.mkdir(parents=True)
    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open(a.out / "per_event_recenter.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n", restval=""); w.writeheader(); w.writerows(rows)
    ev_ = [r for r in rows if "clouds" in r]
    before = Counter(r["status_before"] for r in rows); after = Counter(r["status_after"] for r in rows)
    ood = [r for r in ev_ if r["status_before"] == "OUT_OF_DOMAIN"]
    rec_ = [r for r in ood if r["status_after"] == "IN_DOMAIN"]
    still = [r for r in ood if r["status_after"] != "IN_DOMAIN"]
    chg = lambda rs: sum(r["charge_fC"] for r in rs)
    # charge coverage before/after (per cloud, using the same tests)
    L = ["# Per-primary strip recentering: D anatomy and coverage recovery (non-sealed input metadata only)", "",
         f"kernel {a.kernel}: wires {dom.wires}, r [{dom.rn[0]:.6g}, {dom.rn[-1]:.6g}] um, cloud D support "
         f"[{dom.D_lo:.3f}, {dom.D_hi:.3f}) mm for ±{H} strips; events {len(rows)} (with clouds {len(ev_)})", "",
         "## 1. Why D reaches -14 ... +21 mm", ""]
    sp = spearman([r["D_U_span"] for r in ev_], [r["gap_w_extent_mm"] if r["gap_w_extent_mm"] != "" else np.nan for r in ev_])
    absD = [max(abs(r["D_U_min"]), abs(r["D_U_max"]), abs(r["D_V_min"]), abs(r["D_V_max"])) for r in ev_]
    L += [f"* D is measured from the strip of the muon ENTRY point; |A (u_in - patch_shift)| (entry-offset term) is at most "
          f"{max([abs(r['entry_offset_term_mm']) for r in ev_ if r['entry_offset_term_mm'] == r['entry_offset_term_mm']] or [float('nan')]):.3f} mm, "
          "so large |D| is the clouds' own longitudinal distance from the entry point.",
          f"* Spearman rank correlation of the per-event D_U span with the track's w extent across the 6 mm gap "
          f"(6 |tan theta_w|): {sp:.3f} (n = {len(ev_)}).",
          f"* max |D| quantiles over events: p50 {np.quantile(absD, .5):.2f}, p90 {np.quantile(absD, .9):.2f}, "
          f"p99 {np.quantile(absD, .99):.2f}, max {max(absD):.2f} mm.",
          f"* Events with charge whose local strip is MORE than {H} strips from the entry strip (signal the current "
          f"base-centred 11-strip output cannot contain, in S0 and K2A alike): "
          f"{sum(1 for r in ev_ if r['charge_local_strip_outside_base_window'] > 0)}; their charge fraction p50/max "
          + (lambda v: f"{np.quantile(v, .5):.3g} / {max(v):.3g}" if v else "n/a")(
              [r['charge_local_strip_outside_base_window'] for r in ev_ if r['charge_local_strip_outside_base_window'] > 0]),
          f"* absolute-strip union width per plane (local windows of all primaries): p50 "
          f"{np.quantile([r['union_strips_U'] for r in ev_], .5):.0f}, p95 {np.quantile([r['union_strips_U'] for r in ev_], .95):.0f}, "
          f"max {max(r['union_strips_U'] for r in ev_)} U strips", "",
          "Top 12 events by max |D|:", "",
          "| event | channel | station | D_U [min, max] | D_V [min, max] | theta_w | 6|tan theta_w| | path mm | m_U range | charge outside base ±H |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(ev_, key=lambda r: -max(abs(r["D_U_min"]), abs(r["D_U_max"]), abs(r["D_V_min"]), abs(r["D_V_max"])))[:12]:
        L.append(f"| {r['event']} | {r['channel']} | {r['station']} | [{r['D_U_min']:.2f}, {r['D_U_max']:.2f}] | "
                 f"[{r['D_V_min']:.2f}, {r['D_V_max']:.2f}] | {r['theta_w_deg']:.1f} | "
                 f"{r['gap_w_extent_mm'] if r['gap_w_extent_mm'] == '' else round(r['gap_w_extent_mm'], 2)} | {r['path_length_mm']} | "
                 f"{r['m_U_min']}..{r['m_U_max']} | {r['charge_local_strip_outside_base_window']:.3g} |")
    L += ["", "## 2. Coverage recovery by recentering alone (existing kernel, no new nodes)", "",
          f"status before: {dict(before)}", f"status after:  {dict(after)}", "",
          f"out-of-domain before: {len(ood)} events; recovered by recentering alone: {len(rec_)}; still out: {len(still)}",
          f"charge of out-of-domain events recovered: {chg(rec_):.4g} fC of {chg(ood):.4g} fC", "",
          f"events still out of domain, by remaining cause: {dict(causes_after)}", "",
          "| event | channel | station | causes before | causes after recentering | wires |", "|---|---|---|---|---|---|"]
    for r in still:
        L.append(f"| {r['event']} | {r['channel']} | {r['station']} | {r['causes_before']} | {r['causes_after']} | {r['wires']} |")
    by_ch = defaultdict(lambda: [0, 0, 0])
    for r in ev_:
        c = by_ch[r["channel"]]; c[0] += 1; c[1] += r["status_before"] == "IN_DOMAIN"; c[2] += r["status_after"] == "IN_DOMAIN"
    L += ["", "| channel | events | in domain before | in domain after recentering |", "|---|---|---|---|"]
    L += [f"| {k} | {v[0]} | {v[1]} ({v[1] / v[0]:.2%}) | {v[2]} ({v[2] / v[0]:.2%}) |" for k, v in sorted(by_ch.items())]
    # radial
    step = dom.rn[1] / dom.rn[0]
    inward = dom.rn[0] / step
    L += ["", "## 3. Radial grid near the wire", "",
          f"* innermost node {float(dom.rn[0])!r} um (clearance {dom.rn[0] - R_WIRE_UM:.6g} um), outermost {float(dom.rn[-1])!r} um, "
          f"log step factor {step:.4f}",
          f"* one log step inward would be r = {inward:.4f} um: clearance {inward - R_WIRE_UM:+.4f} um -> "
          + ("INSIDE THE ANODE (invalid)" if inward <= R_WIRE_UM else "still in gas"),
          f"* clouds with r_eff below the grid: {len(cloud_r['low'])}"
          + (f" (r_eff {min(cloud_r['low']):.6g} ... {max(cloud_r['low']):.6g} um)" if cloud_r["low"] else ""),
          f"* clouds with r_eff above the grid: {len(cloud_r['high'])}"
          + (f" (r_eff {min(cloud_r['high']):.6g} ... {max(cloud_r['high']):.6g} um)" if cloud_r["high"] else ""),
          "* proposed inner treatment: nodes geometric in CLEARANCE c = r - R (c_k = c_min (c_0/c_min)^(k/n), c_min = half the "
          "smallest observed clearance, c_0 = existing innermost clearance), never on/inside the wire; outer extension in the "
          "existing log step (no wire constraint). See RECENTERING_DESIGN_REVIEW.md section 4.", ""]
    L += ["## 4. K256 reservoir fraction (stored rows / exact ion births), from event.json", ""]
    for c, v in sorted(sf.items()):
        if v:
            L.append(f"* {c}: n {len(v)}; p5 {np.quantile(v, .05):.4g}, p50 {np.quantile(v, .5):.4g}, p95 {np.quantile(v, .95):.4g}, "
                     f"min {min(v):.4g}, max {max(v):.4g}; events with every birth stored (fraction 1): {sum(1 for x in v if x >= 1 - 1e-12)}")
    if not any(sf.values()):
        L.append("* event.json lacks ion_rows_saved: converter older than 47ac8ef")
    gs = a.pf / "SOURCE_GIT_SHA.txt"
    cov = (a.pf / "domain" / "coverage.md").read_text() if (a.pf / "domain" / "coverage.md").exists() else ""
    L += ["", "## 5. Provenance of the preflight", "",
          f"* SOURCE_GIT_SHA.txt: {gs.read_text().strip() if gs.exists() else 'missing'} (written by run_k2a_campaign.sh at the domain stage)",
          f"* coverage.md has the K256 reservoir line (added in 78e6f1f): {'yes' if 'K256 reservoir' in cov else 'no -> written by the 47ac8ef domain check'}",
          f"* coverage.md has the patch-edge line (47ac8ef and later): {'yes' if 'outermost wires' in cov else 'no'}", ""]
    (a.out / "recovery.md").write_text("\n".join(L) + "\n")
    (a.out / "recovery.json").write_text(json.dumps({
        "before": dict(before), "after": dict(after), "ood_before": len(ood), "recovered": len(rec_), "still_out": len(still),
        "still_out_causes": dict(causes_after), "recovered_ids": [r["event"] for r in rec_], "still_out_ids": [r["event"] for r in still],
        "inward_log_step_node_um": float(inward), "r_grid_um": [float(dom.rn[0]), float(dom.rn[-1])], "log_step": float(step),
        "reservoir_fraction": {c: {"n": len(v), "p50": float(np.median(v)) if v else None, "min": min(v) if v else None} for c, v in sf.items()},
        "spearman_Dspan_vs_gap_w_extent": sp}, indent=1) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fast-only operational domain scan of the FROZEN compiled K2A over every permitted non-sealed archived muon.
No S0, no ion transport, no avalanche, no archive read: the compact moments come from the non-sealed preflight
(PF/moments_nonsealed, written by p3_birth_moments --batch); sealed moments are never read.

  k2a_domain_scan.py NEW_DIR --pf PF --kernel KERNEL_DIR --kernel-id ID --eval k2a_strip_eval
                     [--roles ROLE,ROLE] [--slots 0001,0002] [--prune-waveforms]
                     [--check-moments VALIDATION_RUN ...]

Events: every slot of PF/inventory/slots.csv whose role is not sealed (or only --roles), status USABLE, with preflight
moments. Per event: k2a_strip_eval with its production grouping rule (exactly as the primary K2A of the validation; kernel id
verified on the first and the last event, mapped read-only in between). Recorded: execution (EVALUATED / BLOCKED /
ERROR), status names, clouds, groups, primaries, cloud charge, blocked charge and fraction, blocked clouds by cause
(count, charge, r_eff vs the kernel radial grid, local D vs the D table), informational flags (T_AFTER_WINDOW,
SHIFT_NEAR_TIE), max |m|, strips, evaluation time. Cross-checks: the preflight's own domain status and the
recentering recovery's status_after (agreement table); --check-moments: the preflight moments of every validation
muon must equal that run's archive moments (s_compare_moments at 1e-9), so the scan input is the production input.
Aggregate: statuses, blocked-fraction distribution and thresholds 1e-6 / 1e-5 / 1e-4 / 1e-3 (INFORMATION for a
future, separately frozen policy; nothing is reclassified here), causes, r_eff of blocked clouds, by channel/station.
--prune-waveforms removes only the waveform CSVs this scan wrote, after their sha256 is recorded (they are
reproducible bit-identically); default: keep everything (~1-2 MB per event).
Writes NEW_DIR/scan.csv, blocked_clouds.csv, scan.json, scan.md, provenance.txt.
"""
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCM = HERE.parents[1] / "spatial-ablation" / "s_compare_moments.py"
INFO = {"T_AFTER_WINDOW", "SHIFT_NEAR_TIE", "OK"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--pf", type=Path, required=True)
    ap.add_argument("--kernel", type=Path, required=True); ap.add_argument("--kernel-id", required=True)
    ap.add_argument("--eval", type=Path, required=True); ap.add_argument("--roles", default="", help="default: every non-sealed role")
    ap.add_argument("--slots", default=""); ap.add_argument("--prune-waveforms", action="store_true")
    ap.add_argument("--check-moments", type=Path, nargs="*", default=[])
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    slots = {r["slot"]: r for r in csv.DictReader(open(a.pf / "inventory" / "slots.csv"))}
    sealed = {s for s, r in slots.items() if "sealed" in r["role"]}
    dom = {r["event"]: r for r in csv.DictReader(open(a.pf / "domain" / "per_event.csv"))} if (a.pf / "domain" / "per_event.csv").is_file() else {}
    recs = sorted(a.pf.glob("recenter_*/per_event_recenter.csv"))
    rec = {r["event"]: r for r in csv.DictReader(open(recs[-1]))} if recs else {}
    mroot = a.pf / "moments_nonsealed"
    mfiles = {p.name.replace(".moments.csv", "")[5:]: p for p in mroot.rglob("muon_*.moments.csv")}
    if any(s in sealed for s in mfiles):
        raise SystemExit("refusing: a sealed slot has moments under moments_nonsealed")
    roles = set(a.roles.split(",")) if a.roles else {r["role"] for r in slots.values() if "sealed" not in r["role"]}
    want = {f"{int(x):04d}" for x in a.slots.split(",") if x.strip()}
    if want & sealed:
        raise SystemExit(f"refusing: sealed slot(s) requested {sorted(want & sealed)}")
    todo, skipped = [], Counter()
    for s, r in sorted(slots.items()):
        if s in sealed:
            continue
        if want and s not in want:
            continue
        if r["role"] not in roles:
            skipped[f"role {r['role']}"] += 1; continue
        if r.get("status", "USABLE") != "USABLE":
            skipped[f"status {r.get('status')}"] += 1; continue
        if s not in mfiles:
            skipped["no preflight moments (not converted / no ions)"] += 1; continue
        todo.append(s)
    a.out.mkdir(parents=True); (a.out / "events").mkdir()
    kj = json.load(open(a.kernel / "kernel.json"))
    ax = kj.get("green_axes", {}); rr = ax.get("r_um") or [float("nan")]; DD = ax.get("D_mm") or [float("nan")]
    prov = {"tool": "k2a_domain_scan.py", "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "pf": str(a.pf),
            "recovery": str(recs[-1]) if recs else None, "kernel": str(a.kernel), "kernel_id": a.kernel_id,
            "eval": str(a.eval), "eval_sha256": sha(a.eval), "roles": sorted(roles), "events": len(todo), "skipped": dict(skipped),
            "grouping": "production rule (no imposed shifts)", "kernel_r_um": [rr[0], rr[-1]], "kernel_D_mm": [DD[0], DD[-1]],
            "kernel_wires": ax.get("wires")}
    (a.out / "provenance.txt").write_text("".join(f"{k}={v}\n" for k, v in prov.items()))
    rows, cl_rows, pruned = [], [], []
    t_start = time.time()
    for i, s in enumerate(todo):
        d = a.out / "events" / f"muon_{s}"
        verify = i == 0 or i == len(todo) - 1
        cmd = [str(a.eval), "--kernel", str(a.kernel)] + (["--kernel-id", a.kernel_id] if verify else ["--no-verify"]) + \
              ["--moments", str(mfiles[s]), "--out", str(d)]
        t0 = time.time(); p = subprocess.run(cmd, capture_output=True, text=True); wall = time.time() - t0
        r = slots[s]; e = f"muon_{s}"; dm = dom.get(e, {}); rc_ = rec.get(e, {})
        row = {"slot": s, "role": r["role"], "channel": dm.get("channel", r.get("channel", "")), "station": dm.get("station", ""),
               "exit": p.returncode, "kernel_verified": verify, "process_wall_s": round(wall, 4),
               "preflight_status": dm.get("status", ""), "recovery_status_after": rc_.get("status_after", ""),
               "recovery_causes_after": rc_.get("causes_after", "")}
        res = None
        try:
            res = json.load(open(d / "result.json"))
        except (OSError, ValueError):
            pass
        if res is None or p.returncode not in (0, 4):
            row.update(execution="ERROR", status_names=(p.stderr or p.stdout).strip().splitlines()[-1][:200] if (p.stderr or p.stdout) else "")
            rows.append(row); continue
        q = res.get("cloud_charge_fC") or 0.0; qb = res.get("blocked_charge_fC") or 0.0
        row.update(execution="EVALUATED" if p.returncode == 0 else "BLOCKED", status_names=res["status_names"],
                   clouds=res["clouds"], groups=res["groups"], primaries=res["primaries"], strips_U=res["strips_U"],
                   strips_V=res["strips_V"], cloud_charge_fC=q, blocked_charge_fC=qb, blocked_frac=qb / q if q else None,
                   evaluate_s=res["timing"]["evaluate_cold_s"])
        causes, info, nblk = Counter(), Counter(), 0
        for c in csv.DictReader(open(d / "clouds.csv")):
            names = set(c["status_names"].split("|"))
            info.update(n for n in names if n in ("T_AFTER_WINDOW", "SHIFT_NEAR_TIE"))
            blk = sorted(names - INFO)
            if not blk:
                continue
            nblk += 1; causes.update(blk)
            reff = float(c["r_eff_um"])
            cl_rows.append({"slot": s, "primary": c["primary"], "wire": c["wire"], "status_names": c["status_names"],
                            "N_fC": float(c["N_fC"]), "charge_frac": float(c["N_fC"]) / q if q else None, "t_ns": float(c["t_ns"]),
                            "r_eff_um": reff, "r_eff_over_grid_max": reff / rr[-1], "r_eff_over_grid_min": reff / rr[0],
                            "D_U_local_mm": float(c["D_U_local_mm"]), "D_V_local_mm": float(c["D_V_local_mm"])})
        sh = list(csv.DictReader(open(d / "shifts.csv")))
        row.update(blocked_clouds=nblk, blocked_causes=";".join(f"{k}:{v}" for k, v in sorted(causes.items())),
                   t_after_window_clouds=info["T_AFTER_WINDOW"], near_tie_clouds=info["SHIFT_NEAR_TIE"],
                   max_abs_m=max([max(abs(int(x["m_U"])), abs(int(x["m_V"]))) for x in sh] or [0]))
        rows.append(row)
        if a.prune_waveforms:
            for f in ("replay/waveforms.csv", "prompt/prompt_waveforms.csv"):
                if (d / f).is_file():
                    pruned.append(f"{sha(d / f)}  {d / f}"); (d / f).unlink()
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(todo)} events, {time.time() - t_start:.0f} s", flush=True)
    if pruned:
        (a.out / "pruned_waveforms.sha256").write_text("\n".join(pruned) + "\n")
    cols = ["slot", "role", "channel", "station", "execution", "status_names", "blocked_frac", "blocked_charge_fC", "cloud_charge_fC",
            "blocked_clouds", "blocked_causes", "clouds", "groups", "primaries", "max_abs_m", "strips_U", "strips_V",
            "t_after_window_clouds", "near_tie_clouds", "evaluate_s", "process_wall_s", "kernel_verified", "exit",
            "preflight_status", "recovery_status_after", "recovery_causes_after"]
    with open(a.out / "scan.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with open(a.out / "blocked_clouds.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cl_rows[0]) if cl_rows else ["slot"]); w.writeheader(); w.writerows(cl_rows)

    # moments cross-check against validation runs (the scan input must be the production input)
    mcheck = {}
    for run in a.check_moments:
        for V in sorted((run / "muons").glob("muon_*")):
            s = V.name[5:]; pm = V / "prod" / "birth_moments.csv"
            if s in mfiles and pm.is_file():
                p = subprocess.run([sys.executable, str(SCM), str(pm), str(mfiles[s])], capture_output=True, text=True)
                mcheck[s] = "PASS" if p.returncode == 0 else "FAIL"
    ex = Counter(r["execution"] for r in rows)
    ev = [r for r in rows if r["execution"] in ("EVALUATED", "BLOCKED")]
    fr = np.array([r["blocked_frac"] for r in ev if r["execution"] == "BLOCKED" and r["blocked_frac"] is not None], float)
    agree = Counter((r["execution"], r.get("recovery_status_after") or "-") for r in ev)
    by = defaultdict(Counter)
    for r in ev:
        by[f"channel {r['channel']}"][r["execution"]] += 1; by[f"station {r['station']}"][r["execution"]] += 1
    cause_ev = Counter(); cause_q = defaultdict(float)
    for c in cl_rows:
        for n in c["status_names"].split("|"):
            if n not in INFO:
                cause_ev[n] += 1; cause_q[n] += c["N_fC"]
    reff = np.array([c["r_eff_over_grid_max"] for c in cl_rows if "R_EFF_ABOVE_GRID" in c["status_names"]], float)
    summ = {"provenance": prov, "execution_counts": dict(ex), "events_scanned": len(rows),
            "in_domain_frac": ex["EVALUATED"] / len(ev) if ev else None,
            "blocked_frac_quantiles": {q: float(np.quantile(fr, q)) for q in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0)} if len(fr) else None,
            "blocked_events_below": {f"{t:g}": int(np.sum(fr < t)) for t in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)} if len(fr) else None,
            "blocked_clouds_by_cause": dict(cause_ev), "blocked_charge_by_cause_fC": dict(cause_q),
            "r_eff_above_grid_ratio_quantiles": {q: float(np.quantile(reff, q)) for q in (0.0, 0.5, 0.9, 1.0)} if len(reff) else None,
            "agreement_with_recovery": {f"{k[0]} / recovery {k[1]}": v for k, v in sorted(agree.items())},
            "by_channel_station": {k: dict(v) for k, v in sorted(by.items())}, "moments_check": mcheck,
            "evaluate_s": {"median": float(np.median([r["evaluate_s"] for r in ev])), "max": float(np.max([r["evaluate_s"] for r in ev]))} if ev else None,
            "scan_wall_s": time.time() - t_start,
            "note": "the blocked-fraction thresholds are INFORMATION for a future, separately frozen policy; no event is reclassified"}
    (a.out / "scan.json").write_text(json.dumps(summ, indent=1, default=str) + "\n")
    L = ["# Fast-only K2A operational domain scan (frozen kernel, production grouping; no S0)", "",
         f"PF `{a.pf}`; kernel `{a.kernel}` (id verified on the first and last event); {len(rows)} events scanned; skipped: "
         + (", ".join(f"{k} {v}" for k, v in skipped.items()) or "none") + "; sealed slots never read.", "",
         "Execution: " + ", ".join(f"{k} {v}" for k, v in sorted(ex.items())) + f" (in domain {100 * (summ['in_domain_frac'] or 0):.3g} % of evaluated)", ""]
    if summ["blocked_frac_quantiles"]:
        L += ["Blocked charge fraction of BLOCKED events: " + ", ".join(f"q{k:g} {v:.3g}" for k, v in summ["blocked_frac_quantiles"].items()),
              "BLOCKED events with blocked fraction below (information only): " + ", ".join(f"< {k}: {v}" for k, v in summ["blocked_events_below"].items()), ""]
    L += ["Blocked clouds by cause (clouds / charge fC): " + (", ".join(f"{k} {v} / {cause_q[k]:.4g}" for k, v in cause_ev.most_common()) or "none")]
    if summ["r_eff_above_grid_ratio_quantiles"]:
        L.append("R_EFF_ABOVE_GRID clouds, r_eff / r_max: " + ", ".join(f"q{k:g} {v:.3g}" for k, v in summ["r_eff_above_grid_ratio_quantiles"].items()))
    L += ["", "Agreement with the recentering recovery (status_after): " + ", ".join(f"{k}: {v}" for k, v in summ["agreement_with_recovery"].items()),
          "", "By channel / station: " + "; ".join(f"{k} " + ", ".join(f"{a_}:{b}" for a_, b in v.items()) for k, v in summ["by_channel_station"].items()),
          "", "Moments cross-check (preflight == validation archive moments): " + (", ".join(f"{k} {v}" for k, v in mcheck.items()) or "not requested"),
          "", f"Evaluation time per event (cold, kernel mapped): median {summ['evaluate_s']['median'] * 1e3:.3g} ms" if summ["evaluate_s"] else "",
          "", "Blocked events (all):", "", "| slot | channel | station | status | blocked frac | blocked clouds | causes | recovery |", "|---|---|---|---|---|---|---|---|"]
    for r in sorted((r for r in ev if r["execution"] == "BLOCKED"), key=lambda r: -(r["blocked_frac"] or 0)):
        L.append(f"| {r['slot']} | {r['channel']} | {r['station']} | {r['status_names']} | {r['blocked_frac']:.3g} | {r['blocked_clouds']} | "
                 f"{r['blocked_causes']} | {r['recovery_status_after']} |")
    errs = [r for r in rows if r["execution"] == "ERROR"]
    if errs:
        L += ["", "ERRORS:"] + [f"- {r['slot']}: {r['status_names']}" for r in errs]
    (a.out / "scan.md").write_text("\n".join(L) + "\n")
    print(f"domain scan: {len(rows)} events, " + ", ".join(f"{k} {v}" for k, v in sorted(ex.items())) + f"; {a.out / 'scan.md'}")
    return 1 if errs or any(v == "FAIL" for v in mcheck.values()) else 0


if __name__ == "__main__":
    sys.exit(main())

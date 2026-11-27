#!/usr/bin/env python3
"""Read-only diagnostics of finished run_k2a_validation.sh runs: where the K2A-vs-S0 residual lives, what S0 carries
beyond the K2A window, what was blocked, and how the runtime scales. No model, kernel, threshold or reference changes.

  k2a_residual_diagnostics.py NEW_DIR --from RUN_DIR [RUN_DIR ...] [--plots] [--native NAME=MUON:DIR ...]

For every muon of the runs (variants S0T1 = far E/S0_local, S0T0 = V/T0/S0_local, K2A = V/K2A_prod (own grouping),
K2Aimp = V/K2A_imposed; all on the first nb bins of the K2A window, after the same exact time-grid check as the run):

1. RESIDUAL ANATOMY of each component (model S0T1->K2A, total S0T0->K2A, timing S0T0->S0T1, imposed S0T1->K2Aimp),
   ion-only and ion + electron step, per plane, R = cand - ref on common strips (unshaped, fC/ns per 10 ns bin):
     L2 = ||R|| / ||ref||; share of ||R||^2 per time window (0-50-100-200-500-1000-2000-4000-7000 ns) and the window's
     own relative L2; share per strip rank (by reference charge: lead, 2nd, 3rd-5th, rest); global scale removed
     (alpha = <c,r>/<r,r>: L2 of c - alpha r) and per-strip scales removed (pure shape error); best common time shift
     (first order, delta = <R, dr/dt>/<dr/dt, dr/dt>) and the L2 left after removing it; high-frequency share of ||R||^2
     (R minus its 50 ns running mean) and the share that survives the provisional Tp = 100 ns CR-RC^3 shaper.
   S0 NUMERICAL FLOOR: the frozen replay's own time-compression pilot of the chosen late step per group
     (timing_attempts.json: relative L2 of 32 weighted rows vs the native 10 ns clock, per window), max over groups.
   K2A kernel tau grid: relative step at 10, 100, 1000 and 5000 ns.
   --native NAME=MUON:DIR adds an S0 variant (e.g. a native 10 ns S0_T1 built by k2a_s0_native_step.sh) and the
   components S0native->K2A and S0native->S0T1 (the latter = the S0 late-step compression error itself).
2. TRUNCATION of S0 beyond the K2A window (both S0 references, full frozen clock): signed and absolute ion charge
   beyond W per plane and on the lead strip, as a fraction of the full-clock charge; shaped effect (Tp 100 / 200 ns):
   the largest |shaped(full) - shaped(truncated at W)| over all strips and times, relative to the lead shaped peak
   (zero before W by causality; reported for t >= W).
3. BLOCKED CLOUDS (any muon whose primary K2A has a blocking status): primary, wire, charge and its fraction, t_ns,
   r_eff vs the kernel radial grid [r_min, r_max] (ratio to r_max), local D vs the D table, wire vs the kernel wires,
   min_clearance from the moments, other clouds of the same primary.
4. RUNTIME SCALING: per muon K256 rows, clouds, groups vs measured S0 and K2A times; log-log slopes when n >= 4.
Writes NEW_DIR/residual_diagnostics.json, residual_diagnostics.md and (--plots) PNG figures.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True                # never write __pycache__ next to the frozen timing-ablation tools

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "timing-ablation"))
import k2a_time_grid as tg          # noqa: E402
import t_shaped_compare as tsc      # noqa: E402  (frozen; only its shaper functions are used)

WIN = [0, 50, 100, 200, 500, 1000, 2000, 4000, 7000]
COMP = [("model", "S0T1", "K2A"), ("total", "S0T0", "K2A"), ("timing", "S0T0", "S0T1"), ("imposed", "S0T1", "K2Aimp")]


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def load(d):
    lab, t, ion = tsc.load_wave(Path(d) / "replay" / "waveforms.csv")
    lab2, t2, ele = tsc.load_wave(Path(d) / "prompt" / "prompt_waveforms.csv")
    if lab2 != lab or not np.array_equal(t, t2):
        raise ValueError(f"{d}: replay and prompt grids differ")
    return lab, t, ion, ele


def plane(lab, s):
    return [i for i, l in enumerate(lab) if l.startswith(s + "_")]


def shaper(tp, dt=1.0):
    h, _ = tsc.impulse(tp, 3, dt, max(12 * tp, 2000.0), None)
    return h


def shaped(cur, dt_in, h):
    return tsc.shape(cur, dt_in, h, 1.0)


def anatomy(ref, cand, t, dt, h100):
    """ref, cand: (nb, ns) on common strips of one plane."""
    R = cand - ref
    nr = float(np.sum(ref ** 2))
    if nr <= 0:
        return None
    e = float(np.sum(R ** 2)); L2 = math.sqrt(e / nr)
    o = {"L2_rel": L2, "strips": int(ref.shape[1])}
    o["window_share"], o["window_L2_rel"] = {}, {}
    for a, b in zip(WIN[:-1], WIN[1:]):
        m = (t >= a) & (t < b)
        ew = float(np.sum(R[m] ** 2)); nw = float(np.sum(ref[m] ** 2))
        o["window_share"][f"{a}-{b}"] = ew / e if e > 0 else 0.0
        o["window_L2_rel"][f"{a}-{b}"] = math.sqrt(ew / nw) if nw > 0 else None
    q = ref.sum(0); rank = np.argsort(-np.abs(q))
    cls = {"lead": rank[:1], "2nd": rank[1:2], "3rd-5th": rank[2:5], "rest": rank[5:]}
    o["strip_rank_share"] = {k: (float(np.sum(R[:, v] ** 2)) / e if e > 0 and len(v) else 0.0) for k, v in cls.items()}
    lead = rank[0]
    o["lead_strip_L2_rel"] = float(np.linalg.norm(R[:, lead]) / np.linalg.norm(ref[:, lead])) if np.any(ref[:, lead]) else None
    alpha = float(np.sum(cand * ref) / nr)
    o["global_scale"] = alpha
    o["L2_after_global_scale"] = float(np.linalg.norm(cand - alpha * ref) / math.sqrt(nr))
    ps = np.array([np.sum(cand[:, j] * ref[:, j]) / np.sum(ref[:, j] ** 2) if np.any(ref[:, j]) else 0.0 for j in range(ref.shape[1])])
    o["L2_after_per_strip_scale"] = float(np.linalg.norm(cand - ref * ps[None, :]) / math.sqrt(nr))
    o["lead_strip_scale"] = float(ps[lead])
    dr = np.gradient(ref, dt, axis=0)
    dd = float(np.sum(dr * dr))
    delta = float(np.sum(R * dr) / dd) if dd > 0 else 0.0
    o["best_time_shift_ns"] = delta
    o["L2_after_time_shift"] = float(np.linalg.norm(R - delta * dr) / math.sqrt(nr))
    k = max(1, int(round(50.0 / dt)))
    ker = np.ones(k) / k
    sm = np.apply_along_axis(lambda x: np.convolve(x, ker, mode="same"), 0, R)
    o["high_freq_share_50ns"] = float(np.sum((R - sm) ** 2)) / e if e > 0 else 0.0
    sR, sS = shaped(R, dt, h100), shaped(ref, dt, h100)
    o["shaped_Tp100_L2_rel"] = float(np.linalg.norm(sR) / np.linalg.norm(sS))
    o["shaped_over_unshaped_L2"] = o["shaped_Tp100_L2_rel"] / L2 if L2 > 0 else None
    return o


def s0_floor(sg):
    out = {"groups": 0, "late_step_ns": {}, "max_pilot_L2": {}}
    for f in sorted(Path(sg).glob("g_*/timing_attempts.json")) if Path(sg).is_dir() else []:
        j = jload(f) or {}
        ch = [x for x in j.get("attempts", []) if x.get("late_step_ns") == j.get("chosen_late_step_ns")]
        if not ch:
            continue
        out["groups"] += 1; out["late_step_ns"][f.parent.name] = j.get("chosen_late_step_ns")
        for w, v in (ch[0].get("pilot_relative_L2") or {}).items():
            if v is not None:
                out["max_pilot_L2"][w] = max(out["max_pilot_L2"].get(w, 0.0), float(v))
    return out


def truncation(d, W, hs):
    lab, t, ion, ele = load(d)
    dt = t[1] - t[0]; nb = int(round(W / dt))
    o = {"bins_full": int(len(t)), "t_last_ns": float(t[-1]), "window_ns": W}
    for s in ("U", "V"):
        c = plane(lab, s); I = ion[:, c]
        qf, qw = I.sum() * dt, I[:nb].sum() * dt
        af, aw = np.abs(I).sum() * dt, np.abs(I[:nb]).sum() * dt
        q = I.sum(0); ld = int(np.argmax(np.abs(q)))
        p = {"charge_full_fC": float(qf), "charge_beyond_frac": float((qf - qw) / qf) if qf else None,
             "abs_charge_beyond_frac": float((af - aw) / af) if af else None,
             "lead_strip": lab[c[ld]], "lead_charge_beyond_frac": float((q[ld] - I[:nb, ld].sum()) / q[ld]) if q[ld] else None}
        tot = ion[:, c] + ele[:, c]
        cut = tot.copy(); cut[nb:] = 0.0
        for tp in (100.0, 200.0):
            h = hs[tp]
            yf, yc = shaped(tot, dt, h), shaped(cut, dt, h)
            lead_pk = float(np.max(np.abs(yf)))
            p[f"shaped_Tp{tp:g}_tail_effect_over_lead_peak"] = float(np.max(np.abs(yf - yc))) / lead_pk if lead_pk else None
        o[s] = p
    return o


def blocked(V, KD):
    k = jload(Path(KD) / "kernel.json") or {}
    ax = k.get("green_axes", {})
    rr, Dd, ww = ax.get("r_um") or [], ax.get("D_mm") or [], ax.get("wires") or []
    res = jload(V / "K2A_prod" / "result.json") or {}
    out = {"status_names": res.get("status_names"), "cloud_charge_fC": res.get("cloud_charge_fC"),
           "blocked_charge_fC": res.get("blocked_charge_fC"),
           "kernel_grid": {"r_um": [rr[0], rr[-1]] if rr else None, "r_nodes": len(rr), "D_mm": [Dd[0], Dd[-1]] if Dd else None,
                           "wires": ww, "tau_max_ns": (ax.get("tau_ns") or [None])[-1]}, "clouds": []}
    cf = V / "K2A_prod" / "clouds.csv"
    if not cf.is_file():
        return out
    rows = list(csv.DictReader(open(cf)))
    mom = {}
    mp = V / "prod" / "birth_moments.csv"
    if mp.is_file():
        for r in csv.DictReader(open(mp)):
            mom[(r["primary"], r["wire"])] = r
    q = res.get("cloud_charge_fC") or sum(float(r["N_fC"]) for r in rows)
    for r in rows:
        names = r["status_names"]
        blk = [x for x in names.split("|") if x and x not in ("OK", "T_AFTER_WINDOW", "SHIFT_NEAR_TIE")]
        if not blk:
            continue
        reff = float(r["r_eff_um"]); same = [x for x in rows if x["primary"] == r["primary"] and x is not r]
        m = mom.get((r["primary"], r["wire"]), {})
        out["clouds"].append({"primary": int(r["primary"]), "wire": int(r["wire"]), "status": names, "N_fC": float(r["N_fC"]),
                              "charge_frac": float(r["N_fC"]) / q if q else None, "t_ns": float(r["t_ns"]), "r_eff_um": reff,
                              "r_eff_over_grid_max": reff / rr[-1] if rr else None, "r_eff_over_grid_min": reff / rr[0] if rr else None,
                              "D_U_local_mm": float(r["D_U_local_mm"]), "D_V_local_mm": float(r["D_V_local_mm"]),
                              "wire_in_kernel": int(r["wire"]) in ww if ww else None,
                              "min_clearance_um": float(m["min_clearance_um"]) if m.get("min_clearance_um") else None,
                              "rows": int(float(m["rows"])) if m.get("rows") else None,
                              "A": float(m["A"]) if m.get("A") else None,
                              "same_primary_other_clouds": [{"wire": int(x["wire"]), "N_fC": float(x["N_fC"]), "r_eff_um": float(x["r_eff_um"]),
                                                             "status": x["status_names"]} for x in same]})
    return out


def slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if ok.sum() < 4 or np.ptp(np.log(x[ok])) == 0:
        return None
    return float(np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)[0])


def stat(xs):
    x = np.asarray([v for v in xs if isinstance(v, (int, float)) and math.isfinite(v)], float)
    if not len(x):
        return {"n": 0}
    return {"n": int(len(x)), "median": float(np.median(x)), "min": float(x.min()), "max": float(x.max()), "mean": float(x.mean())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--from", dest="src", type=Path, nargs="+", required=True)
    ap.add_argument("--plots", action="store_true"); ap.add_argument("--native", nargs="*", default=[])
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    natives = {}
    for x in a.native:
        n, rest = x.split("=", 1); m, d = rest.split(":", 1); natives.setdefault(f"{int(m):04d}", {})[n] = Path(d)
    a.out.mkdir(parents=True)
    hs = {tp: shaper(tp) for tp in (100.0, 200.0)}
    res = {"tool": "k2a_residual_diagnostics.py", "sources": [str(s) for s in a.src], "muons": {}, "windows_ns": WIN}
    for src in a.src:
        prov = dict(l.split("=", 1) for l in (src / "provenance.txt").read_text().splitlines() if "=" in l and not l.startswith("tool_sha256"))
        KD = prov.get("kernel"); label = prov.get("label")
        for V in sorted((src / "muons").glob("muon_*")):
            M = V.name[5:]; o = jload(V / "outcome.json") or {}
            E = Path((V / "far_dir.txt").read_text().strip()) if (V / "far_dir.txt").is_file() else None
            rec = {"run": str(src), "label": label, "execution": o.get("outcome"), "accuracy_v1": (o.get("accuracy") or {}).get("status"),
                   "k2a_domain": (o.get("k2a_domain") or {}).get("status")}
            var = {}
            if E and (E / "S0_local" / "replay" / "waveforms.csv").is_file():
                var["S0T1"] = E / "S0_local"
            for n, d in (("S0T0", V / "T0" / "S0_local"), ("K2A", V / "K2A_prod"), ("K2Aimp", V / "K2A_imposed")):
                if (d / "replay" / "waveforms.csv").is_file():
                    var[n] = d
            var.update(natives.get(M, {}))
            W = float((jload(V / "K2A_prod" / "result.json") or {}).get("window_ns") or 7000.0)
            comps = list(COMP) + [(f"{n}_vs_K2A", n, "K2A") for n in natives.get(M, {})] + \
                    [(f"{n}_vs_S0T1", n, "S0T1") for n in natives.get(M, {})]
            if len(var) >= 2:
                ok, grid, nb = tg.check(list(var.items()), W)
                rec["time_grid"] = {"status": grid["status"], "failures": grid["failures"], "bins_compared": nb}
                if ok:
                    data = {n: load(d) for n, d in var.items()}
                    t = data[next(iter(data))][1][:nb]; dt = float(t[1] - t[0])
                    rec["residual"] = {}
                    for cname, rn, cn in comps:
                        if rn not in data or cn not in data:
                            continue
                        lr, _, ir, er = data[rn]; lc, _, ic, ec = data[cn]
                        rec["residual"][cname] = {}
                        for sig in ("ion_only", "ion_plus_electron"):
                            rec["residual"][cname][sig] = {}
                            for s in ("U", "V"):
                                com = [l for l in lr if l.startswith(s + "_") and l in lc]
                                jr = [lr.index(l) for l in com]; jc = [lc.index(l) for l in com]
                                R0 = ir[:nb][:, jr] + (er[:nb][:, jr] if sig == "ion_plus_electron" else 0)
                                C0 = ic[:nb][:, jc] + (ec[:nb][:, jc] if sig == "ion_plus_electron" else 0)
                                rec["residual"][cname][sig][s] = anatomy(R0, C0, t, dt, hs[100.0])
                    if a.plots and "S0T1" in data and "K2A" in data:
                        plot_muon(a.out / f"muon_{M}_residual.png", M, data, nb, rec["residual"].get("model", {}).get("ion_only", {}))
            else:
                rec["time_grid"] = {"status": "NOT_RUN", "failures": [], "bins_compared": None}
            rec["s0_floor"] = {"S0T1": s0_floor(E / "S0groups") if E else None, "S0T0": s0_floor(V / "T0" / "S0groups")}
            rec["truncation"] = {n: truncation(var[n], W, hs) for n in ("S0T1", "S0T0") if n in var}
            if rec["k2a_domain"] == "OUT_OF_DOMAIN" or ((jload(V / "K2A_prod" / "result.json") or {}).get("blocked_charge_fC") or 0) > 0:
                rec["blocked"] = blocked(V, KD)
            rt = o.get("runtime") or {}; S, K = rt.get("S0", {}), rt.get("K2A", {})
            kr = jload(V / "K2A_prod" / "result.json") or {}
            rec["runtime"] = {"rows": K.get("rows"), "clouds": kr.get("clouds"), "groups": kr.get("groups"),
                              "S0_T0_s": S.get("total_T0_single_pass"), "S0_T1_s": S.get("total_T1_single_pass"),
                              "rkf_export_s": S.get("rkf_export_archive_read_plus_ion_drift"),
                              "replay_T0_s": S.get("weighting_replay_T0_chosen_attempt_sum"),
                              "K2A_eval_warm_s": K.get("response_eval_warm_median"), "K2A_moments_s": K.get("moments_construction_min"),
                              "K2A_archive_io_s": K.get("archive_io_read_validate_min"),
                              "speedup_response": (rt.get("speedup_measured") or {}).get("S0_T0/cathode_response_eval")}
            res["muons"][M] = rec
    kj = None
    for src in a.src:
        prov = dict(l.split("=", 1) for l in (src / "provenance.txt").read_text().splitlines() if "=" in l and not l.startswith("tool_sha256"))
        kj = jload(Path(prov.get("kernel", "")) / "kernel.json")
        if kj:
            break
    if kj:
        tau = np.asarray(kj["green_axes"]["tau_ns"], float)
        res["kernel_tau_rel_step"] = {f"{x:g}ns": float(np.interp(x, tau[1:-1], np.diff(tau[1:]) / tau[1:-1])) for x in (10, 100, 1000, 5000)}
    valid = [m for m, r in res["muons"].items() if r["execution"] in ("VALID", "VALID_T1_ONLY")]
    rt = [res["muons"][m]["runtime"] for m in valid]
    res["runtime_scaling"] = {"n": len(rt),
                              "S0_T0_vs_rows_loglog_slope": slope([r["rows"] for r in rt], [r["S0_T0_s"] for r in rt]),
                              "replay_T0_vs_rows_loglog_slope": slope([r["rows"] for r in rt], [r["replay_T0_s"] for r in rt]),
                              "K2A_eval_vs_clouds_loglog_slope": slope([r["clouds"] for r in rt], [r["K2A_eval_warm_s"] for r in rt]),
                              "S0_T0_per_row_s": stat([r["S0_T0_s"] / r["rows"] for r in rt if r["S0_T0_s"] and r["rows"]]),
                              "K2A_eval_per_cloud_s": stat([r["K2A_eval_warm_s"] / r["clouds"] for r in rt if r["K2A_eval_warm_s"] and r["clouds"]]),
                              "speedup_response": stat([r["speedup_response"] for r in rt])}
    (a.out / "residual_diagnostics.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    write_md(a.out, res, valid)
    if a.plots and rt:
        plot_scaling(a.out / "runtime_scaling.png", rt)
    print(f"residual diagnostics: {len(res['muons'])} muons ({len(valid)} valid executions) -> {a.out / 'residual_diagnostics.md'}")


def f(v, k="", pct=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "n/a"
    if pct:
        return f"{100 * v:.3g} %"
    return f"{v:.3g}"


def write_md(out, res, valid):
    L = ["# K2A residual, truncation, blocked-cloud and runtime diagnostics", "",
         f"Sources: {', '.join(res['sources'])}. Labels: {sorted({r['label'] for r in res['muons'].values()})}. Read-only; "
         "nothing re-run. Signal: unshaped 10 ns bins on the K2A window unless stated.", "",
         "## 1. Residual anatomy (model = S0_T1 -> K2A, ion only)", "",
         "| muon | plane | L2 | after global scale | after per-strip scale | after time shift (ns) | high-freq share (<50 ns) | shaped/unshaped | lead strip share | S0 late-step pilot L2 max |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for m in valid:
        r = res["muons"][m]; an = ((r.get("residual") or {}).get("model") or {}).get("ion_only") or {}
        fl = (r.get("s0_floor") or {}).get("S0T1") or {}
        pl = fl.get("max_pilot_L2", {})
        for s in ("U", "V"):
            x = an.get(s)
            if not x:
                continue
            fw = max([v for k, v in pl.items() if k.startswith(s)] or [None], key=lambda v: v or 0)
            L.append(f"| {m} | {s} | {f(x['L2_rel'], pct=True)} | {f(x['L2_after_global_scale'], pct=True)} | {f(x['L2_after_per_strip_scale'], pct=True)} "
                     f"| {f(x['L2_after_time_shift'], pct=True)} ({x['best_time_shift_ns']:+.3g}) | {f(x['high_freq_share_50ns'], pct=True)} "
                     f"| {f(x['shaped_over_unshaped_L2'])} | {f(x['strip_rank_share']['lead'], pct=True)} | {f(fw, pct=True)} |")
    L += ["", "Share of the squared residual per time window (model, ion only):", "",
          "| muon | plane | " + " | ".join(f"{a}-{b}" for a, b in zip(WIN[:-1], WIN[1:])) + " |", "|---" * (len(WIN) + 1) + "|"]
    for m in valid:
        an = ((res["muons"][m].get("residual") or {}).get("model") or {}).get("ion_only") or {}
        for s in ("U", "V"):
            x = an.get(s)
            if x:
                L.append(f"| {m} | {s} | " + " | ".join(f(v, pct=True) for v in x["window_share"].values()) + " |")
    L += ["", "Same residual, other components (ion only, U / V L2; native = K2A vs a native 10 ns S0, S0 step error = native "
          "vs late-step S0, if provided):", "", "| muon | " + " | ".join(c for c, *_ in COMP) + " | native | S0 step error |", "|---" * (len(COMP) + 3) + "|"]
    for m in valid:
        rr = res["muons"][m].get("residual") or {}
        cells = []
        for c, *_ in COMP + [("native", 0, 0), ("s0step", 0, 0)]:
            base = [x for x, *_ in COMP]
            keys = [c] if c not in ("native", "s0step") else [k for k in rr if k not in base and k.endswith("_vs_K2A" if c == "native" else "_vs_S0T1")]
            v = [(rr.get(k) or {}).get("ion_only") or {} for k in keys]
            cells.append(" ; ".join(f"{f((x.get('U') or {}).get('L2_rel'), pct=True)} / {f((x.get('V') or {}).get('L2_rel'), pct=True)}" for x in v) or "n/a")
        L.append(f"| {m} | " + " | ".join(cells) + " |")
    if res.get("kernel_tau_rel_step"):
        L += ["", "K2A kernel tau grid relative step: " + ", ".join(f"{k} {100 * v:.3g} %" for k, v in res["kernel_tau_rel_step"].items())]
    L += ["", "## 2. S0 beyond the K2A window", "",
          "| muon | ref | plane | charge beyond W | |charge| beyond W | lead strip beyond W | shaped tail effect Tp100 | Tp200 |", "|---|---|---|---|---|---|---|---|"]
    for m, r in res["muons"].items():
        for n, tr in (r.get("truncation") or {}).items():
            for s in ("U", "V"):
                p = tr[s]
                L.append(f"| {m} | {n} | {s} | {f(p['charge_beyond_frac'], pct=True)} | {f(p['abs_charge_beyond_frac'], pct=True)} | "
                         f"{f(p['lead_charge_beyond_frac'], pct=True)} | {f(p['shaped_Tp100_tail_effect_over_lead_peak'], pct=True)} | "
                         f"{f(p['shaped_Tp200_tail_effect_over_lead_peak'], pct=True)} |")
    L += ["", "## 3. Blocked clouds of the primary K2A", ""]
    nb = 0
    for m, r in res["muons"].items():
        b = r.get("blocked")
        if not b:
            continue
        g = b["kernel_grid"]
        L.append(f"- muon {m} ({r['execution']}): {b['status_names']}; blocked {b['blocked_charge_fC']} of {b['cloud_charge_fC']} fC; "
                 f"kernel r grid {g['r_um']} um ({g['r_nodes']} nodes), D table {g['D_mm']} mm, wires {g['wires']}")
        for c in b["clouds"]:
            nb += 1
            L.append(f"  - primary {c['primary']} wire {c['wire']}: {c['status']}, {c['N_fC']:.4g} fC ({f(c['charge_frac'])} of the event), "
                     f"t {c['t_ns']:.4g} ns, r_eff {c['r_eff_um']:.4g} um = {f(c['r_eff_over_grid_max'])} x r_max, D_local "
                     f"{c['D_U_local_mm']:.3g}/{c['D_V_local_mm']:.3g} mm, clearance {c['min_clearance_um']} um, K256 rows {c['rows']}, "
                     f"other clouds of this primary: {len(c['same_primary_other_clouds'])}")
    if not nb:
        L.append("none")
    rs = res["runtime_scaling"]
    L += ["", f"## 4. Runtime scaling (valid executions, n = {rs['n']}; measured)", "",
          "| muon | K256 rows | clouds | groups | S0_T0 s | RKF export s | replay s | K2A eval warm s | moments s | archive I/O s | response speedup |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in valid:
        x = res["muons"][m]["runtime"]
        L.append(f"| {m} | {x['rows']} | {x['clouds']} | {x['groups']} | {f(x['S0_T0_s'])} | {f(x['rkf_export_s'])} | {f(x['replay_T0_s'])} | "
                 f"{f(x['K2A_eval_warm_s'])} | {f(x['K2A_moments_s'])} | {f(x['K2A_archive_io_s'])} | {f(x['speedup_response'])} |")
    L += ["", f"log-log slopes: S0_T0 vs rows {f(rs['S0_T0_vs_rows_loglog_slope'])}, replay vs rows {f(rs['replay_T0_vs_rows_loglog_slope'])}, "
          f"K2A eval vs clouds {f(rs['K2A_eval_vs_clouds_loglog_slope'])} (n < 4 -> n/a). Avalanche generation excluded on both sides."]
    (out / "residual_diagnostics.md").write_text("\n".join(L) + "\n")


def plot_muon(path, M, data, nb, an):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ls, ts, i_s, _ = data["S0T1"]; lk, _, i_k, _ = data["K2A"]
    t = ts[:nb]
    fig, ax = plt.subplots(2, 2, figsize=(11, 6.5))
    for k, s in enumerate(("U", "V")):
        c = [j for j, l in enumerate(ls) if l.startswith(s + "_") and l in lk]
        q = i_s[:nb, c].sum(0); ld = c[int(np.argmax(np.abs(q)))]; lab = ls[ld]; jk = lk.index(lab)
        y0, y1 = i_s[:nb, ld], i_k[:nb, jk]
        a0 = ax[0, k]
        a0.plot(t, y0, color="#2b2b29", lw=1.6, label="S0_T1 (reference)")
        a0.plot(t, y1, color="#2a78d6", lw=1.1, ls="--", label="K2A (own grouping)")
        a0.plot(t, 10 * (y1 - y0), color="#eb6834", lw=1.0, label="10 x (K2A - S0_T1)")
        a0.set(xscale="symlog", title=f"muon {M} {s}: lead strip {lab}, ion current", xlabel="time [ns]", ylabel="fC / ns")
        a0.legend(frameon=False, fontsize=7)
        x = an.get(s) or {}
        if x:
            ws = list(x["window_share"].items())
            ax[1, k].bar(range(len(ws)), [100 * v for _, v in ws], color="#2a78d6", width=0.6)
            ax[1, k].set_xticks(range(len(ws)), [w for w, _ in ws], rotation=40, fontsize=7)
            ax[1, k].set(title=f"{s}: share of squared residual per window (L2 {100 * x['L2_rel']:.3g} %)", ylabel="% of ||R||^2")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def plot_scaling(path, rt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    x = [r["rows"] for r in rt]; y = [r["S0_T0_s"] for r in rt]
    ax[0].scatter(x, y, color="#2b2b29", s=30); ax[0].set(xscale="log", yscale="log", xlabel="K256 rows", ylabel="S0_T0 response time [s]",
                                                        title="S0 (RKF drift + weighting replay + electron step)")
    x = [r["clouds"] for r in rt]; y = [r["K2A_eval_warm_s"] for r in rt]
    ax[1].scatter(x, y, color="#2a78d6", s=30); ax[1].set(xscale="log", yscale="log", xlabel="clouds", ylabel="K2A evaluation, warm [s]",
                                                        title="K2A kernel evaluation (moments given)")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()

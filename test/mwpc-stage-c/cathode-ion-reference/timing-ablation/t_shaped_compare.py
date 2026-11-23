#!/usr/bin/env python3
"""Compare shaped signed U/V strip signals between T-ladder variants and T0.

Input layout (written by run_t_ladder.sh):
  ROOT/<muon>/<variant>/replay/waveforms.csv          frozen replay, ion term (fC/ns per 10 ns bin)
  ROOT/<muon>/<variant>/prompt/prompt_waveforms.csv   avalanche-electron term (same grid and channels)
Variant "T0" is the identity export (byte-identical to the frozen export) and is the reference.

Shaper: by default a unit-peak CR-RC^n stand-in, h(t) ~ (t/tau)^n exp(-t/tau), tau = Tp/n (n = 3).
THIS IS NOT THE VMM3a IMPULSE RESPONSE (the VMM3a user guide gives a 3rd-order 1 real + 2 complex-pole shaper without
pole values). Supply --impulse FILE (two columns: t_ns, h; any normalization; renormalized to unit peak) when a
measured or vendor response is available; the output is labelled with the shaper used.

Per muon, plane and strip: main-polarity peak (sign of the plane's largest |extreme|), peak time (parabolic
refinement on a 1 ns grid), opposite-polarity excursion, threshold-crossing times at fractions of the plane's leading
peak and at absolute thresholds in fC-equivalent. Per plane: leading strip, second/leading ratio, multiplicity above
each threshold, peak-weighted centroid of strips above the lowest threshold (strip units), and w, u from U/V
(U = w - A u, V = w + A u on the exporter's unnormalized strip coordinate, pitch 1.7 mm, A = 0.1).
Paired differences vs T0, with bootstrap over muons (percentile 68% and 95% intervals).
"""
import argparse
import json
from pathlib import Path

import numpy as np

PITCH_STRIP = 1.7
A_STEREO = 0.1


def load_wave(path):
    with open(path) as f:
        head = f.readline().strip().split(",")
    d = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return head[1:], d[:, 0], d[:, 1:]


def impulse(tp, n, dt, tmax, table=None):
    t = np.arange(0.0, tmax, dt)
    if table is None:
        tau = tp / n
        h = (t / tau) ** n * np.exp(-t / tau)
        label = f"CR-RC^{n} stand-in, Tp={tp:g} ns (NOT VMM3a poles)"
    else:
        tt, hh = np.loadtxt(table, delimiter=None, unpack=True)
        h = np.interp(t, tt, hh, left=0.0, right=0.0)
        label = f"impulse table {table}"
    return h / h.max(), label


def shape(current, dt_in, h, dt):
    """current: (nbins, nch) fC/ns on dt_in bins (piecewise constant). Returns shaped (unit-peak per fC) on dt grid."""
    k = int(round(dt_in / dt))
    up = np.repeat(current, k, axis=0) * dt          # charge per fine sample, fC
    n = up.shape[0] + len(h)
    nf = 1 << int(np.ceil(np.log2(n)))
    H = np.fft.rfft(h, nf)
    out = np.fft.irfft(np.fft.rfft(up, nf, axis=0) * H[:, None], nf, axis=0)[:up.shape[0]]
    return out


def parabolic(y, i):
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        if den != 0:
            x = 0.5 * (a - c) / den
            return i + x, b - 0.25 * (a - c) * x
    return float(i), y[i]


def first_cross(y, thr, dt):
    idx = np.flatnonzero(y >= thr)
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return 0.0
    return (i - 1 + (thr - y[i - 1]) / (y[i] - y[i - 1])) * dt


def observables(labels, y, dt, fracs, abs_thr):
    res = {}
    planes = {}
    for side in ("U", "V"):
        cols = [j for j, l in enumerate(labels) if l.startswith(side + "_")]
        idx = np.array([int(labels[j].split("_")[1]) for j in cols])
        Y = y[:, cols]
        ext = np.where(np.abs(Y.max(0)) >= np.abs(Y.min(0)), Y.max(0), Y.min(0))
        lead = int(np.argmax(np.abs(ext)))
        sgn = 1.0 if ext[lead] > 0 else -1.0
        S = sgn * Y
        pk = np.empty(len(cols)); tpk = np.empty(len(cols)); opp = S.min(0)
        for c in range(len(cols)):
            i = int(np.argmax(S[:, c]))
            x, v = parabolic(S[:, c], i)
            pk[c], tpk[c] = v, x * dt
        lead_pk = pk[lead]
        thr_t = {}
        mult = {}
        for f in fracs:
            thr = f * lead_pk
            thr_t[f"frac{f:g}"] = [first_cross(S[:, c], thr, dt) for c in range(len(cols))]
            mult[f"frac{f:g}"] = int(np.sum(pk >= thr))
        for q in abs_thr:
            thr_t[f"abs{q:g}fC"] = [first_cross(S[:, c], q, dt) for c in range(len(cols))]
            mult[f"abs{q:g}fC"] = int(np.sum(pk >= q))
        sel = pk >= min(fracs) * lead_pk
        cen = float(np.sum(pk[sel] * idx[sel]) / np.sum(pk[sel]))
        sec = np.sort(pk)[-2] / lead_pk if len(pk) > 1 else np.nan
        planes[side] = dict(sign=sgn, lead_index=int(idx[lead]), lead_peak=float(lead_pk),
                            lead_peak_time=float(tpk[lead]), second_over_lead=float(sec), centroid_strip=cen,
                            mult=mult, thr_time_lead={k: v[lead] for k, v in thr_t.items()},
                            peaks=pk.tolist(), peak_times=tpk.tolist(), opposite=(-opp).tolist() if sgn < 0 else opp.tolist(),
                            strip_index=idx.tolist(), thr_times=thr_t)
    Umm = planes["U"]["centroid_strip"] * PITCH_STRIP
    Vmm = planes["V"]["centroid_strip"] * PITCH_STRIP
    res["planes"] = planes
    res["w_mm"] = 0.5 * (Umm + Vmm)
    res["u_mm"] = (Vmm - Umm) / (2 * A_STEREO)
    return res


def scalar_rows(o):
    r = {"w_mm": o["w_mm"], "u_mm": o["u_mm"]}
    for s, p in o["planes"].items():
        r[f"{s}_lead_peak"] = p["lead_peak"]
        r[f"{s}_lead_peak_time_ns"] = p["lead_peak_time"]
        r[f"{s}_second_over_lead"] = p["second_over_lead"]
        r[f"{s}_centroid_strip"] = p["centroid_strip"]
        for k, v in p["thr_time_lead"].items():
            r[f"{s}_t_{k}_ns"] = v
        for k, v in p["mult"].items():
            r[f"{s}_mult_{k}"] = v
        r[f"{s}_sign"] = p["sign"]
    return r


def boot(x, nb, rng):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return dict(n=0)
    m = rng.integers(0, len(x), (nb, len(x)))
    means = x[m].mean(1)
    rms = np.sqrt((x[m] ** 2).mean(1))
    q = lambda a, p: float(np.percentile(a, p))
    return dict(n=int(len(x)), mean=float(x.mean()), mean_ci68=[q(means, 16), q(means, 84)],
                mean_ci95=[q(means, 2.5), q(means, 97.5)], rms=float(np.sqrt((x ** 2).mean())),
                rms_ci95=[q(rms, 2.5), q(rms, 97.5)], max_abs=float(np.max(np.abs(x))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--variants", default="T0,T1,TW,TM")
    ap.add_argument("--tp", default="100,200")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--impulse", type=Path, default=None)
    ap.add_argument("--fracs", default="0.02,0.05,0.1,0.2")
    ap.add_argument("--abs-thr-fC", default="", help="absolute thresholds on the unit-peak (fC-equivalent) output")
    ap.add_argument("--window-ns", type=float, default=6000.0)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--no-prompt", action="store_true", help="ion term only (diagnostic)")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    variants = a.variants.split(",")
    fracs = [float(x) for x in a.fracs.split(",")]
    abs_thr = [float(x) for x in a.abs_thr_fC.split(",") if x]
    muons = sorted(p.name for p in a.root.iterdir() if p.is_dir() and (p / variants[0] / "replay" / "waveforms.csv").exists())
    if not muons:
        raise SystemExit("no muons found")
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    a.out.mkdir(parents=True)
    rng = np.random.default_rng(12345)
    report = {"muons": muons, "variants": variants, "prompt_included": not a.no_prompt, "shapers": {}}
    lines = []
    for tp in [float(x) for x in a.tp.split(",")]:
        h, hlabel = impulse(tp, a.order, a.dt, max(12 * tp, 2000.0), a.impulse)
        report["shapers"][f"{tp:g}"] = hlabel
        per = {v: [] for v in variants}
        for m in muons:
            for v in variants:
                labels, tc, cur = load_wave(a.root / m / v / "replay" / "waveforms.csv")
                if not a.no_prompt:
                    l2, tc2, pc = load_wave(a.root / m / v / "prompt" / "prompt_waveforms.csv")
                    if l2 != labels or not np.array_equal(tc, tc2):
                        raise SystemExit(f"{m}/{v}: prompt and replay grids differ")
                    cur = cur + pc
                dt_in = tc[1] - tc[0]
                nb = int(round(a.window_ns / dt_in))
                y = shape(cur[:nb], dt_in, h, a.dt)
                o = observables(labels, y, a.dt, fracs, abs_thr)
                o["muon"] = m
                per[v].append(o)
        ref = per[variants[0]]
        with open(a.out / f"per_muon_Tp{tp:g}.jsonl", "w") as f:
            for v in variants:
                for o in per[v]:
                    f.write(json.dumps({"variant": v, "muon": o["muon"], **scalar_rows(o),
                                        "U_peaks": o["planes"]["U"]["peaks"], "V_peaks": o["planes"]["V"]["peaks"]}) + "\n")
        summ = {}
        for v in variants[1:]:
            d = {}
            keys = scalar_rows(ref[0]).keys()
            for k in keys:
                if k.endswith("_sign"):
                    d[k + "_flips"] = int(sum(scalar_rows(x)[k] != scalar_rows(y_)[k] for x, y_ in zip(per[v], ref)))
                    continue
                diff = [scalar_rows(x)[k] - scalar_rows(y_)[k] for x, y_ in zip(per[v], ref)]
                if k.endswith("_lead_peak"):
                    diff = [dd / scalar_rows(y_)[k] for dd, y_ in zip(diff, ref)]
                    k = k + "_rel"
                d[k] = boot(diff, a.bootstrap, rng)
            # all-strip relative peak differences (strips with |T0 peak| >= 5% of lead)
            for s in ("U", "V"):
                rel = []
                for x, y_ in zip(per[v], ref):
                    px, py = np.array(x["planes"][s]["peaks"]), np.array(y_["planes"][s]["peaks"])
                    sel = py >= 0.05 * py.max()
                    rel += ((px[sel] - py[sel]) / py.max()).tolist()
                d[f"{s}_strip_peak_minus_T0_over_T0_lead"] = boot(rel, a.bootstrap, rng)
            summ[v] = d
        report[f"Tp{tp:g}"] = summ
        lines.append(f"== shaper: {hlabel}; muons={len(muons)}; prompt term {'ON' if not a.no_prompt else 'OFF'}")
        for v, d in summ.items():
            lines.append(f"-- {v} minus {variants[0]}")
            for k, s in d.items():
                if isinstance(s, dict) and s.get("n"):
                    lines.append(f"   {k:34s} n={s['n']:3d} mean={s['mean']:+.4g} CI95=[{s['mean_ci95'][0]:+.3g},{s['mean_ci95'][1]:+.3g}]"
                                 f" rms={s['rms']:.4g} max|d|={s['max_abs']:.4g}")
                elif isinstance(s, int):
                    lines.append(f"   {k:34s} {s}")
    (a.out / "t_ladder_summary.json").write_text(json.dumps(report, indent=1) + "\n")
    (a.out / "t_ladder_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

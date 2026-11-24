#!/usr/bin/env python3
"""Compare shaped, signed U/V strip signals of T-ladder variants with T0 (same frozen paths, LUT and replay).

Input layout (written by run_t_ladder.sh):
  ROOT/<muon>/<variant>/replay/waveforms.csv          frozen replay: ion term, fC/ns per 10 ns bin, signed
  ROOT/<muon>/<variant>/prompt/prompt_waveforms.csv   APPROXIMATE electron term: instantaneous step at birth time
Run ion-only (--no-prompt) first. The ion term is the time-dependent Shockley-Ramo signal of the frozen RKF paths; the
electron step is only an approximation of the time-dependent electron signal (see t_prompt_electron_term.py).
The FIRST variant is the reference (T0 for the T-ladder: identity export; S0 for the S-ladder: T1 frozen replay).

Shaping times 100/200 ns are PROVISIONAL until the NA60+ VMM3a configuration is verified. Default shaper: unit-peak
CR-RC^n stand-in, h(t) ~ (t/tau)^n exp(-t/tau), tau = Tp/n, n = 3 -- NOT the VMM3a transfer function. --impulse FILE
(t_ns, h) replaces it. Output unit: "fC-equivalent peak" (a 1 fC delta charge gives peak 1).

Thresholds are ABSOLUTE and fixed from T0: thr_f = f x (T0 leading-strip peak of that plane), and the SAME value is
applied to every variant, so discrete hit changes near threshold are counted, not normalized away.

Per muon, Tp and plane:
  * normalized L2 of the signed shaped waveforms over all strips: ||y_v - y_0|| / ||y_0|| (window 0..--window-ns)
  * leading strip: index change, peak (relative change), peak time, leading-edge times at each threshold
  * all strips with T0 peak >= 5% of the T0 lead: peak changes relative to the T0 lead peak
  * strip multiplicity above each threshold and the number of strips whose above/below status flips
  * plane amplitude (sum of main-polarity peaks over strips above the lowest threshold), U/V amplitude ratio,
    unshaped plane charge in the window (U/V charge)
  * unshaped normalized L2 of the current and of the cumulative charge (numerical closure without the shaper)
  * peak-weighted centroid (strip units) -> stereo w, u; also alternative centroids for diagnostics: membership fixed
    from the reference, lead +-1, lead +-2, all positive strips (u_mm_<def>, w_mm_<def>) with U = w - A u, V = w + A u on the exporter's strip
    coordinate (pitch 1.7 mm, A = 0.1)
Differences vs T0 per muon (jsonl), and over muons: mean, rms, max |d| and a muon bootstrap when >= 3 muons.
--plots writes per-muon PNGs: shaped waveforms of the 3 largest strips per plane, variants over T0, and residuals.
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
        label = f"CR-RC^{n} stand-in, Tp={tp:g} ns (PROVISIONAL; not VMM3a poles)"
    else:
        tt, hh = np.loadtxt(table, unpack=True)
        h = np.interp(t, tt, hh, left=0.0, right=0.0)
        label = f"impulse table {table}, Tp label {tp:g} ns"
    return h / h.max(), label


def shape(current, dt_in, h, dt):
    """current (nbins, nch) in fC/ns, piecewise constant on dt_in bins -> shaped output on a dt grid."""
    k = int(round(dt_in / dt))
    up = np.repeat(current, k, axis=0) * dt
    n = up.shape[0] + len(h)
    nf = 1 << int(np.ceil(np.log2(n)))
    out = np.fft.irfft(np.fft.rfft(up, nf, axis=0) * np.fft.rfft(h, nf)[:, None], nf, axis=0)[:up.shape[0]]
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
    return 0.0 if i == 0 else (i - 1 + (thr - y[i - 1]) / (y[i] - y[i - 1])) * dt


def plane_cols(labels, side):
    cols = [j for j, l in enumerate(labels) if l.startswith(side + "_")]
    return cols, np.array([int(labels[j].split("_")[1]) for j in cols])


def plane_obs(Y, idx, dt, sgn, thr, charge):
    """Y: shaped (nt, nstrip) of one plane; sgn: main polarity (fixed from T0); thr: dict name -> absolute value."""
    S = sgn * Y
    pk = np.empty(S.shape[1]); tpk = np.empty(S.shape[1])
    for c in range(S.shape[1]):
        x, v = parabolic(S[:, c], int(np.argmax(S[:, c])))
        pk[c], tpk[c] = v, x * dt
    lead = int(np.argmax(pk))
    o = dict(lead_index=int(idx[lead]), lead_peak=float(pk[lead]), lead_peak_time=float(tpk[lead]),
             peaks=pk, peak_times=tpk, charge=float(charge))
    o["thr_time"] = {k: first_cross(S[:, lead], v, dt) for k, v in thr.items()}
    o["above"] = {k: pk >= v for k, v in thr.items()}
    o["mult"] = {k: int(np.sum(a)) for k, a in o["above"].items()}
    sel = pk >= min(thr.values())
    o["amplitude"] = float(np.sum(pk[sel]))
    o["centroid_strip"] = float(np.sum(pk[sel] * idx[sel]) / np.sum(pk[sel])) if sel.any() else np.nan
    sk = np.sort(pk)
    o["second_over_lead"] = float(sk[-2] / sk[-1]) if len(sk) > 1 else np.nan
    return o


def analyse_muon(waves, labels, dt, fracs, h, window_ns, dt_in):
    """waves: dict variant -> (nbins, nch) current; first key is T0. Returns per-variant observables and shaped arrays."""
    nb = int(round(window_ns / dt_in))
    shaped = {v: shape(c[:nb], dt_in, h, dt) for v, c in waves.items()}
    charge = {v: c[:nb].sum(0) * dt_in for v, c in waves.items()}
    ref = next(iter(waves))
    res = {v: {} for v in waves}
    for side in ("U", "V"):
        cols, idx = plane_cols(labels, side)
        Y0 = shaped[ref][:, cols]
        ext = np.where(np.abs(Y0.max(0)) >= np.abs(Y0.min(0)), Y0.max(0), Y0.min(0))
        sgn = 1.0 if ext[int(np.argmax(np.abs(ext)))] > 0 else -1.0
        lead0 = float(np.max(sgn * Y0))
        thr = {f"f{f:g}": f * lead0 for f in fracs}
        I0 = waves[ref][:nb][:, cols]
        Q0 = np.cumsum(I0, 0)
        ref_obs = None
        for v in waves:
            Y = shaped[v][:, cols]
            o = plane_obs(Y, idx, dt, sgn, thr, charge[v][cols].sum())
            o["sign"] = sgn
            o["thresholds_fCeq"] = thr
            o["l2_rel_vs_T0"] = float(np.linalg.norm(Y - Y0) / np.linalg.norm(Y0))
            I = waves[v][:nb][:, cols]
            o["l2_unshaped_I"] = float(np.linalg.norm(I - I0) / np.linalg.norm(I0))
            o["l2_unshaped_Q"] = float(np.linalg.norm(np.cumsum(I, 0) - Q0) / np.linalg.norm(Q0))
            if ref_obs is None:
                ref_obs = o
            # alternative centroid definitions (strip units); 'thr' (default) = strips above the lowest threshold
            pk, pk0 = o["peaks"], ref_obs["peaks"]
            memb = pk0 >= min(thr.values())                      # membership fixed from the reference
            l0 = int(np.argmax(pk0))
            def cen(sel):
                return float(np.sum(pk[sel] * idx[sel]) / np.sum(pk[sel])) if np.sum(pk[sel]) > 0 else np.nan
            near = lambda h: (np.abs(np.arange(len(pk)) - l0) <= h)
            o["centroids"] = {"thr": o["centroid_strip"], "fixed_ref_membership": cen(memb),
                              "lead3": cen(near(1) & (pk > 0)), "lead5": cen(near(2) & (pk > 0)), "all_positive": cen(pk > 0)}
            o["membership_changes_at_inclusion"] = int(np.sum((pk >= min(thr.values())) != memb))
            res[v][side] = o
    for v in waves:
        U, V = res[v]["U"], res[v]["V"]
        res[v]["w_mm"] = 0.5 * PITCH_STRIP * (U["centroid_strip"] + V["centroid_strip"])
        res[v]["u_mm"] = PITCH_STRIP * (V["centroid_strip"] - U["centroid_strip"]) / (2 * A_STEREO)
        res[v]["u_mm_by_centroid"] = {k: PITCH_STRIP * (V["centroids"][k] - U["centroids"][k]) / (2 * A_STEREO)
                                      for k in U["centroids"]}
        res[v]["w_mm_by_centroid"] = {k: 0.5 * PITCH_STRIP * (V["centroids"][k] + U["centroids"][k]) for k in U["centroids"]}
        res[v]["UV_amp_ratio"] = U["amplitude"] / V["amplitude"] if V["amplitude"] else np.nan
        res[v]["UV_charge_ratio"] = U["charge"] / V["charge"] if V["charge"] else np.nan
    return res, shaped


def diffs(r, r0, fracs):
    d = {"w_mm": r["w_mm"] - r0["w_mm"], "u_mm": r["u_mm"] - r0["u_mm"],
         "UV_amp_ratio_rel": r["UV_amp_ratio"] / r0["UV_amp_ratio"] - 1,
         "UV_charge_ratio_rel": r["UV_charge_ratio"] / r0["UV_charge_ratio"] - 1}
    for k in r["u_mm_by_centroid"]:
        if k != "thr":
            d[f"u_mm_{k}"] = r["u_mm_by_centroid"][k] - r0["u_mm_by_centroid"][k]
            d[f"w_mm_{k}"] = r["w_mm_by_centroid"][k] - r0["w_mm_by_centroid"][k]
    for s in ("U", "V"):
        a, b = r[s], r0[s]
        d[f"{s}_L2_rel"] = a["l2_rel_vs_T0"]
        d[f"{s}_L2_unshaped_current"] = a["l2_unshaped_I"]
        d[f"{s}_L2_unshaped_charge"] = a["l2_unshaped_Q"]
        d[f"{s}_membership_changes_at_inclusion"] = a["membership_changes_at_inclusion"]
        d[f"{s}_lead_index_changed"] = int(a["lead_index"] != b["lead_index"])
        d[f"{s}_lead_peak_rel"] = a["lead_peak"] / b["lead_peak"] - 1
        d[f"{s}_lead_peak_time_ns"] = a["lead_peak_time"] - b["lead_peak_time"]
        d[f"{s}_charge_rel"] = a["charge"] / b["charge"] - 1
        d[f"{s}_amplitude_rel"] = a["amplitude"] / b["amplitude"] - 1
        d[f"{s}_second_over_lead"] = a["second_over_lead"] - b["second_over_lead"]
        d[f"{s}_centroid_strip"] = a["centroid_strip"] - b["centroid_strip"]
        big = b["peaks"] >= 0.05 * b["lead_peak"]
        d[f"{s}_strip_peak_maxabs_over_T0lead"] = float(np.max(np.abs(a["peaks"][big] - b["peaks"][big])) / b["lead_peak"])
        for f in fracs:
            k = f"f{f:g}"
            d[f"{s}_t_{k}_ns"] = a["thr_time"][k] - b["thr_time"][k]
            d[f"{s}_mult_{k}"] = a["mult"][k] - b["mult"][k]
            d[f"{s}_hitflips_{k}"] = int(np.sum(a["above"][k] != b["above"][k]))
    return d


def summarize(vals, nb, rng):
    x = np.asarray(vals, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return dict(n=0)
    s = dict(n=int(len(x)), mean=float(x.mean()), rms=float(np.sqrt((x ** 2).mean())), max_abs=float(np.abs(x).max()))
    if len(x) >= 3:
        m = rng.integers(0, len(x), (nb, len(x)))
        s["mean_ci95"] = [float(np.percentile(x[m].mean(1), 2.5)), float(np.percentile(x[m].mean(1), 97.5))]
    return s


def plot_muon(path, muon, tp, labels, shaped, res, variants, dt, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#8a5cf6", "#c99a06", "#0fa3b1", "#b5332b", "#6b7c1a", "#d65db1", "#555"]
    cols_v = {variants[0]: "#2b2b29", **{v: pal[i % len(pal)] for i, v in enumerate(variants[1:])}}
    fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    t = np.arange(shaped[variants[0]].shape[0]) * dt
    tmax = min(t[-1], 8 * tp + 400)
    for k, side in enumerate(("U", "V")):
        cols, idx = plane_cols(labels, side)
        pk0 = res[variants[0]][side]["peaks"]
        top = np.argsort(pk0)[::-1][:3]
        for j in top:
            for v in variants:
                y = shaped[v][:, cols[j]]
                ax[0, k].plot(t, y, color=cols_v.get(v, "#888"), lw=1.6 if v == variants[0] else 1.1,
                              ls="-" if v == variants[0] else "--", label=f"{v}" if j == top[0] else None)
                if v != variants[0]:
                    ax[1, k].plot(t, y - shaped[variants[0]][:, cols[j]], color=cols_v.get(v, "#888"), lw=1.0,
                                  label=f"{v}-{variants[0]}, strip {idx[j]}")
        for f, thr in res[variants[0]][side]["thresholds_fCeq"].items():
            ax[0, k].axhline(res[variants[0]][side]["sign"] * thr, color="#bbb", lw=0.6)
        ax[0, k].set(title=f"{side}: 3 largest strips (signed, fC-eq.)", xlim=(0, tmax))
        ax[1, k].set(xlabel="time [ns]", ylabel=f"variant - {variants[0]} [fC-eq.]", xlim=(0, tmax))
        ax[0, k].legend(frameon=False, fontsize=7); ax[1, k].legend(frameon=False, fontsize=6)
    fig.suptitle(f"muon {muon}, Tp={tp:g} ns: {title}", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    global PITCH_STRIP, A_STEREO                 # stereo reconstruction w, u of every variant uses this geometry
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--variants", default="T0,T1,TW")
    ap.add_argument("--tp", default="100,200")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--impulse", type=Path, default=None)
    ap.add_argument("--fracs", default="0.02,0.05,0.1,0.2,0.5")
    ap.add_argument("--window-ns", type=float, default=6000.0)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--no-prompt", action="store_true", help="ion term only")
    ap.add_argument("--plots", action="store_true")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--strip-pitch-mm", type=float, default=PITCH_STRIP, help="readout geometry of ALL variants (default frozen 1.7)")
    ap.add_argument("--stereo-tan", type=float, default=A_STEREO, help="readout geometry of ALL variants (default frozen 0.1)")
    a = ap.parse_args()
    PITCH_STRIP, A_STEREO = a.strip_pitch_mm, a.stereo_tan
    variants = a.variants.split(",")
    REF = variants[0]   # reference variant: T0 for the T-ladder, S0 (= T1 frozen replay) for the S-ladder
    fracs = [float(x) for x in a.fracs.split(",")]
    muons = sorted(p.name for p in a.root.iterdir()
                   if p.is_dir() and all((p / v / "replay" / "waveforms.csv").exists() for v in variants))
    if not muons:
        raise SystemExit("no muon directory has replay outputs for all variants " + a.variants)
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    a.out.mkdir(parents=True)
    rng = np.random.default_rng(12345)
    signal = "ION TERM ONLY (time-dependent Shockley-Ramo, frozen RKF paths)" if a.no_prompt else \
        "ion term + INSTANTANEOUS electron-step APPROXIMATION"
    report = {"muons": muons, "variants": variants, "reference": REF, "signal": signal, "fracs_of_reference_lead": fracs,
              "shaping_status": "PROVISIONAL (NA60+ VMM3a peaking time and transfer function unverified)", "Tp": {},
              "readout_geometry": {"strip_pitch_mm": PITCH_STRIP, "stereo_tan": A_STEREO}}
    lines = []
    for tp in [float(x) for x in a.tp.split(",")]:
        h, hlabel = impulse(tp, a.order, a.dt, max(12 * tp, 2000.0), a.impulse)
        per = {v: [] for v in variants[1:]}
        with open(a.out / f"per_muon_Tp{tp:g}.jsonl", "w") as fo:
            for m in muons:
                waves, labels, dt_in = {}, None, None
                for v in variants:
                    lab, tc, cur = load_wave(a.root / m / v / "replay" / "waveforms.csv")
                    if not a.no_prompt:
                        l2, tc2, pc = load_wave(a.root / m / v / "prompt" / "prompt_waveforms.csv")
                        if l2 != lab or not np.array_equal(tc, tc2):
                            raise SystemExit(f"{m}/{v}: prompt and replay grids differ")
                        cur = cur + pc
                    if labels is not None and lab != labels:
                        raise SystemExit(f"{m}/{v}: channel set differs from T0")
                    labels, dt_in = lab, tc[1] - tc[0]
                    waves[v] = cur
                res, shaped = analyse_muon(waves, labels, a.dt, fracs, h, a.window_ns, dt_in)
                for v in variants[1:]:
                    d = diffs(res[v], res[REF], fracs)
                    per[v].append(d)
                    fo.write(json.dumps({"muon": m, "variant": v, "Tp_ns": tp, **d,
                                         "U_peaks": res[v]["U"]["peaks"].tolist(), "V_peaks": res[v]["V"]["peaks"].tolist()}) + "\n")
                fo.write(json.dumps({"muon": m, "variant": REF, "Tp_ns": tp,
                                     **{f"{s}_{k}": res[REF][s][k] for s in ("U", "V")
                                        for k in ("lead_index", "lead_peak", "lead_peak_time", "charge", "amplitude",
                                                  "centroid_strip", "sign")},
                                     "w_mm": res[REF]["w_mm"], "u_mm": res[REF]["u_mm"],
                                     "U_peaks": res[REF]["U"]["peaks"].tolist(), "V_peaks": res[REF]["V"]["peaks"].tolist(),
                                     "U_strip_index": plane_cols(labels, "U")[1].tolist(), "V_strip_index": plane_cols(labels, "V")[1].tolist(),
                                     "thresholds_fCeq": {s: res[REF][s]["thresholds_fCeq"] for s in ("U", "V")}}) + "\n")
                if a.plots:
                    plot_muon(a.out / f"{m}_Tp{tp:g}.png", m, tp, labels, shaped, res, variants, a.dt,
                              f"{signal}; {hlabel}")
        summ = {v: {k: summarize([d[k] for d in per[v]], a.bootstrap, rng) for k in per[v][0]} for v in variants[1:]}
        report["Tp"][f"{tp:g}"] = {"shaper": hlabel, "reference": REF, "differences_vs_reference": summ}
        lines.append(f"== Tp={tp:g} ns  {hlabel}\n   signal: {signal}; muons={len(muons)}; thresholds = f x {REF} lead peak, fixed from {REF}")
        for v in variants[1:]:
            lines.append(f"-- {v} minus {REF}")
            for k, s in summ[v].items():
                if s.get("n"):
                    ci = f" CI95=[{s['mean_ci95'][0]:+.3g},{s['mean_ci95'][1]:+.3g}]" if "mean_ci95" in s else ""
                    lines.append(f"   {k:34s} n={s['n']:2d} mean={s['mean']:+.4g}{ci} rms={s['rms']:.4g} max|d|={s['max_abs']:.4g}")
    (a.out / "t_ladder_summary.json").write_text(json.dumps(report, indent=1) + "\n")
    (a.out / "t_ladder_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

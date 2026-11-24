#!/usr/bin/env python3
"""Focused S-ladder report: every variant against S0 (total error) AND against G (spatial-model error), per event.

  s_ladder_report.py --vs-s0 COMPARE_DIR_S0 --vs-g COMPARE_DIR_G --tag TAG --out REPORT.md [--tp 100,200]

Both directories are t_shaped_compare.py outputs over the SAME muons and variants; the first variant of
--vs-s0 is S0, the first of --vs-g is G_TAG. Every number is read from the per-muon jsonl, i.e. computed from the
actual waveforms of the two signals being compared. Nothing is derived by subtracting L2 norms. Per event and over
events (mean, max |.|) the report lists: shaped L2 per plane, leading-peak change, leading-peak time, leading-edge
times (2%, 5%, 10%), hit flips summed over thresholds, unshaped window charge per plane and the U/V charge ratio,
second/lead ratio (charge sharing), u and w (default centroid and fixed reference membership).
Section "K2 charge-ratio decomposition" lists the window-charge errors vs G of the single compressions (angular S3,
radial SR, longitudinal SW), their pairs (S3R, S3W, SRW), K2, and the hemisphere-preserving K2A, plus the additivity
check K2 - (S3 + SR + SW).
"""
import argparse
import json
from pathlib import Path

import numpy as np

KEYS = [("U_L2_rel", "L2 U", 1e2, "%"), ("V_L2_rel", "L2 V", 1e2, "%"),
        ("U_lead_peak_rel", "peak U", 1e2, "%"), ("V_lead_peak_rel", "peak V", 1e2, "%"),
        ("U_lead_peak_time_ns", "tpk U", 1, "ns"), ("V_lead_peak_time_ns", "tpk V", 1, "ns"),
        ("U_t_f0.02_ns", "t2% U", 1, "ns"), ("V_t_f0.02_ns", "t2% V", 1, "ns"),
        ("U_t_f0.05_ns", "t5% U", 1, "ns"), ("U_t_f0.1_ns", "t10% U", 1, "ns"),
        ("hitflips", "flips", 1, ""), ("U_charge_rel", "Q U", 1e2, "%"), ("V_charge_rel", "Q V", 1e2, "%"),
        ("UV_charge_ratio_rel", "Q U/V", 1e2, "%"), ("U_second_over_lead", "2nd/lead U", 1, ""),
        ("u_mm", "u", 1e3, "um"), ("u_mm_fixed_ref_membership", "u fixed", 1e3, "um"), ("w_mm", "w", 1e3, "um")]


def load(d, tp):
    rows = [json.loads(l) for l in open(Path(d) / f"per_muon_Tp{tp}.jsonl")]
    ref = next(r["variant"] for r in rows if "thresholds_fCeq" in r)
    out = {}
    for r in rows:
        if r["variant"] == ref:
            continue
        r["hitflips"] = sum(v for k, v in r.items() if "_hitflips_" in k)
        out.setdefault(r["variant"], {})[r["muon"]] = r
    return ref, out


def fmt(x, sc):
    return "n/a" if x is None or not np.isfinite(x) else f"{x * sc:+.3g}"


def table(L, data, variants, title):
    muons = sorted({m for v in variants if v in data for m in data[v]})
    L += [f"### {title}", "", "| variant | " + " | ".join(f"{lab} [{u}]" if u else lab for _, lab, _, u in KEYS) + " |",
          "|---|" + "---|" * len(KEYS)]
    for v in variants:
        if v not in data:
            continue
        cells = []
        for k, _, sc, _ in KEYS:
            x = np.array([data[v][m].get(k, np.nan) for m in muons if m in data[v]], float)
            x = x[np.isfinite(x)]
            cells.append("n/a" if len(x) == 0 else (f"{x.mean() * sc:+.3g} / {np.max(np.abs(x)) * sc:.3g}"
                                                   if len(x) > 1 else f"{x[0] * sc:+.3g}"))
        L.append(f"| {v} | " + " | ".join(cells) + " |")
    L += ["", f"{len(muons)} event(s). With several events each cell is mean / max |.| over events.", ""]
    if len(muons) > 1:
        L += [f"Per event, L2 U / L2 V [%] and |du| [um]:", "", "| variant | " + " | ".join(muons) + " |",
              "|---|" + "---|" * len(muons)]
        for v in variants:
            if v in data:
                L.append(f"| {v} | " + " | ".join(
                    f"{data[v][m]['U_L2_rel'] * 1e2:.3g} / {data[v][m]['V_L2_rel'] * 1e2:.3g} / {abs(data[v][m]['u_mm']) * 1e3:.3g}"
                    if m in data[v] else "-" for m in muons) + " |")
        L.append("")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vs-s0", type=Path, required=True)
    ap.add_argument("--vs-g", type=Path, required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--tp", default="100,200")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    T = a.tag
    order = ["G", "S1", "S1A", "S2", "S2A", "S3", "S3A", "SR", "S2R", "S3R", "SW", "SWW", "SRW", "S3W", "K1", "K1A",
             "K2", "K2A", "K1kern", "K2kern", "K1Akern", "K2Akern"]
    sig = json.loads((a.vs_s0 / "t_ladder_summary.json").read_text()).get("signal", "?")
    L = [f"# S-ladder report, tag {T}", "", f"Signal: {sig}. Shaping PROVISIONAL (CR-RC^3 stand-in).", "",
         f"Sources: `{a.vs_s0}` (reference S0) and `{a.vs_g}` (reference G_{T}); all differences are computed from "
         "the two waveforms being compared (no subtraction of norms).", "",
         "Interpretation: **vs S0** = total error of the description; **vs G** = spatial-model error alone (the Green "
         "table interpolation is common to both and cancels to first order); **G vs S0** = numerical closure.", ""]
    for tp in a.tp.split(","):
        r0, d0 = load(a.vs_s0, tp)
        rg, dg = load(a.vs_g, tp)
        vs = [f"{v}_{T}" for v in order]
        L += [f"## Tp = {tp} ns", ""]
        table(L, d0, vs, f"vs {r0} (total error)")
        table(L, dg, [v for v in vs if v != f"G_{T}"], f"vs {rg} (spatial-model error)")
        dec = [f"{v}_{T}" for v in ("S3", "SR", "SW", "S3R", "S3W", "SRW", "K2", "K2A")]
        muons = sorted({m for v in dec if v in dg for m in dg[v]})
        L += ["### K2 charge-ratio decomposition (vs G, unshaped window charge)", "",
              "| variant | " + " | ".join(f"{m}: Q U / Q V / U/V [%]" for m in muons) + " |", "|---|" + "---|" * len(muons)]
        for v in dec:
            if v in dg:
                L.append(f"| {v} | " + " | ".join(
                    f"{dg[v][m]['U_charge_rel'] * 1e2:+.3f} / {dg[v][m]['V_charge_rel'] * 1e2:+.3f} / {dg[v][m]['UV_charge_ratio_rel'] * 1e2:+.3f}"
                    if m in dg[v] else "-" for m in muons) + " |")
        if all(f"{v}_{T}" in dg for v in ("S3", "SR", "SW", "K2")):
            L.append("| K2 - (S3+SR+SW) | " + " | ".join(
                f"{(dg[f'K2_{T}'][m]['UV_charge_ratio_rel'] - sum(dg[f'{v}_{T}'][m]['UV_charge_ratio_rel'] for v in ('S3', 'SR', 'SW'))) * 1e2:+.3f} (U/V)"
                for m in muons) + " |")
        L += ["", "Reading: the compression whose single-effect row carries most of the K2 U/V error is its origin; a "
              "large K2 - sum row means the effects interact. K2A keeps the exact up/down (hemisphere) split of each cloud.", ""]
    a.out.write_text("\n".join(L) + "\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()

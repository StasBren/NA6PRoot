#!/usr/bin/env python3
"""
Analyze +u/-u reflection symmetry in the high-stat C2 boundary scan.

For the present 2.5+2.5 mm, B=0, w0=0 geometry the continuum model has the
reflection symmetry

    u -> -u,   cathode minus <-> cathode plus.

The same random seed at +u and -u is *not* expected to give an event-by-event
mirror: diffusion/avalanche random kicks are replayed with the same signs, not
with reflected signs.  Symmetry should therefore be tested primarily at the
ensemble/distribution level.

Outputs:
  c2_edge_symmetry_summary.csv
  c2_edge_symmetry_pairwise.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row, key):
    return float(row[key])


def i(row, key):
    return int(float(row[key]))


def tag(u):
    s = f"{u:.1f}".replace("-", "m").replace(".", "p")
    return f"u_{s}"


def load_point(scan_dir: Path, u):
    d = scan_dir / "events" / tag(u)
    ev_path = d / "c2_u_dense_event_summary.csv"
    st_path = d / "c2_u_dense_strip_summary.csv"
    if not ev_path.is_file() or not st_path.is_file():
        raise SystemExit(f"Missing completed point files for u0={u}: {d}")

    events = read_csv(ev_path)
    strips = read_csv(st_path)

    ids = sorted({i(r, "strip") for r in strips})
    by_event = defaultdict(lambda: {"minus": {}, "plus": {}})
    for r in strips:
        iev = i(r, "event")
        side = r["side"].strip()
        k = i(r, "strip")
        a = f(r, "A_fC")
        if not math.isfinite(a):
            a = abs(f(r, "Q_fC"))
        by_event[iev][side][k] = a

    out = {}
    for r in events:
        if i(r, "readout_valid") != 1:
            continue
        iev = i(r, "event")
        sides = {}
        ok = True
        for side in ("minus", "plus"):
            a = np.array([by_event[iev][side].get(k, math.nan) for k in ids])
            if not np.all(np.isfinite(a)) or np.sum(a) <= 0:
                ok = False
                break
            sides[side] = a / np.sum(a)
        if not ok:
            continue

        seed = i(r, "random_seed")
        out[seed] = {
            "u_raw": f(r, "u_stereo_mm"),
            "w_raw": f(r, "w_stereo_mm"),
            "active_wires": i(r, "active_wires"),
            "wire_cog_u": f(r, "wire_cog_u_mm"),
            "minus_coverage": f(r, "minus_coverage"),
            "plus_coverage": f(r, "plus_coverage"),
            "minus": sides["minus"],
            "plus": sides["plus"],
        }

    return out, tuple(ids)


def mean_std(x):
    a = np.asarray(x, dtype=float)
    return float(np.mean(a)), float(np.std(a, ddof=1))


def write_rows(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", type=Path, required=True)
    args = ap.parse_args()

    pair_rows = []
    summary_rows = []

    for absu in (1.7, 1.9):
        neg, ids_n = load_point(args.scan_dir, -absu)
        pos, ids_p = load_point(args.scan_dir, +absu)
        if ids_n != ids_p:
            raise SystemExit(f"Strip-id mismatch at |u|={absu}")
        ids = ids_n

        common = sorted(set(neg) & set(pos))
        if not common:
            raise SystemExit(f"No common seeds for |u|={absu}")

        for seed in common:
            rn, rp = neg[seed], pos[seed]

            # Exact reflection of the *readout geometry* swaps cathodes at the
            # same strip index for central strip coordinates at w0=0.
            d1 = rp["minus"] - rn["plus"]
            d2 = rp["plus"] - rn["minus"]
            mirrored = np.concatenate([d1, d2])

            pair_rows.append({
                "abs_u_mm": absu,
                "random_seed": seed,
                "u_raw_minus_mm": rn["u_raw"],
                "u_raw_plus_mm": rp["u_raw"],
                "u_pair_sum_mm": rn["u_raw"] + rp["u_raw"],
                "w_raw_minus_mm": rn["w_raw"],
                "w_raw_plus_mm": rp["w_raw"],
                "wire_cog_pair_sum_mm": rn["wire_cog_u"] + rp["wire_cog_u"],
                "active_wires_minus": rn["active_wires"],
                "active_wires_plus": rp["active_wires"],
                "mirrored_strip_rms_fraction": float(np.sqrt(np.mean(mirrored**2))),
                "mirrored_strip_l1_fraction": float(np.mean(np.abs(mirrored))),
            })

        nvals = [neg[s] for s in common]
        pvals = [pos[s] for s in common]

        u_nm, u_ns = mean_std([r["u_raw"] for r in nvals])
        u_pm, u_ps = mean_std([r["u_raw"] for r in pvals])
        w_nm, w_ns = mean_std([r["w_raw"] for r in nvals])
        w_pm, w_ps = mean_std([r["w_raw"] for r in pvals])

        minus_neg = np.mean(np.vstack([r["minus"] for r in nvals]), axis=0)
        plus_neg = np.mean(np.vstack([r["plus"] for r in nvals]), axis=0)
        minus_pos = np.mean(np.vstack([r["minus"] for r in pvals]), axis=0)
        plus_pos = np.mean(np.vstack([r["plus"] for r in pvals]), axis=0)

        mean_mirror_diff = np.concatenate([
            minus_pos - plus_neg,
            plus_pos - minus_neg,
        ])

        pair_sum = [r["u_pair_sum_mm"] for r in pair_rows if r["abs_u_mm"] == absu]
        wire_sum = [r["wire_cog_pair_sum_mm"] for r in pair_rows if r["abs_u_mm"] == absu]
        strip_rms = [r["mirrored_strip_rms_fraction"] for r in pair_rows if r["abs_u_mm"] == absu]

        summary_rows.append({
            "abs_u_mm": absu,
            "n_pairs": len(common),
            "u_minus_mean_mm": u_nm,
            "u_minus_std_mm": u_ns,
            "u_plus_mean_mm": u_pm,
            "u_plus_std_mm": u_ps,
            "u_mean_pair_sum_mm": u_nm + u_pm,
            "u_std_ratio_plus_over_minus": u_ps / u_ns if u_ns > 0 else math.nan,
            "paired_u_sum_mean_mm": float(np.mean(pair_sum)),
            "paired_u_sum_std_mm": float(np.std(pair_sum, ddof=1)),
            "w_minus_mean_mm": w_nm,
            "w_minus_std_mm": w_ns,
            "w_plus_mean_mm": w_pm,
            "w_plus_std_mm": w_ps,
            "multiwire_fraction_minus": float(np.mean([r["active_wires"] > 1 for r in nvals])),
            "multiwire_fraction_plus": float(np.mean([r["active_wires"] > 1 for r in pvals])),
            "wire_cog_pair_sum_mean_mm": float(np.mean(wire_sum)),
            "wire_cog_pair_sum_std_mm": float(np.std(wire_sum, ddof=1)),
            "minus_coverage_mean_minus_u": float(np.mean([r["minus_coverage"] for r in nvals])),
            "plus_coverage_mean_minus_u": float(np.mean([r["plus_coverage"] for r in nvals])),
            "minus_coverage_mean_plus_u": float(np.mean([r["minus_coverage"] for r in pvals])),
            "plus_coverage_mean_plus_u": float(np.mean([r["plus_coverage"] for r in pvals])),
            "ensemble_mirrored_strip_rms_fraction": float(np.sqrt(np.mean(mean_mirror_diff**2))),
            "ensemble_mirrored_strip_l1_fraction": float(np.mean(np.abs(mean_mirror_diff))),
            "paired_mirrored_strip_rms_mean_fraction": float(np.mean(strip_rms)),
        })

    pair_file = args.scan_dir / "c2_edge_symmetry_pairwise.csv"
    summary_file = args.scan_dir / "c2_edge_symmetry_summary.csv"

    write_rows(pair_file, pair_rows, list(pair_rows[0]))
    write_rows(summary_file, summary_rows, list(summary_rows[0]))

    print("C2 high-stat edge symmetry")
    for r in summary_rows:
        print(
            f"  |u|={r['abs_u_mm']:.1f} mm: "
            f"<u->={r['u_minus_mean_mm']:+.4f}, "
            f"<u+>={r['u_plus_mean_mm']:+.4f} mm, "
            f"sum={r['u_mean_pair_sum_mm']:+.4f} mm; "
            f"multiwire(-,+)=({r['multiwire_fraction_minus']:.3f}, "
            f"{r['multiwire_fraction_plus']:.3f}); "
            f"mean-pattern mirror RMS={r['ensemble_mirrored_strip_rms_fraction']:.5f}"
        )
    print()
    print(f"Wrote: {summary_file}")
    print(f"       {pair_file}")


if __name__ == "__main__":
    main()

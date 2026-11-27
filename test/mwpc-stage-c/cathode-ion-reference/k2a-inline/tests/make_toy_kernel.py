#!/usr/bin/env python3
"""SYNTHETIC toy kernel directory + toy cloud moments for the compiled strip-evaluator tests (NOT physics).

  make_toy_kernel.py OUTDIR [--seed 5]
  make_toy_kernel.py --single-cloud SRC DST | --gappy SRC DST | --corrupt WAVEFORMS.csv   (test-file helpers)

OUTDIR/kern/{K.npy, KPHI0.npy, kernel.json}   s_kernel.py layout: basis K0 Kc1 Ks1 Kc2 Ks2 KA1 KA2, K (basis, side,
                                              wire, r, D, tau) float32, KPHI0 (basis, side, wire, r, D) float64; smooth
                                              in r, D and tau so interpolation errors would show
OUTDIR/moments_*.csv                          toy events (BirthMoments column schema):
  multi     12 primaries on wires -2..2, primary-mean D_U / D_V from -13 to +13 mm (shifts m = -7..+7, both signs, both
            planes, several primaries per group, multi-wire primaries)
  near      every cloud near the entry window (all m = 0)
  tie       one primary with Dbar_U / P + 1/2 exactly an integer (shift tie)
  late      one cloud after the 7000 ns window (informational) besides normal clouds
  bad_wire / bad_r_low / bad_r_high / bad_spread / bad_t / bad_q / bad_nan / bad_dnan   one blocking case each
"""
import argparse
import json
from pathlib import Path

import numpy as np

COLS = ["primary", "wire", "t_ns", "t_rms_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A",
        "min_clearance_um", "base_U", "base_V", "t_raw_rms_ns", "rows"]


def kernel(out, rng):
    wires = [-2, -1, 0, 1, 2]
    r_um = list(np.geomspace(18.0, 1500.0, 8))
    M = 16; theta = list((np.arange(M) + 0.5) * 360.0 / M)
    P, dsub = 1.7, 10; dD = P / dsub
    nD = int(np.floor(16.0 / dD)); D = np.arange(-nD, nD + 1) * dD
    tau = np.r_[0.0, np.geomspace(0.01, 8000.0, 39)]
    names = ["K0", "Kc1", "Ks1", "Kc2", "Ks2", "KA1", "KA2"]
    shape = (len(names), 2, len(wires), len(r_um), len(D), len(tau))
    b = np.arange(len(names))[:, None, None, None, None, None]
    s = np.arange(2)[None, :, None, None, None, None]
    k = np.array(wires)[None, None, :, None, None, None]
    lr = np.log(np.array(r_um))[None, None, None, :, None, None]
    d = D[None, None, None, None, :, None]
    t = tau[None, None, None, None, None, :]
    amp = (1.0 / (1 + b)) * (1 + 0.3 * s) * (1 + 0.1 * k) * (1 + 0.05 * lr)
    prof = np.exp(-((d - 0.4 * k - 0.2 * b) / (1.3 + 0.1 * lr)) ** 2)
    rise = 1 - np.exp(-t / (30.0 * (1 + 0.2 * b) * (1 + 0.1 * lr)))
    K = (amp * prof * rise * (1 + 0.02 * np.sin(0.37 * t ** 0.5 + d))).astype(np.float32)
    KP = (0.01 * amp[..., 0] * prof[..., 0]) * np.ones(shape[:-1])
    out.mkdir(parents=True)
    np.save(out / "K.npy", K); np.save(out / "KPHI0.npy", KP)
    ax = {"wires": wires, "r_um": r_um, "theta_deg": theta, "D_mm": D.tolist(), "dD_mm": dD, "d_sub": dsub,
          "tau_ns": tau.tolist(), "pitch_mm": P, "stereo_tan": 0.1, "shape": [2, len(wires), len(r_um), M, len(D), len(tau)],
          "lut": "SYNTHETIC", "lut_sha256": "0" * 64, "grid_vertices_sha256": "1" * 64,
          "grid_seal": "SYNTHETIC toy kernel\nline 2 \"quoted\"\n", "convention": "SYNTHETIC (not physics)"}
    (out / "kernel.json").write_text(json.dumps({"basis": names, "green": "SYNTHETIC", "green_axes": ax, "shape": list(shape),
                                                 "bytes": int(K.nbytes)}, indent=1) + "\n")
    return P


def cloud(rng, p, w, t, N, D_U, D_V, r=None):
    th = rng.uniform(0, 2 * np.pi, 400)
    return {"primary": p, "wire": w, "t_ns": t, "t_rms_ns": 0.0, "N_fC": N, "r_eff_um": r if r is not None else float(rng.uniform(25, 1200)),
            "D_U_mm": D_U, "D_V_mm": D_V, "C1": float(np.cos(th).mean()), "S1": float(np.sin(th).mean()),
            "C2": float(np.cos(2 * th).mean()), "S2": float(np.sin(2 * th).mean()), "A": float(np.sign(np.sin(th)).mean()),
            "min_clearance_um": 1.0, "base_U": 37, "base_V": -12, "t_raw_rms_ns": 0.1, "rows": 400}


def write(path, rows):
    rows = sorted(rows, key=lambda r: (r["primary"], r["wire"]))
    with open(path, "w") as f:
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(str(int(r[c])) if c in ("primary", "wire", "base_U", "base_V", "rows") else repr(float(r[c])) for c in COLS) + "\n")


def helpers(argv):
    """Small file transformations for the shell tests (no heredocs needed there)."""
    import csv
    if argv[0] == "--single-cloud":       # SRC DST: first wire-0 cloud only (one wire, m = 0)
        r = list(csv.DictReader(open(argv[1])))
        w = csv.DictWriter(open(argv[2], "w", newline=""), fieldnames=list(r[0])); w.writeheader()
        w.writerows([x for x in r if x["wire"] == "0"][:1])
    elif argv[0] == "--gappy":            # SRC DST: source primaries 0, 1 renumbered to ORIGINAL indices 2, 5
        r = list(csv.DictReader(open(argv[1])))
        out = [dict(x, primary={"0": "2", "1": "5"}[x["primary"]]) for x in r if x["primary"] in ("0", "1")]
        w = csv.DictWriter(open(argv[2], "w", newline=""), fieldnames=list(r[0])); w.writeheader(); w.writerows(out)
    elif argv[0] == "--corrupt":          # WAVEFORMS.csv: scale one sample by 1.001
        p = argv[1]; h = open(p).readline().strip(); a = np.loadtxt(p, delimiter=",", skiprows=1); a[50, 2] *= 1.001
        np.savetxt(p, a, delimiter=",", header=h, comments="")
    else:
        raise SystemExit(f"unknown helper {argv[0]}")


def main():
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ("--single-cloud", "--gappy", "--corrupt"):
        return helpers(sys.argv[1:])
    ap = argparse.ArgumentParser(); ap.add_argument("out", type=Path); ap.add_argument("--seed", type=int, default=5)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    P = kernel(a.out / "kern", rng)
    multi = []
    means = np.linspace(-13, 13, 12)
    for p, mu in enumerate(means):
        mv = -0.7 * mu + 1.1
        for w in sorted(rng.choice([-2, -1, 0, 1, 2], size=int(rng.integers(1, 4)), replace=False)):
            multi.append(cloud(rng, p, int(w), float(rng.uniform(5, 300)), float(rng.uniform(0.5, 40)),
                               float(mu + rng.uniform(-0.6, 0.6)), float(mv + rng.uniform(-0.6, 0.6))))
    write(a.out / "moments_multi.csv", multi)
    near = [cloud(rng, p, w, float(rng.uniform(5, 300)), float(rng.uniform(1, 20)), float(rng.uniform(-0.5, 0.5)),
                  float(rng.uniform(-0.5, 0.5))) for p in range(4) for w in (-1, 0, 1)]
    write(a.out / "moments_near.csv", near)
    write(a.out / "moments_tie.csv", [cloud(rng, 0, 0, 20.0, 2.0, 1.5 * P, 0.2), cloud(rng, 1, 1, 30.0, 3.0, 0.1, 0.1)])
    write(a.out / "moments_late.csv", near[:3] + [cloud(rng, 9, 0, 7500.0, 2.0, 0.2, 0.1)])
    bad = {"bad_wire": cloud(rng, 5, 3, 20.0, 2.0, 0.1, 0.1), "bad_r_low": cloud(rng, 5, 0, 20.0, 2.0, 0.1, 0.1, r=10.0),
           "bad_r_high": cloud(rng, 5, 0, 20.0, 2.0, 0.1, 0.1, r=5000.0), "bad_t": cloud(rng, 5, 0, -1.0, 2.0, 0.1, 0.1),
           "bad_q": cloud(rng, 5, 0, 20.0, 0.0, 0.1, 0.1), "bad_nan": cloud(rng, 5, 0, float("nan"), 2.0, 0.1, 0.1),
           "bad_dnan": cloud(rng, 5, 0, 20.0, 2.0, float("nan"), 0.1)}
    for name, c in bad.items():
        write(a.out / f"moments_{name}.csv", near[:3] + [c])
    # one primary whose two clouds straddle +-12 mm: mean 0 -> m = 0, local D' = +-12 mm outside the +-5-strip table
    write(a.out / "moments_bad_spread.csv", near[:3] + [cloud(rng, 5, 0, 20.0, 2.0, 12.0, 0.1), cloud(rng, 5, 1, 20.0, 2.0, -12.0, 0.1)])
    print(f"SYNTHETIC toy kernel and moments in {a.out}")


if __name__ == "__main__":
    main()

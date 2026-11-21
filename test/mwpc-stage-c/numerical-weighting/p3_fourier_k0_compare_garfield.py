#!/usr/bin/env python3
"""Compare the normalized k=0 MFS mode with independent Garfield++ data.
This is a whole-plane benchmark, not a stereo-strip validation.
No input files are changed. Garfield's analytic wire treatment can itself
have a finite-radius approximation; discrepancies are not assigned to MFS
without further analysis.
"""
import argparse
import csv
import numpy as np
from p3_fourier_mfs_mode_pilot import fit, ring

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("reference_csv")
    args=ap.parse_args()
    with open(args.reference_csv,newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=1080:
        raise ValueError(f"Expected 1080 reference samples, got {len(rows)}")
    x=np.array([[float(r["u_mm"]),float(r["v_mm"])] for r in rows])
    ref=np.array([[float(r["phi_plus"]),float(r["phi_minus"])] for r in rows])
    if not np.all(np.isfinite(ref)):
        raise ValueError("Nonfinite reference values")
    f24,c24=fit(0.,24,0)
    f32,c32=fit(0.,32,0)
    def predict(fn):
        mirrored=x.copy()
        mirrored[:,1]*=-1
        return np.column_stack((fn(x).real,fn(mirrored).real))
    p24,p32=predict(f24),predict(f32)
    def report(label,sel):
        d=(p32-ref)[sel]
        aa=ref[sel]
        md=(p32-p24)[sel]
        bb=p32[sel]
        e=np.linalg.norm(d)/max(np.linalg.norm(aa),1e-300)
        em=np.linalg.norm(md)/max(np.linalg.norm(bb),1e-300)
        print(f"{label:18s} N={len(aa):4d}  ref_relL2={e*100:.6f}% "
              f"max_abs={np.max(np.abs(d)):.8g} "
              f"MFS24_32_relL2={em*100:.6f}%")
    print("Both cathodes combined; k=0 mode, nine wires, mm.")
    report("ALL",np.ones(len(rows),dtype=bool))
    ratio=np.array([float(r["r_over_R"]) for r in rows])
    for q in (1.01,1.1,1.5,2.,5.,20.):
        report(f"r/R={q:g}",ratio==q)
    report("bulk",ratio==0.)
    wire=np.max(np.abs(f32(ring(97,.015,.193))))
    print(f"MFS32 independent wire residual={wire:.9g}; cond={c32:.6g}")
    print(f"MFS32 sampled range=[{p32.min():.9g}, {p32.max():.9g}]")
    print("Comparison only: inspect near-wire discrepancies before accepting.")
    print("No inverse Fourier reconstruction, strip waveforms or production run.")

if __name__=="__main__":
    main()

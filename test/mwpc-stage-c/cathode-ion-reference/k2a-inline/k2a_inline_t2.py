#!/usr/bin/env python3
"""T2 (DIAGNOSTIC ONLY): K256 reservoir moments vs the in-memory full-birth moments of the SAME avalanches.

  k2a_inline_t2.py FULL_PREFIX K256_MOMENTS.csv OUTDIR

FULL_PREFIX     inline generator outputs (the full-population moments: the truth for this comparison)
K256_MOMENTS    p3_birth_moments on the weighted K=256 archive written by the same generator run family
Reports, per cloud, the K256 estimate minus the full value (relative for N_fC and r_eff, absolute otherwise), clouds the
reservoir missed (and their share of the full charge), and charge-weighted summaries. It never fails on the size of the
differences (the reservoir is a statistical estimate, not a reference); exit 3 only if an input is missing.
"""
import csv
import json
import math
import sys
from pathlib import Path

REL = ["N_fC", "r_eff_um"]
ABS = ["D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A", "t_ns", "t_raw_rms_ns", "min_clearance_um"]


def load(p):
    with open(p, newline="") as f:
        return {(int(r["primary"]), int(r["wire"])): r for r in csv.DictReader(f)}


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    full_p, k_p, out = sys.argv[1] + ".moments.csv", Path(sys.argv[2]), Path(sys.argv[3])
    if out.exists():
        raise SystemExit(f"refusing: {out} exists")
    out.mkdir(parents=True)
    if not (Path(full_p).is_file() and k_p.is_file()):
        (out / "t2.json").write_text(json.dumps({"status": "INCONCLUSIVE", "missing": [full_p, str(k_p)]}) + "\n")
        print("INCONCLUSIVE T2: input missing"); return 3
    F, Kc = load(full_p), load(k_p)
    Q = sum(float(r["N_fC"]) for r in F.values())
    missed = sorted(set(F) - set(Kc)); extra = sorted(set(Kc) - set(F))
    rows, agg = [], {f: [0.0, 0.0] for f in REL + ABS}
    for key in sorted(set(F) & set(Kc)):
        f, k = F[key], Kc[key]
        q = float(f["N_fC"])
        row = {"primary": key[0], "wire": key[1], "births_full": f["rows"], "rows_k256": k["rows"], "N_full_fC": f["N_fC"]}
        for x in REL:
            d = float(k[x]) / float(f[x]) - 1; row[f"rel_{x}"] = f"{d:.6e}"; agg[x][0] += q * abs(d); agg[x][1] = max(agg[x][1], abs(d))
        for x in ABS:
            d = float(k[x]) - float(f[x]); row[f"d_{x}"] = f"{d:.6e}"; agg[x][0] += q * abs(d); agg[x][1] = max(agg[x][1], abs(d))
        rows.append(row)
    with open(out / "t2_clouds.csv", "w", newline="") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    qm = sum(float(F[k]["N_fC"]) for k in missed)
    rep = {"status": "DIAGNOSTIC", "note": "K256 is a reservoir estimate; the full-birth moments are the truth here",
           "clouds_full": len(F), "clouds_k256": len(Kc), "clouds_missed_by_k256": [list(k) for k in missed],
           "missed_charge_fraction": qm / Q if Q else 0.0, "clouds_only_in_k256": [list(k) for k in extra],
           "charge_weighted_mean_abs": {x: (agg[x][0] / Q if Q else math.nan) for x in agg},
           "max_abs": {x: agg[x][1] for x in agg},
           "units": {"N_fC": "relative", "r_eff_um": "relative", "D_U_mm": "mm", "D_V_mm": "mm", "t_ns": "ns", "t_raw_rms_ns": "ns",
                     "min_clearance_um": "um", "C1": "1", "S1": "1", "C2": "1", "S2": "1", "A": "1"}}
    (out / "t2.json").write_text(json.dumps(rep, indent=1) + "\n")
    print(f"DIAGNOSTIC T2 K256 vs full-birth: {len(rows)} common clouds, {len(missed)} missed by K256 "
          f"({rep['missed_charge_fraction']:.3g} of the charge); charge-weighted |rel N| {rep['charge_weighted_mean_abs']['N_fC']:.3g}, "
          f"|rel r_eff| {rep['charge_weighted_mean_abs']['r_eff_um']:.3g}, |d D_U| {rep['charge_weighted_mean_abs']['D_U_mm']:.3g} mm, "
          f"|d A| {rep['charge_weighted_mean_abs']['A']:.3g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

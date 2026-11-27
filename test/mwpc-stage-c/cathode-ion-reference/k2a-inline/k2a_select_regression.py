#!/usr/bin/env python3
"""Choose a small REAL avalanche for the multiwire / displaced-primary regression (read-only; no simulation).

  k2a_select_regression.py --moments-dir PF/moments_nonsealed --slots PF/inventory/slots.csv --kernel KERNEL_DIR
                           --out CANDIDATES.csv [--max-births 2000000] [--top 15]

Input: the preflight's K256 batch moments of the NON-sealed archives (p3_birth_moments --batch; one CSV per event).
The sealed moments directory is refused, and so is any event whose slot is sealed or unknown in slots.csv.
For every event and every prefix N of Heed primaries (the run_t1_inline.sh --max-heed-primaries option, which runs the
ORIGINAL indices 0..N-1), using the K256 estimates of the clouds of primaries 0..N-1. Zero-ion primaries have no row in
the moments CSV, so a prefix ending at primary p is N = p + 1 (never the count of non-empty primaries); budget and domain
are evaluated over exactly that prefix. The CSV primaries are checked against the archive metadata <moments>.event.json
(n_primaries, n_primaries_zero_ions): an inconsistent event is skipped, a missing one ranks the candidate last
("census unverified"); run_t1_inline.sh MWPC_T1_PREVIEW=1 then reads the frozen census of the same prefix.
  births      sum of N_fC / e (the reservoir weights restore each primary's exact census)      <= --max-births
  multiwire   >= 2 distinct wires, each carrying >= 5 % of the selected charge
  displaced   >= 1 primary with m_U != 0 or m_V != 0 under the production rule floor(Dbar/P + 1/2), with its Dbar/P + 1/2
              at least 0.05 from an integer (so the full-birth moments are expected to give the same m)
  domain      every selected cloud inside the kernel: wire tabulated, r_eff within the radial grid by >= 5 % in log r,
              local D' = D - m P with all +-5 strips inside the D table with >= 2 grid steps of margin
Ranked: both criteria first, then a primary spanning >= 2 wires, then fewer births. This uses estimates only: the
selected run is then validated in full by run_t1_inline.sh (full-birth moments) and run_m2_adapter.sh (evaluation).
"""
import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

EFC = 1.602176634e-4


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--moments-dir", type=Path, required=True); ap.add_argument("--slots", type=Path, required=True)
    ap.add_argument("--kernel", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-births", type=float, default=2e6); ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--half-strips", type=int, default=5)
    a = ap.parse_args()
    if a.moments_dir.resolve().name == "moments_sealed" or "SEALED" in str(a.moments_dir).upper().replace("NONSEALED", ""):
        raise SystemExit("refusing: this is the sealed moments directory")
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    roles = {r["slot"]: r["role"] for r in csv.DictReader(open(a.slots, newline=""))}
    ax = json.loads((a.kernel / "kernel.json").read_text())["green_axes"]
    P, D0, dD, nD, dsub = ax["pitch_mm"], ax["D_mm"][0], ax["dD_mm"], len(ax["D_mm"]), ax["d_sub"]
    wires, lr0, lr1 = set(int(w) for w in ax["wires"]), math.log(ax["r_um"][0]), math.log(ax["r_um"][-1])
    hs = a.half_strips
    cands, scanned, skipped, inconsistent = [], 0, 0, 0
    for f in sorted(a.moments_dir.rglob("*.moments.csv")):
        m = re.findall(r"(\d{4})", f.name)
        slot = m[-1] if m else None
        if slot is None or roles.get(slot) in (None, "sealed_validation"):
            skipped += 1; continue
        rows = list(csv.DictReader(open(f, newline="")))
        if not rows:
            continue
        scanned += 1
        byp = defaultdict(list)
        for r in rows:
            byp[int(r["primary"])].append(r)
        prims = sorted(byp)
        # the moments CSV holds only primaries WITH ions: an absent index is a zero-ion primary. --max-heed-primaries N
        # runs ORIGINAL indices 0..N-1, so a prefix ending at primary p is N = p + 1 and includes every absent index < p
        # (zero births, no cloud). Cross-check against the archive metadata (<moments>.event.json) when present.
        meta_p = Path(str(f) + ".event.json")
        census = "unverified (no event.json)"
        if meta_p.is_file():
            meta = json.loads(meta_p.read_text())
            n_all, n_zero = int(meta["n_primaries"]), int(meta["n_primaries_zero_ions"])
            if prims[0] < 0 or prims[-1] >= n_all or n_all - len(prims) != n_zero:
                inconsistent += 1
                print(f"  skip {f.name}: primaries in CSV {len(prims)} (max index {prims[-1]}) inconsistent with archive "
                      f"n_primaries {n_all}, zero-ion {n_zero}")
                continue
            census = "verified"
        elif prims[0] < 0:
            inconsistent += 1; continue
        births = 0.0; clouds = []; shifts = {}; ok_domain = True
        for p in prims:
            cl = byp[p]
            N = [float(r["N_fC"]) for r in cl]; Ns = sum(N)
            births += Ns / EFC
            dU = sum(x * float(r["D_U_mm"]) for x, r in zip(N, cl)) / Ns; dV = sum(x * float(r["D_V_mm"]) for x, r in zip(N, cl)) / Ns
            xu, xv = dU / P + .5, dV / P + .5
            mu, mv = math.floor(xu), math.floor(xv)
            margin = min(abs(xu - round(xu)), abs(xv - round(xv)))
            shifts[p] = (mu, mv, margin, len({int(r["wire"]) for r in cl}))
            for r in cl:
                w = int(r["wire"]); lr = math.log(float(r["r_eff_um"]))
                inside = w in wires and lr >= lr0 + 0.05 * (lr1 - lr0) and lr <= lr1 - 0.05 * (lr1 - lr0)
                for D, mm in ((float(r["D_U_mm"]), mu), (float(r["D_V_mm"]), mv)):
                    jd = math.floor((D - P * mm - D0) / dD)
                    inside &= jd - dsub * hs >= 2 and jd + 1 + dsub * hs < nD - 2
                ok_domain &= inside
                clouds.append((w, float(r["N_fC"])))
            if births > a.max_births:
                break
            if not ok_domain:
                break
            Q = sum(q for _, q in clouds)
            wq = defaultdict(float)
            for w, q in clouds:
                wq[w] += q
            multi = sum(1 for q in wq.values() if q >= 0.05 * Q) >= 2
            disp = [pp for pp, s in shifts.items() if (s[0] or s[1]) and s[2] >= 0.05]
            mwp = any(s[3] >= 2 for s in shifts.values())
            if multi or disp:
                cands.append({"slot": slot, "role": roles[slot], "N": p + 1, "last_primary_index": p,
                              "zero_ion_primaries_in_prefix": p + 1 - sum(1 for q in prims if q <= p), "census": census,
                              "expected_births": int(round(births)),
                              "multiwire": int(multi), "displaced_primaries": " ".join(map(str, disp)), "multiwire_primary": int(mwp),
                              "wires": " ".join(f"{w}:{q / Q:.2f}" for w, q in sorted(wq.items())),
                              "shifts": " ".join(f"{pp}:{s[0]},{s[1]}" for pp, s in sorted(shifts.items())),
                              "min_shift_margin": min(s[2] for s in shifts.values()), "source": str(f)})
    cands.sort(key=lambda c: (c["census"] != "verified", -(c["multiwire"] and bool(c["displaced_primaries"])), -c["multiwire_primary"],
                              c["expected_births"]))
    with open(a.out, "w", newline="") as fh:
        if cands:
            w = csv.DictWriter(fh, fieldnames=list(cands[0])); w.writeheader(); w.writerows(cands[:max(a.top, 1) * 20])
    print(f"scanned {scanned} non-sealed events ({skipped} files skipped: sealed/unknown slot; {inconsistent} skipped: census "
          f"inconsistent); {len(cands)} candidate prefixes")
    for c in cands[:a.top]:
        print(f"  slot {c['slot']} N={c['N']} (primaries 0..{c['last_primary_index']}, {c['zero_ion_primaries_in_prefix']} zero-ion, "
              f"census {c['census']}): ~{c['expected_births']} births, multiwire={c['multiwire']}, displaced primaries "
              f"[{c['displaced_primaries']}], wires {c['wires']}, shifts {c['shifts']}")
    return 0 if cands else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""SYNTHETIC test of the read-only preflight (k2a_inventory.py, k2a_domain_check.py). Not NA60+ data.

Builds a fake 30-slot archive (seal files + SHA256SUMS over a dummy ion_births.root), a fake selection directory
(inventory, excluded slots, development/sealed lists, frozen-source checksums), a fake kernel.json (axes only), and
converter-style outputs with one event of every kind: in domain, wire outside, D outside, r_eff outside, late cloud,
converter error, no ions, metadata mismatch, never converted; plus sealed events in a separate directory.
Asserts every status/cause, the sealed blinding, the separation guard and the extension estimate.
  test_preflight_synthetic.py WORKDIR
"""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent
W = Path(sys.argv[1]); W.mkdir(parents=True, exist_ok=False)
A, S, K, M, MS = W / "archive", W / "selection", W / "kernel", W / "moments", W / "sealed_moments"
for d in (A, S, K, M, MS):
    d.mkdir()
NS = 30
excluded = [7, 19, 23]
sealed, dev = [3, 4, 11], [0, 1, 2, 5]
chan = {i: ("jpsi", "phi", "omega")[i % 3] for i in range(NS)}
for i in range(NS):
    d = A / f"muon_{i:04d}"; d.mkdir()
    if i == 7:
        (d / "failed.txt").write_text("worker exited 137 (synthetic)\n"); continue
    for f in ("ion_births.root", "microscopic_wire_muon.csv", "root_audit.log"):
        (d / f).write_bytes(f"synthetic {i} {f}\n".encode())
    if i in (19, 23):
        continue                                     # unsealed: no complete.txt / SHA256SUMS
    (d / "complete.txt").write_text("ok\n")
    sums = "".join(f"{hashlib.sha256((d / f).read_bytes()).hexdigest()}  {f}\n" for f in ("ion_births.root", "microscopic_wire_muon.csv"))
    (d / "SHA256SUMS").write_text(sums)
    if i == 25:                                      # corrupted after sealing, but listed as usable
        (d / "ion_births.root").write_bytes(b"corrupted\n")
with open(S / "archive_inventory.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["index", "channel", "pdg", "station", "det_id", "heed_primaries", "ion_births_exact"])
    for i in range(NS):
        if i not in excluded:
            w.writerow([i, chan[i], 13, i % 6, 2 * (i % 6), 10, 1000])
(S / "excluded_slots.txt").write_text("".join(f"{i:04d}\n" for i in excluded))
(S / "previously_inspected_not_sealed.txt").write_text("0006\n")
for name, ids in (("development20.csv", dev), ("sealed_validation20.csv", sealed)):
    (S / name).write_text("role,index,channel\n" + "".join(f"x,{i},{chan[i]}\n" for i in ids))
(S / "FROZEN_SOURCE_SHA256.txt").write_text("".join(
    f"{hashlib.sha256((S / n).read_bytes()).hexdigest()}  {n}\n" for n in
    ("archive_inventory.csv", "excluded_slots.txt", "development20.csv", "sealed_validation20.csv")))
# reference run: muons 0 and 1 have S0
for i in (0, 1):
    (W / "ref" / f"muon_{i:04d}" / "S0" / "replay").mkdir(parents=True)
    (W / "ref" / f"muon_{i:04d}" / "S0" / "replay" / "waveforms.csv").write_text("t\n")
inv = subprocess.run([sys.executable, str(HERE / "k2a_inventory.py"), "--archive", str(A), "--selection", str(S),
                      "--out", str(W / "inv"), "--slots", str(NS), "--reference-run", str(W / "ref")], capture_output=True, text=True)
print("(expected FAIL below: synthetic corrupted slot 0025)\n" + inv.stdout[-300:])
assert inv.returncode == 1 and "slot 0025" in inv.stdout, "corrupted inventory slot 25 must be reported INCONSISTENT"
slots = {r["slot"]: r for r in csv.DictReader(open(W / "inv" / "slots.csv"))}
assert slots["0007"]["status"] == slots["0019"]["status"] == slots["0023"]["status"] == "EXCLUDED"
assert slots["0025"]["status"] == "INCONSISTENT" and slots["0025"]["sha_check"].startswith("FAIL")
assert slots["0000"]["role"] == "dev_reference" and slots["0002"]["role"] == "development" and slots["0003"]["role"] == "sealed_validation"
assert slots["0006"]["role"] == "previously_inspected" and slots["0008"]["sha_check"] == "ok"
assert "failed.txt" in (W / "inv" / "excluded_slots_detail.md").read_text()
nonsealed = (W / "inv" / "events_nonsealed.txt").read_text().split("\n")
assert not any(l.startswith(f"muon_{i:04d}") for i in sealed for l in nonsealed)
print("PASS inventory: excluded/unsealed/corrupted slots, roles, SHA check, sealed lists separated")
# fake kernel axes
ax = {"wires": [-1, 0, 1], "r_um": list(np.geomspace(15.01, 400, 14)), "D_mm": list(np.arange(-9.9, 9.91, 0.17)),
      "tau_ns": list(np.linspace(0, 8000, 50)), "theta_deg": list(range(0, 360, 5)), "dD_mm": 0.17, "d_sub": 10,
      "build_s": 100.0, "pitch_mm": 1.7, "stereo_tan": 0.1, "grid_vertices_sha256": "t" * 64, "lut_sha256": "w" * 64}
ax["D_mm"] = list(np.round(np.arange(-129, 130) * 0.017, 12)); ax["dD_mm"] = 0.017   # d_sub 10: dD = P/10
(K / "kernel.json").write_text(json.dumps({"basis": ["K0"], "green_axes": ax, "shape": [7, 2, 3, 14, 259, 50],
                                           "bytes": 7 * 2 * 3 * 14 * 259 * 50 * 4, "build_s": 50.0}))
COLS = "primary,wire,t_ns,t_rms_ns,N_fC,r_eff_um,D_U_mm,D_V_mm,C1,S1,C2,S2,A,min_clearance_um,base_U,base_V,t_raw_rms_ns,rows"
rng = np.random.default_rng(3)


def write(dirp, eid, slot, n=8, wire=None, D=None, r=None, t=None, meta_station=None, empty=False):
    p = dirp / f"{eid}.moments.csv"
    rows = [COLS]
    for c in range(0 if empty else n):
        rows.append(",".join(map(str, [c, (wire if wire is not None and c == 0 else int(rng.integers(-1, 2))),
                                       (t if t is not None and c == 0 else rng.uniform(5, 80)), 3, rng.uniform(1, 50),
                                       (r if r is not None and c == 0 else rng.uniform(20, 60)),
                                       (D if D is not None and c == 0 else rng.uniform(-.8, .8)), rng.uniform(-.8, .8),
                                       .1, .2, .0, .0, .3, 1.0, 5, 6, 3, 100])))
    p.write_text("\n".join(rows) + "\n")
    ev = {"source_channel": {"jpsi": "Jpsi_mumu", "phi": "phi_mumu", "omega": "omega_mumu"}[chan[slot]],
          "station": meta_station if meta_station is not None else slot % 6, "det_id": 2 * (slot % 6), "pdg": 13,
          "n_primaries": 10, "ion_births_exact": 1000, "ion_rows_saved": 300, "rows_read": 100 * n, "theta_u_deg": 1.0, "theta_w_deg": 2.0,
          "momentum_GeV": 10.0, "n_primaries_zero_ions": 0}
    Path(str(p) + ".event.json").write_text(json.dumps(ev))
    Path(str(p) + ".frame.json").write_text("{}")


usable = [i for i in range(NS) if slots[f"{i:04d}"]["status"] == "USABLE"]
special = {8: dict(wire=3), 9: dict(D=9.5), 10: dict(r=12.0e3), 12: dict(t=7500.0), 14: dict(meta_station=5), 15: dict(empty=True)}
for i in usable:
    if i in sealed:
        write(MS, f"muon_{i:04d}", i, **({"wire": 4} if i == 11 else {}))
    elif i == 13:
        (M / f"muon_{i:04d}.error.txt").write_text("weighted ion census not conserved\n")
    elif i == 16:
        continue                                     # never converted
    else:
        write(M, f"muon_{i:04d}", i, **special.get(i, {}))
dc = [sys.executable, str(HERE / "k2a_domain_check.py"), "--kernel", str(K), "--inventory", str(W / "inv"), "--selection", str(S)]
r = subprocess.run(dc + ["--moments", str(M), "--sealed-moments", str(MS), "--out", str(W / "dom")], capture_output=True, text=True)
assert r.returncode == 0, r.stderr
pe = {x["event"]: x for x in csv.DictReader(open(W / "dom" / "per_event.csv"))}
exp = {8: ("OUT_OF_DOMAIN", "WIRE"), 9: ("OUT_OF_DOMAIN", "D_U"), 10: ("OUT_OF_DOMAIN", "R_EFF_HIGH"),
       12: ("OUT_OF_DOMAIN", "T_WINDOW"), 13: ("READ_FAILED", "census"), 14: ("META_MISMATCH", "station"),
       15: ("NO_IONS", "NO_IONS"), 16: ("NOT_CONVERTED", "NOT_CONVERTED"), 0: ("IN_DOMAIN", "")}
for i, (st, cause) in exp.items():
    e = pe[f"muon_{i:04d}"]
    assert e["status"] == st and cause in e["causes"], (i, e["status"], e["causes"])
assert not any(f"muon_{i:04d}" in pe for i in sealed), "sealed event in the main per-event table"
cov = json.loads((W / "dom" / "coverage.json").read_text())
assert cov["sealed_blinded"]["events"] == 3 and "OUT_OF_DOMAIN" in cov["sealed_blinded"]["status_counts"]
md = (W / "dom" / "coverage.md").read_text()
assert "muon_0011" not in md and "muon_0003" not in md, "sealed ids leaked into the report"
ext = json.loads((W / "dom" / "extension_estimate.json").read_text())
need = ext["requirement_non_sealed"]
assert need["missing_wires"] == [3], need["missing_wires"]            # sealed wire 4 must NOT drive the design
assert need["D_nodes_to_add"]["high"] > 0 and need["radial_nodes_to_add"]["outer"] > 0 and need["window_extension_needed"]
print(f"PASS domain check: statuses {sorted({x['status'] for x in pe.values()})}; sealed blinded "
      f"{cov['sealed_blinded']['status_counts']}; extension from non-sealed only (missing wires {need['missing_wires']}, "
      f"table factor {ext['cost']['factor_table_and_kernel_size']:.3g})")
# separation guard: a sealed event among the non-sealed moments must stop the tool
(M / "muon_0011.moments.csv").write_text((MS / "muon_0011.moments.csv").read_text())
r = subprocess.run(dc + ["--moments", str(M), "--out", str(W / "dom2")], capture_output=True, text=True)
assert r.returncode != 0 and "separation violated" in (r.stderr + r.stdout)
print("PASS sealed/non-sealed separation guard")

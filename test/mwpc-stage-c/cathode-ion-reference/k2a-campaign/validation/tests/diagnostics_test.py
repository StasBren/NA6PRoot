#!/usr/bin/env python3
"""SYNTHETIC test of k2a_residual_diagnostics.py (part J of test_validation_synthetic.sh), on the run's outputs and on a
COPY in which one cloud is marked R_EFF_ABOVE_GRID (blocked-cloud path). NOT physics.

  diagnostics_test.py VALIDATION_DIR VAL_OUTDIR NEW_WORK_DIR
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

H, O, W = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3])
W.mkdir()
D = str(H / "k2a_residual_diagnostics.py")
p = subprocess.run([sys.executable, D, str(W / "diag"), "--from", str(O), "--plots"], capture_output=True, text=True)
assert p.returncode == 0, p.stderr[-800:]
j = json.load(open(W / "diag/residual_diagnostics.json")); m = j["muons"]["0002"]
assert m["time_grid"]["status"] == "PASS" and m["time_grid"]["bins_compared"] == 700
an = m["residual"]["model"]["ion_only"]["U"]
cmp_ = json.load(open(O / "muons/muon_0002/cmp_S0T1_vs_K2A.json"))["terms"]["ion_only"]["U"]["unshaped_L2_rel"]
assert abs(an["L2_rel"] - cmp_) <= 1e-9 * max(1.0, cmp_), (an["L2_rel"], cmp_)       # same number as the run's comparison
assert abs(sum(an["window_share"].values()) - 1) < 1e-9 and abs(sum(an["strip_rank_share"].values()) - 1) < 1e-9
assert an["L2_after_global_scale"] <= an["L2_rel"] + 1e-15 and an["L2_after_per_strip_scale"] <= an["L2_after_global_scale"] + 1e-15
assert an["L2_after_time_shift"] <= an["L2_rel"] + 1e-15
assert set(m["residual"]) == {"model", "total", "timing", "imposed"}
assert m["residual"]["imposed"]["ion_only"]["U"]["L2_rel"] == an["L2_rel"]          # agreeing grouping: identical
fl = m["s0_floor"]["S0T1"]; assert fl["groups"] == 2 and fl["max_pilot_L2"]["U_10000ns"] > 0
tr = m["truncation"]["S0T1"]["U"]
assert tr["charge_beyond_frac"] is not None and tr["shaped_Tp100_tail_effect_over_lead_peak"] >= 0
assert "blocked" not in m and (W / "diag/muon_0002_residual.png").is_file() and (W / "diag/runtime_scaling.png").is_file()
md = (W / "diag/residual_diagnostics.md").read_text()
for x in ("Residual anatomy", "beyond the K2A window", "Blocked clouds", "Runtime scaling"):
    assert x in md, x
# blocked-cloud path on a copy
R = W / "run_ood"; (R / "muons").mkdir(parents=True)
shutil.copy(O / "provenance.txt", R)
V = R / "muons/muon_0002"; shutil.copytree(O / "muons/muon_0002", V, symlinks=True)
rows = list(csv.DictReader(open(V / "K2A_prod/clouds.csv")))
rows[3]["status"] = "8"; rows[3]["status_names"] = "R_EFF_ABOVE_GRID"
(V / "K2A_prod/clouds.csv").unlink()
with open(V / "K2A_prod/clouds.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
res = json.load(open(V / "K2A_prod/result.json")); res.update(usable=False, status_names="R_EFF_ABOVE_GRID", blocked_charge_fC=float(rows[3]["N_fC"]))
(V / "K2A_prod/result.json").unlink(); (V / "K2A_prod/result.json").write_text(json.dumps(res))
p = subprocess.run([sys.executable, D, str(W / "diag_ood"), "--from", str(R)], capture_output=True, text=True)
assert p.returncode == 0, p.stderr[-800:]
b = json.load(open(W / "diag_ood/residual_diagnostics.json"))["muons"]["0002"]["blocked"]
c = b["clouds"]
assert len(c) == 1 and c[0]["primary"] == int(rows[3]["primary"]) and c[0]["status"] == "R_EFF_ABOVE_GRID" \
    and abs(c[0]["charge_frac"] - float(rows[3]["N_fC"]) / res["cloud_charge_fC"]) < 1e-12 and c[0]["r_eff_over_grid_max"] > 0 \
    and b["kernel_grid"]["r_nodes"] > 0 and c[0]["min_clearance_um"] is not None, b
print("PASS residual diagnostics (SYNTHETIC): model L2 equals the run's comparison; window and strip-rank shares sum to 1; "
      "scale / shift removal never increases L2; S0 late-step pilot floor, truncation and runtime tables present; plots written; "
      "a blocked cloud is listed with primary, charge fraction, r_eff vs the kernel radial grid and clearance")

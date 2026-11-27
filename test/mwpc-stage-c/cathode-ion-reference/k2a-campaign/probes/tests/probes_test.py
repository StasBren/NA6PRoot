#!/usr/bin/env python3
"""SYNTHETIC test of the read-only probes (part K of validation/tests/test_validation_synthetic.sh). NOT physics.

  probes_test.py PROBES_DIR TOY_RUN TAG SCAN_WORK_DIR VAL_OUTDIR NEW_WORK_DIR

TOY_RUN: the S-ladder toy run (muon_0001 with S0, G_TAG ... K2Akern_TAG; kernel_TAG; green_TAG/axes.json -> grid);
SCAN_WORK_DIR: the work directory of tests/domain_scan_test.py (pf/, toy/kern, scan/); VAL_OUTDIR: a pruned validation run.
"""
import json
import subprocess
import sys
from pathlib import Path

P, RUN, TAG, SW, VAL, W = (Path(sys.argv[1]).resolve(), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), Path(sys.argv[5]), Path(sys.argv[6]))
W.mkdir()
py = lambda *a: subprocess.run([sys.executable] + [str(x) for x in a], capture_output=True, text=True)

# closure anatomy on the existing S-ladder variants
p = py(P / "k2a_closure_anatomy.py", W / "ca", "--sladder", RUN, "--tag", TAG)
assert p.returncode == 0, p.stderr[-600:]
j = json.load(open(W / "ca/closure_anatomy.json"))["muons"]["muon_0001"]
assert j["time_grid"] == "PASS" and {"G_vs_S0", "SR_vs_G", "S3A_vs_G", "SW_vs_G", "K2A_vs_G", "K2Akern_vs_K2A"} <= set(j["pairs"]), sorted(j["pairs"])
for s in "UV":
    assert j["pairs"]["K2Akern_vs_K2A"][s]["L2_rel"] < 1e-6                # harmonic kernel == its synthesis (identity)
    assert abs(sum(j["pairs"]["SR_vs_G"][s]["window_share"].values()) - 1) < 1e-9
    assert list(j["pairs"]["G_vs_S0"][s]["window_share"])[:3] == ["0-10", "10-20", "20-50"]

# kernel probe (tau / spatial interpolation, wire translation, time horizon)
grid = Path(json.load(open(RUN / f"green_{TAG}/axes.json"))["grid"])
p = py(P / "k2a_kernel_probe.py", W / "kp", "--kernel", RUN / f"kernel_{TAG}", "--green-paths", grid)
assert p.returncode == 0, p.stderr[-600:]
k = json.load(open(W / "kp/kernel_probe.json"))
assert k["tau_interp_first10bins_rel_to_max_bin"]["U"]["max"] < 0.05
k0 = [r for r in k["wire_translation"] if r["basis"] == "K0"]
assert k0 and all(r["predicted_sign_early_0_50ns_rel_L2"] < r["opposite_sign_early_0_50ns_rel_L2"] for r in k0), k0[:2]
th = k["time_horizon"]
assert th["start_points"] > 0 and th["drift_time_ns_quantiles"]["q1"] > 0
assert th["all_cathode"] and th["classes"]["CATHODE"] == th["start_points"] and th["endpoint_dv_mm"]["max_abs"] <= 0.01, th["classes"]
# fault injection on a COPY of the grid: one path stops short of the cathode, one ends at the lateral area edge, one
# has a vertex record inconsistent with its endpoint; the probe must classify each, not pass them as CATHODE
import csv, shutil
bad = W / "bad_grid"; bad.mkdir()
for f in ("ions.csv", "complete.txt"):
    shutil.copy(grid / f, bad / f)
(bad / "vertices.f64").symlink_to((grid / "vertices.f64").resolve())
rows = list(csv.DictReader(open(grid / "grid.csv"))); cols = list(rows[0])
rows[0]["end_v_mm"] = str(float(rows[0]["end_v_mm"].replace("np.float64(", "").rstrip(")")) * 0.9)
rows[1]["end_u_mm"] = "13.5"
rows[2]["end_ns"] = str(float(rows[2]["end_ns"].replace("np.float64(", "").rstrip(")")) + 50.0)
with open(bad / "grid.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols); wr.writeheader(); wr.writerows(rows)
p = py(P / "k2a_kernel_probe.py", W / "kp_bad", "--kernel", RUN / f"kernel_{TAG}", "--green-paths", bad)
assert p.returncode == 0, p.stderr[-600:]
tb = json.load(open(W / "kp_bad/kernel_probe.json"))["time_horizon"]
cl = {e["row"]: e["class"] for e in tb["examples_not_cathode"]}
assert not tb["all_cathode"] and tb["classes"]["CATHODE"] == th["start_points"] - 3, tb["classes"]
assert cl[int(rows[1]["row"])] == "LATERAL" and cl[int(rows[2]["row"])] == "INCONSISTENT" and cl[int(rows[0]["row"])] in ("INCONSISTENT", "OTHER"), cl
assert "NOT every path" in (W / "kp_bad/kernel_probe.md").read_text()
(bad / "vertices.f64").unlink()
p = py(P / "k2a_kernel_probe.py", W / "kp_novert", "--kernel", RUN / f"kernel_{TAG}", "--green-paths", bad)
assert p.returncode == 0 and "NOT DONE: missing vertices.f64" in json.load(open(W / "kp_novert/kernel_probe.json"))["time_horizon"]["cathode_check"]

# wire probe on the toy preflight (slot 0015: a cloud on a wire outside the toy kernel)
p = py(P / "k2a_wire_probe.py", W / "wp", "--pf", SW / "pf", "--kernel", SW / "toy/kern", "--slots", "15,10", "--scan", SW / "scan")
assert p.returncode == 0, p.stderr[-600:]
w = json.load(open(W / "wp/wire_probe.json"))["events"]
assert w["0015"]["charge_frac_outside_kernel_wires"] > 0 and not w["0015"]["wires"]["3"]["in_kernel"], w["0015"]
assert w["0010"]["charge_frac_outside_kernel_wires"] == 0 and w["0015"]["scan"]["execution"] == "BLOCKED"
# a --scan directory without scan.csv is reported and skipped, never required
p = py(P / "k2a_wire_probe.py", W / "wp_noscan", "--pf", SW / "pf", "--kernel", SW / "toy/kern", "--slots", "15", "--scan", W / "no_such_scan")
assert p.returncode == 0 and "absent" in p.stdout, p.stderr[-600:]
wn = json.load(open(W / "wp_noscan/wire_probe.json"))
assert "absent" in wn["scan"] and wn["events"]["0015"]["scan"] is None
p = py(P / "k2a_wire_probe.py", W / "wp_sealed", "--pf", SW / "pf", "--kernel", SW / "toy/kern", "--slots", "18")
assert p.returncode != 0 and "sealed" in p.stdout + p.stderr

# group synthesis refuses pruned exports
p = subprocess.run(["bash", str(P / "k2a_group_synth.sh"), str(VAL / "muons/muon_0002"), str(RUN), TAG, str(W / "gs")], capture_output=True, text=True)
assert p.returncode != 0 and "pruned" in p.stderr, p.stderr[-300:]
print("PASS probes (SYNTHETIC): closure anatomy on existing S-ladder variants (harmonic kernel == synthesis, early windows); "
      "kernel probe (tau and spatial interpolation estimates, wire-translation sign convention, drift-time horizon, Green-path "
      "termination on the physical cathode with injected short / lateral / inconsistent paths classified); wire probe "
      "(charge on wires outside the kernel, missing scan.csv tolerated, sealed refused); group synthesis refuses pruned exports")

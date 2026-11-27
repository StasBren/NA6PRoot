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
assert k["time_horizon"]["start_points"] > 0 and k["time_horizon"]["drift_time_ns_quantiles"][str(1.0)] > 0

# wire probe on the toy preflight (slot 0015: a cloud on a wire outside the toy kernel)
p = py(P / "k2a_wire_probe.py", W / "wp", "--pf", SW / "pf", "--kernel", SW / "toy/kern", "--slots", "15,10", "--scan", SW / "scan")
assert p.returncode == 0, p.stderr[-600:]
w = json.load(open(W / "wp/wire_probe.json"))["events"]
assert w["0015"]["charge_frac_outside_kernel_wires"] > 0 and not w["0015"]["wires"]["3"]["in_kernel"], w["0015"]
assert w["0010"]["charge_frac_outside_kernel_wires"] == 0 and w["0015"]["scan"]["execution"] == "BLOCKED"
p = py(P / "k2a_wire_probe.py", W / "wp_sealed", "--pf", SW / "pf", "--kernel", SW / "toy/kern", "--slots", "18")
assert p.returncode != 0 and "sealed" in p.stdout + p.stderr

# group synthesis refuses pruned exports
p = subprocess.run(["bash", str(P / "k2a_group_synth.sh"), str(VAL / "muons/muon_0002"), str(RUN), TAG, str(W / "gs")], capture_output=True, text=True)
assert p.returncode != 0 and "pruned" in p.stderr, p.stderr[-300:]
print("PASS probes (SYNTHETIC): closure anatomy on existing S-ladder variants (harmonic kernel == synthesis, early windows); "
      "kernel probe (tau and spatial interpolation estimates, wire-translation sign convention, drift-time horizon); wire probe "
      "(charge on wires outside the kernel, sealed refused); group synthesis refuses pruned exports")

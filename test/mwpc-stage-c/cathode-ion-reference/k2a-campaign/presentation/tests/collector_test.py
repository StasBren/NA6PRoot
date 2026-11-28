#!/usr/bin/env python3
"""SYNTHETIC test of collect_presentation_inputs.py (no real data).

  python3 collector_test.py WORK_DIR     (WORK_DIR must not exist)

Builds a fake validation tree and checks: both preflight naming conventions with relative paths preserved (no
filename collisions); sealed slots never copied and a report with a sealed slot refused before anything is written;
inventory missing or without sealed rows refused; strict kernel slice (correct basis -> slice with axes; missing basis
or wrong shape -> reported, nothing substituted); missing essential metadata -> INCOMPLETE (exit 2) and listed;
in-progress slots not required; upload folder inside an input refused; --also-report and run-dir discovery without
combined_from.txt; and every input file unchanged (sha256 and mtime) after all runs.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

TOOL = Path(__file__).resolve().parents[1] / "collect_presentation_inputs.py"
W = Path(sys.argv[1]).resolve()
assert not W.exists()
SRC = W / "src"


def w(p, t):
    p = SRC / p; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(t)


def per_muon(path, rows):
    p = SRC / path; p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["slot", "execution", "accuracy_v1_primary", "model_U_shaped_L2_Tp100"])
        wr.writeheader(); wr.writerows(rows)


def run(*args):
    return subprocess.run([sys.executable, str(TOOL), *map(str, args)], capture_output=True, text=True)


VAL = SRC / "val"
slots = {"0110": ("VALID", "PASS", "0.0017"), "0637": ("VALID", "PASS", "0.0026"), "0548": ("VALID", "PASS", "0.0018"),
         "0790": ("FAIL", "", ""), "0801": ("INCONCLUSIVE", "", "")}
for s, (ex, acc, l2) in slots.items():
    R = f"val/overnight_800/muon_{s}"
    w(f"{R}/slots.txt", s + "\n"); w(f"{R}/provenance.txt", "git_commit=abc\nkernel=/k/kernel_r12t72\n")
    if ex == "INCONCLUSIVE":
        continue
    V = f"{R}/muons/muon_{s}"
    w(f"{V}/outcome.json", "{}"); w(f"{V}/timing.tsv", "label\twall_s\n")
    w(f"{V}/steps.tsv", "births_moments\tOK\tx\t\n" + ("moments_archive_vs_export\tFAIL\tt_rms\tp\n" if s == "0790" else ""))
    w(f"{V}/cmp_S0T1_vs_K2A.json", "{}"); w(f"{V}/shaped_refS0T1_ionelec/per_muon_Tp100.jsonl", "{}\n")
    w(f"{V}/prod/compare_moments.md", "| t_rms_ns | abs | 6.24e-07 |\n"); w(f"{V}/far_dir.txt", str(SRC / f"far/{s}") + "\n")
    w(f"{V}/shaped_layout/muon_{s}/S0T1/replay/waveforms.csv", "t,U_1\n0,1\n")
    w(f"{V}/K2A_prod/replay/waveforms.csv", "t,U_1\n0,1\n")
    w(f"far/{s}/groups/groups.json", "{}"); w(f"far/{s}/groups/g_0/vertices.f64", "x" * 64)
    w(f"far/{s}/S0groups/g_0/timing_attempts.json", "{}"); w(f"far/{s}/S0_local/replay/waveforms.csv", "t,U_1\n0,1\n")
per_muon("val/analysis_new/per_muon.csv", [dict(slot=s, execution=e, accuracy_v1_primary=a, model_U_shaped_L2_Tp100=l)
                                           for s, (e, a, l) in slots.items()])
w("val/analysis_new/validation_report.md", "# new\n")
per_muon("val/analysis_snapshot/per_muon.csv", [dict(slot="0110", execution="VALID", accuracy_v1_primary="PASS",
                                                     model_U_shaped_L2_Tp100="0.0017")])
w("val/kernel_probe/kernel_probe.md", "k\n"); w("val/overnight_800/queue.tsv", "slot\tstate\n0110\tdone\n")
inv = "slot,role,channel,station\n" + "".join(f"{s},pool,phi,1\n" for s in slots) + "0999,sealed_validation,phi,2\n"
w("pf/inventory/slots.csv", inv); w("pf/domain/per_event.csv", "slot\n0110\n")
# naming A (p3_birth_moments): <id>.moments.csv(.event.json/.frame.json), two batch dirs with the SAME file name
for d, ss in (("b1", ("0110", "0637", "0790")), ("b2", ("0110", "0999", "0801"))):
    for s in ss:
        w(f"pf/moments_nonsealed/{d}/muon_{s}.moments.csv", f"{d}\n")
        w(f"pf/moments_nonsealed/{d}/muon_{s}.moments.csv.event.json", "{}")
        w(f"pf/moments_nonsealed/{d}/muon_{s}.moments.csv.frame.json", "{}")
    w(f"pf/moments_nonsealed/{d}/batch_timing.json", "{}")
# naming B: <id>.event.json / <id>.frame.json
w("pf/moments_nonsealed/b3/0548.moments.csv", "b3\n"); w("pf/moments_nonsealed/b3/0548.event.json", "{}")
w("pf/moments_nonsealed/b3/0548.frame.json", "{}")
axes = {"wires": [-1, 0, 1], "r_um": [20, 40, 80, 160], "D_mm": [-0.17, 0, 0.17, 0.34, 0.51], "tau_ns": list(range(6))}
basis = ["K0", "Kc1", "Ks1", "Kc2", "Ks2", "KA1", "KA2"]
for name, b, shape in (("kern", basis, (7, 2, 3, 4, 5, 6)), ("kern_badbasis", ["B0", "B1", "B2"], (3, 2, 3, 4, 5, 6)),
                       ("kern_badshape", basis, (7, 2, 3, 4, 4, 6))):
    k = SRC / name; k.mkdir(parents=True)
    K = np.arange(np.prod(shape), dtype=np.float32).reshape(shape); np.save(k / "K.npy", K)
    (k / "kernel.json").write_text(json.dumps({"basis": b, "green_axes": axes}))


def snapshot():
    return {str(p): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns) for p in SRC.rglob("*") if p.is_file()}


before = snapshot()
U = W / "up"
common = ["--analysis", VAL / "analysis_new", "--pf", SRC / "pf", "--val-root", VAL, "--queue", VAL / "overnight_800/queue.tsv",
          "--also-report", VAL / "analysis_snapshot"]
p = run(U, *common, "--kernel", SRC / "kern")
assert p.returncode == 0, p.stdout + p.stderr
summ = (U / "SUMMARY.txt").read_text()
assert summ.startswith("status COMPLETE") and "INCONCLUSIVE 1" in summ, summ
man = (U / "MANIFEST.tsv").read_text()
assert "0999" not in man.replace("[1 sealed rows removed]", ""), "sealed slot copied"
for d in ("b1", "b2"):
    assert (U / f"preflight/moments_nonsealed/{d}/muon_0110.moments.csv").read_text() == f"{d}\n"   # no collision
assert (U / "preflight/moments_nonsealed/b3/0548.event.json").is_file() and (U / "preflight/moments_nonsealed/b3/0548.frame.json").is_file()
assert (U / "preflight/moments_nonsealed/b1/muon_0637.moments.csv.frame.json").is_file()
assert "0999" not in (U / "preflight/inventory_slots_nonsealed.csv").read_text()
z = np.load(U / "kernel/K_slice_centralD.npz")
K = np.load(SRC / "kern/K.npy")
assert list(z["basis"]) == ["K0", "Kc1", "KA2"] and np.array_equal(z["K"], K[[0, 1, 6]][:, :, :, :, 1, :]) and float(z["D_mm"]) == 0.0
assert not list((U / "main/muons/muon_0110/far").rglob("*.f64")) and (U / "main/muons/muon_0110/far/groups.json").is_file()
assert (U / "reports/analysis_snapshot/per_muon.csv").is_file() and (U / "queue/queue.tsv").is_file()
assert (U / "main/muons/muon_0790/prod/compare_moments.md").is_file()
ck = (U / "CHECKS.tsv").read_text()
assert "moments_archive_vs_export" in ck and "PASS\t" in ck
assert set(os.listdir(U / "examples")) == {"muon_0548", "muon_0637"}, os.listdir(U / "examples")
# kernel: wrong basis names / wrong shape -> reported, nothing substituted
for kd, msg in (("kern_badbasis", "lacks one of"), ("kern_badshape", "!= (basis")):
    p = run(W / f"up_{kd}", *common, "--kernel", SRC / kd, "--examples", "none")
    assert p.returncode == 0 and not (W / f"up_{kd}/kernel/K_slice_centralD.npz").exists(), p.stdout
    assert "NOT generated" in (W / f"up_{kd}/CHECKS.tsv").read_text() and msg in (W / f"up_{kd}/CHECKS.tsv").read_text()
# missing essential metadata -> INCOMPLETE, listed (temporarily hide one file by renaming a COPY tree)
import shutil
shutil.copytree(SRC / "pf", W / "pf_missing")
os.remove(W / "pf_missing/moments_nonsealed/b1/muon_0637.moments.csv.frame.json")
p = run(W / "up_missing", *[x if x != SRC / "pf" else W / "pf_missing" for x in common], "--examples", "none")
assert p.returncode == 2 and "status INCOMPLETE" in p.stdout, p.stdout + p.stderr
assert "0637\tpreflight frame file" in (W / "up_missing/MISSING_ESSENTIAL.tsv").read_text()
# refusals: nothing written
per_muon("val/analysis_sealed/per_muon.csv", [dict(slot="0999", execution="VALID", accuracy_v1_primary="PASS",
                                                   model_U_shaped_L2_Tp100="0.001")])
for out, args, msg in ((W / "r1", ["--analysis", VAL / "analysis_sealed", "--pf", SRC / "pf", "--val-root", VAL], "sealed"),
                       (VAL / "upload_inside", common, "inside the input"),
                       (W / "r3", ["--analysis", VAL / "analysis_new", "--pf", W / "nopf", "--val-root", VAL], "not found")):
    p = run(out, *args)
    assert p.returncode == 1 and msg in p.stderr and not out.exists(), (msg, p.stderr)
(W / "pf_nosealed/inventory").mkdir(parents=True)
(W / "pf_nosealed/inventory/slots.csv").write_text("slot,role\n0110,pool\n")
p = run(W / "r4", "--analysis", VAL / "analysis_new", "--pf", W / "pf_nosealed", "--val-root", VAL)
assert p.returncode == 1 and "no sealed rows" in p.stderr and not (W / "r4").exists()
after = snapshot()
changed = [k for k in before if before[k] != after.get(k)]
assert not changed, f"inputs modified: {changed}"
print("PASS collector (SYNTHETIC): both preflight namings with paths preserved; sealed never copied and refused; "
      "inventory checks; strict kernel slice; INCOMPLETE on missing metadata; in-progress slots not required; "
      "upload inside an input refused; --also-report; inputs unchanged (sha256 + mtime)")

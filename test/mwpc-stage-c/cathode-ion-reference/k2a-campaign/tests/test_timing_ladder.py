#!/usr/bin/env python3
"""SYNTHETIC test of k2a_timing_ladder.py (not NA60+ data).

  test_timing_ladder.py WORKDIR NUMW_DIR LUT GROUPS_DIR

A. retry logic against a STUB replay module (exit codes and messages exactly as the frozen replay prints them):
   pilot V L2 = 0.567671 % x (step / 100)  -> 100 ns rejected, 50 ns accepted (chosen 50, two attempts recorded)
   closure failure (RuntimeError, exit 1)   -> FAIL at once, no retry
   phase-range refusal                       -> REFUSED at once, no retry
   pilot rejects every step                  -> FAIL (exit 2), nothing chosen, no replay link
   invalid ladder (25 ns, increasing)        -> refused before running
B. the REAL frozen replay on a toy recentering group (GROUPS_DIR from k2a_group_export.py): ladder "100" -> chosen 100;
   the same group with ladder "50" -> chosen 50; the summed local S0 of {group A at 100, group B at 50} is accepted
   (identical output clock), the late step per group is recorded, and the waveform differs from all-100 only at the
   level of the frozen time compression (< the 0.5 % gate); a replay whose pilot did not pass is refused by the sum.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent.parent / "recenter"
LADDER = HERE / "k2a_timing_ladder.py"
STUB = r'''
import argparse, json, os, sys
from pathlib import Path
ap = argparse.ArgumentParser(); ap.add_argument("table"); ap.add_argument("export"); ap.add_argument("out", type=Path)
ap.add_argument("--half-strips", type=int, default=5); ap.add_argument("--late-step-ns", type=int, default=100)
a = ap.parse_args(); mode = os.environ.get("FAKE_MODE", "normal")
if not 10 <= a.late_step_ns <= 200 or a.late_step_ns % 10:
    raise ValueError("invalid strip or time controls")
if mode == "refuse":
    raise ValueError("query outside declared strip phase range")
a.out.mkdir(parents=True, exist_ok=False)
if mode == "closure":
    raise RuntimeError("charge conservation failure: 1e-3, 0")
v = 0.00567671 * a.late_step_ns / 100 * (10 if mode == "allfail" else 1)
errs = {"U_100ns": 0.0, "U_100000ns": v / 3, "V_10000ns": 0.00162685, "V_100000ns": v}
(a.out / "summary.json").write_text(json.dumps({"rows": 256, "exact_ions": 1000, "late_step_ns": a.late_step_ns, "dt_ns": 10.0,
    "subset_relative_L2": errs, "time_grid_pilot_pass": max(errs.values()) < .005, "wall_s": 4.35, "queries": 6830000,
    "window_closure_fC": 1e-12, "endpoint_closure_fC": 1e-12}))
if max(errs.values()) >= .005:
    print("REVIEW: refine late-step-ns in a new output directory; preserve paths."); sys.exit(2)
print("PASS complete weighted muon replay timing pilot"); sys.exit(0)
'''


def ladder(out, numw, steps, export="E", env=None, lut="L"):
    import os
    e = dict(os.environ, **(env or {}))
    p = subprocess.run([sys.executable, str(LADDER), str(export), str(out), "--lut", str(lut), "--numw", str(numw),
                        "--steps", steps, "--link-replay"], capture_output=True, text=True, env=e)
    rec = json.loads((out / "timing_attempts.json").read_text()) if (out / "timing_attempts.json").exists() else None
    return p.returncode, rec, p.stdout + p.stderr


def main():
    W, NUMW, LUT, GROUPS = (Path(x) for x in sys.argv[1:5])
    W.mkdir(parents=True, exist_ok=False)
    stub = W / "stub_numw"; stub.mkdir(); (stub / "p3_replay_full_spectral_muon.py").write_text(STUB)
    (W / "E").mkdir(); (W / "L").write_text("x")
    rc, rec, _ = ladder(W / "a", stub, "100,50,20,10", W / "E", lut=W / "L")
    assert rc == 0 and rec["chosen_late_step_ns"] == 50 and [x["outcome"] for x in rec["attempts"]] == ["TIME_COMPRESSION_REJECTED", "PASS"], rec
    assert abs(rec["attempts"][0]["pilot_relative_L2"]["V_100000ns"] - 0.00567671) < 1e-12
    assert (W / "a" / "replay").is_symlink() and (W / "a" / "replay" / "summary.json").exists()
    rc, rec, _ = ladder(W / "b", stub, "100,50", W / "E", {"FAKE_MODE": "closure"}, lut=W / "L")
    assert rc == 1 and rec["status"] == "FAIL" and len(rec["attempts"]) == 1, rec
    rc, rec, _ = ladder(W / "c", stub, "100,50", W / "E", {"FAKE_MODE": "refuse"}, lut=W / "L")
    assert rc == 3 and rec["status"] == "REFUSED" and len(rec["attempts"]) == 1, rec
    rc, rec, _ = ladder(W / "d", stub, "100,50,20,10", W / "E", {"FAKE_MODE": "allfail"}, lut=W / "L")
    assert rc == 2 and rec["status"] == "FAIL" and rec["chosen_late_step_ns"] is None and len(rec["attempts"]) == 4
    assert not (W / "d" / "replay").exists()
    for bad in ("100,25", "50,100", "100,100"):
        rc, rec, out = ladder(W / f"x{bad.replace(',', '_')}", stub, bad, W / "E", lut=W / "L")
        assert rc != 0 and rec is None and "steps must be" in out, (bad, rc, out)
    print("PASS ladder logic (stub replay): retry only on the time-compression rejection (100 -> 50 chosen, metrics "
          "recorded); closure failure and refusal stop at once; all-rejected gives FAIL with nothing chosen; invalid "
          "ladders (25 ns, increasing, repeated) refused")
    # B. real frozen replay on toy groups
    gs = sorted(p for p in GROUPS.glob("g_*") if p.is_dir())
    assert len(gs) >= 2, "need a toy grouping with at least two groups"
    TA = Path(__file__).resolve().parents[2] / "timing-ablation"
    S100, MIX, BAD = W / "S100", W / "MIX", W / "BAD"
    for root, steps_for in ((S100, lambda i: "100"), (MIX, lambda i: "100" if i == 0 else "50")):
        for i, g in enumerate(gs):
            rc, rec, out = ladder(root / g.name, NUMW, steps_for(i), g, lut=LUT)
            assert rc == 0, out
            subprocess.run([sys.executable, str(TA / "t_prompt_electron_term.py"), str(LUT), str(g), str(root / g.name / "prompt"),
                            "--replay-dir", str(NUMW), "--half-strips", "5"], check=True, stdout=subprocess.DEVNULL)
    SUM = HERE / "k2a_local_sum.py"
    subprocess.run([sys.executable, str(SUM), str(GROUPS), str(S100), str(W / "sum100")], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, str(SUM), str(GROUPS), str(MIX), str(W / "summix")], check=True, stdout=subprocess.DEVNULL)
    js = json.loads((W / "summix" / "local_sum.json").read_text())
    assert sorted(js["late_step_ns_per_group"].values()) == [50, 100], js["late_step_ns_per_group"]
    a = np.loadtxt(W / "sum100" / "replay" / "waveforms.csv", delimiter=",", skiprows=1)
    b = np.loadtxt(W / "summix" / "replay" / "waveforms.csv", delimiter=",", skiprows=1)
    assert np.array_equal(a[:, 0], b[:, 0])
    rel = float(np.linalg.norm(a[:, 1:] - b[:, 1:]) / np.linalg.norm(a[:, 1:]))
    assert rel < 0.005, rel
    # a replay whose pilot failed must not be summed
    shutil.copytree(MIX, BAD, symlinks=True)
    sj = BAD / gs[1].name / "late50" / "replay" / "summary.json"
    d = json.loads(sj.read_text()); d["time_grid_pilot_pass"] = False; sj.write_text(json.dumps(d))
    p = subprocess.run([sys.executable, str(SUM), str(GROUPS), str(BAD), str(W / "sumbad")], capture_output=True, text=True)
    assert p.returncode != 0 and "time-compression pilot did not pass" in (p.stdout + p.stderr), p.stdout + p.stderr
    print(f"PASS real frozen replay (toy groups): ladder chose 100 and 50 ns; groups at different late steps summed on the "
          f"identical 10 ns output clock (late steps {js['late_step_ns_per_group']}); all-100 vs mixed rel L2 {rel:.2e} "
          "(< 0.5 % gate); a pilot-rejected replay is refused by the sum")


if __name__ == "__main__":
    main()

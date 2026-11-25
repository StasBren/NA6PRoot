#!/usr/bin/env python3
"""Frozen replay of ONE (sub-)export with a late-time-step convergence ladder. Frozen code, LUT and gate unchanged.

  k2a_timing_ladder.py EXPORT OUT_DIR --lut LUT.npz --numw NUMW_DIR [--steps 100,50,20,10] [--half-strips 5]
                       [--link-replay]

The frozen p3_replay_full_spectral_muon.py exposes --late-step-ns (integer ns, multiple of 10, 10..200: between 1000 ns
after birth and 1000 ns before collection it samples the RKF path / LUT every late-step on the fixed 10 ns event clock;
the output waveform grid is ALWAYS the same 10 000 x 10 ns clock). Its own pilot compares that compression with the
native 10 ns clock on 32 weighted rows; if any window's L2 >= 0.5 % it writes its outputs, prints "REVIEW: refine
late-step-ns in a new output directory; preserve paths." and exits 2. This tool does exactly that, nothing more:
  attempt k: frozen replay with --late-step-ns steps[k] into OUT_DIR/late<step>/replay (new directory each time)
  retry ONLY if the attempt is this time-compression rejection: exit 2 + the REVIEW line + summary.json with
        time_grid_pilot_pass false and charge/endpoint closure intact (the frozen code raises on closure failure);
  any other outcome stops: exit 0 + pilot pass -> chosen; window/phase/cross-section refusal -> REFUSED; else FAIL.
--link-replay makes OUT_DIR/replay a relative symlink to the chosen attempt (no data copied or moved).
Writes OUT_DIR/timing_attempts.json: every attempt (step, exit code, pilot metrics per window, pass flag, wall time,
queries, census/closure numbers, outcome) and the final status, chosen step and output grid (dt, bins).
Exit: 0 chosen, 1 FAIL, 2 all steps rejected by the timing pilot (FAIL, nothing chosen), 3 REFUSED.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REFUSAL = re.compile(r"selected channels omit the landing strip|query outside declared strip phase range|"
                     r"query outside validated cross-section")
REVIEW = "REVIEW: refine late-step-ns in a new output directory; preserve paths."


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=Path); ap.add_argument("out", type=Path)
    ap.add_argument("--lut", type=Path, required=True); ap.add_argument("--numw", type=Path, required=True)
    ap.add_argument("--steps", default="100,50,20,10"); ap.add_argument("--half-strips", type=int, default=5)
    ap.add_argument("--link-replay", action="store_true")
    a = ap.parse_args()
    steps = [int(x) for x in a.steps.split(",")]
    if any(s % 10 or not 10 <= s <= 200 for s in steps) or steps != sorted(steps, reverse=True) or len(set(steps)) != len(steps):
        raise SystemExit("steps must be distinct, decreasing multiples of 10 in [10, 200] (the frozen replay's own limits)")
    if (a.out / "timing_attempts.json").exists() or (a.out / "replay").exists():
        raise SystemExit(f"refusing: {a.out} already has a ladder result")
    a.out.mkdir(parents=True, exist_ok=True)
    rec = {"tool": "k2a_timing_ladder.py", "export": str(a.export), "lut": str(a.lut), "numw": str(a.numw),
           "half_strips": a.half_strips, "steps_allowed": steps, "gate": "frozen: every subset window L2 < 0.5 %",
           "attempts": [], "status": None, "chosen_late_step_ns": None}
    code = 1
    for s in steps:
        d = a.out / f"late{s}"
        if d.exists():
            raise SystemExit(f"refusing: {d} exists")
        d.mkdir()
        cmd = [sys.executable, "p3_replay_full_spectral_muon.py", str(a.lut.resolve()), str(a.export.resolve()),
               str((d / "replay").resolve()), "--half-strips", str(a.half_strips), "--late-step-ns", str(s)]
        p = subprocess.run(cmd, cwd=a.numw, capture_output=True, text=True)
        (d / "replay.log").write_text(p.stdout + p.stderr)
        att = {"late_step_ns": s, "exit_code": p.returncode, "dir": str(d / "replay")}
        sj = d / "replay" / "summary.json"
        if sj.exists():
            S = json.loads(sj.read_text())
            att.update(pilot_relative_L2=S.get("subset_relative_L2"), pilot_pass=S.get("time_grid_pilot_pass"),
                       wall_s=S.get("wall_s"), queries=S.get("queries"), rows=S.get("rows"), dt_ns=S.get("dt_ns"),
                       window_closure_fC=S.get("window_closure_fC"), endpoint_closure_fC=S.get("endpoint_closure_fC"),
                       worst_window=max(S["subset_relative_L2"], key=S["subset_relative_L2"].get) if S.get("subset_relative_L2") else None)
        out = p.stdout + p.stderr
        if p.returncode == 0 and att.get("pilot_pass") is True:
            att["outcome"] = "PASS"; rec["attempts"].append(att)
            rec["status"] = "OK"; rec["chosen_late_step_ns"] = s; rec["chosen_dir"] = att["dir"]; code = 0
            break
        if (p.returncode == 2 and REVIEW in out and att.get("pilot_pass") is False
                and att.get("window_closure_fC", 1) <= 1e-7 and att.get("endpoint_closure_fC", 1) <= 1e-7):
            att["outcome"] = "TIME_COMPRESSION_REJECTED"; rec["attempts"].append(att)
            code = 2
            continue                                  # the one failure that is retried, with the next smaller step
        if REFUSAL.search(out):
            att["outcome"] = "REFUSED"; att["reason"] = REFUSAL.search(out).group(0)
            rec["attempts"].append(att); rec["status"] = "REFUSED"; code = 3
            break
        att["outcome"] = "FAIL"; att["log_tail"] = out.strip().splitlines()[-3:] if out.strip() else []
        rec["attempts"].append(att); rec["status"] = "FAIL"; code = 1
        break
    if rec["status"] is None:                         # every step rejected by the timing pilot
        rec["status"] = "FAIL"; rec["reason"] = "every late step rejected by the frozen time-compression pilot"; code = 2
    if code == 0:
        S = json.loads((Path(rec["chosen_dir"]) / "summary.json").read_text())
        rec["output_grid"] = {"dt_ns": S["dt_ns"], "note": "frozen fixed event clock (NBINS x DT), independent of late step"}
        if a.link_replay:
            (a.out / "replay").symlink_to(Path(f"late{rec['chosen_late_step_ns']}") / "replay")
    (a.out / "timing_attempts.json").write_text(json.dumps(rec, indent=1) + "\n")
    tried = ", ".join(f"{x['late_step_ns']} ns {x['outcome']}"
                      + (f" (worst {x['worst_window']} {100 * x['pilot_relative_L2'][x['worst_window']]:.4g} %)" if x.get("worst_window") else "")
                      for x in rec["attempts"])
    print(f"{rec['status']} timing ladder {a.export.name}: {tried}"
          + (f"; chosen late step {rec['chosen_late_step_ns']} ns" if code == 0 else ""))
    return code


if __name__ == "__main__":
    raise SystemExit(main())

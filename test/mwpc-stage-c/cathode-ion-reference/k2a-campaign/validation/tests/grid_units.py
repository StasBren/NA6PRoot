#!/usr/bin/env python3
"""Time-grid regression (part H of test_validation_synthetic.sh): k2a_time_grid.py and the grid check inside
k2a_label_union.py must FAIL explicitly on any cross-variant time-axis mismatch; and the real SYNTHETIC run's variants
(S0 frozen clock 10 000 bins, K2A 700 bins) share one grid on the K2A window.

  grid_units.py VALIDATION_DIR NEW_WORK_DIR [VAL_OUTDIR]
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

H, U = Path(sys.argv[1]), Path(sys.argv[2]); U.mkdir()
TG, LU = str(H / "k2a_time_grid.py"), str(H / "k2a_label_union.py")


def wave(d, times, labels=("U_0", "V_0"), prompt_times=None):
    for f, ts in (("replay/waveforms.csv", times), ("prompt/prompt_waveforms.csv", prompt_times or times)):
        (d / f).parent.mkdir(parents=True, exist_ok=True)
        with open(d / f, "w") as fo:
            fo.write(",".join(["time_center_ns", *labels]) + "\n")
            for k, t in enumerate(ts):
                fo.write(",".join([repr(t)] + [repr(0.1234567890123 * (i + 1) * (k + 1)) for i in range(len(labels))]) + "\n")
    return d


def grid(name, window, *variants):
    out = U / f"{name}.json"
    rc = subprocess.run([sys.executable, TG, str(out), "--window-ns", str(window)] + [f"{n}={d}" for n, d in variants],
                        capture_output=True, text=True).returncode
    return rc, json.load(open(out))


def union(name, window, *variants):
    o = U / f"lu_{name}"
    rc = subprocess.run([sys.executable, LU, str(o), "--window-ns", str(window)] + [f"{n}={d}" for n, d in variants],
                        capture_output=True, text=True).returncode
    return rc, o


ref = wave(U / "ref", [5.0 + 10 * k for k in range(5)])                  # "S0": 5 bins
k2a = wave(U / "k2a", [5.0 + 10 * k for k in range(3)])                  # "K2A": the 30 ns window
rc, j = grid("ok", 30, ("S0", ref), ("K2A", k2a))
assert rc == 0 and j["status"] == "PASS" and j["bins_compared"] == 3 and j["variants"]["S0"]["ion_abs_current_beyond_window_frac"] > 0, j
rc, o = union("ok", 30, ("S0", ref), ("K2A", k2a))
assert rc == 0
rs = list(csv.reader(open(o / "S0/replay/waveforms.csv"))); rk = list(csv.reader(open(o / "K2A/replay/waveforms.csv")))
src = list(csv.reader(open(ref / "replay/waveforms.csv")))
assert len(rs) == len(rk) == 4 and [r[0] for r in rs] == [r[0] for r in rk] and rs == src[:4], (len(rs), len(rk))
assert (o / "K2A").is_symlink() and not (o / "S0").is_symlink()
cases = {
    "dt": wave(U / "dt", [5.0 + 20 * k for k in range(3)]),                                   # 20 ns sampling
    "offset": wave(U / "offset", [6.0 + 10 * k for k in range(3)]),                            # first bin centre 6 ns
    "short": wave(U / "short", [5.0, 15.0]),                                                    # does not cover 30 ns
    "value": wave(U / "value", [5.0, 15.0, 25.000001]),                                         # one bin moved
    "nonuniform": wave(U / "nonuniform", [5.0, 15.0, 25.0, 36.0]),                              # non-uniform spacing
    "prompt_axis": wave(U / "prompt_axis", [5.0, 15.0, 25.0], prompt_times=[5.0, 15.0, 26.0]),  # replay != prompt axis
}
for n, d in cases.items():
    rc, j = grid(n, 30, ("S0", ref), (n, d))
    assert rc == 1 and j["status"] == "FAIL" and j["failures"], (n, j)
    rc, o = union(n, 30, ("S0", ref), (n, d))
    assert rc == 1 and not o.exists(), n                                # nothing written for t_shaped_compare
rc, j = grid("window", 25, ("S0", ref), ("K2A", k2a))                  # 25 ns is not an integer number of 10 ns bins
assert rc == 1 and any("integer number" in f for f in j["failures"]), j
rc, j = grid("missing", 30, ("S0", ref), ("none", U / "absent"))
assert rc == 1 and any("missing" in f for f in j["failures"]), j

if len(sys.argv) > 3:                                                   # the real SYNTHETIC run's variants
    O = Path(sys.argv[3]); V = O / "muons/muon_0002"; E = Path((V / "far_dir.txt").read_text().strip())
    real = [("S0T1", E / "S0_local"), ("S0T0", V / "T0/S0_local"), ("K2A", V / "K2A_prod"), ("K2Aimp", V / "K2A_imposed")]
    rc, j = grid("real", 7000, *real)
    assert rc == 0 and j["bins_compared"] == 700 and j["variants"]["S0T1"]["replay/waveforms.csv"]["bins"] == 10000 \
        and j["variants"]["K2A"]["replay/waveforms.csv"]["bins"] == 700, j
    bad = U / "K2A_shifted"; shutil.copytree(V / "K2A_prod", bad)
    for f in ("replay/waveforms.csv", "prompt/prompt_waveforms.csv"):
        lines = (bad / f).read_text().splitlines()
        out = [lines[0]] + [",".join([repr(float(l.split(",", 1)[0]) + 1.0), l.split(",", 1)[1]]) for l in lines[1:]]
        (bad / f).write_text("\n".join(out) + "\n")
    rc, j = grid("real_shifted", 7000, real[0], ("K2A_shifted", bad))
    assert rc == 1 and any("first bin centre" in f for f in j["failures"]), j
    st = json.load(open(V / "outcome.json"))["steps"]["validation"]["time_grid"]
    assert st["status"] == "PASS", st
print("PASS time grid: longer S0 clock and the K2A window share one grid (S0 rows truncated to the window bit-exactly); "
      "different dt, first-bin offset, a too-short variant, one moved bin, non-uniform spacing, replay != prompt axis, a "
      "non-integer window and a missing file each FAIL explicitly and nothing is written for the shaped comparison"
      + ("; the real SYNTHETIC run (S0 10000 bins, K2A 700) passes and a 1 ns-shifted K2A copy fails" if len(sys.argv) > 3 else ""))

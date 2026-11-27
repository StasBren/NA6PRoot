#!/usr/bin/env python3
"""Unit regression of k2a_shift_assignment.py and k2a_label_union.py (part F of test_validation_synthetic.sh).

  assignment_units.py VALIDATION_DIR NEW_WORK_DIR
"""
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

H, U = Path(sys.argv[1]), Path(sys.argv[2]); U.mkdir()
SA, LU = str(H / "k2a_shift_assignment.py"), str(H / "k2a_label_union.py")


def k2a_dir(name, shifts, imposed=0):
    d = U / name; d.mkdir()
    with open(d / "shifts.csv", "w") as f:
        f.write("primary,m_U,m_V,computed_m_U,computed_m_V,imposed,Dbar_U_mm,Dbar_V_mm,margin_U,margin_V\n")
        for p, (a, b) in shifts.items():
            f.write(f"{p},{a},{b},{a},{b},{imposed},{1.7 * a},{1.7 * b},{0.3 if p != 2 else 1e-8},0.3\n")
    with open(d / "clouds.csv", "w") as f:
        f.write("primary,wire,N_fC\n" + "".join(f"{p},0,{1.0 + p}\n" for p in shifts))
    return d


g = U / "groups.json"; g.write_text(json.dumps({"primary_shift": {"0": [0, 0], "1": [7, 7], "2": [0, 0]}}))


def run(d):
    p = subprocess.run([sys.executable, SA, str(g), str(d), str(d / "sa.json")], capture_output=True, text=True)
    return p.returncode, json.load(open(d / "sa.json"))


rc, j = run(k2a_dir("agree", {0: (0, 0), 1: (7, 7), 2: (0, 0)}))
assert rc == 0 and j["status"] == "AGREE" and j["near_tie_primaries"] == [2], j
rc, j = run(k2a_dir("disagree", {0: (0, 0), 1: (6, 7), 2: (0, 0)}))
dd = j["disagreements"]
assert rc == 2 and j["n_disagree"] == 1 and dd[0]["primary"] == 1 and dd[0]["s0_m"] == [7, 7] and dd[0]["k2a_m"] == [6, 7] \
    and abs(j["disagree_charge_frac"] - 2.0 / 6.0) < 1e-12, j
rc, j = run(k2a_dir("extra_primary", {0: (0, 0), 1: (7, 7), 2: (0, 0), 3: (1, 1)}))
assert rc == 2 and j["only_k2a"] == [3], j
rc, j = run(k2a_dir("imposed_run", {0: (0, 0), 1: (7, 7), 2: (0, 0)}, imposed=1))
assert rc == 1 and j["status"] == "FAIL", j
d = U / "missing"; d.mkdir(); rc, j = run(d)
assert rc == 3 and j["status"] == "INCONCLUSIVE"


def wave(d, labels, base):
    for f in ("replay/waveforms.csv", "prompt/prompt_waveforms.csv"):
        (d / f).parent.mkdir(parents=True, exist_ok=True)
        with open(d / f, "w") as fo:
            fo.write(",".join(["time_center_ns"] + labels) + "\n")
            for t in range(3):
                fo.write(",".join([f"{5 + 10 * t}"] + [repr(base + 0.1234567890123456789 * (i + 1) * (t + 1)) for i in range(len(labels))]) + "\n")


a, b = U / "A", U / "B"
wave(a, ["U_-1", "U_0", "V_0"], 1.0); wave(b, ["U_-1", "U_0", "V_0"], 2.0)
o = U / "lu_same"
assert subprocess.run([sys.executable, LU, str(o), "--window-ns", "30", f"A={a}", f"B={b}"], capture_output=True).returncode == 0
assert os.path.islink(o / "A") and os.path.islink(o / "B") and json.load(open(o / "label_union.json"))["identical"]
c = U / "C"; wave(c, ["U_0", "U_1", "V_-1", "V_0"], 3.0)
o = U / "lu_diff"
assert subprocess.run([sys.executable, LU, str(o), "--window-ns", "30", f"A={a}", f"C={c}"], capture_output=True).returncode == 0
j = json.load(open(o / "label_union.json"))
assert j["union"] == ["U_-1", "U_0", "U_1", "V_-1", "V_0"] and j["variants"]["A"]["padded_labels"] == ["U_1", "V_-1"] \
    and j["variants"]["C"]["padded_labels"] == ["U_-1"], j
ra = list(csv.reader(open(a / "replay/waveforms.csv"))); rp = list(csv.reader(open(o / "A/replay/waveforms.csv")))
assert rp[0][1:] == j["union"]
for k in range(1, 4):
    assert rp[k][1] == ra[k][1] and rp[k][2] == ra[k][2] and rp[k][3] == "0" and rp[k][4] == "0" and rp[k][5] == ra[k][3], (rp[k], ra[k])
bad = U / "BAD"; wave(bad, ["U_0"], 0.0)
(bad / "prompt/prompt_waveforms.csv").write_text("time_center_ns,U_1\n5,0\n")
assert subprocess.run([sys.executable, LU, str(U / "lu_bad"), "--window-ns", "30", f"A={a}", f"BAD={bad}"], capture_output=True).returncode == 1
print("PASS shift assignment / label union units: AGREE (near tie listed); DISAGREE with primary, both (m_U, m_V) and charge "
      "fraction; K2A-only primary -> DISAGREE; imposed run refused (FAIL); missing -> INCONCLUSIVE; identical labels symlinked; "
      "different labels zero-padded on the U-then-V union with values copied bit-exactly; inconsistent variant refused")

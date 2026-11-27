#!/usr/bin/env python3
"""SYNTHETIC test of k2a_domain_scan.py (part I of test_validation_synthetic.sh): toy kernel and toy moments of
k2a-inline/tests/make_toy_kernel.py in a fake preflight; NOT physics.

  domain_scan_test.py VALIDATION_DIR NEW_WORK_DIR K2A_STRIP_EVAL [VALIDATION_RUN_FOR_MOMENT_CHECK]
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

H, W, EV = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), sys.argv[3]
RUN = Path(sys.argv[4]) if len(sys.argv) > 4 else None
CIR = H.parents[1]
W.mkdir()
subprocess.run([sys.executable, str(CIR / "k2a-inline/tests/make_toy_kernel.py"), str(W / "toy")], check=True, capture_output=True)
sys.path.insert(0, str(CIR / "spatial-ablation/engine")); import cathode_engine as ce  # noqa: E402
KD = W / "toy/kern"; KID = ce.KernelTable(str(KD)).kernel_id
pf = W / "pf"; md = pf / "moments_nonsealed" / "chunk_000"; md.mkdir(parents=True); (pf / "inventory").mkdir()
plan = {"0010": ("multi", "pool"), "0011": ("near", "development"), "0012": ("tie", "pool"), "0013": ("late", "pool"),
        "0014": ("bad_r_high", "pool"), "0015": ("bad_wire", "previously_inspected"), "0016": ("near", "dev_reference")}
for s, (m, _) in plan.items():
    shutil.copy(W / "toy" / f"moments_{m}.csv", md / f"muon_{s}.moments.csv")
rows = [{"slot": s, "role": r, "status": "USABLE"} for s, (_, r) in plan.items()]
rows += [{"slot": "0017", "role": "pool", "status": "EXCLUDED"}, {"slot": "0018", "role": "sealed_validation", "status": "USABLE"},
         {"slot": "0019", "role": "pool", "status": "USABLE"}]                                   # 0019: no moments
with open(pf / "inventory/slots.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["slot", "role", "status"]); w.writeheader(); w.writerows(rows)
(pf / "recenter_20261127_000000").mkdir()
with open(pf / "recenter_20261127_000000/per_event_recenter.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["event", "status_after", "causes_after"]); w.writeheader()
    for s, (m, _) in plan.items():
        w.writerow({"event": f"muon_{s}", "status_after": "OUT_OF_DOMAIN" if m.startswith("bad") else "IN_DOMAIN", "causes_after": ""})
SC = str(H / "k2a_domain_scan.py")
args = ["--pf", str(pf), "--kernel", str(KD), "--kernel-id", KID, "--eval", EV, "--roles", "pool,development,previously_inspected"]
p = subprocess.run([sys.executable, SC, str(W / "scan")] + args + (["--check-moments", str(RUN)] if RUN else []), capture_output=True, text=True)
assert p.returncode == 0, p.stdout[-500:] + p.stderr[-500:]
j = json.load(open(W / "scan/scan.json"))
sc = {r["slot"]: r for r in csv.DictReader(open(W / "scan/scan.csv"))}
assert set(sc) == {"0010", "0011", "0012", "0013", "0014", "0015"}, sorted(sc)          # dev_reference role, EXCLUDED, sealed, no moments skipped
assert sc["0010"]["execution"] == "EVALUATED" and sc["0013"]["execution"] == "EVALUATED" and int(sc["0013"]["t_after_window_clouds"]) >= 1
assert sc["0012"]["execution"] == "EVALUATED" and int(sc["0012"]["near_tie_clouds"]) >= 1
assert sc["0014"]["execution"] == "BLOCKED" and "R_EFF_ABOVE_GRID" in sc["0014"]["blocked_causes"] and 0 < float(sc["0014"]["blocked_frac"]) <= 1
assert sc["0015"]["execution"] == "BLOCKED" and "WIRE_NOT_IN_KERNEL" in sc["0015"]["blocked_causes"]
assert sc["0010"]["kernel_verified"] == "True" and any(r["kernel_verified"] == "False" for r in sc.values())
bc = list(csv.DictReader(open(W / "scan/blocked_clouds.csv")))
assert any(r["slot"] == "0014" and float(r["r_eff_over_grid_max"]) > 1 for r in bc), bc
assert j["agreement_with_recovery"].get("BLOCKED / recovery OUT_OF_DOMAIN") == 2, j["agreement_with_recovery"]
assert (W / "scan/events/muon_0010/replay/waveforms.csv").is_file()                      # kept by default
if RUN:
    assert j["moments_check"] == {}, j["moments_check"]                                   # toy slots differ from the run's slot
# default roles = every non-sealed role (adds dev_reference 0016); sealed never
p = subprocess.run([sys.executable, SC, str(W / "scan_all"), "--pf", str(pf), "--kernel", str(KD), "--kernel-id", KID, "--eval", EV], capture_output=True, text=True)
assert p.returncode == 0 and {r["slot"] for r in csv.DictReader(open(W / "scan_all/scan.csv"))} == set(sc) | {"0016"}
# sealed moments under moments_nonsealed -> refused before anything runs
shutil.copy(W / "toy/moments_near.csv", md / "muon_0018.moments.csv")
p = subprocess.run([sys.executable, SC, str(W / "scan_sealed")] + args, capture_output=True, text=True)
assert p.returncode != 0 and "sealed" in (p.stdout + p.stderr) and not (W / "scan_sealed").exists()
(md / "muon_0018.moments.csv").unlink()
p = subprocess.run([sys.executable, SC, str(W / "scan_req_sealed")] + args + ["--slots", "18"], capture_output=True, text=True)
assert p.returncode != 0 and "sealed" in (p.stdout + p.stderr)
# opt-in pruning: only the scan's own waveform files, sha256 recorded
p = subprocess.run([sys.executable, SC, str(W / "scan_prune")] + args + ["--slots", "10,11", "--prune-waveforms"], capture_output=True, text=True)
assert p.returncode == 0, p.stderr[-300:]
pr = (W / "scan_prune/pruned_waveforms.sha256").read_text().splitlines()
assert len(pr) == 4 and not (W / "scan_prune/events/muon_0010/replay/waveforms.csv").exists() and (W / "scan_prune/events/muon_0010/result.json").is_file()
print("PASS domain scan (SYNTHETIC toy kernel): production grouping on preflight moments; in-domain, late-cloud and near-tie events "
      "EVALUATED with flags; r_eff-above-grid and unknown-wire events BLOCKED with cause, charge fraction and r_eff/r_max; kernel id "
      "verified on first/last; sealed / EXCLUDED / no-moments / non-permitted roles skipped; sealed moments or a sealed request refused; "
      "agreement with the recovery tabulated; opt-in pruning records sha256")

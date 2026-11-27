#!/usr/bin/env python3
"""SYNTHETIC test helper: the complete-coverage run_m2_adapter.sh verdict must record every required item."""
import json
import sys

v = json.load(open(sys.argv[1]))
c = v["coverage"]
assert v["verdict"] == "PASS" and v["label"] == "SYNTHETIC", (v["verdict"], v["label"])
assert c["e2b_muon_0560"] and c["e2b_muon_0997"] and c["t1_0000"], c
assert c["regression_validated"] == ["regression_0042"], c
assert not v["missing"] and not v["not_pass"] and not v["failed_steps"], v
bench = [s for s in v["steps"] if s["step"].startswith("bench_")]
assert len(bench) == 3 * 4 and all(s["status"] == "PASS" for s in bench), bench
print("verdict record complete")

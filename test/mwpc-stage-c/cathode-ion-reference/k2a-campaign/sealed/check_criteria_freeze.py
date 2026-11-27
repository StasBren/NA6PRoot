#!/usr/bin/env python3
"""Verify the frozen sealed-set acceptance criteria (exit 0 = frozen files intact, 1 = mismatch).

  check_criteria_freeze.py [--version 1]

Reads FREEZE_v<N>.sha256 next to this script and requires every listed file to have exactly that sha256, and the JSON
to parse with version == N. Any sealed-set runner must call this first and refuse to start on a non-zero exit; it also
records the criteria sha256 in its run-configuration freeze record.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--version", type=int, default=1)
    a = ap.parse_args()
    rec = HERE / f"FREEZE_v{a.version}.sha256"
    if not rec.is_file():
        print(f"FAIL: no freeze record {rec.name}"); return 1
    bad = []
    for line in rec.read_text().splitlines():
        if not line.strip():
            continue
        sha, name = line.split(maxsplit=1)
        p = HERE / name.strip()
        got = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else "missing"
        if got != sha:
            bad.append(f"{name}: recorded {sha[:12]}..., found {got[:12]}...")
    j = json.loads((HERE / f"sealed_acceptance_criteria_v{a.version}.json").read_text())
    if j.get("version") != a.version:
        bad.append(f"JSON version {j.get('version')} != {a.version}")
    for b in bad:
        print("  MISMATCH", b)
    crit = hashlib.sha256((HERE / f"sealed_acceptance_criteria_v{a.version}.json").read_bytes()).hexdigest()
    print(("PASS" if not bad else "FAIL") + f" sealed acceptance criteria v{a.version} freeze; criteria sha256 {crit}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())

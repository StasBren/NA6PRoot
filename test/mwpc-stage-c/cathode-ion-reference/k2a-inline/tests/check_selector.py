#!/usr/bin/env python3
"""SYNTHETIC test helper: selector prefix N and births for an event with zero-ion (absent) primaries.

  check_selector.py CANDIDATES.csv SOURCE_MOMENTS.csv
muon_0004 holds the source primaries 0 and 1 renumbered to ORIGINAL indices 2 and 5 (0, 1, 3, 4, 6 zero-ion).
"""
import csv
import sys

EFC = 1.602176634e-4
c = list(csv.DictReader(open(sys.argv[1])))
assert c[0]["census"] == "verified" and c[0]["slot"] == "0004", c[0]       # only census-verified event ranks first
d4 = [x for x in c if x["slot"] == "0004"]
got = sorted((int(x["N"]), int(x["last_primary_index"]), int(x["zero_ion_primaries_in_prefix"])) for x in d4)
assert got == [(3, 2, 2), (6, 5, 4)], got
m = list(csv.DictReader(open(sys.argv[2])))
b2 = sum(float(x["N_fC"]) for x in m if x["primary"] == "0") / EFC
b5 = b2 + sum(float(x["N_fC"]) for x in m if x["primary"] == "1") / EFC
births = {int(x["N"]): int(x["expected_births"]) for x in d4}
assert births == {3: round(b2), 6: round(b5)}, (births, b2, b5)
assert all(x["census"] == "unverified (no event.json)" for x in c if x["slot"] in ("0001", "0002"))
assert not any(x["slot"] in ("0003", "0005") for x in c)
print(f"prefixes for primaries {{2}} and {{2, 5}}: N = 3 and 6, births {births}")

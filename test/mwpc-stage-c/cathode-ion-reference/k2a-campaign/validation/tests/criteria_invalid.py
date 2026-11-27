#!/usr/bin/env python3
"""criteria_dry_run with missing / invalid leading-strip and hit-flip metrics (part E2 of test_validation_synthetic.sh):
nothing missing or invalid may be read as an ambiguous lead, a marginal flip or a passed criterion.

  criteria_invalid.py VALIDATION_DIR CRITERIA_V1_JSON
"""
import json
import sys

sys.path.insert(0, sys.argv[1])
import k2a_validation_muon as vm  # noqa: E402

crit = json.load(open(sys.argv[2]))


def lines(mod):
    out = {}
    for tp in ("100", "200"):
        ref = {"variant": "S0T1", "U_peaks": [10.0, 9.95, 1.0], "V_peaks": [8.0, 4.0, 0.805],   # U top-2 within 1 %; V strip 2 at 1.006 x thr f0.1
               "thresholds_fCeq": {s: {f"f{f:g}": f * p for f in (0.02, 0.05, 0.1, 0.2, 0.5)} for s, p in (("U", 10.0), ("V", 8.0))}}
        var = {"variant": "K2A", "U_peaks": [10.0, 9.95, 1.0], "V_peaks": [8.0, 4.0, 0.805], "w_mm": 0.0, "u_mm": 0.0,
               "UV_charge_ratio_rel": 0.0}
        for s in ("U", "V"):
            var.update({f"{s}_lead_index_changed": 0, f"{s}_hitflips_f0.05": 0, f"{s}_hitflips_f0.1": 0, f"{s}_t_f0.1_ns": 0.0,
                        f"{s}_charge_rel": 0.0, f"{s}_second_over_lead": 0.0, f"{s}_L2_rel": 0.0})
        mod(ref, var); out[tp] = {"S0T1": ref, "K2A": var}
    return out


def run(mod):
    return vm.criteria_dry_run(crit, lines(mod), "S0T1", "K2A")


r = run(lambda r, v: None); assert r["status"] == "PASS" and r["n_checks"] == 34, r


def fails_on(key, r):
    hit = [c for c in r["checks"] if c["key"] == key]
    return r["status"] == "FAIL" and hit and all(c["status"] == "FAIL" for c in hit) and not any(key in e for e in r["exemptions"])


cases = {
    # leading strip: the U plane is AMBIGUOUS-eligible (S0 top-2 within 1 %), so a lax check would exempt these
    "lead_missing": (lambda r, v: v.pop("U_lead_index_changed"), "U_lead_index_changed"),
    "lead_none": (lambda r, v: v.update(U_lead_index_changed=None), "U_lead_index_changed"),
    "lead_nan": (lambda r, v: v.update(U_lead_index_changed=float("nan")), "U_lead_index_changed"),
    "lead_inf": (lambda r, v: v.update(U_lead_index_changed=float("inf")), "U_lead_index_changed"),
    "lead_string": (lambda r, v: v.update(U_lead_index_changed="0"), "U_lead_index_changed"),
    "lead_bool": (lambda r, v: v.update(U_lead_index_changed=True), "U_lead_index_changed"),
    "lead_negative": (lambda r, v: v.update(U_lead_index_changed=-1), "U_lead_index_changed"),
    "lead_fraction": (lambda r, v: v.update(U_lead_index_changed=0.5), "U_lead_index_changed"),
    "lead_change_no_ref_peaks": (lambda r, v: (v.update(U_lead_index_changed=1), r.pop("U_peaks")), "U_lead_index_changed"),
    "lead_change_nan_peak": (lambda r, v: (v.update(U_lead_index_changed=1), r["U_peaks"].__setitem__(1, float("nan"))), "U_lead_index_changed"),
    "lead_change_zero_peaks": (lambda r, v: (v.update(U_lead_index_changed=1), r.update(U_peaks=[0.0, 0.0, 0.0])), "U_lead_index_changed"),
    # hit flips: the V strip 2 sits within 2 % of the f0.1 threshold (0.8): a lax check would call these MARGINAL
    "flip_missing": (lambda r, v: v.pop("V_hitflips_f0.1"), "V_hitflips_f0.1"),
    "flip_none": (lambda r, v: v.update(**{"V_hitflips_f0.1": None}), "V_hitflips_f0.1"),
    "flip_nan": (lambda r, v: v.update(**{"V_hitflips_f0.1": float("nan")}), "V_hitflips_f0.1"),
    "flip_negative": (lambda r, v: v.update(**{"V_hitflips_f0.1": -1}), "V_hitflips_f0.1"),
    "flip_no_peaks": (lambda r, v: (v.update(**{"V_hitflips_f0.1": 1}), v.pop("V_peaks")), "V_hitflips_f0.1"),
    "flip_peak_length": (lambda r, v: (v.update(**{"V_hitflips_f0.1": 1}), v.update(V_peaks=[8.0, 4.0])), "V_hitflips_f0.1"),
    "flip_nan_peak": (lambda r, v: (v.update(**{"V_hitflips_f0.1": 1}), v.update(V_peaks=[8.0, 4.0, float("nan")])), "V_hitflips_f0.1"),
    "flip_no_threshold": (lambda r, v: (v.update(**{"V_hitflips_f0.1": 1}), v.update(V_peaks=[8.0, 4.0, 0.795]),
                                        r["thresholds_fCeq"]["V"].pop("f0.1")), "V_hitflips_f0.1"),
    "flip_count_inconsistent": (lambda r, v: v.update(**{"V_hitflips_f0.1": 1}), "V_hitflips_f0.1"),   # peaks show no flip
}
for name, (mod, key) in cases.items():
    r = run(mod)
    assert fails_on(key, r), (name, r["status"], [c for c in r["checks"] if c["key"] == key], r["exemptions"])
# the valid exemptions still work: ambiguous lead, marginal flip with a consistent count
r = run(lambda r, v: v.update(U_lead_index_changed=1)); assert r["status"] == "PASS" and any("AMBIGUOUS" in e for e in r["exemptions"]), r
r = run(lambda r, v: (v.update(**{"V_hitflips_f0.1": 1}), v.update(V_peaks=[8.0, 4.0, 0.795])))
assert r["status"] == "PASS" and any("MARGINAL" in e for e in r["exemptions"]), r
print(f"PASS criteria dry run, invalid metrics: {len(cases)} missing / None / NaN / inf / non-integer / bool / negative lead "
      "and hit-flip metrics, missing or invalid peak lists or threshold, and an inconsistent flip count all FAIL (never "
      "AMBIGUOUS, MARGINAL or PASS); valid ambiguous-lead and marginal-flip exemptions unchanged")

#!/usr/bin/env python3
"""Outcome / report negatives on COPIES of a SYNTHETIC validation run (part D of test_validation_synthetic.sh).

  validation_negatives.py VAL_OUTDIR VALIDATION_DIR NEW_WORK_DIR CRITERIA_V1_JSON
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

O, H, N, CRIT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
MU, REP = str(H / "k2a_validation_muon.py"), str(H / "k2a_validation_report.py")
V0 = O / "muons/muon_0002"; E0 = Path((V0 / "far_dir.txt").read_text().strip()); N.mkdir()


def setstep(f, name, status):
    out = []
    for l in f.read_text().splitlines():
        p = l.split("\t")
        if p[0] == name:
            p[1] = status
        out.append("\t".join(p))
    f.write_text("\n".join(out) + "\n")


def block(V, d):
    j = json.load(open(V / d / "result.json")); j.update(usable=False, status_names="R_EFF_ABOVE_GRID", blocked_charge_fC=1.0)
    (V / d / "result.json").write_text(json.dumps(j))


def case(name, edit, expect):
    d = N / name; V = d / "V"; E = d / "E"
    shutil.copytree(V0, V, symlinks=True); shutil.copytree(E0, E, symlinks=True)
    (V / "far_dir.txt").write_text(str(E) + "\n"); (V / "outcome.json").unlink()
    edit(V, E)
    p = subprocess.run([sys.executable, MU, str(V), "--muon", "0002", "--label", "SYNTHETIC", "--independent", "yes",
                        "--t0", "1", "--criteria", CRIT], capture_output=True, text=True)
    o = json.load(open(V / "outcome.json"))
    assert o["outcome"] == expect, (name, o["outcome"], o["failures"], o["inconclusive"], p.stderr[-300:])
    return o


def prod_blocked(V, E):
    setstep(V / "steps.tsv", "K2A_prod", "BLOCKED"); block(V, "K2A_prod")


def all_blocked(V, E):
    prod_blocked(V, E); setstep(V / "steps.tsv", "K2A_imposed", "BLOCKED"); block(V, "K2A_imposed")
    setstep(V / "steps.tsv", "K2A_cc_export", "BLOCKED"); setstep(E / "steps.tsv", "E2_K2A_local", "FAIL")


def imposed_blocked_only(V, E):
    setstep(V / "steps.tsv", "K2A_imposed", "BLOCKED"); block(V, "K2A_imposed")


cases = [
    ("identity_fail", lambda V, E: setstep(V / "steps.tsv", "identity_cc_vs_python", "FAIL"), "FAIL"),
    ("domain_changed_by_moment_path", imposed_blocked_only, "FAIL"),
    ("consistent_refusal", all_blocked, "OUT_OF_DOMAIN"),
    ("rule_blocked_imposed_ok", prod_blocked, "OUT_OF_DOMAIN"),
    ("agree_but_responses_differ", lambda V, E: setstep(V / "steps.tsv", "grouping_identity", "FAIL"), "FAIL"),
    ("assignment_not_compared", lambda V, E: setstep(V / "steps.tsv", "shift_assignment", "INCONCLUSIVE"), "INCONCLUSIVE"),
    ("assignment_malformed", lambda V, E: setstep(V / "steps.tsv", "shift_assignment", "FAIL"), "FAIL"),
    ("s0_replay_refused", lambda V, E: setstep(E / "steps.tsv", "E2_replay_g_0_0", "INCONCLUSIVE"), "S0_FAILED"),
    ("archive_changed", lambda V, E: (V / "archive_sha256_after.txt").write_text("0" * 64 + "  x\n"), "FAIL"),
    ("moments_differ", lambda V, E: setstep(V / "steps.tsv", "moments_archive_vs_export", "FAIL"), "FAIL"),
    ("t0_failed", lambda V, E: setstep(V / "steps.tsv", "T0_replay_g_0_0", "FAIL"), "VALID_T1_ONLY"),
    ("shaped_missing", lambda V, E: shutil.rmtree(V / "shaped_refS0T1_ion"), "INCONCLUSIVE"),
    ("grouping_effect_missing", lambda V, E: shutil.rmtree(V / "shaped_refK2Aimp_ionelec"), "INCONCLUSIVE"),
    ("shaped_failed", lambda V, E: setstep(V / "steps.tsv", "shaped_refS0T1_ionelec", "FAIL"), "FAIL"),
    ("label_union_failed", lambda V, E: setstep(V / "steps.tsv", "label_union", "FAIL"), "FAIL"),
    ("cmp_numerical", lambda V, E: setstep(V / "steps.tsv", "cmp_S0T1_vs_K2A", "FAIL"), "FAIL"),
]
for c in cases:
    case(*c)
o = case("ood_named", prod_blocked, "OUT_OF_DOMAIN")
assert "R_EFF_ABOVE_GRID" in json.dumps(o["k2a_domain"]) and o["accuracy"]["status"] == "NOT_EVALUATED"
assert any("grouping changes the K2A domain" in n for n in o["notes"]), o["notes"]

# report: an intended-OOD slot is listed separately and never enters the accuracy set, even if it executes VALID
R = N / "rep_ood"; (R / "muons").mkdir(parents=True)
shutil.copy(O / "provenance.txt", R); (R / "slots.txt").write_text("0002\n0005\n")
shutil.copytree(V0, R / "muons/muon_0002", symlinks=True)
shutil.copytree(N / "ood_named/V", R / "muons/muon_0005", symlinks=True)
sel = N / "sel_ood.csv"
sel.write_text("slot,reason,station,channel\n0002,out_of_domain,3,phi\n0005,typical,1,omega\n")
assert subprocess.run([sys.executable, REP, str(R), "--selection", str(sel)], capture_output=True).returncode == 0
j = json.load(open(R / "validation_report.json"))
assert j["intended_ood"] == ["0002"] and j["accuracy_set"] == [] and j["execution_counts"] == {"VALID": 1, "OUT_OF_DOMAIN": 1}, j
rows = {r["slot"]: r for r in csv.DictReader(open(R / "per_muon.csv"))}
assert rows["0002"]["accuracy_v1_primary"].startswith("not counted") and rows["0005"]["execution"] == "OUT_OF_DOMAIN"
assert "Intended out-of-domain events" in (R / "validation_report.md").read_text()

# report: a slot without outcome is INCONCLUSIVE (exit 3); a FAIL outcome gives exit 1; combined report refuses duplicates
R = N / "rep"; R.mkdir(); (R / "muons").mkdir()
shutil.copy(O / "provenance.txt", R); (R / "slots.txt").write_text("0002\n0004\n")
shutil.copytree(V0, R / "muons/muon_0002", symlinks=True)
assert subprocess.run([sys.executable, REP, str(R)], capture_output=True).returncode == 3
j = json.load(open(R / "validation_report.json")); assert j["execution_counts"] == {"VALID": 1, "INCONCLUSIVE": 1}, j["execution_counts"]
shutil.copy(N / "identity_fail/V/outcome.json", R / "muons/muon_0002/outcome.json"); (R / "slots.txt").write_text("0002\n")
assert subprocess.run([sys.executable, REP, str(R)], capture_output=True).returncode == 1
assert subprocess.run([sys.executable, REP, str(N / "comb_dup"), "--from", str(O), str(O)], capture_output=True).returncode != 0
assert not (N / "comb_dup").exists()
assert subprocess.run([sys.executable, REP, str(N / "comb"), "--from", str(O)], capture_output=True).returncode == 0
assert json.load(open(N / "comb/validation_report.json"))["execution_counts"] == {"VALID": 1}
print(f"PASS outcome negatives: {len(cases) + 1} edited copies give FAIL / OUT_OF_DOMAIN (status named; the PRIMARY K2A decides) / "
      "S0_FAILED / VALID_T1_ONLY / INCONCLUSIVE as specified; report keeps an intended-OOD slot out of the accuracy set; "
      "missing outcome -> exit 3, FAIL -> exit 1; combined report refuses a duplicated slot")

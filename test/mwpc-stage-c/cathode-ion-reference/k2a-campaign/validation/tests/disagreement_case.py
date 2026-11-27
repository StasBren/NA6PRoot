#!/usr/bin/env python3
"""Disagreement end-to-end case (part G of test_validation_synthetic.sh), on COPIES of a SYNTHETIC validation run.

  disagreement_case.py prepare VAL_OUTDIR NEW_DIR K2A_RECENTERED_PY MAKE_TOY_KERNEL_PY ENGINE_DIR
  disagreement_case.py check NEW_DIR

prepare: copy the far-mode muon directory and the S0_T0 part of the validation muon; in the COPY, move the S0 grouping of
the first (0,0) primary to (-1,-1) (the production rule keeps (0,0)) and regenerate the Python K2A_local with that
grouping, so the equivalence tests still compare like with like. The S0 waveforms keep the original grouping (both
are exact S0 physics up to window truncation): this exercises the disagreement path, it is not a physics check.
A one-strip disagreement away from a tie puts the imposed variant's strip-relative D beyond the coarse toy run kernel's
D table (+-5.5 P), so this case uses the wider SYNTHETIC toy kernel of k2a-inline/tests/make_toy_kernel.py (+-16 mm)
for every K2A variant; NEW_DIR/kernel_dir.txt and kernel_id.txt name it.
check:   after k2a_validation_k2a.sh and k2a_validation_muon.py ran on the copy.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path


def prepare(O, D, REC, MTK, ENG):
    V0 = O / "muons/muon_0002"; E0 = Path((V0 / "far_dir.txt").read_text().strip())
    V, E = D / "muon_0002", D / "E"
    shutil.copytree(E0, E, symlinks=True)
    V.mkdir(parents=True)
    st = (V0 / "steps.tsv").read_text().splitlines()
    (V / "steps.tsv").write_text("\n".join([st[0]] + [l for l in st[1:] if l.startswith(("far_experiment", "T0_"))]) + "\n")
    tt = (V0 / "timing.tsv").read_text().splitlines()
    (V / "timing.tsv").write_text("\n".join([tt[0]] + [l for l in tt[1:] if l.startswith("T0_")]) + "\n")
    shutil.copytree(V0 / "T0", V / "T0", symlinks=True)
    for f in ("archive_sha256_before.txt", "archive_sha256_after.txt"):
        shutil.copy(V0 / f, V / f)
    (V / "far_dir.txt").write_text(str(E) + "\n")
    g = json.load(open(E / "groups/groups.json"))
    p = next(k for k, v in sorted(g["primary_shift"].items(), key=lambda kv: int(kv[0])) if v == [0, 0])
    g["primary_shift"][p] = [-1, -1]
    (E / "groups/groups.json").write_text(json.dumps(g, indent=1))
    (D / "edited_primary.txt").write_text(p + "\n")
    subprocess.run([sys.executable, MTK, str(D / "toy")], check=True, capture_output=True)
    KD = str(D / "toy" / "kern")
    sys.path.insert(0, ENG); import cathode_engine as ce
    KID = ce.KernelTable(KD).kernel_id
    (D / "kernel_dir.txt").write_text(KD + "\n"); (D / "kernel_id.txt").write_text(KID + "\n")
    shutil.rmtree(E / "K2A_local")
    subprocess.run([sys.executable, REC, KD, str(E / "prepared/cloud_moments.csv"), str(E / "K2A_local"), "--groups",
                    str(E / "groups/groups.json"), "--half-strips", "5", "--kernel-id", KID], check=True, capture_output=True)
    print(f"S0 grouping of primary {p} edited (0,0) -> (-1,-1) in a copy")


def check(D):
    V = D / "muon_0002"; p = int((D / "edited_primary.txt").read_text())
    o = json.load(open(V / "outcome.json")); v = o["steps"]["validation"]
    assert o["outcome"] == "VALID", (o["outcome"], o["failures"], o["inconclusive"])
    sa = o["shift_assignment"]; d0 = sa["disagreements"][0]
    assert sa["status"] == "DISAGREE" and sa["n_disagree"] == 1 and d0["primary"] == p and d0["s0_m"] == [-1, -1] \
        and d0["k2a_m"] == [0, 0] and 0 < sa["disagree_charge_frac"] < 1, sa
    assert "grouping_identity" not in v                                  # only checked when the assignments agree
    for k in ("identity_cc_vs_python", "moment_path_archive_vs_export"):
        assert v[k]["status"] == "PASS", (k, v[k])                       # equivalence tests unchanged
    # grouping effect, measured and reported: in the strip-translation-invariant model a one-strip recentering change is
    # an exact relabelling on the strips both windows cover; what changes is WHICH strips are computed (window edge)
    g = o["comparisons"]["K2Aimp_vs_K2A"]
    raw = json.load(open(V / "cmp_K2Aimp_vs_K2A.json"))["terms"]["ion_plus_electron_step"]
    assert raw["U"]["ref_only_strips"] == ["U_-6"] and raw["V"]["ref_only_strips"] == ["V_-6"], (raw["U"]["ref_only_strips"], raw["V"]["ref_only_strips"])
    assert g["ion_plus_electron_step.U.ref_charge_on_ref_only_strips_frac"] >= 0
    assert g["worst_L2_rel"] < 1e-10, g["worst_L2_rel"]
    lu = o["label_union"]
    assert lu["identical"] is False and "U_-6" in lu["variants"]["K2A"]["padded_labels"] \
        and not lu["variants"]["K2Aimp"]["padded_labels"], lu
    assert set(o["shaped"]["refK2Aimp_ionelec"]["metrics"]) == {"K2A"}
    assert set(o["shaped"]["refS0T1_ionelec"]["metrics"]) == {"K2A", "K2Aimp", "S0T0"}
    assert o["criteria_v1_dry_run"]["primary_vs_S0T1"]["candidate"] == "K2A"
    rp = json.load(open(V / "K2A_prod/result.json")); assert rp["shifts_from"].startswith("production rule")
    print(f"PASS disagreement end-to-end (SYNTHETIC): primary {p} S0 (-1,-1) vs production rule (0,0) reported with "
          f"{100 * sa['disagree_charge_frac']:.3g} % of the charge; grouping effect: common strips L2 {g['worst_L2_rel']:.2g} (exact "
          f"relabelling), extra strips U_-6/V_-6 carry {raw['U']['ref_charge_on_ref_only_strips_frac']:.2g} of the charge; "
          "variants placed on one zero-padded strip label set; primary K2A still on its own grouping; equivalence tests PASS")


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])
    elif sys.argv[1] == "check":
        check(Path(sys.argv[2]))
    else:
        raise SystemExit("usage: disagreement_case.py prepare|check ...")

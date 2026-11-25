#!/usr/bin/env bash
# SYNTHETIC test of the recentering experiment chain on a TOY run (toy transport, coarse toy kernel; not NA60+ physics).
#   bash tests/test_recenter_synthetic.sh TOY_RUN TAG NUMW_DIR LUT WORK
# TOY_RUN/muon_0001 needs T1_export, S0/{replay,prompt}, cloud_moments.csv, K2Akern_TAG; TOY_RUN/kernel_TAG, green_TAG
# (geometry-study/tests/test_pipeline_geometry.sh produces such a run in WORK/base). The toy run is not modified: a
# scratch copy of its directory structure (symlinks + a provenance file) is used as the "reference run".
# 1. environment check: a wrong LUT must abort before anything runs
# 2. refs mode: E1 both shifts PASS; 2b. aggregate completeness negatives (13 variants); E2 single group: S0_local == S0 (identity PASS), K2A_local == K2A base (PASS);
#    aggregate PASS; reference files unchanged
# 3. compare negatives: zero common strips -> INCONCLUSIVE, nonfinite -> FAIL, perturbed e1 candidate -> FAIL,
#    missing files -> INCONCLUSIVE; none can exit 0
# 4. FAR toy event (half of the primaries moved 7 strips along w; valid in the w-invariant toy field): frozen
#    base-centred replay and base-centred engine refuse; grouped S0_local conserves the census; K2A_local vs S0_local
#    as accurate as on the unshifted muon
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; RC="$HERE/recenter"; TA="$HERE/../timing-ablation"; SA="$HERE/../spatial-ablation"
SRC="${1:?}"; TAG="${2:?}"; NUMW="${3:?}"; LUT="${4:?}"; W="${5:?}"
[[ ! -e "$W" ]] || { echo "$W exists"; exit 2; }; mkdir -p "$W"
RUN="$W/refrun"; mkdir -p "$RUN/muon_0001"
for f in kernel_$TAG green_$TAG; do ln -s "$SRC/$f" "$RUN/$f"; done
for f in T1_export S0 cloud_moments.csv K2Akern_$TAG; do ln -s "$SRC/muon_0001/$f" "$RUN/muon_0001/$f"; done
{ echo "lut_sha256=$(sha256sum "$LUT" | cut -d' ' -f1)"; echo "replay_sha256=$(sha256sum "$NUMW/p3_replay_full_spectral_muon.py" | cut -d' ' -f1)"; } > "$RUN/provenance_prepare_0001.txt"
python3 - "$SRC/green_$TAG/axes.json" "$LUT" <<'PY'
import json, sys, hashlib
a = json.load(open(sys.argv[1])); h = hashlib.sha256(open(sys.argv[2], "rb").read()).hexdigest()
assert a["lut_sha256"] == h, "toy green table was not built with this LUT"
PY
KID=$(python3 -c "import sys; sys.path.insert(0, '$SA/engine'); import cathode_engine as ce; print(ce.KernelTable('$RUN/kernel_$TAG').kernel_id)")
export MWPC_NUMW_DIR="$NUMW" MWPC_K2A_OUT="$W/out" MWPC_K2A_KERNEL_ID="$KID"
# 1. wrong LUT must abort
cp "$LUT" "$W/other.npz"; printf x >> "$W/other.npz"
if MWPC_SPECTRAL_LUT="$W/other.npz" bash "$RC/k2a_recenter_experiment.sh" refs "$RUN" "$TAG" 0001 > "$W/wrong_env.log" 2>&1; then
  echo "FAIL: a wrong LUT was accepted"; exit 1; fi
grep -q "MWPC_SPECTRAL_LUT sha256" "$W/wrong_env.log" && echo "PASS environment check refuses a LUT that differs from the validated one"
rm -rf "$W/out"
export MWPC_SPECTRAL_LUT="$LUT"
# 2. refs
bash "$RC/k2a_recenter_experiment.sh" refs "$RUN" "$TAG" 0001 > "$W/refs.log" 2>&1 || { tail -20 "$W/refs.log"; echo "FAIL refs"; exit 1; }
X=$(ls -d "$W"/out/recenter_experiment_*)
python3 - "$X" <<'PY'
import json, sys
r = json.load(open(f"{sys.argv[1]}/results.json")); m = r["muons"]["muon_0001"]
assert r["overall"] == "PASS" and m["E1"] == "PASS" and m["reference_integrity"] == "PASS", (r["overall"], m["E1"], m["reference_integrity"])
c = m["comparisons"]
assert c["E2_S0local_vs_S0"]["status"] == "PASS" and c["E2_S0local_vs_S0"]["kind"] == "identity"
assert c["E2_K2Alocal_vs_K2Abase"]["status"] == "PASS"
assert c["E2_K2Alocal_vs_S0local"]["status"] in ("OK", "WARN")      # coarse toy kernel: WARN expected
print(f"PASS refs (toy): E1 {m['E1']}, S0_local==S0 and K2A_local==K2A base (identity PASS), K2A vs S0 "
      f"{c['E2_K2Alocal_vs_S0local']['status']} (coarse toy kernel), references unchanged, overall {r['overall']}")
PY
# 2b. aggregate completeness negatives (copies of the real experiment directory, never the original)
python3 - "$X" "$RC/k2a_recenter_aggregate.py" "$W/agg" <<'PY'
import json, shutil, subprocess, sys
from pathlib import Path
X, AGG, W = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]); W.mkdir()
def variant(name, edit):
    d = W / name; shutil.copytree(X, d, symlinks=True)
    for f in d.rglob("steps.tsv"):                    # step json paths point into X: retarget to the copy
        f.write_text(f.read_text().replace(str(X), str(d)))
    edit(d)
    for f in ("results.json", "results.md"):
        (d / f).unlink(missing_ok=True)
    rc = subprocess.run([sys.executable, AGG, str(d)], capture_output=True, text=True).returncode
    return rc, json.loads((d / "results.json").read_text())["overall"]
def drop(prefixes):
    def e(d):
        t = d / "muon_0001" / "steps.tsv"
        t.write_text("".join(l + "\n" for l in t.read_text().splitlines() if not any(l.startswith(p) for p in prefixes)))
    return e
def setstatus(step, st):
    def e(d):
        t = d / "muon_0001" / "steps.tsv"
        t.write_text("".join((l.split("\t")[0] + "\t" + st + "\t" + "\t".join(l.split("\t")[2:]) if l.startswith(step + "\t") else l) + "\n"
                             for l in t.read_text().splitlines()))
    return e
def envedit(fn):
    def e(d):
        j = json.loads((d / "env_check.json").read_text()); fn(j); (d / "env_check.json").write_text(json.dumps(j))
    return e
cases = [
    ("baseline_copy", lambda d: None, 0, "PASS"),
    ("no_E1", drop(["E1_"]), 3, "INCONCLUSIVE"),
    ("one_E1", drop(["E1_shift_-1"]), 3, "INCONCLUSIVE"),
    ("no_E2_K2Abase", drop(["E2_K2Abase_vs_S0"]), 3, "INCONCLUSIVE"),
    ("no_E2_at_all", drop(["E2_"]), 3, "INCONCLUSIVE"),
    ("empty_steps", lambda d: (d / "muon_0001" / "steps.tsv").write_text(""), 3, "INCONCLUSIVE"),
    ("missing_report", lambda d: (d / "muon_0001" / "cmp_E2_K2Alocal_vs_S0local.json").unlink(), 3, "INCONCLUSIVE"),
    ("status_mismatch", setstatus("E2_K2Alocal_vs_S0local", "OK"), 1, "FAIL"),   # toy report says WARN
    ("missing_integrity", lambda d: (d / "muon_0001" / "reference_sha256_after.txt").unlink(), 3, "INCONCLUSIVE"),
    ("changed_reference", lambda d: (d / "muon_0001" / "reference_sha256_after.txt").write_text("x\n"), 1, "FAIL"),
    ("missing_muon", envedit(lambda j: j["provenance"].update({"0099": "x"})), 3, "INCONCLUSIVE"),
    ("env_not_pass", envedit(lambda j: j.update(status="FAIL")), 1, "FAIL"),
    ("far_with_E1", envedit(lambda j: j.update(mode="far")), 3, "INCONCLUSIVE"),
]
for name, edit, rc_want, st_want in cases:
    rc, st = variant(name, edit)
    assert (rc, st) == (rc_want, st_want), (name, rc, st)
print(f"PASS aggregate completeness: {len(cases)} variants (missing E1 / E2 steps / reports / integrity / muons, status "
      "mismatch, changed reference, wrong mode) give the expected FAIL or INCONCLUSIVE; the unmodified copy PASSes")
PY
# 3. compare negatives
CMP="$RC/k2a_recenter_compare.py"; N="$W/neg"; mkdir -p "$N"
python3 - "$SRC/muon_0001/S0" "$N" <<'PY'
import shutil, sys
from pathlib import Path
import numpy as np
src, n = Path(sys.argv[1]), Path(sys.argv[2])
def make(name, f):
    d = n / name
    for sub, fn in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv")):
        (d / sub).mkdir(parents=True)
        h = open(src / sub / fn).readline().strip().split(",")
        a = np.loadtxt(src / sub / fn, delimiter=",", skiprows=1)
        h, a = f(h, a)
        np.savetxt(d / sub / fn, a, delimiter=",", header=",".join(h), comments="")
make("disjoint", lambda h, a: ([h[0]] + [l.split("_")[0] + "_" + str(int(l.split("_")[1]) + 100) for l in h[1:]], a))
def nanf(h, a):
    a = a.copy(); a[5, 3] = np.nan; return h, a
make("nonfinite", nanf)
def pert(h, a):
    a = a.copy(); a[:, 2] *= 1 + 1e-6; return h, a
make("perturbed", pert)
PY
expect(){ local want="$1" code="$2"; shift 2; local rc=0; python3 "$CMP" "$@" > /dev/null 2>&1 || rc=$?
  local st; st=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['status'])" "${@: -1}" 2>/dev/null || echo none)
  [[ $rc -eq $code ]] || { echo "FAIL: expected exit $code got $rc ($*)"; exit 1; }; }
expect INCONCLUSIVE 3 "$SRC/muon_0001/S0" "$N/disjoint" --kind k2a_vs_s0 --out "$N/c1.json"
expect FAIL 1 "$SRC/muon_0001/S0" "$N/nonfinite" --kind k2a_vs_s0 --out "$N/c2.json"
expect FAIL 1 "$SRC/muon_0001/S0" "$N/perturbed" --kind e1 --out "$N/c3.json"
expect INCONCLUSIVE 3 "$SRC/muon_0001/S0" "$N/missing" --kind e1 --out "$N/c4.json"
for k in 1 2 3 4; do python3 -c "import json,sys; d=json.load(open(sys.argv[1])); assert d['status'] in ('FAIL','INCONCLUSIVE'), d['status']" "$N/c$k.json"; done
echo "PASS compare negatives: zero common strips INCONCLUSIVE, nonfinite FAIL, 1e-6 perturbation FAIL under e1, missing files INCONCLUSIVE"
# far toy event
F="$W/far/muon_0001"; mkdir -p "$F"
python3 - "$RUN/muon_0001/T1_export" "$F/T1_export" <<'PY'
import csv, shutil, sys
from pathlib import Path
import numpy as np
src, dst = Path(sys.argv[1]), Path(sys.argv[2]); dst.mkdir(parents=True)
seal = dict(l.split("=", 1) for l in (src / "complete.txt").read_text().split() if "=" in l)
V = np.fromfile(src / "vertices.f64", dtype="<f8").reshape(int(seal["vertices"]), 5)
rows = list(csv.DictReader(open(src / "ions.csv")))
prim = sorted({int(r["primary"]) for r in rows}); far = set(prim[::2])
for r in rows:
    if int(r["primary"]) in far:
        s, n = int(r["first_vertex"]), int(r["vertices"]); V[s:s + n, 3] += 7 * 1.7; V[s:s + n, 4] += 7 * 1.7
# toy weights are not integer per primary (real archives are: weight = exact/saved per primary); make them so,
# keeping the total census, because the frozen replay checks every (sub-)export's integer census
from collections import defaultdict
tot = defaultdict(float)
for r in rows:
    tot[int(r["primary"])] += float(r["weight"])
target = {p: max(1, round(t)) for p, t in tot.items()}
diff = int(round(sum(tot.values()))) - sum(target.values()); target[prim[0]] += diff
for r in rows:
    r["weight"] = repr(float(r["weight"]) * target[int(r["primary"])] / tot[int(r["primary"])])
with open(dst / "ions.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
V.tofile(dst / "vertices.f64")
ws = sum(float(r["weight"]) for r in rows)
(dst / "complete.txt").write_text((src / "complete.txt").read_text().replace(f"weight_sum={seal['weight_sum']}", f"weight_sum={ws!r}"))
print(f"toy far event: {len(far)} of {len(prim)} primaries moved +7 strips along w")
PY
python3 "$SA/s_kernel.py" moments "$F/T1_export" "$F/cloud_moments.csv" > /dev/null
if (cd "$NUMW" && python3 p3_replay_full_spectral_muon.py "$LUT" "$F/T1_export" "$F/S0_base/replay" --half-strips 5) > "$F/base.log" 2>&1; then
  echo "FAIL: frozen base-centred replay accepted a far event"; exit 1; fi
why=$(grep -o -E "selected channels omit the landing strip[^\"]*|query outside declared strip phase range" "$F/base.log" | head -1)
[[ -n "$why" ]] || { echo "FAIL: frozen replay failed for another reason"; tail -3 "$F/base.log"; exit 1; }
echo "PASS frozen base-centred S0 refuses the far toy event: $why"
python3 - "$RUN/kernel_$TAG" "$F/cloud_moments.csv" "$SA/engine" <<'PY'
import sys
sys.path.insert(0, sys.argv[3]); import cathode_engine as ce
k = ce.KernelTable(sys.argv[1]); e = ce.CathodeResponseEngine(k, "K2A")
try:
    e.evaluate(ce.read_moments_csv(sys.argv[2])); print("FAIL: base-centred engine accepted far clouds"); sys.exit(1)
except ValueError as x:
    print(f"PASS base-centred engine refuses the far toy event: {x}")
PY
E="$W/far/exp"; mkdir -p "$E"
python3 "$RC/k2a_group_export.py" "$F/T1_export" "$E/groups" --mode primary
mkdir -p "$E/S0groups"
for g in "$E"/groups/g_*; do
  (cd "$NUMW" && python3 p3_replay_full_spectral_muon.py "$LUT" "$g" "$E/S0groups/${g##*/}/replay" --half-strips 5) > /dev/null
  python3 "$TA/t_prompt_electron_term.py" "$LUT" "$g" "$E/S0groups/${g##*/}/prompt" --replay-dir "$NUMW" --half-strips 5 > /dev/null
done
python3 "$RC/k2a_local_sum.py" "$E/groups" "$E/S0groups" "$E/S0_local"
python3 "$RC/k2a_recentered_k2a.py" "$RUN/kernel_$TAG" "$F/cloud_moments.csv" "$E/K2A_local" --groups "$E/groups/groups.json"
python3 "$RC/k2a_recenter_compare.py" "$E/S0_local" "$E/K2A_local" --kind k2a_vs_s0 --out "$E/cmp.json" --label "far toy K2A_local vs S0_local" > "$E/cmp.txt" || true   # status asserted below
python3 - "$E/cmp.json" "$X/muon_0001/cmp_E2_K2Alocal_vs_S0local.json" <<'PY'
import json, sys
far, ref = (json.load(open(p)) for p in sys.argv[1:])
nf = far["terms"]["ion_only"]["U"]["common_strips"]
assert nf > 11, nf                                  # union of two separated local windows
assert far["status"] in ("OK", "WARN") and far["worst_L2_rel"] is not None, far["status"]
assert far["worst_L2_rel"] <= 1.5 * ref["worst_L2_rel"] + 1e-3, (far["worst_L2_rel"], ref["worst_L2_rel"])
print(f"PASS far toy: {nf} absolute U strips; K2A_local vs S0_local worst rel L2 {far['worst_L2_rel']:.3g} "
      f"(unshifted toy muon {ref['worst_L2_rel']:.3g}, coarse toy kernel)")
PY
echo "PASS recentering synthetic test; work dir $W"

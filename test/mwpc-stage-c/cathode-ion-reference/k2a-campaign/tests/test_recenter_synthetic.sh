#!/usr/bin/env bash
# SYNTHETIC test of the recentering experiment chain on a TOY run (toy transport, coarse toy kernel; not NA60+ physics).
#   bash tests/test_recenter_synthetic.sh TOY_RUN TAG NUMW_DIR LUT WORK
# TOY_RUN/muon_0001 needs T1_export, S0/{replay,prompt}, cloud_moments.csv, K2Akern_TAG; TOY_RUN/kernel_TAG.
# (geometry-study/tests/test_pipeline_geometry.sh produces such a run in WORK/base.)
# 1. refs mode on the toy muon: E1 symmetry (frozen replay relabelled by +-1 strip == frozen S0 to 1e-10),
#    E2 with one group: S0_local == S0 exactly, K2A_local == base-centred K2A to 1e-12
# 2. FAR toy event: half of the primaries moved 7 strips along w (D + 7P, same indices; the toy field is w-invariant, so
#    this is a valid toy event): the frozen base-centred replay must REFUSE (landing outside the 11 strips), the
#    base-centred engine must refuse (D outside the kernel), the grouped S0_local must conserve the collected-ion
#    census exactly, and K2A_local must agree with S0_local as well as the toy kernel agrees on the unshifted muon.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; RC="$HERE/recenter"; TA="$HERE/../timing-ablation"; SA="$HERE/../spatial-ablation"
RUN="${1:?}"; TAG="${2:?}"; NUMW="${3:?}"; LUT="${4:?}"; W="${5:?}"
[[ ! -e "$W" ]] || { echo "$W exists"; exit 2; }; mkdir -p "$W"
export MWPC_SPECTRAL_LUT="$LUT" MWPC_NUMW_DIR="$NUMW" MWPC_K2A_OUT="$W/out"
bash "$RC/k2a_recenter_experiment.sh" refs "$RUN" "$TAG" 0001 > "$W/refs.log" 2>&1
X=$(ls -d "$W"/out/recenter_experiment_*)
grep -q "PASS E1 shift 1:" "$X/muon_0001/E1.txt" && grep -q "PASS E1 shift -1:" "$X/muon_0001/E1.txt"
python3 - "$X/muon_0001" <<'PY'
import json, sys
d = sys.argv[1]
a = json.load(open(f"{d}/cmp_S0local_vs_S0.json")); b = json.load(open(f"{d}/cmp_K2Alocal_vs_K2Abase.json"))
assert a["worst_L2_rel"] == 0.0, a["worst_L2_rel"]
assert b["worst_L2_rel"] < 1e-12, b["worst_L2_rel"]
print(f"PASS refs (toy): E1 symmetry, S0_local == S0 exactly, K2A_local == K2A base ({b['worst_L2_rel']:.1e})")
PY
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
python3 "$RC/k2a_recenter_compare.py" "$E/S0_local" "$E/K2A_local" --out "$E/cmp.json" --label "far toy K2A_local vs S0_local" > "$E/cmp.txt"
python3 - "$E/cmp.json" "$X/muon_0001/cmp_K2Alocal_vs_S0local.json" <<'PY'
import json, sys
far, ref = (json.load(open(p)) for p in sys.argv[1:])
nf = far["terms"]["ion_only"]["U"]["common_strips"]
assert nf > 11, nf                                  # union of two separated local windows
assert far["worst_L2_rel"] <= 1.5 * ref["worst_L2_rel"] + 1e-3, (far["worst_L2_rel"], ref["worst_L2_rel"])
print(f"PASS far toy: {nf} absolute U strips; K2A_local vs S0_local worst rel L2 {far['worst_L2_rel']:.3g} "
      f"(unshifted toy muon {ref['worst_L2_rel']:.3g}, coarse toy kernel)")
PY
echo "PASS recentering synthetic test; work dir $W"

#!/usr/bin/env bash
# SYNTHETIC tests of the response engine (no ROOT, toy transport): needs the work dir of test_birth_moments.sh
# (births, moments from both paths, a toy kernel and the s_kernel.py K2A outputs).
#   bash test_engine.sh BIRTH_TEST_DIR WORKDIR
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; B="${1:?test_birth_moments.sh work dir}"; W="${2:?work dir}"
mkdir -p "$W"; [[ -z "$(ls -A "$W")" ]] || { echo "workdir not empty" >&2; exit 2; }
g++ -O3 -std=c++17 -Wall -Wextra -shared -fPIC -o "$W/libkernel_eval.so" "$HERE/kernel_eval.cc"
export NA60_KERNEL_EVAL_LIB="$W/libkernel_eval.so"
python3 - "$B" "$HERE" <<'PY'
import csv, sys, numpy as np
B, H = sys.argv[1], sys.argv[2]; sys.path.insert(0, H)
import cathode_engine as ce
# 1. numpy moments_from_births == BirthMoments.hh core (both from the same synthetic births)
rows = list(csv.DictReader(open(f"{B}/syn/births.csv")))
fr = list(csv.DictReader(open(f"{B}/syn/frame.csv")))[0]
g = lambda c: np.array([float(r[c]) for r in rows])
m_np = ce.moments_from_births(g("primary").astype(int), g("u_cm"), g("v_cm"), g("w_cm"), g("t_ns"), g("weight"),
                              float(fr["patch_shift_cm"]), float(fr["u_in_mm"]), float(fr["w_in_mm"]))
m_cc = ce.read_moments_csv(f"{B}/mom_births.csv")
assert (m_np.base_U, m_np.base_V) == (m_cc.base_U, m_cc.base_V), "strip frame differs"
worst = max(float(np.max(np.abs(m_np.cols[c] - m_cc.cols[c]) / np.maximum(1, np.abs(m_cc.cols[c])))) for c in ce.MOMENT_COLS)
assert worst < 1e-12, worst
print(f"PASS numpy moments_from_births == BirthMoments.hh (max scaled diff {worst:.2e})")
# 2. kernel identity: mismatch refused
k = ce.KernelTable(f"{B}/kern", verify=True, preload=True)
try:
    ce.KernelTable(f"{B}/kern", expected_id="0" * 64); raise SystemExit("FAIL: wrong kernel id accepted")
except ValueError:
    print("PASS wrong kernel id refused")
# 3. implementations vs the validated loop and vs s_kernel.py apply output
ref = np.loadtxt(f"{B}/kb/K2Akern/replay/waveforms.csv", delimiter=",", skiprows=1)[:, 1:]
refp = np.loadtxt(f"{B}/kb/K2Akern/prompt/prompt_waveforms.csv", delimiter=",", skiprows=1)[:, 1:]
for model in ("K2A", "K2", "K1A", "K1"):
    e = ce.CathodeResponseEngine(k, model); rl = e.evaluate(m_cc, "loop")
    for impl in ("vector", "compiled"):
        r = e.evaluate(m_cc, impl)
        di = np.linalg.norm(r.ion - rl.ion) / np.linalg.norm(rl.ion); de = np.linalg.norm(r.electron - rl.electron) / np.linalg.norm(rl.electron)
        assert di < 1e-13 and de < 1e-13, (model, impl, di, de)
    if model == "K2A":
        d = np.linalg.norm(rl.ion.reshape(len(rl.ion), -1) - ref) / np.linalg.norm(ref)
        dp = np.linalg.norm(rl.electron.reshape(len(rl.ion), -1) - refp) / np.linalg.norm(refp)
        assert d < 1e-13 and dp < 1e-13, (d, dp)
        print(f"PASS engine loop == s_kernel.py apply (ion {d:.1e}, electron {dp:.1e})")
print("PASS vector and compiled == loop for K1, K2, K1A, K2A (rel L2 < 1e-13)")
# 4. electronics interface == t_shaped_compare stand-in shaper
sys.path.insert(0, H + "/../../timing-ablation"); import t_shaped_compare as tsc
r = ce.CathodeResponseEngine(k, "K2A").evaluate(m_cc)
el = ce.Electronics.crrc_standin(100.0)
h, _ = tsc.impulse(100.0, 3, 1.0, 2000.0)
a1 = el.shape(r.ion, r.dt_ns).reshape(-1, r.ion.shape[1] * r.ion.shape[2]); a2 = tsc.shape(r.ion.reshape(len(r.ion), -1), r.dt_ns, h, 1.0)
assert np.allclose(a1, a2, rtol=0, atol=1e-12 * np.abs(a2).max())
print("PASS Electronics.crrc_standin == t_shaped_compare stand-in")
PY
echo "PASS engine synthetic tests"

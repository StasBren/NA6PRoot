# Phase-C: external 120-parent Geant4→TrackHeed→full Garfield holdout

**Prepared 2026-10-31; NOT physically run on user's WSL yet.**

## Motivation and fixed comparison

Previous original A20/A100, 100k microscopic LUT samples, and 120 true muons showed:
- Conditional-only, 90%-PI gain coverage A20: 368/480=76.7%; guarded spatially pooled A20: 411/480=85.6%.
- Shared global training-cell bootstrap, 5,000 iterations: A20 406/480=84.6%, A100 399/480=83.1%, guarded A20 417/480=86.9%, guarded A100 414/480=86.3%. Mean absolute error of predictive mean 22,391/20,858/18,750/18,843, respectively.
- Improvements are **source conditioned** because prior 120 positions guided A100 cell selection, and pooling radius, weights and guard were investigated before external validation.
- Fixed model **before new muon truth**: radius 2 bins each axis, weight 25% own A20 / own A100, 75% neighboring **old A20** library only; Gaussian neighbor weights `exp(-(Δiu²+Δiv²)/2)`; no crossing v=0, near-wire u sign and 0.35 mm central/field-wire protection; global cell training bootstrap shared across all parents.

## Real new physics versus reuse of preexisting source files

`test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py` (new standalone wrapper) uses **the original SHA-sealed Geant4 J/psi chunk 2**, 50k generated parent events, where previous 120 validation used ONLY 120 separate parent events. It chooses **120 different, previously unused physical parent events**, without looking at their microscopic outcomes, balancing 10 per each of 6 stations × 2 signed-angle strata. It rejects all old 120 parent event IDs (not only track IDs) and requires unique event IDs in the new set. A new deterministic salt fixes source selection. **These are genuinely new physical events for the analysis, but not freshly generated Geant4 events:** they were produced and SHA archived earlier. The reason for this reuse is to get independent original physical muon input without re-running Geant4 needlessly.

The existing, native, previously validated C++ `--muon-cache-physics` is invoked to run **new TrackHeed clusters/electrons per newly selected parent**, then **one genuinely new full microscopic Garfield AvalancheMicroscopic realization** for every accepted Heed primary and each of the four B=0 geometries. Event truth is **never obtained from the LUT**. Source and full reference are SHA checkpointed for each 10 parents (12 total) in a **distinct new output root**. The 84k A20, 16k A100, 120 old full muon references and all V1/V2/V3 RKF artifacts are unchanged.

The wrapper loads the immutable source of `pi_holdout120_A20_from_D500` under a **separate Python module name** using `importlib.util.spec_from_file_location`, and redirects only that isolated module instance's directory globals `ROOT/INPUT/REF/COMP` to the separate root. The canonical `pi_holdout120_A20_from_D500` module remains untouched, so the standard SHA-verified old A20/A100 source loader continues using the original 120-muon reference. It does **not** invoke the original old120 `prepare`, `evaluate` or `publish`. Its `run` delegates to original `--muon-cache-physics` + checked C++ parser+simulator, and original `finish` SHA-seals the new full physics. The common fixed-neighbor / bootstrap math from `pi_spatial_lut_muon120_hier_v1.py` is called unchanged; the external wrapper only removes the *old fixed 5881-Heed-primary count*, replacing it with the verified NEW Heed count, and does a new source/Geant4-truth join.

### Run 1: freeze independent parent list before knowing new full micro truth

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py selftest
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py prepare
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py check
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py status
)
```

**Stop if any step fails.** The `prepare` writes the fixed input CSV SHA and freezes *the old 120 source SHA*, the physical source file SHA, the original full Garfield runner source SHA, and the prior guarded-model code SHAs. The script refuses to overwrite an existing campaign. Selection uses neither the original full 120 muons' physics output nor the future new muon truth. The original old120 ID list is only used for EXCLUSION.

### Run 2: preflight real native input without microscopic avalanches

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py preflight
)
```

This invokes the native B1 C++ input parser and shared TrackHeed validation for the 12 input groups but performs **no full microscopic avalanches**. If it finds a native problem, stop and preserve the frozen source list for forensic review. **Do not re-run `prepare`**.

### Run 3: first real Garfield 10-muon checkpoint

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py run --max-chunks 1 --hours 8
python3 test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py status
)
```

This runs full Garfield **for the first 10 new muon parents** in each of 4 geometries, SHA sealing an `full_garfield_120/chunks/chunk_00/complete.json` and preserving untouched raw physics logs, seed, cluster, sector and birth-sample CSVs. Stop and send output for review if the C++ physics aborts, any caps/unknown nonmetal failure or an incomplete checkpoint is detected. A failed/incomplete chunk is intentionally NOT overwritten and requires manual forensic action.

### Run 4: remaining new physical muons after first successful native checkpoint

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py run --max-chunks 12 --hours 8
python3 test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py status
)
```

If after 8h only part completed, **repeat the exact same Run 4 block**. It verifies and skips already SHA-sealed checkpoints, no overwrite, no new selection, and stops between chunks. On full 12/12, old proven `finish()` writes an immutable full 120-muon real Garfield reference and SHA manifests. No publication/upload is performed.

### Run 5: blind, fixed-hyperparameter external evaluation

Only after 12/12 SHA sealed:

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py evaluate --reps 3000
python3 test/mwpc-stage-c/pi_external_muon120_guarded_lut_v1.py status
)
```

The 3000 global bootstrap iterations use the **original A20/A100 training pools and fixed 5x5 protected-neighbor algorithm**, with the new Heed count. The algorithm never fits weights or revisits selection after seeing new Garfield truth. Evaluate `A20`, `A100`, `guarded_A20`, `guarded_A100`. Each new physical parent supplies four correlated geometries; 95% bootstrap intervals clustered across parents. Reports 90%-nominal PI coverage of full micro ion births, mean absolute error of predicted expectation, mean PI width, ion hemisphere and top/bottom asymmetries. Expected outputs `blind_evaluation_fixed_guard_v1/per_external_muon_4_models.csv`, `summary.json`, `manifest.json`.

## Limitations

- Independent from old 120 **parent event identities**, not independent of the common Geant4 generator/Jpsi energy/configuration, nor a random flux-unbiased experiment (10/station/sign stratification). Existing A100 cell selection was conditioned on old120 source only.
- Only **one** stochastic full microscopic Garfield response per primary per configuration on new Heed. Clustered 120-parent error bars are essential; do not expect each outcome to equal a conditional mean.
- No magnetic field, gas variation, ion drift, or induced U/V strip charge in this phase.
- The original shared-cell bootstrap approximates finite LUT sample uncertainty, not a fully calibrated Bayesian posterior. Absolute coverage may fail nominal 90% even on external muons.
- **Do not run multi-day A100-all-cells** until this new independent holdout is examined; then decide targeted 48k additional micro avalanches or more.

**Validation status:** New Python source `py_compile` PASS, synthetic inherited global 120-parent correlated geometry/ion sector selftest PASS, synthetic variable-Heed-size 120-parent guarded kernel shape PASS. Native WSL Geant4/TrackHeed/full Garfield run **not yet executed or confirmed**.

**Critical isolation verification:** initial version redirected canonical legacy module globals. Before any native new-120 preparation/run, this was corrected in commit `8d75adf` to an **isolated driver module instance**; avoids accidentally loading the new still-unfinished reference during verification of frozen A20/A100 old120 source. Python syntax check and isolated mock-module selftest PASS; native full Geant4/TrackHeed/Garfield WSL remains unrun.

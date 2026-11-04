# NA60+ MWPC p3-only avalanche model decision and full-ion cloud pilot

**Status 2026-11-04:** Two independent 3 mm geometries: `p3_noFW`,
`p3_FW75`. These are separate scientific hypotheses and must never
be combined for error/coverage/CRPS model acceptance. No induced signal,
ion drift, weighting fields or electronics are calculated here.

## Physics decision

1. Keep original immutable A20, A100 and independent old/new120 Geant4→Heed→
   full microscopic Garfield archives. The fixed `guarded_A20`
   (25% local + 75% Gaussian nearby A20 with geometric wire/hemisphere guards)
   is the **current validated count/sector baseline**, NOT a proven spatial
   cloud model. The original 5 ion birth samples per avalanche cannot certify
   eventwise space/time response.
2. Reassess Pólya fairly: `pi_p3_fair_avalanche_model_selection_v1.py` uses
   SAME old A20 training for local20/guarded20/geometry-only zero+Gamma;
   holds out 80 independent new A100 outcomes in each of 50 previously selected
   cells **per geometry**. Reports CRPS and rare tails, with cell-bootstrap CI.
   This is development, not untouched new-muon validation; near-wire A100
   exposure is very sparse.
3. Challenge count models on old, previously-seen **external new120** physically
   independent Geant4 parent muons with shared frozen Heed electrons:
   `pi_p3_muon_count_model_challenge_v1.py`. Gamma fit on old A20 only;
   protected A20 fixed before original external new120 Garfield truth. Both
   sampled with 3000 replicas per parent. Separate parent bootstrap per
   geometry, proper CRPS and q99 pinball. **This is retrospective**, not a
   newly blind hypothesis-selection test, because that reference has already
   been analyzed in previous work.
4. **Critical next physics**: microphysical complete ion BIRTH cloud (all
   Garfield ionization callback positions `u,v,w,t`, not ion drift);
   `pi_p3_full_cloud_convergence_v1.py` uses untouched old D500 initial
   electrons with deterministic sample passes >=100, excluding **all A20
   indices 0..19 and A100 extra indices 20..99**, and intentionally sampling
   distinct spatial cells. The independent Garfield ionization collision
   RNG is generated from original D500 IDs. The two p3 geometries have
   individually balanced anode-near/transition/bulk and plus/minus v samples;
   FW75 additionally gets a *near field-wire* stratum.
5. Pilot compares *exact same Garfield avalanche* full birth cloud to
   repeated 32/128/512-ion WITHOUT-replacement samples. The exact
   per-avalanche all-birth moments, angular sector counts and number of
   births are checked against the C++ full-census output. The output
   compares radius/time moments and percentiles, angle fractions,
   top-bottom asymmetry and 12x6 joint azimuth/time histogram TV.
   This tests storage compression, not geometry-conditioned model
   generalization. No automatic go/no-go threshold is asserted.
6. Based on each geometry's result, select spatial response storage: 128,
   512, or a more explicit cloud representation. The separate C++ target
   `mwpc_phase_b_joint_cloud_lut_v1` now supports
   `--lut-sidecar-joint128`, `--lut-sidecar-joint512`, and
   `--lut-sidecar-full-births`. **Do not launch a uniform multi-day
   new cloud LUT campaign before the exact/full cloud convergence report.**
   That would risk spending compute and disk on an unverified schema.
7. Then make one further SHA-frozen learning-curve campaign for the **two p3
   geometries only**, at increasing independent samples per spatial bin.
   Evaluate trial-level proper scores and multivariate spatial/temporal
   distances against separate untouched D500 positions; finally test
   the new model on reserved genuinely independent Geant4 muons. Spatial
   pooling must be relearned/validated for cloud shapes, not blindly
   inherited from the 5-sector count baseline.

**Not validated:** this work does not provide an induced cathode current,
SPSC chamber readout, detector resolution, or response at new HV/B/gas;
the current analytic-wire gas model is 1800 V, B=0, Ar/CO2 70:30,
cathode spacing ±3 mm. Frozen D500 `w=0` starts rely on translation
invariance along wires, which must be revisited for finite-strip/end effects.

## 0. Update and test Python only (no new Garfield physics)

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" -m py_compile \
  test/mwpc-stage-c/pi_p3_fair_avalanche_model_selection_v1.py \
  test/mwpc-stage-c/pi_p3_muon_count_model_challenge_v1.py \
  test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py
"$PY" test/mwpc-stage-c/pi_p3_fair_avalanche_model_selection_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_muon_count_model_challenge_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py selftest
```

**Optional immediate no-physics model evaluations** (separate new sibling
output directories; each refuses to overwrite):

```bash
"$PY" -u test/mwpc-stage-c/pi_p3_fair_avalanche_model_selection_v1.py audit --bootstrap 4000
"$PY" -u test/mwpc-stage-c/pi_p3_muon_count_model_challenge_v1.py audit --reps 3000 --bootstrap 4000
```

Both compare two p3 geometries in separately reported groups, never as a
pooled four-geometry score. The muon challenge on the already-known new120
sample is a retrospective check, not an untouched final validation.

## 1. Build the SEPARATE full-ion-cloud microscopic research target

The CMake target was added separately from all pinned original physical
binaries. Rebuilding this C++ code has no effect on already-sealed A20,
A100 or 2400-parent muon samples.

```bash
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build
cmake --build test/mwpc-avalanche/build --target mwpc_phase_b_joint_cloud_lut_v1 -j2
```

If CMake needs the Garfield package discovery hint, reuse your original
successful Garfield++ `-DGarfield_DIR=...` setup. **Do not change** the gas,
RNG, field or old source code.

## 2. Pilot with exact full ion birth point clouds, 24 per geometry

**Each of the 48 independent source sites is one frozen old D500 original
ELECTRON START**, never a muon crossing. D500 drift endpoint flags are
provenance only and never used to accept/reject a true new microscopic
avalanche. Default is 24 per geometry (6 native full-cloud 8-start
checkpoints). This is a compression-pilot scale, not a tail
statistics/production model sample.

```bash
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py prepare --per-geometry 24
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py preflight
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py run --max-chunks 1 --hours 2
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py status
```

Inspect first eight true full-cloud avalanches and their `physics.log`
before leaving WSL unattended. If the native process stops midway through
one block, the directory is preserved **unsealed**, and the runner will
fail instead of pretending it is complete or overwriting it. Never delete
or edit a failed physical checkpoint without a recovery plan.

## 3. Finish the remaining five SHA-frozen chunks

```bash
nohup nice -n 15 "$PY" -u \
  test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py \
  run --max-chunks 5 --hours 16 \
  > "$HOME/na6p/p3_full_cloud_pilot_20261104.log" 2>&1 < /dev/null &
echo "P3 full-cloud PID $!"
```

Read-only progress:

```bash
tail -n 25 "$HOME/na6p/p3_full_cloud_pilot_20261104.log"
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py status
```

If the soft 16h limit is reached BETWEEN complete chunks, repeat the same
`run --max-chunks 5 --hours 16` command. Every completed chunk is
cryptographically verified and skipped. Never start two runs at once.

## 4. Compare 32/128/512 points against FULL identical microscopic clouds

Run only after `status` reports `6 / 6` complete chunks:

```bash
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py analyze --reps 64
```

Outputs for p3 noFW and p3 FW75 separately:
`test_runs/mwpc_stage_c/pi_p3_full_cloud_convergence_v1/analysis_complete/`
contains `summary.json`, `per_avalanche.json`, `manifest.json`.
The full particle births remain SHA-sealed in
`.../full_cloud_chunks/chunk_0000/cloud_birth_reservoir.csv` etc.
No Cathode signal is simulated.

The origin and physics configuration, exact 4-D birth points, ion counts,
electron endpoint 64-point reservoir, timing and birth moments are
preserved in each sealed chunk. Full-census comparison validates every
surviving cloud in the pilot, *not* only pooled averages.

**Decision process:** for *each* geometry independently, examine the
full-clock results and all near-wire/bulk strata. If 128 ion samples
miss significant angle-time correlations, radial width, upper time tails
or near-wire distribution, consider 512 or a different compression.
The precise tolerances are engineering precriteria to be agreed before
any induced-signal study; do not call a model 'cathode-ready' based on
these pilot distributions alone. Retain full clouds as ground truth.

### What will and will not be determined

After commands (0)–(4), we can decide whether:
- `guarded_A20` or count-only `geometry_hurdle_gamma` better predicts
  **counts** on selected cells and on the known new120 muons;
- `128` is an adequate *exploratory* compression for **ion birth space/time**
  versus more points, separately near anode and field wire.

We **cannot** yet certify the complete joint spatio-temporal surrogate as
a replacement for Garfield: that requires an independently trained
conditional joint bank, increasing learning-curve sample sizes and
untouched test source electrons/muons. Preserve all physical outputs to
make that later campaign auditably reusable.

## Reproducibility/security

New driver/source target only. The old A20/A100 original and native muon
source/binary remain SHA-frozen. New records refuse overwrite, every
completed original exact-source full-cloud chunk stores source and binary
SHA. New script does not upload physical data externally and has no
network dependency. Native C++ compilation/execution must occur on
user's installed WSL Garfield++; Python mock/synthetic selftests are NOT
proof of native physics operation.

## 2026-11-04: Corrected new120 in-metal electron audit (after first WSL run)

The first `pi_p3_muon_count_model_challenge_v1.py audit` halted during
the **archived legacy validator**, before training/evaluating any muon
predictions. The recorded failing source is

```text
parent=jpsi40-bore|2|8904|5|2|5|266
geometry=p3_noFW, primary=7,
seed_u_mm=0.00237593484574722, seed_v_mm=0.0107657602328344,
successful=0, capped=0, garfield_ne=1, garfield_ni=0,
ion_births=0, collected_e=1, attached=0, other=0
```

This start has radial distance about **0.0110 mm** from the center of
a 30um-diameter anode: STRICTLY INSIDE its physical metal radius
0.015 mm. It is not a gas transport failure; it is a physically
nontransportable start. A historical simple endpoint radius heuristic
nevertheless classified the one endpoint as collected. The old
`validate_reference` recognized *field-wire* metal starts only and
rejected this anode-metal zero. **Do not change it or the original C++,
CSV or SHA manifests retroactively.** The raw `collected_e=1` here is
not a valid physical collected electron for any future signal study.

New helper `pi_p3_external_new120_metal_provenance_v1.py` checks:
- unchanged original complete Geant4 selection, full SHA-locked 12
  Garfield chunks, SHA of every aggregated original file;
- all 120 actual parents, original TrackHeed cluster→electron membership,
  exact four-geometry Garfield per-electron count and sector sums;
- `successful=0` only for a geometrically **interior** 30um anode or
  pre-existing allowed 75um field-wire metal position, with rigorously
  constrained one-electron/no-birth/no-cap physical outcomes. Unknown
  failures abort. The helper does not modify or overwrite any inputs.

The patched **p3-only** muon challenge uses this helper instead of
`legacy.verify_full()`, retaining zero ion births exactly and reporting
a separate parent-bootstrap **sensitivity exclusion** for any affected
whole Geant4 parents, for each 3mm geometry separately. It continues
to make no claim of signal validation or new blind holdout independence.

### Run after the observed first full cloud pilot 1 / 6

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python

"$PY" test/mwpc-stage-c/pi_p3_external_new120_metal_provenance_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_external_new120_metal_provenance_v1.py audit

# Only after the audit says PASS, retry the model challenge.
"$PY" -u test/mwpc-stage-c/pi_p3_muon_count_model_challenge_v1.py audit --reps 3000 --bootstrap 4000

# Already selected, preflighted, native compiled, and SHA-sealed 8/48 starts:
# DO NOT re-prepare, rebuild the C++ executable, or delete original chunk_0000.
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py run --max-chunks 5 --hours 2
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py status
# Only when status shows 6 / 6 verified:
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py analyze --reps 64
```

These commands operate on already frozen independently generated D500
electron starts and existing new120 muon truth; no signal-response.
The 48-avalanche full-cloud pilot is a **compression fidelity pilot**
and cannot establish rare-tail or whole-muon joint-LUT accuracy.

## 2026-11-04: 48/48 full clouds PASS, analysis strict JSON NaN fix

WSL completed all six immutable full Garfield++ chunks: **48 true
independent exact-start D500 electron avalanches**, **99,388 exact
ionization-birth callbacks**, separately `p3_noFW` and `p3_FW75`.
The follow-up `analyze --reps 64` reached report serialization but
refused a nonfinite float because some nonzero avalanches have no ions
in either the top or bottom ±45° cathode-facing angular sector:
`A_TB=(top-bottom)/(top+bottom)=0/0` is *undefined*, not zero.
The physics data and full-ion moment/sector verification were already
successful; the issue was a REPORT ONLY missing-value serialization.

The Python analyzer is now patched to represent such exact asymmetries
as **JSON null**, exclude these physically undefined values from
conditional asymmetry-error statistics and explicitly count how many
nonzero avalanches have undefined `A_TB` for each geometry. It asserts
every output is finite and strictly JSON compliant BEFORE creating an
output directory. An extra synthetic selftest demonstrates the
side-only birth distribution and undefined asymmetry. The failed v1
`analysis_complete/` directory is retained exactly as found:
**never delete/overwrite it or any physical checkpoints**. The
corrected script creates separate `analysis_complete_v2/` containing
SHA-sealed `summary.json`, `per_avalanche.json`, `manifest.json`.

Run from updated source (NO Garfield or C++ rebuild):

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py status
"$PY" -u test/mwpc-stage-c/pi_p3_full_cloud_convergence_v1.py analyze --reps 64
```

The existing independent external new120 p3 count challenge completed
with 120 muons/geometry and parent-bootstrap confidence intervals:
`gamma - guarded` CRPS = **-109.2, 95% CI [-360.4,+133.2]** on
p3_noFW and **-419.4, CI [-1538.1,+671.4]** on p3_FW75;
neither CI excludes zero. A separate strict audit of 24,772
original microscopic trials found exactly four no-transport
non-ionizing metal starts across all four old reference geometries,
including exactly one 30um-anode-metal start in each of the p3
geometries, with one affected parent event shared by these two.
The sensitivity check excludes *whole affected parents* and retains
119 clean parent events per p3 geometry. The old physical source and
legacy verifier remain immutable.

**No physical signal response is included.** Forty-eight total
avalanches are useful as a storage-compression pilot, not a
proof of 99% tails or cathode-ready LUT accuracy.

## 2026-11-04: full-cloud v2 analysis PASS; independent 32/128/512 readout

The original WSL `analyze --reps 64` completed with sealed
`analysis_complete_v2/`. Final source numbers:
- `p3_noFW`: 24 D500 electron-start avalanches, 18 nonzero, 6 zero,
  **10,634** exact ion birth points. Mean 128-point
  coarse 12 azimuth × 6 time joint total-variation (TV) distances:
  near anode 0.0508, transition 0.0715, bulk 0.0569.
- `p3_FW75`: 24 D500 start avalanches, 19 nonzero, 5 zero,
  **88,754** exact ion births. Mean 128-point joint TV:
  near anode 0.0763, near field wire 0.0709,
  transition 0.0746, bulk 0.0703.

**These are compression diagnostics, NOT 128-point bank acceptance:**
the number of independent nonzero avalanches is 18/19; 64 Monte Carlo
subsamplings per avalanche do NOT count as 64 new independent avalanches.
Some avalanches can have <=128 ions so the 128-point sample is the full
cloud; all such exact cases inflate apparent compression quality and
must be shown separately. The report bins angle/time from exact
full-cloud per-avalanche values, not from an independent spatial model.

New read-only report
`pi_p3_cloud_compression_readout_v1.py` displays the **32, 128, 512**
comparisons for each of the two geometries and each near-wire/bulk
stratum independently; verifies both the final analysis output SHA and
all six FULL microscopic physics chunk SHA seals, including each raw
cloud file. Reports unconditional and genuinely compressed
(`N_ion>cap`) results, ion-birth top/bottom missingness,
angular/time joint TV, absolute angular errors, relative radius
error, time q90 error, and paired improvement of 512 vs 128 with
bootstrap over independent full-Garfield *avalanches*.

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_p3_cloud_compression_readout_v1.py selftest
"$PY" -u test/mwpc-stage-c/pi_p3_cloud_compression_readout_v1.py report --bootstrap 3000
```

No Garfield job, no changes to `analysis_complete_v2`, and no
new output files. The report is intentionally exhaustive so the
independent (p=3mm) geometric decisions can be made after inspecting
the real precision, rather than declaring that TV≈0.06 is sufficient.

**Next experiment decision**: validate compression further with
larger independent complete birth-cloud sets in each p3 geometry;
then train a *joint* conditional generator and test on truly
unseen exact D500 cells and independent real Geant4/TrackHeed muons.
The selected-cell and muon-count model comparisons show near-equal
Guarded A20/Gamma count scores; no count-only Pólya fit by itself
provides real ion-birth space/time positions. Don't pool errors
from noFW/FW75 or claim signal readiness.

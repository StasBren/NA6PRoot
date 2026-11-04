# P3 extended full-ion-birth cloud: compression decision gate

**Date:** 2026-11-04. This analyzes the SHA-sealed extended P3 campaign already completed on the user's WSL. No new Garfield or cathode signal has been run.

## Verified microscopic sample

- 960 distinct frozen original D500 initial electrons: 480 **p3_noFW** and 480 **p3_FW75**.
- 120/120 sealed physical 8-trial chunks; 2,981,533 exactly recorded ionization births.
- **p3_noFW:** 417 nonzero / 63 zero avalanches, 284,272 total births.
- **p3_FW75:** 401 nonzero / 79 zero avalanches, 2,697,261 total births.
- This is a deliberately spatially/sign-stratified electron source, **not** a Geant4 muon/Heed electron incidence PDF. Total birth ratio is NOT a general physical muon-gain ratio.
- Statistical uncertainties are grouped by D500 spatial cell, not ions or 32 repeated reservoir subsamples. Geometries must not be pooled.

## Original five-cap result: 12x6 joint angle/birth-time TV error

| Cap | noFW TV | FW75 TV | noFW genuinely compressed | FW75 genuinely compressed |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 0.0590458 | 0.0733142 | 370/417 | 396/401 |
| 256 | 0.0312440 | 0.0508037 | 309/417 | 391/401 |
| 512 | 0.0125024 | 0.0337481 | 213/417 | 383/401 |
| 1024 | 0.00333627 | 0.0212616 | 92/417 | 367/401 |

1024-point TV 95% cell-bootstrap CIs:
- noFW: [0.0025979, 0.0041237], 278 spatial cells.
- FW75: [0.0202651, 0.0221938], 242 spatial cells.

At 1024, **most noFW clouds are no longer compressed**, while **367/401 FW75 clouds still are**. FW75 ion-birth-weighted 1024 TV is 0.02504. Thus it is premature to select a shared 512-point storage format, or even to certify 1024 for FW75 without a detector-observable-based accuracy threshold. Original compression error is NOT new conditional-LUT generalization error.

## Next safe experiment: high-cap analysis with NO new Garfield

New file: test/mwpc-stage-c/pi_p3_highcap_compression_v1.py

It rechecks original chunk and analysis SHA, streams the exact birth CSV, confirms exact C++ moment/sector counts, and estimates storage compression for 2048 and 4096 reservoir points. It reports each geometry/regime separately, all avalanches, genuinely compressed only, ion-birth-weighted sensitivity, and paired improvements using spatial-cell bootstrap. It **does not** write or delete frozen data.

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python

"$PY" -m py_compile test/mwpc-stage-c/pi_p3_highcap_compression_v1.py
"$PY" test/mwpc-stage-c/pi_p3_highcap_compression_v1.py selftest
set -o pipefail
"$PY" -u test/mwpc-stage-c/pi_p3_highcap_compression_v1.py report \
  --caps 2048 4096 --reps 32 --bootstrap 3000 \
  | tee "$HOME/na6p/p3_highcap_compression_v1.log"
~~~

Optional cap 8192 is available; choose it only after inspecting 2048/4096. The new script has not yet been executed against the actual WSL Garfield archive.

## Next decision, only after the high-cap report

1. Review genuinely compressed and large-N tails, angular/time TV, mean radius, upper birth-time quantile, and one-sided undefined asymmetries. Freeze an accuracy requirement before adopting a storage cap.
2. Choose separate cap or adaptive/hybrid representation for each P3 geometry; evaluate disk/runtime/storage tradeoffs. Preserve exact count, source/time coordinates, zero-ion outcomes, and count-cloud covariance per physical trial.
3. Plan a **new** spatially conditional joint-LUT learning curve with tentative 20/50/100 independent avalanches per 0.1 x 0.2 mm D500 cell/geometry (~36k/90k/180k total for the two geometries), distinct spatial-cell holdouts and independent real Geant4/TrackHeed parent holdouts. These are SCALE OPTIONS, not approved Garfield jobs.
4. Only after held-out joint count/shape/time validity can we separately study ion drift, Shockley–Ramo U/V weighting, electronics, and detector-level hit precision.

**No assertion of validated cathode signal or single-hit resolution follows from these data.**

## 2026-11-04: WSL 2048/4096 high-cap report PASS

The user's WSL executed the new high-cap script on the existing exact
960-avalanche archive and independently SHA-audited all 120 native chunks.
The job completed with `PASS high-cap read-only study; frozen files untouched`.
No new Garfield physical samples and no cathode simulation occurred.

**Across-regime, equal-spatial-cell 12x6 angular/time TV** (mean, 95% spatial-cell bootstrap CI):

| Cap | p3_noFW mean [95% CI] | p3_FW75 mean [95% CI] |
|---:|---|---|
|1024|0.00333627 [0.00256563, 0.00412207]|0.0212616 [0.02025997, 0.02222811]|
|2048|0.000310342 [0.000146759, 0.000485737]|0.0116908 [0.01083416, 0.01249914]|
|4096|0.000003467 [0, 0.000010401]|0.00531799 [0.00471979, 0.00590224]|

**How often compression actually occurs:**
- noFW: >1024: 92/417; >2048: 18/417; >4096: only **1/417** nonzero clouds.
- FW75: >1024: 367/401; >2048: 319/401; >4096: **233/401** nonzero clouds.
- Conditional on actually compressed clouds, the 4096-point FW75
  12x6 TV is **0.00879–0.00928** across four regimes; the low global
  noFW TV at 4096 primarily reflects exact storage for 416/417 clouds.
- FW75 **ion-birth-weighted** TV = 0.0250369/0.0158979/0.00889384
  at 1024/2048/4096. This sensitivity analysis emphasizes high-gain clouds.
- At 4096 FW75, mean radius relative absolute error = 0.00112916;
  q90 birth-time absolute error = 0.000419749 ns (0.420 ps), and
  top-bottom asymmetry absolute error = 0.00380857 **only where defined**.
  Undefined sector asymmetries are excluded, not silently set to zero.
- Paired 2048 -> 4096 FW75 TV improvement = 0.00637279,
  95% cell CI [0.00588062, 0.00686369].

**Decision:** 2048 noFW and 4096 FW75 are candidates for *storage
compression diagnostics*, not final cathode-LUT acceptance. Choosing a
cap requires a predefined tolerance and scrutiny of tail clouds. In
particular, a 12x6 angular/time TV is deliberately coarse, uses each
exact cloud's true time range when scoring, and does **not** verify
longitudinal `w` structure, strip-weighting response, induced current,
spatial generalization or independent Geant4-muon event response.

**Next conservative read-only gate:** compare 4096 with 8192 in the same
exact archived FW75 avalanches using the already supported high-cap
script. This reprocesses the noFW archive too (no extra physics) but
the FW75 result is decision-critical. All original data stay untouched.

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_p3_highcap_compression_v1.py selftest
set -o pipefail
"$PY" -u test/mwpc-stage-c/pi_p3_highcap_compression_v1.py \
  report --caps 4096 8192 --reps 32 --bootstrap 3000 \
  | tee "$HOME/na6p/p3_highcap_8192_v1.log"
~~~

The 8192 comparison is *not yet run*. Avoid commissioning a massive
conditional-LUT production job until compression selection and the
held-out joint-response validation protocol are explicitly agreed.

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

# NA60+ MWPC — blind external new120 geometry & wire diagnostics

Date: 2026-11-04. Results **pasted by user from native WSL**. Read-only diagnostics `pi_external_new120_readonly_diagnostics_v1.py` verified SHA seals on both reference truth and all 480 matched paired-model source rows. No detector physics recomputation or model retuning.

## Main results: per geometry

| Geometry | A20 90%-PI coverage | Guarded A20 90%-PI coverage | A20 MAE ions | Guarded A20 MAE ions | A20 PI width | Guarded A20 PI width |
|---|---:|---:|---:|---:|---:|---:|
| p3_noFW | 102/120=85.000% | 105/120=87.500% | 4206.12 | 3443.84 | 14485.21 | 12983.23 |
| p3_FW75 | 108/120=90.000% | 114/120=95.000% | 34940.90 | 29645.61 | 153374.97 | 136672.72 |
| p4_noFW | 105/120=87.500% | 109/120=90.833% | 10589.79 | 9264.90 | 38806.44 | 35335.91 |
| p4_FW75 | 105/120=87.500% | 104/120=86.667% | 37249.05 | 31565.97 | 142439.81 | 130220.30 |
| All four paired per 120 physical parents | 420/480=87.500% | 432/480=90.000% | 21746.46 | 18480.08 | 87276.61 | 78803.04 |

**Critical finding:** Guarded A20 lower MAE *and* narrower 90%-PI width in **all four geometries**, but the coverage improvement is heterogeneous. `p4_FW75` coverage falls by exactly one case/120, from 87.50% to 86.67%, so 90.00% aggregate does not mean each geometry has 90% population calibration. `p3_FW75` coverage 95.00% could reflect some overcoverage; at n=120 one must compute confidence intervals before calling it overcalibrated. Do NOT retune pooling on this new external holdout (it would cease being an independent test).

## Paired 120-parent clustered uncertainty from actual native audit

- Guarded A20 minus A20 coverage **+2.500 pp**, 95% parent cluster bootstrap percentile **[+0.2083,+4.7917]** pp.
- Guarded A20 versus A20 predicted-mean absolute error reduction **3266.3822 ion births**, interval **[1917.0006,4705.6842]**.
- Guarded A20 narrower 90%-PI by **8473.5692 ions**, interval **[5971.1802,11381.2913]**.
- Guarded A100 versus Guarded A20: **+0.4167 pp** coverage, percentile CI **[0,+1.0417]**; MAE reduction **12.6302 ions**, CI **[-83.9901,+113.1811]** (includes zero); nominal PI width reduced 134.9355 ions, CI [58.0109,217.8596]. In particular the extra targeted A100 16k trials have **no demonstrated improvement of predictive-mean accuracy over guarded A20**, despite slightly narrower intervals.

## Source exposure

- Among 480 correlated parent×geometry rows, **264 contain at least one primary in an A100-targeted cell**, 216 none. In the latter 216, baseline A20 and A100 results are numerically identical as required by CRN; guarded A20 and guarded A100 likewise identical.
- In the 264 affected rows, guarded A20 238/264=90.15% and guarded A100 240/264=90.91% coverage. MAE 17184.96 vs 17162.00 ions. This verifies that targeted A100 yields extremely modest marginal benefit in this source-selected ensemble.

## The near-wire comparison needs a better statistical grouping

Posthoc analyzer `pi_external_new120_readonly_diagnostics_v1.py` printed true Heed seed fractions within 0.35 mm of the nearest wire:
- p3_noFW **2.1153%** of 6,193 primary×geometry starts; physically in-metal seeds **1**.
- p3_FW75 **3.0841%**, in-metal **1**.
- p4_noFW **2.4382%**, in-metal **0**.
- p4_FW75 **3.6493%**, in-metal **2**.

**The median PER-PARENT fraction is ZERO in all geometries.** Its diagnostic split into a sorted first 60 vs second 60 cannot be interpreted as a clean comparison of 'near wire' versus 'far from wire': many parents in both halves can have no primary in the defined 0.35mm zone. Sorting tie-breaks by parent identity, not detector physics. DO NOT draw a wire-specific performance claim from those half-sample percentages. Next test should explicitly classify `n_near_wire>0` vs `=0`, report parent counts, ion-birth tail cases and sparse sample uncertainties. As only ~2–4% of individual primaries are near wire, these results may be underpowered to test metal/large local gradients.

## Physics decisions

1. Retain **Guarded A20** as the best-supported fast avalanche **ion count/sector** emulator on this balanced J/psi Geant4 source, with original A20 and Guarded A100 immutable controls. Do not spend several days filling entire A100 LUT based on these data.
2. Do not claim exact 90% calibration in all geoms, source distributions or near-wire physical phase space. Four geometries per parent are correlated, training uncertainty is approximated by shared global cell-bootstrap, and one true Garfield realization per new Heed seed remains stochastic.
3. Near-wire explicit nonzero exposure and high-ion tails need new read-only diagnostics; any retuning after using external truth would require another independent holdout.
4. Future NA60+ MWPC model requires spatial/time ion birth distribution and subsequent ion drift/weighting fields/induced U/V cathode signals with electronics response. The current `ion_up/down/top/bottom` counts are not U/V strip charge or waveforms.

Native sealed sources are under `test_runs/mwpc_stage_c/pi_external_new120_guarded_lut_blind_v1/{new_full_garfield_120,blind_evaluation_fixed_guard_v1}`, original physics untouched.

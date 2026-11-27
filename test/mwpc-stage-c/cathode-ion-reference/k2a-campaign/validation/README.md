# K2A versus detailed S0 on independent non-sealed muons

Goal: K2A predicts the detailed cathode response from compact avalanche moments. It is checked against the frozen S0 workflow on the same archived 2100 V K256 ion births, with quantified errors, domain coverage and a measured speedup. No avalanche is regenerated. The frozen tools, kernel, LUT, archive and sealed set are not changed.

## What is compared

| name | role |
|---|---|
| **S0_T1** | Established E2b reference: frozen RKF paths, per-primary T1 birth time, frozen weighting replay, plus the **approximate** instantaneous electron step. |
| **S0_T0** | Same groups, true per-ion birth times (uncollapsed export), same replay and the same **approximate** electron step. It is not a microscopic electron-current simulation. |
| **K2A** | **Primary prediction.** Archive → `p3_birth_moments` → compiled `k2a_strip_eval` using its **own** production grouping/recentering rule. It receives nothing from S0. |
| **K2Aimp** | **Diagnostic only.** The same moments with the S0 grouping imposed (`--shifts groups.json`). It isolates kernel and interpolation error at identical windows and is never evidence that the prediction is accurate. |

The per-primary `(m_U, m_V)` assignments of K2A and S0 are compared (`k2a_shift_assignment.py`):
- **AGREE:** K2A and K2Aimp must be bit-identical.
- **DISAGREE:** the primaries and their charge are listed, and the effect is measured (`cmp_K2Aimp_vs_K2A`, `shaped_refK2Aimp_ionelec`). All variants are first put on one zero-padded absolute-strip label set (`k2a_label_union.py`).

Equivalence tests:
- compiled K2A on the export-path moments with the S0 grouping == Python `K2A_local` (1e-12);
- archive moments == export moments (`s_compare_moments`, 1e-9);
- the responses on those two moment sets agree (1e-7).

Error components:
- **model** = S0_T1 vs K2A;
- **total** = S0_T0 vs K2A;
- **timing** = S0_T0 vs S0_T1;
- **model_imposed** = S0_T1 vs K2Aimp (diagnostic);
- **grouping** = K2Aimp vs K2A.

Shaping uses the provisional CR-RC^3 stand-in, not VMM3a.

## Execution is not accuracy

**Execution outcome** per muon:
- VALID or VALID_T1_ONLY;
- OUT_OF_DOMAIN: the primary K2A blocks; the statuses are named and S0 is still computed;
- S0_FAILED;
- FAIL;
- INCONCLUSIVE.

VALID means the references and the prediction were computed and are consistent. It does not mean the prediction is accurate.

**Accuracy** is a separate PASS/FAIL: the frozen sealed criteria v1, unchanged, applied to the primary K2A against S0_T1. It is diagnostic on these non-sealed events. Slots selected as the out-of-domain stratum are listed separately and never enter the accuracy set.

## Runtime

Both sides start from the same archived Garfield++ avalanche. Avalanche generation is excluded and no acceleration of it is claimed. The report separates:
- `cathode_response_eval`: S0 response with the archive read subtracted approximately, over the K2A warm kernel evaluation alone;
- the K2A moment construction;
- archive I/O (`end_to_end_from_archive_incl_io`);
- the one-time kernel open.

## Files

| file | role |
|---|---|
| `k2a_validation_select.py` | Stratified pilot (5) and full (16) selection from the non-sealed preflight. |
| `k2a_validation_env.sh` | Resolves and checks the WSL inputs (run, preflight, recovery output, `p3_birth_moments`, inline build, archive) into an env file for review. |
| `run_k2a_validation.sh` | Per muon: S0_T1 (far-mode experiment), S0_T0, then `k2a_validation_k2a.sh`, the outcome and an optional prune. |
| `k2a_validation_k2a.sh` | K2A stage: moments, primary K2A, diagnostic K2Aimp, assignment comparison, equivalence tests, comparisons, shaped observables. |
| `k2a_validation_muon.py`, `k2a_validation_report.py` | Per-muon outcome, accuracy and runtime; aggregate report (`--from` combines runs). |
| `k2a_shift_assignment.py`, `k2a_label_union.py`, `k2a_validation_lib.sh` | Helpers. |
| `tests/test_validation_synthetic.sh` (+ `validation_negatives.py`, `assignment_units.py`, `disagreement_case.py`) | SYNTHETIC tests, including agreement and disagreement of the production and imposed shifts. |

## Disk

The run keeps four ion-path copies per muon. `MWPC_K2A_VAL_PRUNE_VERTICES=1` removes only the three derived copies inside OUTDIR once a muon's outcome is written. Their sha256 is recorded, and they are rebuilt bit-identically by rerunning.

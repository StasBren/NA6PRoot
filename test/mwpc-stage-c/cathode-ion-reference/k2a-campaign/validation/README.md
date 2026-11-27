# K2A versus detailed S0 on independent non-sealed muons

Goal: K2A predicts the detailed cathode response from compact avalanche moments. It is checked against the frozen S0 workflow on the same archived 2100 V K256 ion births, with quantified errors, domain coverage and a measured speedup. No avalanche is regenerated. The frozen tools, kernel, LUT, archive and sealed set are not changed.

| file | role |
|---|---|
| `k2a_validation_select.py` | Deterministic, stratified choice of a 3–5 muon pilot and a 10–20 muon set from the preflight. It reads the non-sealed moments only. It excludes the sealed set, the kernel muons of RUN, the development-decision slots, previously inspected and EXCLUDED slots. |
| `run_k2a_validation.sh` | Runs per muon: S0_T1 (far-mode experiment, instrumented), S0_T0 (true birth times, same groups), production K2A (`p3_birth_moments`, then `k2a_strip_eval`), the consistency checks, the comparisons, `t_shaped_compare`, the outcome and an optional prune. |
| `k2a_validation_muon.py` | Gives each muon an explicit outcome: VALID, VALID_T1_ONLY, OUT_OF_DOMAIN, S0_FAILED, FAIL or INCONCLUSIVE. It also writes the metrics, the runtime split and a DIAGNOSTIC dry run of the frozen sealed criteria v1. |
| `k2a_validation_report.py` | Aggregates the muons: outcome table, accuracy distributions for model / total / timing, feature dependence, coverage, runtime and measured speedup. |
| `tests/test_validation_synthetic.sh` | SYNTHETIC end-to-end test (toy run, fake exporter, stand-in moments reader), including negatives. |

Error components:
- **model** = S0_T1 vs K2A: kernel, recentering and spatial error at the same T1 timing.
- **timing** = S0_T0 vs S0_T1: the T1 time collapse alone.
- **total** = S0_T0 vs K2A.

Shaping is the provisional CR-RC^3 stand-in, not VMM3a.

Runtime:
- The RKF export is the frozen exporter, which does the archive read and the Garfield RKF drift in one process. The archive read alone is measured by `p3_birth_moments`.
- "Single pass" counts only the chosen late-step attempt.
- The K2A per-event time is moments plus a warm kernel evaluation.
- Only measured times are reported.

Cost per muon:
- one RKF export;
- one frozen replay plus electron step per shift group for T1, and again for T0;
- for in-window muons, the base-centred S0 replay of the far-mode experiment, as well.

Disk: four ion-path copies per muon (T1 export, kept frozen export, T1 groups, T0 groups). `MWPC_K2A_VAL_PRUNE_VERTICES=1` deletes only the three derived copies inside OUTDIR once the muon's outcome is written. Their sha256 is recorded, and they are rebuilt bit-identically by rerunning.

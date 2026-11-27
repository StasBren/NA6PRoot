# Sealed-set acceptance criteria, version 1 (FROZEN)

Frozen on 2026-11-27, before any sealed event was opened. The machine-readable definition is
`sealed_acceptance_criteria_v1.json`; this page explains it. Both files are hashed in `FREEZE_v1.sha256`, and
`check_criteria_freeze.py` verifies them. If the two ever disagree, the JSON is authoritative.

Changing anything requires a version 2. A version created after the sealed set has been opened means the set is no
longer an independent validation.

## What is tested

**Comparison.** K2A against S0, per event, on the 20 sealed events (10 φ, 10 ω):
- K2A is the compiled absolute-strip adapter with the frozen kernel and per-primary recentering.
- S0 is the frozen ion-only replay on the frozen export, as local per-group S0 summed on absolute strips (as in E2b),
  plus the approximate prompt electron term.

**Same inputs on both sides.** Both are computed from the same frozen K256 archive of each event, so the birth sample
is identical. No avalanche is generated for this validation. This isolates the K2A model error, which comprises:
- the kernel compression;
- the moment representation;
- the recentering;
- the electron-step approximation.

The full-birth in-memory path is validated separately, on development events (T1). The K256-vs-full sampling difference
is a separate diagnostic (T2).

**Signal and electronics.** The decisive signal is the ion term plus the electron step. Ion-only is reported.
Electronics are the provisional stand-in, unchanged:
- shaper: CR-RC³, Tp = 100 and 200 ns;
- thresholds: fixed from the S0 lead peak at 2, 5, 10, 20 and 50 %;
- tool: `t_shaped_compare.py`, sha256 `d4b1a73e…`.

## Criteria

Each criterion must hold for **every in-domain sealed event**, in **both planes** where it is defined per plane, and at
**both Tp**.

| observable | limit (K2A − S0) | exemption (listed, not a failure) |
|---|---|---|
| lead-strip identity | no change | S0 top-two peaks within 1 % (ambiguous) |
| hit flips at 5 % and 10 % thresholds | 0 | every flipping strip has its S0 peak within ±2 % of the threshold (marginal) |
| threshold time at 10 % | \|Δt\| ≤ 1 ns | — |
| window charge per plane | \|ΔQ\|/Q ≤ 0.5 % | — |
| U/V charge ratio | \|Δ\| ≤ 1 % relative | — |
| second/lead | \|Δ\| ≤ 0.005 | — |
| w centroid | \|Δw\| ≤ 10 µm | — |
| u centroid | \|Δu\| ≤ 50 µm | — |
| shaped rel L2 per plane | ≤ 1 % (ceiling) | — |

These are the draft values of `K2A_CAMPAIGN_PLAN.md` section 5, unchanged. Each maps one to one onto an output key of
`t_shaped_compare.py`, listed in the JSON. Both exemptions are defined from the S0 reference only, so they cannot be
tuned on K2A.

## Domain and coverage

- An event is **in domain** if the adapter evaluates every cloud with no blocking status: frozen kernel, no table
  enlargement, no clamping.
- An out-of-domain event is listed with its cause and charge. It is not compared against the criteria.
- Coverage requirement: **at least 17 of 20 events in domain**. Development coverage after recentering is 929/977 =
  95.1 %, and at that rate fewer than 17 of 20 would occur with probability 1.5 %.

## Integrity (not model criteria)

These checks must hold:
- archive SHA256SUMS;
- exact census (charge relative difference ≤ 1e-9, integer group census);
- the frozen replay's time-compression pilot, using the frozen late-step ladder;
- kernel, transport and weighting ids;
- identical label sets and output clock on both sides;
- frozen inputs unchanged.

A failure makes the run **INVALID**. That is not a model result:
- The infrastructure is repaired on non-sealed data only. Criteria, model, kernel and thresholds are not touched.
- The whole set is then rerun.
- Every attempt is kept and reported.

## Outcome

Each event is PASS, FAIL (with values), OUT_OF_DOMAIN (with cause) or INVALID.

- **ACCEPTED:** no INVALID event, ≥ 17/20 in domain, and every in-domain event PASS.
- **NOT ACCEPTED:** any in-domain FAIL, or fewer than 17 in domain.

Results are reported per channel (φ, ω) and for all 20. Every event appears in the table. Diagnostics are reported but
never decide:
- unshaped L2;
- ion-only;
- multiplicities;
- all threshold fractions;
- alternative centroids.

## Procedure

1. Validate the criteria evaluator on development (non-sealed) events.
2. Write the run-configuration freeze record, which must cite this criteria file's sha256 unchanged. It records:
   - kernel, transport and weighting ids;
   - the adapter configuration (K2A, ±5 strips, 7000 ns, 10 ns);
   - the electronics;
   - the git commit and the sha256 of every tool;
   - the sha256 of `sealed_validation20.csv`.
3. The sealed runner refuses to start unless the freeze record matches this file and the checked-out tools.
4. Run all 20 events. Nothing is excluded, re-selected or rerun selectively.
5. No retuning afterwards. A later model, kernel, recentering or electronics change spends the sealed set.

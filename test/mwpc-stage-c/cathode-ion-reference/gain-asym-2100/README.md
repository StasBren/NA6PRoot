# 2100 V / 997-muon exact ion-birth hemispheric asymmetry

**Scope:** original NA60+ Garfield++ 2100 V microscopic archives. ROOT/C++ only, read-only; no Geant4, Heed, Garfield++, ion drift, FEM, strip response, 1800 V comparison or simulation reruns.

## Find the original archived run

Expected WSL directory:

```
$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/
```

Each sealed `muon_NNNN` folder typically contains:
- `ion_births.root` (trees `muons`, `clusters`, `primaries`, `ions`; K256 sampled **ion birth positions and creation times**, exact per-primary counters).
- `microscopic_seed.csv` (original per-Heed-primary ion and collected-electron census).
- `microscopic_wire_seed.csv` and `microscopic_wire_muon.csv` (actual collecting-wire electrons and exact nearest-ion-birth-wire sectors, including `A_ion_hemi`).
- `complete.txt`, `root_audit.log`, `SHA256SUMS` (archival integrity controls).

The November 21–23 development recorded **997/1000** such verified ROOT archives, plus one failed and two unsealed slots; actual current count must be checked on your WSL machine. The source archive **does not contain RKF ion drift trajectories or U/V cathode waveforms for all 997**. Those are separate, later replay outputs (notably muon 0301). `locate` checks the archive and likely later response directories without scanning the huge original microscopic ion data.

## Commands

From a NEW worktree/clone of this branch on the same WSL machine, run:

```bash
RUNNER=test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_gain_asym_2100.sh
bash "$RUNNER" locate
bash "$RUNNER" prepare
bash "$RUNNER" pilot
bash "$RUNNER" population
```

- `locate`: check all original 1000 slots for sealed ROOT + audit log and inspect likely response-output paths; no modification.
- `prepare`: compile the **existing** vetted `mwpc_v0_archive_stats` ROOT/C++ reader.
- `pilot`: process five SHA-rechecked events as a software smoke test.
- `population`: process **all 997 seal/audit-pass archives** into a fresh, timestamped output directory. If the count differs, it stops rather than silently changing the sample.

The runner does not rehash all 997 files (that can be expensive); for an independent complete SHA audit, use the original campaign status/audit tools. The V0 reader checks ROOT schema, census and physical-wire counters and will abort on inconsistencies.

Outputs:
- `.../census/v0_primaries.csv`, `v0_muons.csv`, `v0_fragments.csv`, `v0_wires.csv` and other existing V0 products.
- `.../figures/01_*.png/.pdf` through `09_*.png/.pdf`: ROOT/ATLAS-style distributions, 2D correlations, signed and absolute A profiles for individual primaries and entire muons.
- `.../figures/gain_asym_2100.root`: ROOT histograms and `TGraphErrors` objects.
- `.../figures/gain_asym_profiles.csv`: bin occupancies and muon-block bootstrap uncertainty.
- `.../figures/summary.txt`: counts and **global mean** ionization gain as ratio of sums.

## Physics definitions and limitations

Garfield's original sector convention is fixed:

```
phi = atan2(v, u - nearest_anode_u)
up: 0 < phi < pi ; down: -pi < phi < 0
A_UD = (N_up-N_down)/(N_up+N_down)
```

Exact dividing-axis births are excluded **only** from the A denominator; they remain in the total ionization gain. `A_UD` is ion **BIRTH** hemispheric asymmetry, not collected cathode charge, and not induced U/V signals.

- Per-Heed-primary `G_ion` is its exact avalanche ion-birth count.
- Per-muon `G_ion` is total exact ion births / number of Heed primaries.
- Dataset-wide mean is total exact ion births / total Heed primaries.
- Zero-gain primaries remain in the gain distribution, but have undefined A.
- Both asymmetry levels use **exact original census** from `microscopic_wire_seed.csv` via V0; the K256 reservoir is **not** used to infer exact asymmetry.
- Individual avalanches share a parent muon. Binned errors are 250 muon-level cluster-bootstrap repetitions. Bins with <20 entries or <4 muons are not profiled.
- The full 997 descriptive sample **includes sealed model-validation slots**, so its outputs are explicitly **NOT FOR CALIBRATION**.

The preserved detector model is an idealized 2100 V P3 no-field-wire geometry, 30 µm anode diameter, 3 mm wire pitch, ±3 mm cathodes, Ar/CO2 70:30 and Penning r=0.57, B=0. This is not the as-built Prototype-3.

The existing V0 reader also reports true physical collecting-anode electrons, separate from the nearest ion-birth anode. The present branch focuses on reliable gain–asymmetry plots for the presentation; the full 27-observable report additionally needs angular and longitudinal reservoir-sampling analysis and charge-sharing plots. No ion drift is necessary to obtain the ion-birth/gain results.

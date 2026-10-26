# Garfield 500-trial p=3 mm / field-wire convergence experiment

Isolated microscopic Garfield LUT study for exactly **two** geometries:
p=3.0 mm, 30 μm anodes, with and without a grounded 75 μm field wire.

- 30 u cells (0.1 mm each) × 30 v cells (0.2 mm each): **900 spatial cells/configuration**.
- 500 independent Garfield microscopic single-electron drift trials per box.
- 450,000 Garfield trajectories per configuration, 900,000 across both.
- Frozen physical settings: 3+3 mm gas gap, Ar/CO₂ 70:30, +1800 V anodes, B=0.
- Dedicated RNG seed **260500**, entirely distinct output namespace
  `test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/pilot/`.
- The original 20-trial map, the separately running 50-trial all-six suite,
  and the 720 frozen 6M-parent Heed+LUT folds are **never changed**.

## Start / monitor / verify

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh check
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh start

# At any moment (read only):
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh status
tail -n 20 test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/runner.log

# When both configurations have finished:
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh verify
```

The runner uses `nice -n 15` and complete immutable u-column checkpoints.
There are **60 u-column checkpoints**, each containing 15,000 independent
Garfield electron trials. Safe interruption does not invalidate completed
columns: repeating the `start` command validates and skips saved columns.
Do not edit the `frozen_config.json` or remove files from the running output.

`status` only counts present columns and checks immutable physics identifiers.
`verify` audits each column's contents and SHA-linked timing sidecars plus the
two assembled map outputs and electron conservation. All generated data are
**microscopic single-electron Garfield LUT samples**, NOT 900,000 Geant4 events.


## First PI question: quantitative LUT angular fluctuation test

After `run_pi_pitch_fieldwire_500.sh verify` passes, use the isolated,
**read-only** 20-vs-500 report. It does not need Garfield, Geant4 or the old
720-fold pipeline and does not modify frozen checkpoints:

```bash
python3 test/mwpc-stage-c/test_pi_lut_20_vs_500_convergence.py
python3 test/mwpc-stage-c/pi_lut_20_vs_500_convergence.py --station 0 --chunks 2
explorer.exe "$(wslpath -w test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_v1/compare.html)"

# After the quick first-two-chunk validation, full existing 2M J/psi exposure:
python3 test/mwpc-stage-c/pi_lut_20_vs_500_convergence.py --station 0 --chunks 40
```

Outputs: `compare.html`, `phi_comparison.csv`, `comparison.json`. The
comparison uses exact deterministic Heed **birth positions** from the same
native J/psi source/chunks and folds their spatial cell occupancies over
**all** stored Garfield outcomes at each LUT cell. This avoids an expensive
electron-by-electron repeated random LUT lookup. Finite physical anode
boundaries and relative Garfield wire migration are retained, including
Heed births just beyond the outermost nearest physical anodes.

The error estimate is conditioned on fixed Heed spatial occupancies and
comes from K independent microscopic Garfield trials *per spatial cell*,
with reuse of each trial at different physical wire locations treated as
correlated. It is **not** a Geant4-event statistical error and not a
substitute for independent repeat LUT seeds. Differences between 20 and
500 should be compared with the displayed Garfield uncertainty, not
claimed to be proof of total physical convergence.

The 20-vs-500 integrated PDF has the same expected electron response as
the frozen old LUT resampling approach but will not yield identical old
integer counts: the new method analytically averages over the finite LUT
sample rather than performing per-electron pseudo-random re-draws.
Likewise, the `--chunks 2` exposure has different Heed cell weights from
the full 40-chunk dashboard; use `--chunks 40` for source-matched plots.

## Existing 720-shard production progress (independent)

```bash
bash test/mwpc-stage-c/run_pi_weekend_2m.sh status
python3 test/mwpc-stage-c/pi_2m_lightweight_progress_snapshot.py
```

Expected progress: 120 Geant4 original chunks, 120 ROOT native exports and
720 Heed/six-Garfield-response station folds (3 sources × 40 source chunks
× 6 spectrometer stations). The status command counts markers; lightweight
snapshot requires both checkpoint manifest and gzip. Neither is a full
SHA256 verification of all 720 gzip files.

## What this does NOT do

This **does not** regenerate or improve the already-produced J/ψ 2M
polar histograms. Those are folded using a frozen 20-trial map. To compare
J/ψ muon-conditioned φ_hit distributions against the new 500-trial map,
perform a new **isolated** Heed+LUT fold using the same Geant4 source/Heed
seeds into a *new* output directory. No Geant4 rerun is necessary, but the
re-fold consumes CPU and disk. Never overwrite the original 720 checkpoints.

An independent LUT study first isolates the hypothesis that the 20
microtrajectories/cell cause visible per-bin noise; it does not prove it
until a controlled convergence test with independently seeded LUTs and
uncertainties has been performed.

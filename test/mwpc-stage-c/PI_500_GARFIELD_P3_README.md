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

Use the **fast analytically expected Heed occupancy** by default. The former
same-RNG sampled-Heed mode is left unchanged in the original script.
The fast mode integrates the **mean Heed electron density along each actual
straight native Geant4 gas chord** and partitions its length EXACTLY at every
0.1mm-u and 0.2mm-v LUT grid boundary, with periodic physical-anode
indexing and finite detector edges. This removes millions of per-event
NumPy RNG, Poisson, mark, and electron-list operations. Geant4 gas chords
and Heed cluster rate and multiplicity mean remain unchanged. This mode
answers **LUT-only angular convergence at fixed Geant4 muon exposure**;
it does not reproduce a particular sampled Heed experiment, and must not
be used for event-by-event asymmetry or Heed fluctuation calculations.

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/test_pi_lut_20_vs_500_convergence.py
python3 test/mwpc-stage-c/test_pi_lut_20_vs_500_expected.py

# Quick speed check, 2 existing 50k-parent Geant4 chunks:
python3 test/mwpc-stage-c/pi_lut_20_vs_500_expected.py --station 0 --chunks 2

# HTML output (distinct from original sampled-Heed mode):
explorer.exe "$(wslpath -w test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_expected_v1/compare.html)"

# All 2M already generated J/psi parents, still without new Geant4/Garfield:
python3 test/mwpc-stage-c/pi_lut_20_vs_500_expected.py --station 0 --chunks 40

# Optional exact original Heed RNG realization (much more CPU intensive):
python3 test/mwpc-stage-c/pi_lut_20_vs_500_convergence.py --station 0 --chunks 2 --max-hits 2000
```

Each mode creates `compare.html`, `phi_comparison.csv`,
`comparison.json` in **distinct** output directories.
The uncertainty is from independent microscopic Garfield drift trials
per LUT cell, weighted by the fixed Heed spatial occupancy; physically
distinct wire locations reusing the same LUT trial are treated as
correlated. The fast mode further conditions on the *mean* Heed
occupancy, not sampled Heed shot noise. These are NOT independent
Geant4 parent-event uncertainties, and they do not quantify systematic
grid/coarse-graining or detector-model bias. Same event/chunk Geant4
track ensemble and same expected Heed physics are used for 20 and 500;
different microscopic LUT seeds are independent.

The original 720 native production folds, original 20-trial LUT and both
500-trial Garfield input archives are read-only. The new script does
not start Garfield or ROOT. It also does not overwrite the old dashboard.
This expectation-based first-PI-question diagnostic is intentionally
**different** from the second PI question about per-event up/down
asymmetry, which needs seeded individual electron responses and event IDs.

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


## PI single-slide ROOT-style presentation (one configuration)

The script reads **only existing** measured phi_comparison.csv from the
mean-Heed first-two-chunk convergence test. It draws two vertically stacked
ROOT TH1D histograms in memory (LUT20 above, LUT500 below) with matching axes,
1-sigma Garfield-LUT sampling uncertainties, and the nominal factor-5
scaling law from K=20 to 500 independent microscopic trajectories per cell.

The script writes **only PNG and PPTX**, not .root, PDF or an extra JSON
provenance file. It does not rerun Garfield, Geant4, or Heed. python-pptx
must be available in the local WSL interpreter for the PPTX output.

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 -m py_compile test/mwpc-stage-c/pi_lut_root_slide.py
python3 test/mwpc-stage-c/pi_lut_root_slide.py --validate-only
python3 test/mwpc-stage-c/pi_lut_root_slide.py

# Display the PNG with Windows:
explorer.exe "$(wslpath -w test_runs/mwpc_stage_c/pi_lut_root_slide_v1/p3p0_field75p0um_slide.png)"

# If python-pptx is absent in WSL, install it there:
# sudo apt update && sudo apt install python3-pptx

# Optional alternate electrode geometry without FW:
python3 test/mwpc-stage-c/pi_lut_root_slide.py --geometry p3p0_anode_only
```

PNG and PPTX outputs are in
`test_runs/mwpc_stage_c/pi_lut_root_slide_v1/`:
`p3p0_field75p0um_slide.png` and
`p3p0_field75p0um_slide.pptx`.
The PowerPoint slide embeds the ROOT-rendered image; original numerical
CSV and Garfield LUTs remain in their read-only sources.

**The fivefold reduction is the expected standard-error scaling** for
independent microscopic trials per LUT cell, not an exact numerical
fivefold reduction of observed neighbouring-bin differences. This is
a fixed-Geant4, analytically expected-Heed diagnostic, not the
event-by-event asymmetry calculation.

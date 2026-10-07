# MWPC Stage C — chamber response integration

This directory is the **canonical working area for Stage C**.

The branch `feature/mwpc-stage-c` was created from
`feature/mwpc-avalanche-visualization`, so the mature Stage-A and Stage-B
work is already present on the same branch under `test/mwpc-avalanche/`.
Do not reimplement those pieces inside Stage C.

## One-branch workflow

Use only this branch for new Stage-C work:

```bash
git switch feature/mwpc-stage-c
git pull --ff-only
```

Older branches are historical development branches:

- `feature/mwpc-avalanche-phase-a`
- `feature/mwpc-avalanche-visuals`
- `feature/mwpc-avalanche-visualization`
- `feature/mwpc-avalanche`

They are useful for history, but new Stage-C commits should not be split across
them.

## Reuse from Stage A — do not duplicate

Canonical code remains in `test/mwpc-avalanche/`.

Important reusable pieces:

- `MWPCCoordinateFrame.h`
  - local `u`: across anode wires
  - local `v`: chamber normal
  - local `w`: along anode wires
  - Garfield coordinates are `(x_G,y_G,z_G)=(u,v,w)`
  - nominal mapping is global `Y -> u`, global `Z -> v`, global `X -> w`
- `phaseA_coordinate_test.cxx`: round-trip and basis checks
- `phaseA_heed_muon.cxx`: Heed clusters and conduction-electron handoff
- `phaseA_heed_avalanche.cxx`: full Heed -> microscopic avalanche chain
- `phaseA_magnetic_drift.cxx`: magnetic-drift response
- gain/collection/statistics utilities already used in the Phase-A studies

Stage C must call or factor these implementations rather than create a second
coordinate transform, Heed model, or avalanche model.

## Reuse from Stage B — do not duplicate

The same `test/mwpc-avalanche/` directory already contains:

- strip weighting-potential/weighting-field tests;
- electron and ion Shockley-Ramo signal calculations;
- real microscopic avalanche -> strip signal;
- finite-time signal integration;
- strip sharing and CoG reconstruction;
- strip pitch/width scans;
- tilted-strip and two-family reconstruction;
- opposite-cathode studies;
- avalanche azimuth / cathode-sharing diagnostics.

The Stage-C task is therefore **integration**, not another strip-response model.

## Stage-C chain

The intended event-level chain is

```text
saved/full Geant4 muon crossing
  -> global position/momentum
  -> existing ChamberFrame global -> local transform
  -> Heed full muon ionisation
  -> many conduction electrons
  -> Garfield drift + many microscopic avalanches
  -> existing Stage-B strip-response kernel for each avalanche
  -> sum strip currents/charges over the full muon event
  -> strip cluster / CoG / timing observables
```

No electronics/noise/threshold/VMM model is included yet.

## C0 — accepted incident phase space

`analyze_saved_bore_angles.py` reads the persistent high-statistics MNP33-bore
production generated in the separate `na60-dice-mwpc-simulation` repository.
It does **not** rerun Geant4.

The production contains 200k parent dimuons per source for J/psi, Phi and
Omega, split into two chunks/source, with transported momentum recorded at
MS0..MS5 and at MNP33 reference planes z = 365, 430, 495 cm.

Selections:

1. `all`
2. `geometry_bore` — analytic MNP33 opening
3. `useful_bore` — 245 x 240 cm useful aperture proxy
4. `useful_bore_working_area` — useful bore plus finite station working area

The analysis reports the Stage-C local incidence angles explicitly:

```text
theta_u = atan2(p_u, p_v) = atan2(pY_global, pZ_global)
theta_w = atan2(p_w, p_v) = atan2(pX_global, pZ_global)
```

for the present nominal station orientation.

This distinction matters: the large downstream magnetic-bending tail is mainly
in `theta_w`, i.e. **along the anode wires**, while the across-wire angle
`theta_u` is much more compact.

## Reference configuration

`stage_c_reference.json` contains only the compact Stage-C analysis defaults:
useful MNP33 aperture, approximate working rectangles, source labels and the
nominal local/global axis mapping.

It is not a replacement for the detector geometry parameters in NA6PRoot.

## Running C0

The helper script accepts an explicit production suite:

```bash
bash test/mwpc-stage-c/run_c0_angles.sh \
  "$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/mnp33-bore-runs/20260912-202112-0wrxhdg_"
```

If no argument is supplied, it tries the known high-statistics suite under
`$HOME/na60-dice-mwpc-simulation`.

Runtime products go under

```text
test_runs/mwpc_stage_c/c0_angles/
```

and are intentionally not committed.

Small, reviewed reference summaries that we want to preserve for the project
belong under `test/mwpc-stage-c/results/`.

## Next implementation step: C1

C1 should add one executable that consumes a **single local muon crossing** and
runs the already validated Stage-A + Stage-B physics chain end-to-end.

Before writing new physics code, first factor reusable Stage-A/B routines out of
the existing executables where necessary. The goal is one implementation of
each physical operation, shared by the old validation programs and Stage C.

The first C1 observables should be:

- number of Heed clusters and primary electrons;
- active wire multiplicity;
- avalanche charge per wire;
- per-strip integrated charge at a chosen observation time;
- reconstructed local coordinate(s);
- residual/bias versus the true local crossing;
- dependence on `(p, theta_u, theta_w)`.

## Scope boundary

Stage C does not yet claim final detector resolution. Absolute gas gain is not
prototype-calibrated, and electronics/noise/threshold/shaping are still absent.

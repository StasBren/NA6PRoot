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

## C1a — first full-muon event response

### Current Stage-C chamber baseline

For new C1 response work the default geometry is the current preferred
two-sided MuTr baseline rather than the asymmetric Prototype-3 sandbox:

```text
gas gap       : 2.5 + 2.5 mm
strip width   : 1.7 mm
tan(alpha)    : 0.1
readout       : two-sided stereo cathode strips
```

The historical Prototype-3 2+4 mm geometry remains available through explicit
runtime arguments for regression studies. The Stage-A/B 1.8 kV operating
voltage is retained for integration continuity only; gain/HV for the symmetric
gap is not yet calibrated.



C1a is now implemented in

```text
c1_full_muon_response.cxx
```

and built as

```text
mwpc_stage_c1_full_muon
```

It starts from one **local muon crossing** defined at the wire plane `v=0`
by `(u0,w0,p,theta_u,theta_w)`.

The chain is:

```text
one muon
  -> Heed clusters
  -> all conduction electrons
  -> one microscopic Garfield avalanche per conduction electron
  -> microscopic electron endpoints + avalanche-ion birth points
  -> finite-time ion transport to Tobs
  -> Shockley-Ramo segment sums
  -> Stage-B hybrid opposite-cathode stereo sharing
  -> event-level strip amplitudes
  -> two projected CoGs
  -> reconstructed (u,w)
```

The low-level finite-time ion propagation and endpoint Shockley-Ramo evaluation
were factored into `../mwpc-avalanche/MWPCSignalUtils.h`.  The Stage-B3c.2
opposite-cathode stereo model was factored into
`../mwpc-avalanche/MWPCStereoReadout.h`.  C1 uses these shared helpers rather
than carrying a second independent strip-response implementation.

The forward track convention is `+v`, matching nominal `+global Z`.
Projected angles are exactly the C0 quantities:

```text
du/dv = tan(theta_u)
dw/dv = tan(theta_w)
```

### Build and smoke test

Load the usual NA6P/Garfield environment, switch to the canonical branch, then:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

source ~/na6p/env.sh
source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh

bash test/mwpc-stage-c/run_c1_smoke.sh
```

The smoke run intentionally limits the full Heed event to eight conduction
electrons.  It is a compile/integration check, not a physics result.

It uses a C0-motivated oblique test point:

```text
p = 5 GeV/c
theta_u = 10 deg
theta_w = 25 deg
u0 = 1 mm
w0 = 0
Tobs = 100 ns
```

Outputs are written under

```text
test_runs/mwpc_stage_c/c1_smoke/
```

with event-, seed- and strip-level CSVs.

### Full event after the smoke test

Run the same executable with `--max-seeds 0` (the default) to process every
Heed conduction electron.  Start with one event because the microscopic
avalanche + explicit ion transport is intentionally expensive.

The convenience runner takes `p, theta_u, theta_w, seed`:

```bash
bash test/mwpc-stage-c/run_c1_full.sh 5 10 25 130001
```

The defaults are the same values shown above.  The output directory name
records the requested phase-space point and random seed.

The event summary reports:

- Heed cluster and primary-electron counts;
- number of processed/zero-gain avalanches;
- total avalanche electron/ion counts;
- active wire multiplicity;
- finite-time cathode charges;
- stereo projected CoGs;
- strip-derived `u,w`;
- avalanche-ion-birth mean `u,w` as a charge-weighted microscopic truth
  diagnostic;
- wire-charge CoG in `u`;
- differences of the readout estimators both from avalanche truth and from the
  requested track intercept `(u0,w0)`;
- numerical diagnostics such as ion-transport failures.

For truncated smoke events, comparisons to `(u0,w0)` are not detector
residuals: the processed seeds may occupy only one part of the physical track.
The avalanche-ion-birth mean is the useful diagnostic in that case.

The strip summary contains the signed and absolute charge on every exported
strip of both stereo families.  The wire summary separately records collected
avalanche electrons/charge per anode wire, so wire multiplication and cathode
strip induction remain distinct observables.

The default Stage-C segmented weighting mode is `full-cathode-gap`: the ideal
planar strip weighting problem spans the two physical cathodes, so it is defined
throughout the gas volume.  `--sharing-domain legacy` retains the historical
Stage-B wire-plane-to-cathode half-gap construction for regression tests.

Both are still hybrid approximations: full cathode coupling is wire-aware,
while segmentation uses Garfield's ideal planar-strip solution rather than a
full numerical 3-D weighting field with discrete wires and PCB strips.

## Next C1 checks

Once the smoke/full normal-incidence event is numerically clean:

1. compare `theta_w=0, 25, 40, 50 deg` at fixed `theta_u`;
2. compare `theta_u=0, 10, 15 deg` at fixed `theta_w`;
3. repeat over several momenta representative of the C0 distributions;
4. only then build an ensemble/lookup-table scan.

The >40--50 deg points are tail validation, not the dense core grid.

## Scope boundary

Stage C does not yet claim final detector resolution. Absolute gas gain is not
prototype-calibrated, and electronics/noise/threshold/shaping are still absent.


## C2 — numerical strip weighting

The C1 full-event tests exposed the analytic hybrid strip-sharing model as the
remaining weak point.  Stage C2 therefore replaces only the weighting-field
piece while keeping the validated microscopic transport chain unchanged.

The first numerical prototype lives in:

```text
test/mwpc-stage-c/numerical-weighting/
```

It uses a parameterised local 3-D Gmsh/Elmer model with explicit cylindrical
anode wires and two physical cathodes.  For each weighting solve, the selected
strip is set to 1 V while all other conductors, including every wire, are set to
0 V.  The resulting map is imported through Garfield++ `ComponentElmer`.

Geometry is controlled by `numerical-weighting/weighting_config.json`; wire
pitch/diameter, both cathode gaps, strip pitch/width/stereo angle, patch size,
mesh refinement and strip indices are all runtime parameters.

Initial map smoke test:

```bash
bash test/mwpc-stage-c/run_c2_weighting_smoke.sh
```

Do not connect C2 maps to the full C1 event accumulator until the numerical map
passes its boundary and mesh/patch convergence checks.


## C2 3-D event visualization

A first ParaView-oriented event viewer is available through:

```bash
bash test/mwpc-stage-c/run_c2_event_3d.sh 120001 2
```

The runner reuses the default C2 FEM map when possible, dumps one canonical
Heed/Garfield event and exports a set of VTK layers:

- transparent cathode planes;
- explicit anode-wire axes;
- the Gmsh/Elmer tetrahedral chamber mesh;
- the muon track;
- all Heed ionisation clusters;
- all primary Heed electrons;
- representative microscopic avalanche electron drift paths for the selected
  detailed seeds;
- a sampled avalanche ion-birth cloud;
- a numerical u-v field slice at w=0 containing E, |E|, electrostatic
  potential and the central minus/plus strip weighting potentials.

The visualization Python environment is expected at
`~/na6p/venvs/mwpc-viz` and needs `meshio`.

Output is written below:

```text
test_runs/mwpc_stage_c/c2_event_3d/seed_<seed>/vtk/
```

The exporter also writes `open_scene.py`.  With ParaView installed:

```bash
paraview --script="test_runs/mwpc_stage_c/c2_event_3d/seed_120001/vtk/open_scene.py"
```

The FEM mesh starts hidden because it is visually dense; enable it from the
Pipeline Browser when inspecting the mesh refinement around wires.


### Animated C2 scene

The 3-D visualization runner now also writes a time-resolved ParaView scene.

```bash
bash test/mwpc-stage-c/run_c2_event_3d.sh 120001 2
```

Static scene:

```bash
paraview --script="$HOME/na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/c2_event_3d/seed_120001/vtk/open_scene.py"
```

Animated scene:

```bash
paraview --script="$HOME/na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/c2_event_3d/seed_120001/vtk/open_animation.py"
```

The time-series uses the Garfield/Heed timestamps already stored in the event
dump.  As time advances:

- Heed clusters appear;
- primary electrons appear;
- sampled microscopic electron trajectories grow point-by-point;
- avalanche ion-birth points accumulate.

Cathodes, wires, FEM mesh and the field slice are static reference layers.
The animation exporter uses cumulative PVD time series, so ParaView's native
Animation toolbar and time slider control the process.

# MWPC Stage C — chamber response integration

This directory is the **canonical working area for Stage C**.

The branch `feature/mwpc-stage-c` was created from
`feature/mwpc-avalanche-visualization`, so the mature Stage-A and Stage-B
work is already present on the same branch under `test/mwpc-avalanche/`.
Do not reimplement those pieces inside Stage C.

## One-branch workflow

**Remote code-edit workflow:** The assistant changes and commits source files
on GitHub directly; the local WSL user only executes `git pull` and bash
commands. Do not ask the user to paste patches or manually copy source files.

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


### Single-wire animation

For studying the local detector physics around one anode, the animation exporter
also creates a focused scene:

```text
open_animation_single_wire.py
```

Default local region:

```text
u = +/-2.0 mm
v = full configured gas gap
w = +/-3.0 mm
```

Only the following layers are loaded:

- the central anode wire;
- small local cathode patches;
- the local segment of the muon track;
- time-resolved Heed clusters;
- time-resolved microscopic electron/avalanche paths.

The chamber-wide wire bank, FEM mesh, field slice, primary-electron point cloud
and ion-birth cloud are deliberately omitted to keep playback responsive.

The ROI remains parameterised:

```bash
python test/mwpc-stage-c/visualization/export_event_animation.py \
  --event-dir <event-dir> \
  --vtk-dir <vtk-dir> \
  --config test/mwpc-stage-c/numerical-weighting/weighting_config_c2b.json \
  --single-wire-u-half-mm 2.0 \
  --single-wire-w-half-mm 3.0
```


## Long dense across-wire scan

For unattended detector-response characterization, use the resumable dense
normal-incidence scan:

```bash
bash test/mwpc-stage-c/run_c2_u_cell_scan_long.sh 30 140001 0.2
```

Default range is `u0=-1.9..+1.9 mm`.  With a 0.2-mm step this gives 20
positions; 30 events per point therefore gives 600 full microscopic events.

The same seed block is reused at every position for paired comparisons.  The
FEM map is built only once and reused.  Completed positions are checkpointed;
rerunning the same command skips them, while an interrupted point is rerun from
scratch.

Useful outputs:

```text
test_runs/mwpc_stage_c/c2_u_cell_scan_dense/
  c2_u_dense_event_summary.csv
  c2_u_dense_point_summary.csv
  scan_progress.log
  events/u_*/...
```

The compact point summary reports raw stereo means/spreads, gain-weighted truth,
wire assignment fractions, multi-wire fraction and strip-charge coverage.

The range can be overridden without editing the script:

```bash
C2_U_MIN_MM=-1.9 C2_U_MAX_MM=1.9 \
bash test/mwpc-stage-c/run_c2_u_cell_scan_long.sh 50 140001 0.1
```


## Full strip-pattern estimator

The dense C2 scan already writes one strip-summary CSV per u0 point, so the
existing microscopic sample can be reused to test reconstruction improvements
without rerunning Garfield or Elmer.

The comparison script uses the normalized two-cathode strip amplitudes

```text
a(side,k) = |Q(side,k)| / sum_j |Q(side,j)|
```

and compares a full strip-pattern template estimator against the current
nonlinear calibration of the stereo CoG.

Run it with:

```bash
~/na6p/venvs/mwpc-viz/bin/python \
  test/mwpc-stage-c/analyze_c2_strip_patterns.py \
  --scan-dir test_runs/mwpc_stage_c/c2_u_cell_scan_dense
```

The comparison is cross-validated by random seed.  Because the dense scan
reuses the same seed block at every u0, all events from a held-out seed are kept
out of the template used to reconstruct that seed.

Outputs:

```text
c2_strip_estimator_predictions.csv
c2_strip_estimator_metrics.csv
c2_strip_estimator_by_u.csv
c2_strip_templates.csv
```

The strip-template estimator deliberately uses no truth coordinates, avalanche
birth positions, gain-weighted seed coordinates, or wire assignment.  It tests
whether information discarded by the two-CoG compression is recoverable from
the cathode strip pattern itself.


## Regularized strip reconstruction v2

The first full-strip nearest-template test is intentionally simple and can be
noise-dominated because the strip fractions change only weakly with sub-wire
position.  A second reconstruction study therefore compares four seed-grouped
cross-validated estimators:

```text
nonlinear_cog_calibration
pooled_strip_template
cog_plus_full_strip_ridge
cog_plus_strip_moment_ridge
```

The pooled template uses one regularized covariance matrix for strip-pattern
residuals.  The ridge estimators use nested seed-grouped cross-validation to
select the regularization strength.

Run:

```bash
~/na6p/venvs/mwpc-viz/bin/python \
  test/mwpc-stage-c/analyze_c2_strip_reco_v2.py \
  --scan-dir test_runs/mwpc_stage_c/c2_u_cell_scan_dense
```

This is reconstruction-only and does not rerun Garfield or Elmer.


## Wire-only Stage-C response

For the current detector-development step, the cathode-strip reconstruction can
be factored out and the coordinate associated with the anode-wire plane can be
studied directly.

The dedicated executable

```text
mwpc_stage_c_wire_response
```

keeps

```text
muon -> Heed clusters -> primary electrons -> microscopic avalanche
     -> collected avalanche electrons on the anode wires
```

but intentionally omits Elmer weighting maps, cathode-strip induction, ion
transport, electronics, thresholds and shaping.

For every event it records the active wires, the dominant wire, its collected
charge fraction, and the charge-weighted wire coordinate.  The binary estimator
is simply the coordinate of the dominant wire.

A focused scan maps the transition between neighbouring 4-mm wire cells:

```bash
bash test/mwpc-stage-c/run_c_wire_boundary_scan.sh 100 160001
```

The default points are

```text
u0 = +/- (1.5, 1.6, ..., 2.5) mm
```

so the geometric boundaries at `u = +/-2 mm` are crossed explicitly.  The run
is resumable point-by-point and is much lighter than the full C2 cathode scan
because no FEM weighting or ion signal is evaluated.

Outputs:

```text
test_runs/mwpc_stage_c/c_wire_boundary_scan/
  c_wire_event_summary.csv
  c_wire_wire_summary.csv
  c_wire_point_summary.csv
  c_wire_transition_summary.csv
  scan_progress.log
```

The point summary contains dominant-wire probabilities, multi-wire fractions,
binary-wire and charge-weighted-wire residuals.  The transition summary
estimates the 16%, 50%, and 84% neighbour-wire transition positions and the
16--84% transition width separately at the positive and negative boundaries.


## Electron transport box map

To separate the gas-transport problem from the later GEANT ionisation pattern
and cathode readout, Stage C now has a reusable single-electron response map.

The executable

```text
mwpc_stage_c_transport_box_map
```

divides one 4-mm wire cell and the full 5-mm gas gap into `(u,v)` boxes.
From every box it launches many electrons with
`AvalancheMicroscopic::DriftElectron` and stores

```text
starting box
  -> destination wire
  -> drift time
  -> hit position around the wire
  -> last microscopic direction
  -> incidence angle with respect to the inward wire normal
```

The default grid is

```text
u = [-2,+2] mm
v = [-2.5,+2.5] mm
du = dv = 0.2 mm
20 x 25 = 500 boxes
```

and the default run uses 50 electrons per box.

Run:

```bash
bash test/mwpc-stage-c/run_c_transport_box_map.sh 50 170001
```

The shell scan is resumable by complete u-column. It uses the current
5-mm, 4-mm-pitch, 30-um-wire, 1800-V, Ar/CO2 70:30 baseline and starts at
B=0.

Outputs:

```text
test_runs/mwpc_stage_c/c_transport_box_map/
  c_transport_box_electrons.csv
  c_transport_box_summary.csv
  scan_progress.log
  columns/
```

The raw electron file contains the endpoint wire, drift time,
`phi_hit`, `cos(phi_hit)`, the last microscopic direction, and the
incidence angle.  The box summary contains

```text
P(wire=-1), P(wire=0), P(wire=+1)
mean/std/median drift time
A_LR = (N(+u hit side)-N(-u hit side))/(N(+u)+N(-u))
<cos(phi_hit)>, <sin(phi_hit)>
mean/std incidence angle
```

This map is intentionally independent of a particular muon track.  A later
digitisation step will weight the boxes with the GEANT cluster / energy-deposit
distribution.  The exact GEANT adapter should be connected to the actual
NA6PRoot cluster or step output rather than inventing a surrogate format.


## PI track + Poisson-cluster workflow

The large historical MNP33 bore scans can be reused for realistic muon
crossing kinematics without rerunning the absorber/magnet transport. Their
`bore_muon_crossings.csv` files contain station position and momentum.

Stage-C coordinate mapping:

```text
global X -> local w (along wires)
global Y -> local u (across wires)
global Z -> local v (chamber normal)

theta_w = atan(px/pz) = theta_x
theta_u = atan(py/pz) = theta_y
```

Prepare compact per-channel/per-station track samples and summaries:

```bash
python3 test/mwpc-stage-c/prepare_c_pi_tracks.py \
  --base <mnp33-bore-run-dir> \
  --output-dir test_runs/mwpc_stage_c/c_pi_tracks \
  --sample-per-channel-station 2000
```

The old global `y_cm` is **not** automatically folded into a 4-mm wire
cell, because the exact per-chamber wire-origin phase is not frozen in those
historical crossing files.

A first-order ionisation model can then generate clusters along each realistic
track:

```bash
python3 test/mwpc-stage-c/generate_c_pi_poisson_clusters.py \
  --tracks test_runs/mwpc_stage_c/c_pi_tracks/c_pi_track_sample.csv \
  --output test_runs/mwpc_stage_c/c_pi_tracks/c_pi_poisson_clusters.csv
```

Default cluster density is `3.4752/mm`, estimated from the current Stage-C
Heed baseline (`<Nclusters>=17.376` in a 5-mm normal-incidence crossing).
The default optional cluster-size approximation is a shifted Poisson with mean
`2.3704` primary electrons/cluster, matching
`<Nprimary>/<Nclusters> = 41.188/17.376` in the same sample.

By default the wire phase at `v=0` is sampled uniformly in one pitch. This is
appropriate for a broad acceptance sample until the exact old-run chamber/wire
phase is tied to the current geometry. The generated cluster coordinates are
folded back into one wire cell and store the integer wire offset separately.

The next stage folds these synthetic clusters with the Garfield transport map:
`P(wire|u,v)`, drift time, and conditional wire-arrival angle.


## PI checkerboard occupancy + incidence angles (18 + 6 figures)

The historic high-statistics MNP33 bore production supplies transported muon
states at the station reference planes. The checkerboard projector maps those
states onto the current 298 staggered sensitive gas boxes without rerunning Geant4.

**Step 1 — ensure full checkerboard projection CSVs exist:**

```bash
bash test/mwpc-stage-c/run_c_checkerboard_projection.sh full
```

This generates six `crossings_<source>_ms2_<z>.csv` files plus `modules.csv`
and `stations.csv` under `test_runs/mwpc_stage_c/checkerboard_full/`.

**Step 2 — make PI-ready figures without rerunning the projection:**

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh
```

The runner first checks geometry and multi-chamber counting with
`test_checkerboard_projection.py` and `test_pi_checkerboard_figures.py`.
Default results: `test_runs/mwpc_stage_c/pi_checkerboard_figures/`.

- `occupancy/`: **18** separate PNG images (3 channels x MS0..MS5).
  Each is a two-panel figure: global XY hit density (counts/cm²) on the left,
  and raw counts per chamber annotated with A/B/C/D checkerboard layer on the
  right. **Every gas-volume crossing gets a spatial hit**, so one muon
  crossing two staggered chambers increments both modules and produces two
  spatial hits.
- `angles/`: **6** images (3 channels x theta_u/theta_w), each with six
  station panels, PNG only. **One angle entry per contributing
  (job,event,track,station)**, irrespective of chamber crossing multiplicity.
  Common angular range/normalization across stations within each figure;
  panels show mean, standard deviation and P95 of |theta|.
- `pi_station_crossing_summary.csv`: unique muons, chamber crossings, and
  counts for 1,2,3,4,5+ chamber intersections.
- `pi_muon_angles_summary.csv`: angular statistics per source/station.
- `READ_ME_FIRST.txt`: settings, approximate wire count, and caveats.

The full figure producer checks **each individual chamber count** against
`modules.csv` and unique-muon/station totals against `stations.csv`.
It fails on disagreement.

To compare the historic reference position of MS2=550 cm while preserving
the current code baseline of MS2=530 cm:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh \
  --ms2-z 550 \
  --out test_runs/mwpc_stage_c/pi_checkerboard_figures_ms2_550
```

These are geometric straight-line projections from the historical transported
states, **not** actual Geant4 gas hits in the current 298-chamber material
geometry, calibrated wire signals, or physical Pb-Pb event occupancy rates.
Do not interpret preselected muon reference states as unbiased efficiency.

The 4-mm reference anode pitch along global Y gives about **128 wire
intervals per standard** chamber (51.20 cm active Y), and **77 per MS0**
narrow chamber (30.72 cm active Y, rounded). A provisional system estimate
is 32 x 77 + 266 x 128 = **36,512** anode wires; physical edge/wire-origin
counting must be confirmed with the PI.


### Occupancy color contrast

The 18 PI occupancy figures use sqrt color normalization by default, with
an upper color limit at the 99.5th percentile of nonempty bins. Bins
with zero hits are white; bins above the limit use the maximum color
with a colorbar extension. Raw hit counts are unchanged and colorbars
show physical values rather than transformed intensities.

No new Geant4 transport or checkerboard projection is needed to update
the plots. Re-use the full checkerboard CSVs and run:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh
```

To show faint occupancy with a logarithmic scale in a separate folder:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh \
  --occupancy-scale log \
  --out test_runs/mwpc_stage_c/pi_checkerboard_figures_log
```

To retain the full color range with no high-end saturation:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh \
  --occupancy-scale sqrt --occupancy-clip-quantile 1 \
  --out test_runs/mwpc_stage_c/pi_checkerboard_figures_sqrt_noclip
```

The old linear color scale is also available with
`--occupancy-scale linear`. The mode and quantile are labelled on each map.


### Shared XY frames and PNG-only PI figures

By default both halves of each occupancy figure use the **same global
physical X/Y range for all six stations**, based on the largest active
gas footprint across the spectrometer: X = [-287.52, +287.52] cm and
Y = [-242.50, +242.50] cm. This shows the smaller MS0/MS1 stations with
empty outer space instead of artificially filling the whole frame.

The left-hand spatial-density map uses identical XY bin edges and
physical bin areas across stations. The right-hand module-count map
is now also drawn in global X/Y centimetres rather than incomparable
row/column indices. Module tiles are centred on the true installed
gas XY positions. Each tile visualizes a module's crossing total;
the true 3-cm gas overlap is preserved in the counts but cannot be
depicted as disjoint rectangular colored tiles.

The color norm remains independently determined per figure (sqrt
with 99.5th-percentile high-end saturation by default); hence equal
colors do NOT automatically mean equal numeric occupancy between
different station figures. Read colorbar values when comparing absolute
chamber hit rates.

Output is **PNG only**. When rerun in an existing result directory,
the script removes older PDFs matching only its own previously generated
occupancy and angle filename patterns, leaving unrelated files alone.

Default shared-frame run:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh
```

To compare to the old station-fitted layout in a separate directory:

```bash
bash test/mwpc-stage-c/run_pi_checkerboard_figures.sh \
  --xy-range station \
  --out test_runs/mwpc_stage_c/pi_checkerboard_figures_station_range
```


## Reconstruct original dimuon momentum/rapidity distributions from saved ROOT

The older `test/mwpc-dimuon-validation/analyzeMWPCDimuon.C` implemented
`kinematics_{Jpsi,Phi,Omega}.png` with generated parent momentum, parent
`pT`, parent rapidity, daughter-muon momentum and daughter-muon `pT`.
The previous 20k-event image files need not be present: the historical
large MNP33 bore production has six `MCKine.root` files (two chunks x
200k generated parent events/channel in aggregate). Rebuild from those
existing ROOT files with the independent macro
`recover_bore_dimuon_kinematics.C`.

```bash
bash test/mwpc-stage-c/run_recover_bore_dimuon_kinematics.sh smoke
bash test/mwpc-stage-c/run_recover_bore_dimuon_kinematics.sh full
```

These commands do **not** rerun Geant4 or Garfield. The smoke mode
processes the first 200 ROOT events **per chunk** for compilation and
schema checks, and should not be used for physical interpretation.
The full mode processes all stored events. The runner sources the
pre-existing `~/na6p/env.sh` only if `root` is not on PATH.

Results: `test_runs/mwpc_stage_c/recovered_dimuon_kinematics/`
(or `_smoke/` for the test):

- `kinematics_Jpsi.png`, `kinematics_Phi.png`, `kinematics_Omega.png`:
  generated **parent p/pT/y**, daughter muon **p/pT/y**, and
  **daughter theta_u/theta_w at production** (8 panels plus statistics).
- `kinematics_truth_summary.csv`: number of events, generated parents,
  direct muon pairs, missing pairs, and maximum four-momentum closure.

**Scope:** These are generator-level parent/daughter distributions
read back from the historical transport files. No downstream bore or
chamber selection is applied. They are *not* the same as Stage-C
accepted muon angular distributions, which select muons that
geometrically cross sensitive gas. Generated muon angles can
be compared with MS0 crossing angles to diagnose whether
J/psi's two-peaked angle spectrum is due to decay kinematics
or detector/acceptance selection.

This restores the analysis using an **independent saved sample**, not
the exact lost PNGs or necessarily the exact old 20k production. ROOT
is required; a short smoke run should be inspected before full processing.


## PI extended 2.5-hour study: collected wire multiplicity PER MODULE

This run processes **every individual 3D projected chamber crossing** from
`checkerboard_full/crossings_<channel>_ms2_530.csv`. If a muon intersects
two checkerboard modules, it contributes two separate gas-response crossings,
one to each `det_id`. Data are not averaged across stations for the PI
deliverable. Baseline station Z uses MS2=530 cm.

The fast model samples a Poisson number of ionisation clusters along the
actual gas path using the existing 2+4 mm Heed density; cluster electron
multiplicities use the saved empirical Heed distribution. Each primary
electron receives a sampled single-electron Garfield wire, collection
flag and drift time from the 20x30 transport-box map. By default the
script automatically chooses the enhanced hybrid high-stat boundary map
if it already exists; otherwise uses the available complete 50-electron
map (it will NEVER start a new Garfield scan on its own).

Before running on high-stat inputs, a fully synthetic regression test
checks per-chamber IDs, the 600-box map, random sampling, checkpoint
restart determinism, and 298-module probability normalization.

```bash
bash test/mwpc-stage-c/run_pi_per_chamber_response.sh smoke
bash test/mwpc-stage-c/run_pi_per_chamber_response.sh full
```

`smoke` processes just the first 100 source crossings in a separate output
folder. `full` runs for at most ~2.5 hours per invocation with checkpoints
every 10k processed chamber crossings, then can be resumed by repeating
the same command. Runtime depends on the WSL hardware and saved
transport map; it may finish sooner or need multiple sessions.

Results: `test_runs/mwpc_stage_c/pi_per_chamber_response/`:

- `wire_crossings_<source>.csv`: sampled per-chamber crossing response,
  including geometric wire cells, Heed clusters, electrons collected per
  crossing, wires with >=1/2/5/10 collected primary electrons and wires
  with electron arrivals by 100/200 ns.
- `pi_wire_multiplicity_per_chamber.csv`: **298 modules x 3 channels**
  with per-module probabilities `P_nwire_0/1/2/3/4/5plus`, >=2,
  >=3, rough electron-count threshold sensitivities, timing and means.
  Includes `processed_chamber_crossings` and `expected_all_crossings`,
  so incompletely processed stations cannot be mistaken for final results.
- `checkpoint_<source>.json`: durable source progress with RNG state;
  rerunning the same command resumes deterministically.
- `complete_<source>.ok`: written only after this source's full
  298-module crossing counts agree exactly with `modules.csv`.
- `settings.json`: immutable input paths, random seed, geometry conventions
  and explicit limitations.

**Physics caveats:** These are STILL straight-line projected chamber
crossings, not updated Geant4 hit/Edep production. The reference Garfield
map is B=0 with asymmetric [-2,+4] mm gas gap. An explicit +1 mm
coordinate offset aligns the six-mm symmetric projected gas box with
that asymmetric simulation; this does not create the same actual
electric field in a symmetric chamber. Wire phase is randomized per
crossing because exact wire origins are not available, and boundary
losses are not simulated. This is a *primary-electron collection*
distribution, NOT measured wire efficiency: no gain fluctuation,
electronic threshold, ion drift pulse shaping or noise is included.
The >=2/5/10 electron thresholds are only a sensitivity diagnostic.


## PI fixed-grid wire occupancy atlas (18 + 894 + 18 PNGs)

PI's chamber-by-chamber physical consistency request is addressed by a
**fixed, finite anode grid**, so the nominal physical wire number is stable
across crossings in the same chamber. Every checkerboard gas crossing
contributes to its own `det_id`, even when the same muon also intersects
another staggered chamber.

Anode pitch is 4 mm, with the nominal wire positions symmetrically
centred inside each module's sensitive Y span. The working assumptions
are **77 wires on each of the 32 narrow MS0 chambers**, and **128 on
each of the remaining 266 chambers**, or 36,512 in total. Exact
mechanical wire origin and wire count still require PI confirmation.
Three physical-grid phases are simulated: -1, 0, +1 mm, using matched
underlying ionisation and microscopic response random quantiles.

**Step 1: smoke test (separate directory; skips 894 per-chamber details)**

```bash
bash test/mwpc-stage-c/run_pi_fixed_wire_atlas.sh smoke
```

**Step 2: full production on the already saved checkerboard crossings,
Heed library and Garfield transport map**

```bash
bash test/mwpc-stage-c/run_pi_fixed_wire_atlas.sh full
```

The runner first tests synthetic physical-wire indexing, finite edge
handling and independent per-source chamber counters. The full run
recomputes the fast *sampling/folding* of Heed+Garfield response on all
~2.49M projected gas-volume crossings. It does **not** rerun Geant4,
Heed microscopic tracks, or Garfield transport. Runtime depends on
CPU and especially on writing 894 detailed PNGs.

Primary output directory:
`test_runs/mwpc_stage_c/pi_fixed_wire_maps/`.

- `station_overview/`: 18 PNG station maps, three conditions on the
  same global physical XY axes with the actual checkerboard gas outlines.
  For readability, each colored short horizontal stripe is the
  **mean firing probability of up to eight neighboring anodes**,
  not a single anode.
- `per_chamber/<source>/MS<n>/wire_DET<id>.png`: **894 PNGs**,
  three true per-anode profiles for each source/module, including modules
  with no accepted muons (marked as such).
- `contact_sheets/`: 18 contact sheets, showing the three criteria
  and the fine per-anode profiles simultaneously across all modules
  in a station.
- `physics_summary/`: 3 channel-wise plots each of station
  multiplicity, P(>=2 anodes) vs absolute incident theta_u,
  and -1/0/+1 mm wire-grid phase sensitivity (9 PNGs).
- `pi_fixed_wire_counts.csv`: 109,536 nominal per-anode rows,
  with actual `station,det_id,row,col,layer,wire_index,y_global_cm`
  and three counts plus fractions per chamber crossing.
- `pi_fixed_wire_station_multiplicity.csv`: 0/1/2/3/4/5+
  active-anode distribution for geometric, >=1e, >=5e and
  first-100ns interpretations at each station, including all phases.
- `pi_fixed_wire_angle_response.csv`: P(>=2 active anodes) vs
  |theta_u| by source, station, condition and wire-grid phase.
- `pi_fixed_wire_settings.json`: all inputs, seed, wire-grid definition
  and physics approximations.

The three signal definitions **are not electronics discriminators**:
`>=1e` = at least one collected primary electron on the wire;
`>=5e` = at least five collected primary electrons;
`<=100ns` = at least one collected primary electron reaching it
within 100 ns in the saved Garfield response. We do not include
avalanche gain fluctuations, ion motion, capacitive cathode signal,
front-end threshold, mechanical edge-field distortion, or noise.

The prior Stage-C statistical module multiplicity sample randomly
shifted wire phase for **each crossing**, which was adequate for
global multiplicity but cannot label individual physical wires.
The atlas is a new fold of the same reference gas projections with
a **stationary module-local grid**. The finite grid excludes Garfield
responses reaching nonexistent edge wires; these lost electron counts
are separately reported in the wire CSV.

Color normalization is shared across the **six stations within each
channel** and all three conditions, with the 99.5th percentile
of nonzero wire probabilities as the color ceiling; colorbars have
an extension when saturated. Do not compare raw colors between
different dimuon channels without reading their numeric colorbars.

## Recover the Stage-C Garfield drift-time (u,v) grid figure

Previously computed timing information is in the baseline
`test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/`:

- `c_transport_box_electrons.csv`: individual launch boxes and drift times
- `c_transport_box_summary.csv`: 600 box means, standard deviations,
  medians and 16%/84% drift time quantiles.

The advanced **hybrid** transport map contains all 600 grid boxes,
with 2000 electron trials in each of 120 boundary boxes and 50 trials
in each of the other 480 boxes:
`c_transport_box_map_hybrid_boundary_2p4/c_transport_box_electrons.csv`.

To reconstruct four 2D heat maps (mean/median collected-electron
drift time; collection fraction; fraction collected within 100 ns)
from whichever existing complete map is available:

```bash
bash test/mwpc-stage-c/run_c_drift_time_grid.sh
```

PNG and 600-row CSV:
`test_runs/mwpc_stage_c/c_drift_grid_figures/`.

The printed image title explicitly labels the conditional vs
unconditional probabilities: mean and median use collected electrons
only; 100-ns fraction uses *all launched electrons* as denominator.


## Physics-first PI wire atlas V2 (after visual review)

The original 18 station wire images and 894 chamber details are kept.
The improved **V2** avoids interpreting a row of horizontal wire stripes
as spatial muon XY density; an anode wire intrinsically integrates along X.

V2 reads the complete already-computed `pi_fixed_wire_counts.csv`
**without resampling Heed/Garfield**, plus the original projected
`checkerboard_full/crossings_<source>_ms2_530.csv`. A deterministic
geometry-control pass counts the true finite anode cells intersected
by each projected gas track and constructs real XY crossing histograms
from gas-segment midpoints, preserving gas-volume det_id.

Preview (all real data, just 18 station + 18 contact figures):

```bash
bash test/mwpc-stage-c/run_pi_wire_validation_v2.sh preview
```

After reviewing, make 894 corrected per-chamber profile sheets as well:

```bash
bash test/mwpc-stage-c/run_pi_wire_validation_v2.sh full
```

Outputs remain **separate from V1**, under:

`test_runs/mwpc_stage_c/pi_fixed_wire_maps/physical_validation/`.

- `pi_geom_wire_counts.csv`: 109,536 rows of physical anode cell
  geometric-hit counts (source,station,det_id,wire_index), with a
  chamber crossing denominator.
- `pi_crossing_density_xy.npz`: 18 real XY crossing-density histograms
  on a physical 5-cm bin grid, with their own bin edges.
- `atlas_v2/station_validation/`: 18 **six-panel** figures, for each
  source/station: real XY track crossing density; geometric cells
  crossed; Garfield anodes collecting >=1 or >=5 primary electrons;
  electrons collected within 100ns; and response-minus-geometry
  on a common physical wire coordinate.
- `atlas_v2/contact_sheets/`: 18 station-wide grids of sparklines
  overlaying gray geometry, orange collected >=1e, and blue >=5e;
  chamber IDs and crossing statistics are labeled; red rims explicitly
  warn when a chamber has fewer than 100 projected crossings.
- `atlas_v2/chamber_profiles/`: 894 detailed PNGs each comparing
  the same anode indices across geometry, >=1e, >=5e, and <=100ns;
  low-stat chambers are flagged instead of silently overstated.

**Improved semantics.** In station-level wire maps, the color is
anode activation count per **10,000 chamber gas-volume crossings in
the same channel/station**. This is an *overall occupancy/rate*
map, not a conditional module efficiency. On the individual chamber
sheet, the vertical axis is activation probability per **crossing
of that specific chamber**. The two denominators must not be
conflated. Each wire runs along X so no apparent X localization
should be inferred from long colored lines. The true XY information
is shown in a dedicated panel. The physical view zooms to each
station, with the SAME axis range across the three channels for a
given station. No physically meaningful single-anode X resolution
is claimed.

**Validation scope.** The geometric and anode-response comparisons
still use the saved GEANT **reference-plane states** extrapolated
into the staggered 298 gas boxes, not native Geant4 MWPC gas Edep.
The 77/128 anodes, mechanical Y origin, and use of the reference
B=0 and asymmetric Garfield 2+4mm transport map remain nominal.


## PI project status and literally numbered physical wires

A task-by-task audit of the PI's requested simulation outputs, exact
interpretation of wire-rate denominators, physical caveats and next modern
Geant4 runs is maintained in
**test/mwpc-stage-c/PI_STATUS_AND_TALK.md**.

A reproducible schematic of an entire station with each chamber's real
Stage-C det_id and stagger layer is now available. Example MS0 DET17
includes a second image with **all 77 nominal anodes literally numbered
0..76**, each with three measured-in-model collection probabilities.

    bash test/mwpc-stage-c/run_pi_station_numbered_anodes.sh

Optional other physical detector and channel:

    bash test/mwpc-stage-c/run_pi_station_numbered_anodes.sh 52 phi40-bore

Output: test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_explanation/.

These are nominal, centred 4-mm-pitch anode grids: MS0 wire indices 0..76;
MS1..5 wire indices 0..127. Increasing indices run along global Y (across
wires); each wire extends along global X. Chamber XY outlines are staggered
in Z according to A/B/C/D. **Exact first wire position and readout thresholds
are not mechanically/electronically confirmed.**


## Interactive standalone station dashboard (offline browser application)

For a PI walkthrough, the existing **full data**, not synthetic demonstration
curves, can be explored in a standalone HTML file that needs **no website,
server, npm, ROOT or external CDNs**.

Run once after completing the fixed-wire atlas and geometric controls:

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_station_dashboard.sh

The runner executes a synthetic end-to-end regression (all 298 modules,
3 source types, 36,512 nominal wires per source), checks syntax and builds
a single file, opening it with Windows Explorer if WSL interop is available:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_station_interactive.html

Features:
- source switch J/psi, phi, omega; station switch MS0..MS5;
- **Muon spatial mode**: 5x5-cm XY crossing footprint OR chamber counts
  overlaid with actual nominal A/B/C/D chamber outlines;
- **Wire mode**: physical horizontal anode lines coloured by activity;
  selectable >=1 collected primary electron / >=5 / first 100ns / geometric;
  group by 1, 4 or 8 adjacent real anodes for station overview;
- selectable normalization **per actual chamber gas crossing**,
  **per transported incoming station-reference muon**, or
  **per 10,000 station chamber crossings**. The choice is explicit;
- click a real detector ID; inspect its 77/128 distinct wire indices and
  4 aligned response/geometric profiles; hover wire groups for numerics;
- zoom/pan the station view, export the selected diagram as SVG and the
  selected detector's exact per-wire counts as CSV;
- warnings for zero-exposure and low-statistic modules.

Input files (also validated against modules.csv/stations.csv and each other):
- checkerboard_full/modules.csv and stations.csv
- pi_fixed_wire_maps/pi_fixed_wire_counts.csv
- pi_fixed_wire_maps/physical_validation/pi_geom_wire_counts.csv
- pi_fixed_wire_maps/physical_validation/pi_crossing_density_xy.npz

This is a **review UI**, not a replacement for full Geant4 hit/Edep production.
The spatial hit map uses midpoint of projected gas-box intersections, not
final digitized hits. The wire coloured stripes integrate over all X along
the wire; do not treat their apparent XY extent as reconstructed muon X.
The current nominal fixed wire grid is not a verified mechanical blueprint.


## Geometry-change-resistant validation strategy and native Geant4 benchmark

**Two distinct datasets, with different validity domains.**

1. Historical high-stat **transport/reference-state library** (200k generated
   dimuon events per channel) is already in the separate MNP33 archive; see
   LOCAL_DATA_INVENTORY.md. Its LocalTrajectoryAudit.root and MCKine.root
   can be reprojected cheaply to alternate candidate chamber placements.
   These projected crossings do NOT reproduce new chamber material scattering
   or gas ionization/Edep, and do not replace true Geant4 chamber hits.
   The historical transport had its OWN upstream apparatus/material settings.
2. New **native Geant4 current-MWPC benchmark** runs the actual 298-chamber
   geometry with source dimuons and preserves native ROOT truth and gas hits.
   It freezes/chunks output for reproducibility, but it is NOT transport
   geometry-independent when the chamber material/layout changes. It is a
   baseline for validating faster alternative layouts and for calculating
   Edep/angle/individual gas-step distributions. The detailed wire pitch,
   zero position, threshold and pulse model may be varied offline using these
   saved actual gas-hit segments and deposited-charge diagnostics.

A safety-first, resumable run is now available; **run a 3-events/channel
smoke test first** (also compiles the ROOT generator/hook before parallel runs):

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh smoke

Then launch longer 40-AGeV native Geant4 production; arguments are
total generated parent dimuon events/channel, events per atomic chunk,
workers in parallel and approximate wall-clock budget in hours:

    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh full 30000 3000 4 2.5

The default is 30k J/psi + 30k phi + 30k omega generated events, stored
in 3000-event chunks. The runner stops *between* chunks when the time
budget expires, may overshoot that time by one chunk, and resumes via the
EXACT SAME command. Never rename the suite directory or change parameters
mid-campaign. For different physics/geometries use a fresh separate suite.

Optional detached launch after smoke success:

    mkdir -p test_runs/mwpc_stage_c
    nohup bash test/mwpc-stage-c/run_pi_geant4_campaign.sh full 30000 3000 4 2.5 \
      > test_runs/mwpc_stage_c/pi_geant4_campaign.log 2>&1 &
    tail -f test_runs/mwpc_stage_c/pi_geant4_campaign.log

Campaign lives in:

    test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1/

It contains frozen manifest campaign.json, source/channel/chunk directories,
full command and simulator logs, per-chunk ROOT checks, SHA256 output
digests and completed.json atomic chunk checkpoints.

**Files preserved in each completed chunk:**
- MCKine.root: original generator parents and retained daughters, including
  the kept signal muon ancestors enforced by mwpcDimuonHooks.C;
- HitsMuonSpec.root: native TTree hitsMuonSpec and gasStepsMuonSpec, the
  latter enabled by NA6P_MWPC_STORE_GAS_STEPS=1;
- geometry.root: exact TGeo geometry for interpreting gas local coordinates;
- na6pLayout.ini: the exact layout/configuration saved by that Geant4 run.

Modern NA6PBaseHit ALREADY stores three-momentum at entry AND at exit.
The updated export_c_geant_muon_crossings.C can export local theta_u/w
at entry/exit, their difference, true length and hit Edep, as well as
positions and status bits. The real ROOT records stay the authoritative
unreduced source.

**Important join caveat:** The optional gasStepsMuonSpec diagnostic is filled
during Geant4 stepping, before final NA6PMC track remapping. Its trackID is
the original stepping/stack ID. MuonSpec hit vectors and MCKine retained IDs
may use remapped IDs. Do NOT join these trees using those raw IDs alone.
The step tree can be analysed separately spatially and by detector; for
track-wise linkage a remapping-aware export is required.

**Validation:** The campaign checks every chunk has all four mandatory files,
equal truth/hit TTree event counts, tracks and MuonSpec branches, and the
gasStepsMuonSpec tree with spatial and Edep branches. A missing diagnostic
tree likely means the locally installed libsimLib was not rebuilt from this
Stage-C branch; STOP and fix the build before any long run. Failed or
interrupted chunks are never reported as complete. They are archived under
_interrupted_* before retry, not deleted.

The saved step-level data allow new wire-grid or readout variants to be
compared without replaying Geant4 WHEN THE GAS AND PASSIVE MATERIAL
GEOMETRY REMAIN THE SAME. Moving/replacing MWPC material changes multiple
scattering and energy loss: the corresponding true detector-hit/Edep
distributions must be rerun or approximated and explicitly labeled.


### NEW: reference-plane library written by the modern Geant4 campaign

To make the **new** production useful for changing the *placement* and
acceptance of the MWPC modules, the branch now has an optional passive
muon-trajectory auditor in NA6PMC and NA6PMuonReferenceAudit.h. It is
enabled by NA6P_MWPC_STORE_REFERENCE_PLANES=1 and writes the fifth mandatory
file for each chunk:

    MuonReferencePlanes.root

The Planes tree contains global position and momentum of each transported
muon at virtual Z slices **280..900 cm every 10 cm** plus the magnet
reference planes 365/430/495 cm. Each row stores event/track/PDG,
plane index, direction and the endpoints of the actual Geant4 step used
for interpolation. The states are **linearly interpolated between
actual Geant4 steps**. The observer introduces no material, no step
limiter and no particle cuts. The track ID is the *original* stack ID,
prior to MCKine track ID remapping. Reusing these data requires correct
event/track mapping, and is valid for prospective layouts only as an
approximation when the chamber material budget changes.

Unlike the separate historical archive, these virtual-plane tracks and
native gas hits come from the **SAME NEW Geant4 campaign**; event and
chunk identifiers can be retained together for cross checks.

Since this is new native compiled C++ code, **rebuild/install the current
feature branch before running the smoke test**:

    bash test/mwpc-stage-c/rebuild_pi_geant4_stagec.sh
    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh smoke

The rebuild helper requires an existing cached CMake build tree for
THIS exact source path and refuses to guess an unrelated build. Supported
paths include the checkout's build/, the documented NA6PBUILDDIR,
and customary ~/na6p/build variants. If none matches, it stops safely
without starting a multi-hour job; review the error and locate the
original CMake build instead of creating duplicate installations.

The smoke test now requires and validates BOTH auxiliary outputs:
- gasStepsMuonSpec (in HitsMuonSpec.root) from NA6P_MWPC_STORE_GAS_STEPS=1
- Planes (in MuonReferencePlanes.root) from NA6P_MWPC_STORE_REFERENCE_PLANES=1.

The preserved five mandatory chunk outputs are MCKine.root,
HitsMuonSpec.root, geometry.root, MuonReferencePlanes.root, and
na6pLayout.ini. The simulator binary AND loaded libsimLib.so hashes are
frozen in campaign.json to ensure resumed chunks do not mix versions.

To analyze actual module entry/exit momentum and local angles after a
successful run, the existing
test/mwpc-stage-c/export_c_geant_muon_crossings.C
now additionally exports px/py/pz-out, theta_u/w-in/out and delta-theta;
the native saved ROOT hit object remains authoritative.


**Worker provenance and disk space.** In addition to the per-chunk
merged ROOT files, the campaign now preserves original per-worker ROOT
subruns under each chunk's raw_workers/ directory. This is deliberate:
the diagnostic gasStepsMuonSpec event/track labels are written before
MCKine remapping and may be worker-local. Keeping raw workers guarantees
we can reconstruct the worker-specific event namespace and align native
steps with parent truth later. Expect more disk use than the merged-only
ROOT output (roughly up to twice the raw data, geometry/log dependent).
Do not delete raw_workers/ for PI studies relying on particle ancestry.


### Diagnosing missing native Geant4 hit and reference-plane files

If smoke reports missing HitsMuonSpec.root and MuonReferencePlanes.root
despite na6psim_parallel reporting exit code zero, this is NOT a valid run.
A stale installed libsimLib.so or alternative active-detector layout is
the main possibility. Git pull does NOT rebuild installed C++ code.

The updated campaign verifies that the *actual dynamically linked*
libsimLib.so contains both native Stage-C diagnostic markers BEFORE
starting Geant4. If missing, rebuild this branch with the existing
CMake build tree, then retry smoke:

    git pull --ff-only
    bash test/mwpc-stage-c/rebuild_pi_geant4_stagec.sh
    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh smoke

A failed smoke with old frozen binary hashes is safely archived under
a sibling smoke_failed_timestamp directory if no chunks completed.
Completed runs are NEVER silently mixed across binary/config changes.
If a post-run file check fails, the exception now lists actual ROOT
files and extracts relevant lines from preserved worker logs.

For manual confirmation, if needed:

    tail -n 90 test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1_smoke/jpsi/chunk0001/raw_workers/worker000/na6psim.log

Do not launch full production until smoke shows PI_GEANT_CHUNK_PASS
for all three channels.


### Old vs MWPC installed prefixes (confirmed in user build log)

The default login environment may still resolve **the old**:

    ~/na6p/install/NA6PRoot/bin/na6psim

even after a fully successful build and install of **this branch**:

    ~/na6p/install/NA6PRoot-mwpc/bin/na6psim

The campaign runner now sources activate_pi_stagec_env.sh and uses the
single CMake build whose CMAKE_HOME_DIRECTORY is the active Git checkout.
It prepends THAT install prefix to PATH, LD_LIBRARY_PATH and
ROOT_INCLUDE_PATH, and sets NA6PROOT_ROOT for the subprocess only.
The Python preflight verifies that the executable, wrapper and ACTUALLY
loaded simLib are from the exact same prefix. Old NA6PRoot remains intact.

After a previously successful compile/install, do **not** rebuild again:

    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh smoke

The failed old smoke run is archived automatically when its library hashes
differ; no manual cleanup. Do not run full until all 3 channel smoke
chunks pass their ROOT object and tree checks.


### Required four-worker parallel smoke before long production

The one-worker, three-events/channel smoke test covers native Geant4 ROOT
serialization and the reference-plane observer. It does NOT exercise merging
and event offsets from four independent workers. A separate fixed
12-events/channel, four-worker regression now checks exactly that pathway:

    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh parallel-smoke

This is kept strictly separate under

    test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1_parallel_smoke/

The ordinary full campaign output remains

    test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1/

Require three PASS lines and CAMPAIGN_COMPLETE before starting 30k/channel
with four workers. The test checks the merged mckine/hitsMuonSpec ROOT entry
counts, the nonempty mandatory ROOT files, the Geant4 gas-step and
MuonReferencePlanes trees, and original per-worker output directories.

Long campaign:

    bash test/mwpc-stage-c/run_pi_geant4_campaign.sh full 30000 3000 4 2.5

The final argument is the approximate *per-invocation* wall-clock budget,
checked at chunk boundaries; another invocation with the exact same command
resumes from completed checkpoints. It can overshoot by the running chunk.


### Per-chamber wire multiplicity inside the interactive station app

The app now shows the full (P(N_{​anodes}=0,1,2,3,4,5+)) for the
**selected physical chamber**, conditional on a projected gas crossing.
A new `spatialLayer = p2` colors every module by its local probability
of **at least two collected anodes** using a fixed absolute 0–50% color
scale to compare source channels and stations.

**Immediate result with NO microscopic or Geant4 recalculation:**
the earlier completed
`test_runs/mwpc_stage_c/pi_per_chamber_response/pi_wire_multiplicity_per_chamber.csv`
contains all 3*298 old per-module PDFs. The builder checks
`completed_source=1`, that observed/existing exposure exactly matches
the nominal module histogram, and that (P0+P1+...+P5plus=1). It
automatically embeds these results, clearly marked as the **HISTORICAL
RANDOM-PHASE** model. That model randomized wire-grid origin per gas
crossing, and did not clip out-of-gas anode destinations to a finite
77/128-element grid. Hence these histograms must NOT be interpreted
as the distribution of the already visualized fixed numbered wires.

You can generate the upgraded HTML while the native Geant4 campaign is
running; this does not rerun Heed, Garfield or Geant4 and reads only
existing small summary and fixed-wire count files:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_station_dashboard.sh

Then in the UI select: "Пересечения мюонов" -> "Кратность: P(≥2
активных анодов)"; click any module to view the complete six-bin
probability histogram and its binomial 1-sigma statistical uncertainty
on P(≥2). For modules with zero simulated exposure the histogram is
correctly undefined, rather than artificially zero.

**Later**, after the long Geant4 production has finished, derive
exact per-module (N0/N1/.../N5plus) for the same fixed nominal wire
grid as the per-anode maps:

    bash test/mwpc-stage-c/run_pi_module_multiplicity.sh

That script does a NEW slow Heed/Garfield empirical *fold* over the
2.49 million existing historical projected gas crossings, but does
NOT rerun native Geant4 or microscopic Garfield electron trajectories.
It creates `pi_fixed_wire_module_multiplicity.csv` for all 3 channels,
298 det IDs, 3 mechanical phases and 4 definitions (geometric, ≥1e,
≥5e, ≤100 ns). The builder prioritizes its **nominal 0-mm phase**
when available, so its six-bin bar graph will then be consistent with
the fixed-wire heatmap. The old PDFs remain a cross-check, never
silently spliced with the fixed-wire probabilities.



### Muon Track Explorer: browse individual muons and candidate collected anodes

A separate new **offline event display** is implemented in
pi_muon_track_explorer_template.html and
build_pi_muon_track_explorer.py. The existing station statistics dashboard
links to it. The builder consumes *existing historical 2.49M projected gas
crossing CSV rows* but deterministically selects an EVENT-level subset,
preserving (source, job, event, track, PDG) identity and therefore the
mu+ and mu- candidates in the same event when both are present.

**While native Geant4 is busy**, inspect trajectories/geometrical anode
cells only without loading Heed/Garfield:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh geometry-only

**To see individual charge-collection outcomes** from the ALREADY saved
Heed cluster library and Garfield transport lookup (no new Garfield or
native Geant4), run after sufficient CPU/RAM is free:

    bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh

Both commands write and open:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_muon_track_explorer.html

The new interactive display includes:
- source selector J/psi / phi / omega, paired mu+/mu- selection by
  genuine same source/job/event, next/previous sampled events;
- two projected trajectories X(Z), Y(Z), as polylines through ACTUAL
  saved-reference projected gas intersection midpoints (the lines
  between stations are NOT actual Geant4 step paths);
- six station buttons, with actual gas modules traversed by the
  selected projected track highlighted in XY;
- selection of a crossed detector module; visualization of ALL 77/128
  fixed nominal anodes and individual collected-primary-electron count;
- geometric wire cells, >=1 collected electron, >=5 electrons, and
  >=1 electron within 100 ns as distinct selectors.

The Heed+stored-Garfield charge table is generated ONLY for displayed
crossings and is one reproducible stochastic realization, individually
seeded using SHA256(source,job,event,track,station,det). It does not
pretend that the old aggregate occupancy map recorded individual
per-event wire identities, nor does it treat fixed-nominal anode counts
as a HARDWARE fired-channel simulation. Model intentionally excludes
avalanche gain, Shockley-Ramo signal and shaping, electronics noise,
threshold, and cross talk.

Importantly, all display event choices are deterministically curated
for interesting crossings, not representative of the source population.
No physics probabilities should be inferred by counting events shown.
No gas crossing in MSn is NOT evidence the track did not cross that
reference plane, since historical input only includes GAS intersection
records. The 6 station trajectory joins are from reference-state
projections, not native current-geometry Geant4 gas hits.

The display is intentionally separate from native 200k/channel new Geant4
campaign outputs; later the viewer can import the new actual MWPC gas
hits/Edep after track-ID and parallel-worker ancestry reconciliation.

### Muon Track Explorer visual validation and anode zoom

The browser display was verified by a user-provided full-page screenshot
in Chrome/Windows+WSL. With two valid MWPC gas crossings at the final
station, the viewer previously rendered all 128 wires at once; this was
accurate but made one active anode almost impossible to inspect.

The new viewer defaults to **Вблизи активных (±10)**, displaying a
CONTIGUOUS window of real numbered 4-mm-pitch anodes around the selected
active wires, intersected geometric cells and projected gas midpoint.
Select **Все 77/128 анодов** to revert to the full detector anode span.
In the focused view the annotations show collected electrons per anode,
or within 100 ns when the corresponding timing criterion is selected.
An amber dashed guide marks the projected Y position of the gas midpoint,
not a measured reconstructed coordinate.

On the station XY geometry view a dot drawn over the chamber rectangles
now locates each projected gas crossing; selecting either point opens the
respective module. Symbol size is VISUAL, not a physical spot radius.

To refresh the stored HTML without stopping the long Geant4 production,
rerun the previously validated builder (does not modify Geant4 outputs):

    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh

If memory/CPU resources are constrained, `geometry-only` remains
available. The charge view simulates PRIMARY-electron collection; a
real electronics-fired-channel model remains a future stage requiring
avalanche charge gain, induced current/weighting fields, shaping,
threshold, noise and cross talk.



## Anode-only electronics development after the 200k/channel Geant4 production

The user requested the strongest possible **anode-only** event response model
with no cathode U/V strip simulation. This is now divided into three auditable
layers and **does not overwrite the frozen Geant4 campaign**.

### First: verify what "long run finished" means

Because the native production was launched with a 2.5-hour *per-invocation*
budget, ending is NOT sufficient proof that all 200k events per source ran.
Run a fast **READ-ONLY** metadata/root-validator audit:

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    python3 test/mwpc-stage-c/audit_pi_geant4_campaign.py

The audit prints completion counters Jpsi/Phi/Omega, verified Geant4 Edep
step entries and interpolated muon reference-plane states from existing ROOT
validator logs, and requires a separate campaign completion certificate.
A complete run displays FINAL CAMPAIGN STATUS: COMPLETE. If it says
INCOMPLETE, inspect the missing chunks; the runner can be resumed with
the **exact original command**, including 200000/3000/4/2.5. Do NOT
switch this campaign's original configuration to another total.

### Second: native Geant4 gas hits in the Track Explorer

Use the FIRST **validated merged** 3000-event chunk of EACH of the three
sources to inspect actual gas detector response (read only). This prototype
avoids huge all-chunk CSV exports before the native dataset has been checked:

    bash test/mwpc-stage-c/run_pi_native_anode_first_chunk.sh

This script:
- runs no new Geant4 or Garfield microscopic jobs;
- activates the correct installed NA6PRoot-mwpc ROOT library environment;
- exports true Geant4 MWPC gas hits from HitsMuonSpec.root and
  geometry.root with the preexisting TGeo local-frame converter;
- **JOINS TParticle PDG** by remapped hit trackID to the matching
  MCKine.root tracks branch, within that merged chunk;
- excludes all nonmuon hits and preserves exact Geant4 actual Edep,
  entry/exit and local u/v, times and 3-momentum for that specific hit;
- samples Heed electron clusters from the existing Heed library
  *independently*, NOT by converting Edep into a fake exact electron count;
- selects up to 32 interesting native same-event groups/channel,
  then samples existing Garfield electron transport endpoints, including
  individual drift times, for a fixed nominal 77/128-wire 4-mm-pitch grid;
- applies the explicitly preliminary anode-only gain/electronics surrogate;
- saves/opens a separate page:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_muon_track_explorer_native.html

Native export logs and first-chunk true muon-hit CSVs are preserved under:

    test_runs/mwpc_stage_c/pi_native_anode_first_chunk/

This is a new dataset. The historical projected-crossing event explorer
remains unchanged at pi_muon_track_explorer.html. The native page explicitly
labels true Geant4 gas hits and independent Edep. The links connecting
true chamber hit midpoints between stations are VISUAL, not step trajectories.

### Third: anode-only preliminary analog/digital response

pi_anode_digitizer.py consumes individual saved Garfield endpoint times and
the fixed anode channel number. It draws an independent positive avalanche
gain for every primary electron; zero-gain avalanches remain zero. It
can use the original **unconditioned** Phase-A gain-fluctuation CSV with
collected_on_wire values INCLUDING ZERO, via the builder option
--anode-gain-csv /path/to/gain_fluctuations.csv; otherwise it uses a
clearly labelled **UNCALIBRATED** gamma/Polya distribution.

The stored configurable fast response then approximates the fast electron
and much slower ion-induced charge, convolves these *SURROGATE* currents
with a normalized nominal CR-RC shaper, converts to hypothetical mV
and adds reproducibly seeded sample noise. It saves an mV waveform for
each wire that collected at least one primary electron (not 128 traces
per every chamber, so noise-only hits are not yet in the efficiency).
The browser redraws the waveform and firing pattern INSTANTLY when the
user changes the threshold.

Defaults are ONLY illustrative (not Stage A/B calibrated detector truth):
mean gain 10000, shape 1.5, zero-avalanche probability 5%, prompt
fraction 20%, fast current 5 ns, ion current 1200 ns, front-end
peaking 80 ns, conversion 0.35 mV/fC, RMS Gaussian sampled noise
0.35 mV and nominal threshold 5 mV. These values must be replaced
by a real **anode Shockley-Ramo weighting field/current impulse** and
measured ASIC timing/noise/threshold before quoting predictions.

The default historical charge viewer also now contains this anode-only
study-level voltage waveform:

    bash test/mwpc-stage-c/run_pi_muon_track_explorer.sh

For faster/no-anode runs, the underlying Python builder supports
--no-anode-digitizer. The geometry-only runner bypasses Heed/Garfield
entirely. The anode model does **not** calculate cathode strip waveforms,
U/V stereo coordinates, strip sharing, strip thresholds or cathode CoG.

### Physical next validation before final numerical claims

1. Validate long-run per-source gas hit rates, Edep/path/angular
   distributions and geometric bias, especially low-exposure chambers.
2. Obtain or calculate an UNCONDITIONED electron-gain ensemble for the
   actual 30-um wire, 4-mm pitch, gas mix, field and HV; do not reuse
   Stage-B conditioned-nonzero avalanche visualization as gain statistics.
3. Generate a genuine anode weighting-field Shockley-Ramo
   **single avalanche** current library (electron AND ion signals);
   compare with direct microscopic Garfield on controlled few-electron
   crossings, then replace the surrogate impulse.
4. Obtain preamplifier transfer function/shaping, equivalent noise charge
   and discriminator/TOT settings from PI/electronics colleagues.
5. Simulate noise-only channels and any electronic cross talk before
   estimating final single-wire efficiency or fake hit multiplicity.

The new anode-specific code is separated from the cathode strip Stage-B
algorithms. Reusable genuine Stage-A gain statistics are consumed by
the importer rather than duplicating avalanche physics. The user's
large Geant4 output is always treated as an immutable transport library.


### ALL 201 real-Geant4 chunks: PI chamber Edep, entrance/exit angles, overlap QA

The completed full campaign has 67 chunks/source and 200,000 events/source
(600,000 total). Its audited 59,046,706 passively interpolated muon Z-plane
states and 7,671,338 positive-Edep Geant4 gas steps are immutable raw
transport evidence. The successful new native first-chunk exporter selected
29,488 / 23,301 / 21,533 true **muon** gas-hit rows for J/psi/phi/omega
out of 3,000 generated events/source.

New read-only, resumable \`pi_native_full_campaign_qa.py\` covers **all native
chamber gas-hit segments**, not the old approximate geometry projections.
It uses \`export_c_geant_muon_crossings.C\` with the ACTUAL source/chunk-specific
MCKine truth PDG join and TGeo transformed local entry/exit momenta.
The analysis runs at most ONE temporary CSV at a time (no huge full-run CSV)
and saves per-chunk QA checkpoints. Re-running skips verified completed
checkpoint chunks. Its input checks include frozen campaign completion and
source ROOT file sizes/mtime. The checkpoint epoch binds both exporter and
analyzer source revisions.

Run one first chunk per channel with synthetic conservation tests:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_full_qa.sh first

Check for \`PASS: PI FULL native gas Edep...\`. Then run all 67 chunks/channel:

    bash test/mwpc-stage-c/run_pi_native_full_qa.sh all

Outputs:

    test_runs/mwpc_stage_c/pi_native_full_quality/
      pi_native_all_chamber_quality.csv    # every 298 det ID * all 3 sources
      pi_native_all_station_quality.csv    # six stations * 3 sources
      pi_native_all_histograms.csv         # per module Edep/angles/path bins
      pi_native_all_QA_report.json
      checkpoints/                         # resumes after power/interruption

Each chamber reports true muon gas-hit frequency, energy and path length,
fraction with Edep>0, approximate Pearson Edep-vs-length correlation and
physics flags. The histogram table contains conditional distributions
for eight variables: Edep, gas path, theta_u/w **in and out**, and
wrapped theta_out-minus-theta_in. Each bin retains its exact edges and
counts. Station summaries distinguish:
- distinct remapped muon tracks with gas hits in the station;
- the number visiting **two or more separate physical chambers**;
- extra distinct chamber crossings beyond the first per track;
- multiple hit segments in the **same** chamber (do not conflate with
  physical detector overlap).

This is **NOT the complete PI deliverable yet**: the analyzer contains
no Garfield electrons, drift-time acceptance, avalanche gain, hardware
trigger, or per-wire 0/1/2/3/4+ probabilities. Those are the separately
defined anode-response layers. In addition, a muon with a hit in active gas
is not identical to every muon incident on a station reference plane.
The true "per station incoming muon" denominator must be extracted
separately from the new MuonReferencePlanes.root snapshots and must
not be silently replaced with gas-crossing tracks or generated events.

The full run computes with existing ROOT files ONLY and performs NO
Geant4 transportation, Garfield electron avalanche or cathode simulation.
It may take a while to read all 201 saved chunks, but will resume safely.


### Physics readout after the completed full native QA

After running \`run_pi_native_full_qa.sh all\` through Omega chunk0067
(the final chunk contains only **2,000**, not 3,000, generated events),
run the independent **read-only** physics and conservation summary:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_qa_readout.sh

It checks completion provenance and 894 chamber/18 station coverage,
all 8 conditional histograms per chamber (7,152 independent spectra,
including consistent zero-exposure handling), per-detector bin
conservation, native muon/nonmuon count conservation, and the sums
of chamber gas hits within each station. It separately reports unique
gas-hitting muon tracks, real multi-module overlaps and repeated
same-module hit segments. It writes the human-readable:

    test_runs/mwpc_stage_c/pi_native_full_quality/pi_native_QA_readout.txt

Share that small text file with the analysis assistant to interpret
physical outliers and decide whether electron-drift-time and incidence
angle binning matter. The readout is the next PHYSICS inspection, not
yet a per-anode firing result. An Edep histogram by itself does not
predict how many physical anodes collect Heed electrons, and a generated
event is not necessarily an incoming muon at a specific station.


### Where is the arrival-time distribution plot? (PI timing question)

The full 201-chunk Geant4 QA report alone contains Edep and entering/exiting
momentum angles, but **not** Garfield electron drift-time distributions.
Geant4 \`NA6PMuonSpecHit::getTime()\` is a particle crossing time, not a
primary-electron arrival time and not electronic threshold crossing.

An unbiased **first-chunk pilot** now resamples Heed on EACH *native* muon
gas-hit segment (3,000 generated events per Jpsi/Phi/Omega, not 32 curated
Track Explorer events) and uses already saved Garfield electron endpoint
time outcomes, without running microscopic Garfield, Geant4, or cathodes:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_drift_times.sh

The synthetic regression checks 0-electron hits, true times >2µs, unknown
time handling, finite-grid anode counts and time-window dependence.
Outputs:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_drift_times.html

    test_runs/mwpc_stage_c/pi_native_drift_time_first_chunk/
      pi_native_drift_time_histograms.csv
      pi_native_time_window_multiplicity.csv
      pi_native_drift_time_provenance.json

The HTML is an interactive **true native Geant4 gas-segment / Heed+saved
Garfield single-electron** timing study, switchable across Jpsi/Phi/Omega,
all 6 stations, every one of 298 modules, time ranges 0–500/1000/2000 ns.
It plots 20-ns-bin distributions of collected primary-electron drift times,
together with \`P(N_active_anodes >= 2)\` as a function of individual
electron arrival-time window (50, 100, 200, 400, 800 ns or all collected).
A full 0/1/2/3/4/5+ conditional multiplicity table appears alongside.
The station and native Track Explorer dashboards provide direct links.

**Critical interpretation:** this is not the already-finished full
200,000/source native Geant4 response simulation. The Geant4 geometric
truth IS real for the first 3,000 events/source; sampled Heed clusters
and Garfield transport endpoints are a response model (B=0, approximate
asymmetric gap mapping, finite nominal 77/128 anodes). The denominator
is one **gas-hit segment** (including cases with 0 collected electrons),
not one incoming station muon, because module overlaps and repeated
segments exist. Times refer to PRIMARY electron arrival and NOT to
analog shaped pulses, avalanches, actual ASIC or hardware-fired wires.
The full 201-chunk anode timing fold is a subsequent CPU-intensive job;
do not present the pilot PDF as full-campaign statistics.


### PI angular acceptance vs chamber XY + centred-anode timing work (2026-10-19)

Full technical specification (data fields, coordinate definitions, two
different meanings of coordinate-by-time plots, acceptance criteria):

    test/mwpc-stage-c/PI_CENTERED_ANODE_ANGLE_TIMING_SPEC.md

**Fast first output — full 201-chunk Geant4 angular position dashboard:**

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_angle_position.sh

This is pure Python and JavaScript; reads the already completed
pi_native_all_chamber_quality.csv, pi_native_all_station_quality.csv and
pi_native_all_histograms.csv and uses NO ROOT/Geant4/Garfield or transport.
It validates 894 source×chamber records, native station-vs-module exposure,
and all four signed TGeo-local entrance/exit u/w angular PDFs per module.

Opens:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_angle_position.html

Browse one of six stations, source Jpsi/phi/omega, signed theta_u/theta_w
entrance or exit and one of these EXACT bin-edge probabilities:
P(theta>=0deg), P(abs(theta)>=10deg), P(abs(theta)<=2deg).
Select any one of 298 chambers for its original nonuniform-bin signed
angle density in percent/degree. Nothing is falsely presented as an exact
mean or percentile reconstructed from the binned QA.

**Independent new microscopic response (different physics computation):**

    bash test/mwpc-stage-c/run_pi_centered_transport_box_map.sh 50 180019

This runs an existing STAGE-C GARFIELD drift generator, but with symmetrical
cathode gaps 3.0+3.0 mm and anode v=0. It does NOT rerun the completed
600k-event Geant4 transport. It DOES take time to compute 30,000 fresh
Garfield electron paths (600 boxes × 50 per box, B=0). It uses its own
versioned output directory:

    test_runs/mwpc_stage_c/c_transport_box_map_gm3p0_gp3p0/

The wrapper freezes seed/geometry/sample count, refuses incompatible
partial maps, and validates all new table rows/box centres using
pi_validate_centered_box_map.py. The output manifest is
pi_centered_transport_manifest.json.

**Do not apply earlier +1-mm response-coordinate shift to this new map.**
Existing 2+4-mm Garfield maps and their plots are legacy approximation.
Running the centred Garfield producer alone does NOT automatically switch
existing timelike pi_native_drift_times or Muon Track Explorer to the
new physically correct field: a separate centred native Heed adapter,
pilot, genuine Garfield trajectory-snapshot exporter and full fold must
be implemented and validated according to the specification.

The full PI suite must distinguish:
(i) a time series of actual drifting-electron positions P(u(t)), P(v(t)),
P(w(t)), which requires Garfield intermediate electron trajectories,
(ii) conditional P(t_arrival | u_birth,v_birth,w_birth) using endpoint times,
(iii) Geant4 muon entry/exit time, (iv) actual electronics comparator time,
not conflated. There is no cathode-strip signal work in this scope.

### Centered 3+3mm native Geant4 muon timing pilot (2026-10-19)

The user-provided centered Garfield response CSV SHA256 matches its
validated manifest. It contains 30,000 single-electron drift trials,
24,291 collected (80.97%), median collected time 29.62ns, p95 67.45ns,
and 98.711% collected by 100ns. These are uniformly sampled response
map entries, NOT an efficiency estimate for real muon detector hits.
Negative and positive normal-v sides yield very similar collection and
timing; the response should still be folded with real Heed birth positions.

First run (synthetic test + 100 true native Geant4 muon gas hits/channel):

    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_centered_native_timing_pilot.sh smoke

Then full FIRST CHUNK from 3000 generated events/source:

    bash test/mwpc-stage-c/run_pi_centered_native_timing_pilot.sh all

Requires an already completed centered Garfield gm3p0_gp3p0 run and
pi_centered_transport_manifest.json. Validates actual SHA256 response,
20x30 geometry and local v=-3..+3 mm. Never maps v+1 (legacy 2+4mm).
Real MCKine-tagged TGeo native gas-hit coordinates seed independent Heed
cluster positions, and saved centered Garfield electron collection/time
draws determine finite 77/128 physical-anode multiplicities.

Output interactive graph:
    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_native_timing.html
Output machine-readable binned first-chunk results:
    test_runs/mwpc_stage_c/pi_centered_timing_native_first_chunk/
      pi_centered_native_first_chunk_angle_timing.json

Graph: 10ns electron arrival-time PDF by signed muon theta_u and theta_w,
1D histograms of *electron birth positions* u0/v0/w0 in different ARRIVAL
time slices, P(two or more anodes) versus timing window and full N0..N5+.
These are NOT intermediate u(t)/v(t)/w(t) electron drift trajectories.
Those require a new microscopic Garfield trajectory step recorder.
Do not claim full 201-chunk native timing, chamber detection efficiency,
or a real electronic discriminator threshold from this pilot.


### NEW centered 3+3mm native Geant4/Heed/Garfield time-angle pilot

Prerequisites: finished native first chunk exports and VERIFIED centered
Garfield B=0 map and pi_centered_transport_manifest.json.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_centered_timing_pilot.sh

The pilot consumes the RAW MCKine-joined true native gas-hit exporter
CSVs at pi_native_anode_first_chunk/native_exports/ (all local u/v/w
coordinates), not the older compatibility track-viewer format without
the w component. It loads the NEW centered 3+3 mm response with
manifest-pinned SHA256 and checks original electric-field provenance.
The historical 2+4mm lookup with the legacy +1mm offset is NOT used.
The seed is deterministic per (source,chunk,event,hit_index,track).

The primary Heed clusters are sampled independently of Geant4 Edep
along actual TGeo-local true gas-hit line segments. No analogue anode
weighting current, cathode strips, ASIC threshold or hardware firing
is claimed. Exactly 77 MS0 or 128 other physical anodes are retained,
with overflow electrons recorded as losses.

Outputs:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_timing_pilot.html
    test_runs/mwpc_stage_c/pi_centered_anode_timing_pilot/
      pi_centered_native_hit_multiplicities.csv
      pi_centered_time_angle_multiplicity.csv
      pi_centered_timing_manifest.json
      pi_centered_PI_pilot_review.zip

Open the interactive HTML and choose Jpsi/Phi/Omega, MS station,
one physical chamber or all, theta_u or theta_w and a signed angular
bin. The endpoint Garfield t_drift PDFs use 10-ns bins. Heatmaps show
initial u_cell, u, v, or w position versus FINAL electron arrival
time for all chambers of a station at the selected angle. Probability
rows P(N=0,1,2,3,4,5+) are gas-hit-inclusive and show 50/100/200/
400/800-ns drift gates and all collected endpoints. Each N is a count
of DISTINCT finite anodes with >=1 collected primary electron.

Scope is first 3000 generated events/source, NOT the full 200000;
native GAS HIT SEGMENTS are the denominator, NOT incident station
muons. Independent Heed ionisation is not forced to match Geant4 Edep.
Birth-coordinate vs arrival time is NOT a true intermediate
u(t)/v(t)/w(t) electron cloud; next stage needs microscopic Garfield
drift trajectory snapshots. Attach the review ZIP and terminal output
to the analysis assistant after a successful local run.


### PI centered 3+3mm dashboard V2: click-on-XY instead of chamber dropdown

After successful run of run_pi_native_centered_timing_pilot.sh, the pilot
contains a FULL self-contained set of the actual initial-coordinate/time
histograms and per-module wire multiplicities in the generated HTML.
That entire 4.8MB original physics payload can be reused without
running Geant4, Heed or Garfield again:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_centered_timing_dashboard_v2.sh

The builder first validates 298 physical chamber positions and that all
source/station/module/angle/time-bin gas hit counts are conserved.
It then only changes visualization, writing:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/
      pi_centered_timing_dashboard_v2.html

No original source CSV/HTML, primary electron response, truth events or
model manifest is modified. The old PI viewer is preserved separately.

Changes explicitly requested by the PI/researcher:

- Select a physical DET by CLICKING its actual XY chamber rectangle
  (visible global X/Y in cm, axes/ticks and module IDs). A separate
  "All chambers" button restores the station-total view; each camera is
  colored by the 100ns P(N>=2) for the current muon angle selection.
- Signed entrance muon theta_u or theta_w selector now lists ONLY bins
  with >0 native gas-hit exposures in the selected source, station and
  module, including their sample sizes. Selection of an empty bin resets
  to a known available all-angles view.
- Initial electron coordinates clearly defined: u_cell relative to
  nearest wire, full local u across active camera, local v from -3mm
  cathode to +3mm cathode with anode at v=0, and w along wire.
- Instead of an unexplained RAW COUNT heatmap, show independently
  normalized one-dimensional P(birth_coordinate_bin | Garfield
  electron ARRIVAL time interval), with 5 colored 1D lines for
  t_drift=0-20,20-40,40-60,60-100,>=100ns. Each curve integrates to 100%
  over the plotted coordinate bins for that time band (conditional on
  collected electrons with known drift-time). THIS IS NOT P(u(t)):
  the true moving-cloud position snapshots still need a Garfield
  trajectory recorder.
- Existing source only stored electron coordinate/time correlations
  aggregated over the STATION and signed MUON ANGLE, not each module.
  Therefore the V2 coordinate plot explicitly says the DET selection
  DOES NOT change it. Module selection DOES change the time distribution
  and anode-multiplicity charts. This deliberate visibility prevents
  the previous silent mixing of selection scopes.
- Active-anode multiplicity chart has auto-zoomed Y axis (% of real
  native gas-hit segments) and the physical formula: distinct finite
  physical anodes receiving >=1 collected primary electron before T,
  never primary-electron count nor ASIC threshold crossing. Zero-hit
  tracks/hits are retained. Table rows display complete mutually
  exclusive P0,P1,P2,P3,P4,P5+ and sum P>=2.
- ALL SVG plots now label X and Y axes with units. Screen captures
  remain reproducible in original global XY frame and angular binnings.

The second requirement, the future p=3,3.5,4mm anode *pitch* sweep, is
separate from this visualization fix and still requires a new
Garfield electric-field calculation at each pitch. Do not fake it from
this single p=4mm response.


### NEW: one native muon and TWO independent physical anode signals

PI question: how do TWO physical anode channels respond to a single native
Geant4 muon crossing, including electron arrival-time offsets, pulse
shapes and (illustrative only) arithmetic sum of the two channels?
Explore real Geant4 gas-hit segments from the already verified first
3000 generated events/source (Jpsi/phi/omega); select chamber on global
XY, select one muon event by signed theta_u, then inspect independent
anode A/B voltage waveforms, Garfield endpoint electron arrival rasters,
wire numbers and separation, and an artificial A+B overlay.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_two_anode_signal_explorer.sh

Output:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/
      pi_two_anode_signal_explorer.html
      pi_two_anode_signal_explorer_manifest.json

The runner reads frozen first-chunk true native Geant4 gas export CSVs and
the completed 3+3mm centered PI pilot per-hit CSV. It filters exactly
N_all=2 independent physical wire IDs and retains up to seven
angle-spread examples per source×module, including late second-wire
collection if available. Most camera/source combinations contain such
events, but some have none. Selected examples are deliberately curated;
NEVER use this set for unconditional detector probabilities or efficiencies.

The Heed+Garfield electron sampling replays the EXACT per-native-hit
SHA256 seed used in pi_native_centered_timing_pilot.py, and validates
N_all=2 and N_t100 against the existing event table. Both wires are
simulated independently by the current pi_anode_digitizer.py:
UNCONDITIONED Stage-A gain file optionally configured, otherwise
an EXPLICITLY UNCALIBRATED Polya+zero gain model, fast electron and slow
ion surrogate current, ideal nominal CR-RC, fixed independent noise and
user-settable display discriminator threshold (default 5 mV).

The display plots TWO distinct per-wire voltage channels on the same
time axis and a THIRD DASHED hypothetical pointwise A+B VOLTAGE SUM.
The real MWPC has two separate anode readouts; this arithmetic overlay
does NOT represent a real analogue summing circuit or an electrode's
induced charge. The interface flags when TWO actual independent
surrogate channels cross the requested voltage threshold, versus
merely collecting at least one primary electron each. Zero threshold
crossings can occur on one or both wires despite N_all=2.
The amplifier output time is NOT equal to Garfield electron first-arrival
time. Time origin here is approximate local electron birth and omits
finite muon-ionisation creation-time spread. This digitizer is NOT a
calibrated Shockley-Ramo anode weighting field or ASIC response.
Further realism requires computing anode weighting potentials, actual
electron+ion avalanche waveforms, cross talk and verified ASIC settings.

Gallery quick picks try near-normal (<3deg), moderate (3..10deg), and
oblique (>=10deg) native muon theta_u; only available real events are
shown. All charts and station XY map have labeled axes and physical
units. No cathode strip output is modeled.

No new native Geant4, Garfield drift or ROOT export is performed. Run
the tiny native/synthetic regression before production; the wrapper
does this automatically, then generates the HTML. The old 2+4mm response
and earlier Track Explorer are untouched.


### PI: actual Geant4 muon u-v and true Garfield microscopic cloud u(t),v(t),w(t)

After the previously completed centered 3+3mm Heed/Garfield
two-anode pilot is available, run:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_two_anode_electron_explorer.sh

This builds pi_two_anode_electron_explorer.html in
test_runs/mwpc_stage_c/pi_fixed_wire_maps/. It reuses collected primary
electron arrival timestamps from the old 3+3mm per-event HTML, but NEVER
the arbitrary model mV voltages. It joins each selected muon to the TRUE
TGeo-local u/v/w Geant4 gas entrance/exit from existing
pi_native_anode_first_chunk/native_exports CSVs (strict source/event/
hit/track geometry matching), and draws the actual u-v segment between
cathodes at v=±3mm and two highlighted finite anode wires at v=0.
Time plot: separate transparent histograms of primary electron counts
on both physical anodes, 5/10/20ns bin selector, optional dashed arithmetic
A+B primary electron COUNT curve, no electronics/analog circuit sum,
no ASIC threshold, no cathode strips. Every axis is labeled.
Original analog viewer remains unchanged in a different HTML.

For genuinely time-evolving microscopic electron positions (not
initial position vs final arrival time), run the NEW, intentionally
limited simulation:

    bash test/mwpc-stage-c/run_pi_centered_true_1d_cloud.sh

This command first tests snapshot conservation on synthetic traces;
next selects independent Heed initial electron locations along ACTUAL
native first-chunk muon gas-hit segments in four |theta_u| strata
(0–2, 2–5, 5–10, >=10 degrees) for every source and MS station.
It rebuilds the EXISTING Stage-C Garfield executable and invokes only
its new --trace-starts mode: symmetric 3+3mm, 4mm pitch, B=0, real
AvalancheMicroscopic::DriftElectron with EnableDriftLines() and 20-collision
sample stride. The original 600-box map, completed Geant4 600k run and
old 2+4mm reference data stay untouched.

The saved real microscopic drift points are evaluated at fixed times
0,5,10,20,40,60,80,100,150ns by interpolation ONLY between adjacent
stored Garfield microscopic samples (NOT the straight line from birth
to final endpoint). Once an electron's trace ends it leaves the
still-drifting population; each coordinate PDF is normalized to the
number of electrons STILL drifting at its time. The survival plot
shows the fraction of launched microscopic electrons not yet stopped.
The standalone pi_centered_true_1d_cloud.html in
test_runs/mwpc_stage_c/pi_fixed_wire_maps/ can display independent
P(u_relative_to_initial_wire(t)), P(v(t)), P(delta_w(t)), and local
P(w(t)), for angle-conditioned real native muon samples.

This is a limited STRATIFIED pilot: default 8 native muon-hit gas
segments per channel/station/angle bin and at most 12 Heed primary
electrons per hit. Therefore these curves do not give absolute muon
acceptance, electronics efficiency or full 201-chunk response statistics.
The new Garfield trajectory mode and HTML require first local execution
in WSL; do not claim the scientific outputs are finished merely because
their code exists in GitHub.


### Geant4 native angular XY dashboard: selected-chamber local u-w + 1D projections

The existing full Geant4 **global XY angular** station dashboard is
extended with local gas-hit midpoint maps of whichever chamber is
selected by clicking its global rectangle (or selecting its DET ID).
The local panel is *u-w*, NOT u-v:
- u across wires, in mm (full gas half-width 153.6mm in MS0, otherwise
  256mm); w along wires (full gas half-width 332.8mm);
- explicit axes/ticks and physical gas boundaries; central (0,0),
  physical 4mm wire-pitch reference marks (every 8th wire drawn);
- modes: nominal local chamber geometry, all native muon gas-hit
  midpoints, and gas-hit midpoints with the existing active angular
  selector (signed theta_u/ theta_w in/out and positive / |theta|>=10 /
  |theta|<2 cut);
- 16x24 two-dimensional *native gas hit segment count per local bin*
  histogram, with its u and w one-dimensional projections in native
  gas-hit segments per coordinate bin (not primary electrons, not
  detector efficiency, not unique throughgoing muon count).

The global XY, complete 200k-event Geant4 angular dashboard remains
unchanged in physics. The local maps require extra actual TGeo-local
hit entrance/exit u,w that the existing full QA histograms do NOT
retain. Their input is exclusively MCKine-linked native GEANT4 gas
hit export, with u_mid=(u_in+u_out)/2 and w_mid=(w_in+w_out)/2.
Their validity definitions are identical to saved QA: finite
in/out u/w signed angles, Edep>=0, gas length>0, muon PDG=+/-13.
Overlaps and repeated same-module segments count separately.

**Fast preview (no ROOT re-export, 3000 events/source)**:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_local_uw.sh first

This uses existing native_exports/*.csv. The page's full-source
angular map STILL uses 200k generated events/channel; the local
u-w occupancy and u/w projections explicitly say PILOT 3000
generated events/channel. Mixing full and pilot scopes under an
undisclosed legend is prohibited.

**Full source (read-only 201 native Geant4 ROOT chunks, resumable):**

    bash test/mwpc-stage-c/run_pi_native_local_uw.sh all

This DOES NOT rerun Geant4 or Garfield. It runs the known working
ROOT native hit exporter on each saved chunk; independent per-chunk
immutable-code+file-stat checkpoints are written in
test_runs/mwpc_stage_c/pi_native_local_uw_density/checkpoints.
The combined per-detector valid muon-hit n and all 4 signed
angle histograms are required to match the existing independent
complete pi_native_full_quality QA exactly, before the dashboard
can claim FULL local source coverage. Interruptions can resume via
the same command; a changed exporter/aggregator or input chunk
requires a clean new --out directory.

Output updated in place:

    test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_angle_position.html

and intermediate local grid data:

    test_runs/mwpc_stage_c/pi_native_local_uw_density/
      pi_native_local_uw_density.json

There is no Geant4 gas re-simulation, no Heed, no Garfield/charge
drift/electronics or cathode-strip model in either local mode.


#### ROOT executable discovery for full local u-w processing

The Stage-C activator enables na6psim and libsimLib but does not
necessarily put ROOT's root executable on PATH. Therefore
run_pi_native_local_uw.sh all now auto-resolves ROOT from an explicitly
configured PI_ROOT_EXECUTABLE, active root in PATH, ROOTSYS,
or ROOT_DIR/CMAKE_PREFIX_PATH/ROOT_LIBRARY_DIR entries from the
EXISTING CMakeCache of NA6PRoot, including ROOT_DIR paths ending in
lib/cmake/ROOT. It populates ROOTSYS, PATH and LD_LIBRARY_PATH for
the ROOT macro exporter before starting the first chunk.
No simulation or ROOT rebuild occurs.

If ROOT is installed outside recorded paths, the program fails
before re-exporting any chunks and provides directions for selecting
the installed bin/root with PI_ROOT_EXECUTABLE. A synthetic regression
covers the WSL missing-PATH scenario using a fake ROOT installation in
a CMake directory and rejects invalid explicitly configured binaries.
Previously completed per-chunk checkpoints remain resumable; all native
input ROOT chunks and their full Geant4 QA remain unchanged.


### Current PI dashboard navigation (no historical/legacy links)

The newest independently useful dashboards are linked bidirectionally
using a common offline navigation header and a self-contained start page:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_latest_dashboards.sh

Open pi_latest_dashboards.html in the pi_fixed_wire_maps directory.

Current pages:
- pi_native_angle_position.html: verified full native Geant4 200k
  generated events/source, global XY plus local u-w positions/1D
  occupancy for pilot or full (as labeled);
- pi_centered_timing_dashboard_v2.html: FIRST-CHUNK centered 3+3mm
  Heed/Garfield electron arrival and time-gated 0,1,2,3,4,5+ wires;
- pi_two_anode_electron_explorer.html: curated native first-chunk
  two-collected-wire event with TRUE local u-v muon segment and
  independent primary electron count-versus-time histograms (NO mV);
- pi_centered_true_1d_cloud.html: micro Garfield P(u(t)),P(v(t)),P(w(t))
  for native angle-stratified electron cloud IF the separate
  microscopic trajectory calculation successfully produced the page.

The common offline navigation tool only links HTML files ACTUALLY on
this workstation. A missing optional microscopic cloud is shown as
not-yet-created, not a dead link. It preserves the exact embedded
source JSON data bytes and does not rerun, overwrite or process ROOT,
Geant4, Heed or Garfield data. It is idempotent, and checks that the
original data script text stays unchanged. Historical legacy 2+4mm,
old projected station estimates, provisional centered native timing,
and parametric-voltage explorations remain on disk for reproducibility
but are deliberately excluded from all new menus.

Modern dashboard builders automatically refresh the latest navigation
when they finish. If only restoring navigation across already-built
pages, run the one-shot command above. Old bookmarked station/track
pages will also acquire links forward to these current analyses;
they are not endorsed as updated physics.


### PI October 20 preflight: quadrants, runtime, LUT precision (no new physics run)

The PI now requests 2 million generated dimuon events **per** Jpsi/Phi/Omega
channel, anode pitches 3/3.5/4 mm with/without grounded field wires,
arrival angle conditional on SIGNED muon incidence, and **mean + centered
RMS + statistical errors** (not median). Details and non-negotiable
physics controls are in PI_CENTERED_ANODE_ANGLE_TIMING_SPEC.md.

Run the fast read-only baseline first:

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_stagec_preflight.sh

This executes synthetic tests and reads the already completed native
201-chunk QA, when present, to write:

    test_runs/mwpc_stage_c/pi_quadrant_preflight/
      pi_quadrant_geometry_preflight.json
      pi_quadrant_native_counts.csv

It also reads existing completed.json run times where available, without
fabricating timings for the old campaign, and produces per-box precision +
vectorized lookup performance for an existing validated centered p=4 mm
Garfield table:

    test_runs/mwpc_stage_c/pi_centered_lut_benchmark/
      pi_centered_lut_box_quality.csv
      pi_centered_lut_benchmark.json

New native Geant4 chunks record perf_counter wall durations automatically
in completed.json. Use a NEW frozen suite for 2M/source, not the existing
200k suite; do not start it until quadrant feasibility and performance
are reviewed. **Do not infer exact quadrant transport symmetry from similar
hit maps:** ABCD z-staggered layers are not fourfold equivalent. The
nearest-anode geometric speed test is not a Garfield physics validation.
These scripts generate no new Geant4, Heed or Garfield electron trajectories.


### Next physics stage: six independent B=0 Garfield wire-pitch/field-wire pilots

After the completed read-only symmetry/LUT preflight, run a tiny, genuinely
NEW microscopic Garfield simulation for the PI's six electrode arrangements:

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_pitch_fieldwire_pilot.sh smoke 75 2

Here **75 um is the selected field-wire DIAMETER for this study**, not
a verified mechanical design requirement. Confirm actual field-wire diameter
and HV with the PI before committing to physically calibrated production.
Field wires are halfway between adjacent +1800 V anodes at potential 0 V;
the cathodes remain at v=+-3 mm and 0 V. v=0 is the wire plane.
The calculation uses Ar/CO2 70:30, B=0 and 30-um anode diameter.

The smoke mode runs ONLY the central u-column for each of the 6 setups,
2 electrons per (u,v) box. This validates new Garfield components, electrode
collection classification and output without producing a physics comparison.
It creates its own immutable namespace:

    test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/smoke/
      all_six_configurations.json
      p3p0_anode_only/
      p3p0_field75p0um/
      p3p5_anode_only/
      p3p5_field75p0um/
      p4p0_anode_only/
      p4p0_field75p0um/

Only after smoke succeeds, a moderate all-unit-cell pilot can be launched:

    bash test/mwpc-stage-c/run_pi_pitch_fieldwire_pilot.sh pilot 75 20

This uses the FULL wire pitch in u, -3...+3 mm cathode gap in v,
du=0.1 mm, dv=0.2 mm; 30/35/40 u bins for p=3/3.5/4mm.
Each individual u-column is SHA-pinned, validated and resumable.
Separate anode/field-wire/not-on-electrode outcomes are recorded.
The pilot outputs arithmetic mean, centered RMS and bootstrap statistical
uncertainties of the anode-collected electron drift time and last-step
electron incidence-normal angle. It does **NOT** reuse or overwrite the
existing verified centered p=4mm map.

**Limitations:** launch electrons uniformly in each cell, independently of
GEANT muon tracks; this first pass cannot show P(theta_e | theta_mu), since
it has NO muon-angle input. Different pitches do not start from identical
global electron positions in this unit-cell scan: use fixed native Geant4/Heed
birth coordinates and common random numbers in the subsequent paired response
comparison. The last-step incidence angle is a microscopic transport proxy,
not a clean field-line angle nor a reconstructed ASIC pulse angle. Results are
not absolute NA60+ response until wire diameters, voltages, gas, B, finite
anode aperture and source muon acceptance are confirmed.

Loader diagnostic for native C++ Garfield runtime: this wrapper now sources
the local Garfield setup, includes Garfield and ROOT library directories in
LD_LIBRARY_PATH, and runs ldd before launching electrons. Failure code 127
indicates a runtime loader or shell-execution problem, not a Garfield physics
result. The Python driver prints the last 45 lines of each failed column log.
If a prior failed attempt left a frozen manifest but no completed CSV columns,
the old manifest is archived automatically on retry after a driver update.
Completed microscopic output is never silently reused across incompatible
physics configurations or deleted.

### October 20: validated six-geometry Garfield smoke (75 um field wires)

The WSL workstation completed ALL six 60-electron Garfield central-u-column
smoke checks: p=3.0/3.5/4.0 mm with and without 75 um grounded midpoint field
wires, 360 tracked electrons total, 0 classified as collected on the field
wires in these CENTRAL-ANODE columns. That is **not** evidence that a field
wire has no effect: the smoke sampled near the anode, not the full unit cell.
The six measured anode-conditional mean drift times are about 18-20 ns with
~1.7-1.8 ns bootstrap errors, and are not meaningful pitch-difference tests.

NEXT new Garfield pilot:

    bash test/mwpc-stage-c/run_pi_pitch_fieldwire_pilot.sh pilot 75 20

This samples all 30/35/40 columns for the full gas cell and both electrode
conditions, 126,000 Garfield trajectories total. The micro-simulation can
take considerable wall time, depending on the host: do NOT infer runtime
from the tiny smoke test. If interrupted between/within columns, run the
same command again. Each completed column is checksum-validated and its
measured subprocess wall time recorded in an adjacent .timing.json file.
Summaries include Garfield subprocess seconds/electron mean, centered RMS
and column-bootstrap uncertainty; missing historical measurements are
represented as unknown, never inferred from timestamps. The sampler
retains exact frozen CPP/binary/geometry/gas/HV/seed identity even when
only its reporting and timing bookkeeping source changes. Not a
muon-angle-conditioned response yet.


### 2026-10-21: full 126k Garfield pitch pilot completed, paired follow-up

WSL pilot run PASSED all six geometries:
- p=3.0mm: 18,000 trials/variant, 14,875 collected without the field wire,
  14,258 with 75um field wire. Mean drift 27.152/26.614 ns,
  centered RMS 14.949/14.727 ns.
- p=3.5mm: 21,000 trials/variant, 17,124 / 16,629 anode collected.
  Mean drift 29.407/28.448 ns, RMS 17.151/16.407 ns.
- p=4.0mm: 24,000 trials/variant, 19,383 / 19,041 anode collected.
  Mean drift 31.964/30.901 ns, RMS 20.247/18.651 ns.
- Exactly 126,000 total Garfield electron trajectories in this initial
  *uniform local-cell* scan. Field-wire endpoints classified as collected
  were zero in all six runs; this is compatible with 0V conductors
  repelling electrons compared with +1800V anodes, but does not establish
  that the classification, field, or endpoint-loss mechanism is validated.

Important: these mean/RMS differences are CONDITIONAL on successful ANODE
collection and may be influenced by differential acceptance of slow
electrons. To distinguish within-electron transport from censoring/losses,
do the read-only paired analysis on the already completed CSVs:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_pitch_fieldwire_paired_analysis.sh

This does NOT launch Garfield/Heed/Geant4 or rerun the 126k trajectories.
It verifies exact matching of each Garfield electron's initial u/v and
random seed across with/without field wire for EACH given pitch; checks
frozen parameters and reports:
- paired 2x2 collected/not-collected transition table;
- unconditional collection fractions with bootstrap errors;
- mean and centered RMS of anode-collected drift times and last-step
  electron incidence-to-normal angles, with errors;
- matched-electron drift-time changes CONDITIONED ON COLLECTION IN BOTH;
- separate time-cut collection probabilities and mean/RMS u-column
  spatial profiles, allowing field-wire edge effects to be localized;
- full endpoint-status counts for every field condition;
- offline interactive local HTML plus machine-readable CSV/JSON.

Files are generated in test_runs/mwpc_stage_c/pi_pitch_fieldwire_paired_v1/.
Open pi_pitch_fieldwire_paired_dashboard.html in a browser. Both conditions
are paired with a stratified bootstrap within fixed (u,v) cells. This is
*finite Garfield trial uncertainty*, not correlated GEANT event error;
different pitches have different u intervals and are NOT exactly paired
to the same native Geant4/Heed electron birth ensemble. A native angle-
conditioned Heed fold must follow for the PI's requested dependence
P(electron incidence angle | signed muon entrance angle).


### Paired Gaussian-drift follow-up: metal birth-point relocation

On the already completed full-pitch 75um pilot, the read-only paired
analysis originally stopped on 'Initial Garfield electron not paired at
14/2'. This is consistent with an *explicit C++ start-inside-electrode
safeguard* in c_transport_box_map.cxx: when a uniform sampled birth lies
inside an installed metallic conductor, the source launch coordinate is
moved to the corresponding (u,v) box center before drift. A field wire at
u=-p/2 or +p/2 can trigger this *only in the field-wire-on geometry*,
so matching Garfield sample ID + seed is insufficient for identical
physical starts. The first observed mismatch is in p=3 mm u-bin 0,
v-bin 14 near the field-wire boundary at v~0.

The corrected pi_pitch_fieldwire_paired_analysis.py permits an exception
ONLY when the bare-run coordinate lies inside the theoretical 75um field
wire and the field-run coordinate is exactly the recorded box center;
all other coordinate or seed mismatches remain fatal. These rare pairs
are tabulated in the JSON as excluded_birth_relocations and excluded from
the paired transition, time and fraction denominators. Raw six Garfield
CSV sets are not modified; the separate ORIGINAL marginal pilot results
still include relocated initial electrons and should be regarded as an
approximate cell-uniform model in metal-overlap regions. There is no
claim that this is a complete physical treatment of primary ionization
inside conductor volume. Changes to the Garfield physics producer should
be versioned as a new independent pilot, not silently patch existing CSVs.

Run without recomputing Garfield:

    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_pitch_fieldwire_paired_analysis.sh

Read printed excluded relocation counts per pitch and confirm that every
excluded launch has a physically justified metal-geometry reason.


### 2026-10-21: first real-native six-geometry 10-ns timing spectra

The NEXT analysis after the successful 126k Garfield six-field pilot and
paired comparison is now implemented as a separate, NON-DESTRUCTIVE
Geant4 native gas-hit + Heed + LUT fold:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_six_geometry_native_timing.sh smoke

This smoke invokes 200 TRUE Geant4 gas-hit segments from each of the three
already exported first chunks. One Heed ionization birth list per unique
native muon gas segment is reused identically for all six Garfield
response configurations. The 20 microscopic trial/cell six pitch-specific
LUTs are unchanged. Outputs and offline HTML are generated under

    test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/smoke/
       pi_six_geometry_native_timing.json
       pi_six_geometry_native_timing.html

On success, launch the FULL existing native FIRST chunk (3,000 generated
parent dimuon events per channel, not 2 million or 3,000 gas segments):

    bash test/mwpc-stage-c/run_pi_six_geometry_native_timing.sh full

Full first-chunk inputs are already frozen MCKine-confirmed native TGeo-local
gas crossings for Jpsi (29488 hits), phi (23301 hits), omega (21533 hits).
No new Geant4 or Garfield microscopic transport is launched. Native muon
angles retain signed theta_u/theta_w; Heed clusters are sampled along
the real chamber-local 3-D segment. Same per-hit seed and same random birth
points/LUT trial quantiles are used across six maps.

Physics dimensions:
- source: each Jpsi/Phi/Omega separately and their joined "all" sample;
- station MS0..MS5;
- individual detector DET (one of 298 physical installed chambers);
- physical anode index in each pitch-dependent finite bank, not just local
  relative wire offset;
- signed muon entrance angular bins, for station and chamber spectra;
- p=3.0/3.5/4.0 mm, no field wire vs 75um grounded midpoint wire.

Existing physical gas transverse width is HELD FIXED: 307.2mm in MS0,
512mm elsewhere. Centered anode-bank wire count N=round(W/p)
(positive floor convention), giving MS0 N=102/88/77 and each other
station N=171/146/128 at p=3/3.5/4mm. At p=4mm this reproduces the
previous 77/128 centered physical-wire phases. A candidate hit on a
nearest infinite-grid anode lying outside this finite bank is lost.

10-ns time bins cover [0,10), ... [290,300) ns; a 31st bin tracks
times >=300ns. The PDF denominator includes the overflow, so plots
never renormalize a truncated tail without noting it. Per group the
raw electron counts, known/unknown collected times, chamber edge losses,
mean, centered RMS, and last-step electron incidence-normal-angle mean
and RMS are stored. The dashboard can select the station, module, anode,
and source, and can overlay all six spectra at matching selections.

CURRENT ACCEPTANCE: only after passing the synthetic regression and the
real 200-segment/channel smoke, run the full first-chunk fold. The latter
has NOT been executed by the assistant's remote tools. Mean and RMS
currently have NO parent-event bootstrap error bars. Because muon gas
segments and electrons are correlated inside a generated dimuon event,
a parent-event bootstrap (or weighted independent parent-event jackknife)
is required before the PI's uncertainties are publishable. Anode-level
"segments" counts segments with >=1 electron collected on that anode,
not the number of incident reference-plane muons. B=0 and no amplifier
or discriminator physics are still strict limitations. Fixed Heed
calibration is a sampled model, not inferred directly from native Edep.

Following this first-chunk check, extend to the PI's newly requested
2,000,000 generated parent events PER channel only with measured runtime,
checkpointed parent-event bootstrap, and independent full Geant4 production
planning. Do NOT label the first chunk or old 200k/source set as 2m/source.

The same first-chunk fold also SAVES real Garfield electron-angle
histograms with 10-degree bins: 0..180deg incidence-to-inward-surface-normal
from the final microscopic transport step; -180..180deg signed phi of the
electron's geometric impact point around the wire. They are each conditioned
on the signed muon theta_u or theta_w bin at station/chamber level, and
on the selected physical anode for angle-integrated native events. The
interactive HTML has an observable selector for drift time, incidence angle
and wire-hit phi, plus a toggle for six-LUT overlay. Angular denominators
are all collected electrons with known angle and may differ from known-time
counts; these are NOT normalized as if they were the same ensemble.


### 2026-10-21: TRUE complete native first-chunk 6-field distribution PASS

The user completed:

    bash test/mwpc-stage-c/run_pi_six_geometry_native_timing.sh full

The synthetic regressions and real fold PASSED. Exactly:
- Jpsi 29,488 MCKine-confirmed native gas segments;
- Phi 23,301;
- Omega 21,533;
- total 74,322 muon gas-hit segments from 3,000 generated dimuon
  parent events/channel (=9,000 total parents);
- 496,809 stored grouped distributions, including native source/station/
  chamber/physical anode plus signed muon angles and Garfield electron
  drift times and two distinct arrival-angle observables;
- all physical-anode --> chamber --> station collected-electron and 10ns
  time-histogram conservation assertions passed.

The original comprehensive JSON and HTML are under

    test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/

The total HTML is potentially huge: 496k groups must NOT all be loaded
into a single browser page to review chamber/anode effects. New READ-ONLY,
standalone local postprocessing:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_six_geometry_postprocess.sh pages

This parses the original huge JSON *incrementally* in stdlib Python
(no ijson dependency) to generate 18 much smaller selectable source x
station offline HTMLs and a compact index, retaining the original JSON.

    .../pi_six_geometry_native_timing_v1/pages/index.html

The overview lists p=4mm anode-only means/RMS for source x station.
Each station page retains time (10ns), approach angle (10deg), signed
wire-impact phi (10deg), plus config/chamber/finite physical anode and
signed muon-angle selections. The splitter writes a small independent
station reference JSON for verifying the exact full-fold results:

    .../pages/pi_six_geometry_station_reference.json

After pages has PASSED, separately bootstrap by the *generated parent*
dimuon events (not independent gas segments or microscopic electrons):

    bash test/mwpc-stage-c/run_pi_six_geometry_postprocess.sh bootstrap

This is a complete re-fold of native first-chunk original CSVs using
identical deterministic per-segment Heed seeds and pre-existing 6 Garfield
LUTs, but launches NO C++ microscopic Garfield, Geant4, or ROOT jobs.
It includes all 3000 generated parent event IDs PER source even if no
true muon gas-hit segments appeared from that event. The bootstrap
draw weights the complete set of hit segments, primary electrons and
six detector-response alternatives of each parent together.

Full source x station x configuration counts, 31-bin drift-time
histograms, 18-bin electron approach-angle histograms and 36-bin signed
wire-impact-angle histograms MUST match the original full-fold station
reference EXACTLY. Independently computed means and centered RMS must
match to a strict floating tolerance; otherwise abort without releasing
any uncertainties.

A first 160-replica parent-bootstrap gives std.dev. statistical errors
for collection fraction, arithmetic mean, centered RMS, and each time
or angle histogram fraction; this does NOT give Garfield LUT finite-trial,
electronics, model, or transport uncertainty and does NOT yet resolve
event-level errors by individual chamber/anode or signed muon-angle bin.
Writes CSV + JSON + offline errorbar dashboard under:

    .../pi_six_geometry_native_timing_v1/parent_bootstrap/

The POSTPROCESS tools themselves are new and must pass the WSL synthetic
and real-data runs; their availability does not imply they already ran.
Do not claim first 3000 parents/source equals the proposed new 2M
generated parent events/source.


### 2026-10-21: original PI V2-style XY map for six Garfield geometries

The user requested that the clickable physical-detector XY map (left) and
10-ns electron drift-time histograms (right) from the ORIGINAL
pi_centered_timing_dashboard_v2.html be extended across all six pitch/
field-wire variants. Its source was located and reused as the visual and
interaction reference: pi_centered_timing_dashboard_v2_template.html and
build_pi_centered_timing_dashboard_v2.py.

The EXISTING frozen six-field first-chunk distribution data and current
detailed per-anode pages remain unchanged. Re-running the read-only
pages postprocessor now ALSO generates 18 new visual V3 HTML pages:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_six_geometry_postprocess.sh pages

Open this NEW V2-style dashboard launcher:

    test_runs/mwpc_stage_c/pi_six_geometry_native_timing_v1/pages/
      pi_pitch_fieldwire_timing_dashboard_v3.html

Each Jpsi/Phi/Omega x MS0..MS5 page has the familiar physical XY clickable
chamber rectangles and camera highlight, the selected physical-chamber
10ns electron drift-time PDF, a toggle overlay of all 6 Garfield
configurations, a selector for any physical anode, a signed muon-angle
selector (chamber/station only), electron arrival-angle PDFs, and a
six-geometry mean/RMS comparison. Same native 3k generated parents/source
and six previously calculated Garfield LUTs; no new Geant4 or Garfield.

IMPORTANT: the full six-field JSON contains single-electron times,
collection and final-step arrival-angle histograms, but NOT the original
V2's 2D primary electron birth-coordinate/time distributions and per-gas-
segment N(T)>=2 distinct-anode multiplicity. The NEW V3 heatmap colors
show P(PRIMARY electron collected on an anode within first 100ns), per
chamber and signed muon-angle interval, with ALL native primary electrons
in that group's denominator. This is NOT the original V2's P(gas hit
with >=2 distinct anodes at T=100ns). The two lower panels therefore
show REAL available Garfield electron arrival-angle spectra and six-field
mean/RMS comparisons rather than fabricating nonexistent physical data.
Anode-specific arrival-angle spectra currently exist ONLY integrated over
native muon entrance angle; when an anode is selected, the signed muon
angle selector is disabled. Parent-event statistical errors by specific
chamber/angle are not yet available.

Original pi_centered_timing_dashboard_v2.html and the original 18 full
per-anode six-geometry pages remain available. The source page
pages/index.html links to the new V3. A Python-created STATIC timing
spectrum remains visible if browser JavaScript aborts, and the page
displays version NA60-SIX-GEOMETRY-XY-V3-20261021.


### 2026-10-21: logarithmic timing PDFs and safe 2M parents/source runner

The original-layout V3 dashboard and the detailed per-anode pages both
now offer a LINEAR/LOGARITHMIC Y-axis for their 10ns drift histograms
(and detailed electron-angle views). The log scale shows actual POSITIVE
percentages, omits zero bins (log(0) is undefined), and keeps the
original normalization INCLUDING the >=300ns overflow. It defaults to
linear. Refresh both page families with:

    bash test/mwpc-stage-c/run_pi_six_geometry_postprocess.sh pages

The PI requested 2,000,000 GENERATED DIMUON PARENTS PER SOURCE for
Jpsi/Phi/Omega (6M total). The existing native Geant4 campaign supports
a MAXIMUM chunk size of 50,000, not 100,000. Thus the production plan
is 40 x 50,000/source, with the first TWO chunks PER source giving
the requested 100,000/source benchmark. Both modes share the SAME
2M/source campaign manifest and all six real physical-chamber layers,
so benchmark checkpoints are safely REUSED by the full continuation.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh plan
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh benchmark
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh status
    # ONLY AFTER actual time + memory + ROOT disk-usage review:
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh full

Defaults: PI_2M_HOURS=2.5 (invocation budget between COMPLETE chunks);
PI_2M_WORKERS=4; PI_2M_MIN_FREE_GIB=20. Changing WORKERS, simulator
binary, ROOT/Geant4 installed libraries, generator, fieldmaps or
detector settings halfway through triggers the existing frozen-manifest
check and MUST NOT be silently accepted. The hourly budget is NOT a
per-chunk timeout and cannot interrupt a running 50k-event simulation
before completion. All original worker-specific ROOT outputs are kept
for unambiguous track ancestry; inspect observed disk use before
starting the full 2M/source run. The disk guard is checked before
each newly started chunk but cannot promise a 50k chunk fits.

This is the HIGH-STATISTICS NATIVE GEANT4 CAMPAIGN, not a completed
full six-field analysis. After Geant4 ROOT truth has been produced,
the existing event-to-muon gas-hit TGeo exporter needs an explicitly
chunk-aware workflow followed by a Heed + frozen 6 Garfield LUT fold.
Only after those results are merged into original-layout histograms,
with source/event/chunk provenance and parent-bootstrap statistics,
may dashboards be labelled 100k/source or 2M/source.


### 2026-10-21: ROOT/hadd missing from PATH during 2M benchmark launch

The real WSL benchmark reported:
- 2M GENERATED parents/source target, 50k/chunk, plan PASS;
  865.52 GiB free disk at startup.
- Synthetic capped-run regression PASS; its SKIP completed Jpsi/Phi/
  Omega lines refer to dummy /tmp fixtures, NOT actual ROOT checkpoints.
- Native Stage-C simulator prefix detected correctly under
  /home/stanm/na6p/install/NA6PRoot-mwpc.
- REAL campaign setup STOPPED before any Geant4 event was generated:
  CAMPAIGN_ERROR: Required hadd missing from PATH.

ROOT is installed independently of NA6PRoot. The Stage-C environment
activator now chooses the installed ROOT prefix, prioritizing explicit
NA6P_ROOT_PREFIX, then ROOTSYS, then /home/stanm/na6p/install/root,
then root-config. It demands executable bin/root and bin/hadd plus
libCore, libRIO and libTree from the SAME installation, sets PATH,
ROOTSYS, LD_LIBRARY_PATH and ROOT_INCLUDE_PATH and checks linked
libraries using ldd. No ROOT/Geant4 rebuild warranted just for PATH.

Run only the read-only executable/environment preflight first:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh env-check

Only after env-check PASS, retry:

    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh benchmark

The env-check never creates a campaign manifest or starts Geant4.
The real benchmark preserves the frozen 2M/source identity and all
completed 50k chunks. Review measured ROOT disk and CPU wall time
before FULL production.


### 2026-10-21: full native 100k per-channel Geant4 benchmark PASSED

User's local WSL native production finished six genuine 50k generated
parent-event chunks, two per Jpsi/Phi/Omega, all with validated ROOT.
Observed six-chunk campaign elapsed 21.3 minutes. Measured per-event
GEANT4 time: Jpsi 0.0048374 seconds, Phi 0.00416779 seconds,
Omega 0.00371075 seconds. NOT a measurement of Heed/Garfield.
Combined 6M parent-events Geant4 would take about 7–7.5 hours
IF no slowdown; do not present as guaranteed. Startup free disk
865.52 GiB, but remaining disk after new ROOT checkpoints unknown.

The old six-geometry native fold uses PREVIOUS campaign first 3k
parents/source and hardcodes chunk0001 in its Heed RNG seeds. It
MUST NOT process two new 50k chunks as if the same source/event pair
were unique across chunks or as if high-stat PDFs already existed.

A separate TGeo gas-crossing export stage has been added for the
actual new 100k/source Geant4 ROOT files WITHOUT new Geant4/Garfield:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_2m_generated_parent_campaign.sh status
    bash test/mwpc-stage-c/run_pi_2m_native_export.sh preflight
    bash test/mwpc-stage-c/run_pi_2m_native_export.sh benchmark

PRE-FLIGHT reads six real completed.json checkpoints and compares
unchanged ROOT input SHA256, verified gas/reference trees, source,
chunk index, generated event count and RNG seed (50k/source/chunk).
BENCHMARK invokes the existing ROOT TGeo/MCKine joined native muon
gas-crossing exporter on each of the six chunks, with local event
indexes verified in [0,50000). Separate CSV and immutable manifest
are saved under:

    test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1/
      jpsi/chunk0001/geant_muon_crossings.csv
      jpsi/chunk0001/crossings_manifest.json
      ... phi and omega chunk0001 and chunk0002 ...

The new outputs never overwrite the existing first 3k/source
pi_six_geometry_native_timing_v1 or Garfield LUT maps.

This is ONLY the native Geant4 crossing export, NOT YET larger-stat
six-configuration Heed/LUT timing or arrival-angle histograms. To
finish the physics deliverable, fold each (source,chunk,local_event)
with shared Heed births and six existing Garfield LUTs using
chunk-aware RNG seeds, preserve 10ns bins and chamber/anode
conservation, merge with parent IDs and compute bootstrap statistics.
Do not claim 100k/2M drift PDF or bootstrap until those stages pass.


### 2026-10-23: verified 100k/channel real gas-crossing exports + six LUT fold

WSL user confirmed all six true native ROOT-to-TGeo exports passed with
SHA-verified ROOT checkpoints and independent per-chunk manifests:

  Jpsi chunk0001 681951 gas intersections, 494627 muon, 346.1 MiB
  Jpsi chunk0002 685118 gas intersections, 495068 muon, 347.7 MiB
  Phi  chunk0001 516974 gas intersections, 391297 muon, 263.3 MiB
  Phi  chunk0002 517062 gas intersections, 392164 muon, 263.3 MiB
  Omega chunk0001 470225 gas intersections, 359247 muon, 239.5 MiB
  Omega chunk0002 475192 gas intersections, 361307 muon, 242.1 MiB

  TOTAL 3,346,522 Geant4 gas intersections, 2,493,710 muon crossings,
        100,000 GENERATED PARENTS/SOURCE, ~1,702 MiB CSV.

The previous 3k-parent/source six-field code originally held all
groups and the duplicate source=all groups in memory, and Heed RNG
seeded by constant chunk0001. For the new 100k/source first stage, a
separate station-sharded, source-and-chunk-aware fold has been added:

  test/mwpc-stage-c/pi_2m_native_six_geometry_fold.py
  test/mwpc-stage-c/test_pi_2m_native_six_geometry_fold.py
  test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh

It uses existing measured native Geant4 CSV, independent frozen six
Garfield LUTs and original empirically parameterized Heed ionization.
It does NOT rerun microscopic Garfield or Geant4. Different chunks
have different shared-Heed random seeds even when local event indices
are equal. Within ONE event and ONE gas hit, all six physical fields
receive the SAME sampled primary electron positions and quantiles.
Each bounded-RAM shard is for one source, one source-specific 50k
Geant4 chunk and one native station, including selected chamber,
anode and signed entrance-angle groups. No copied source=all groups.

Each shard publishes compressed exact bins, known/unknown collection,
edge counters, approach/impact angles, and raw sum(t), sum(t^2)
and other moments (to combine centered RMS EXACTLY across chunks,
not average RMS). Every shard has a separate SHA and a frozen manifest
of SOURCE+CHUNK+STATION, both Heed libraries, six Garfield map SHAs,
input CSV SHA and source code. It fails closed on mismatches.
Full 100k per source = 2 x 50k x 3 sources x 6 stations = 36 shards,
then merges 2 ROOT chunks at a time for each station/source into
NEW 18 V2-style clickable map/time pages with 10 ns histograms, log Y,
angles, mean and centered RMS, using the SAME physical active area.
18 separate detailed per-anode pages also retained in the NEW output.

Start with the SMALL genuine benchmark, NOT all 2.49m crossings:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh smoke

This reads exactly the FIRST 200 real muon gas hits from Jpsi
chunk0001 MS0, checks group conservation, verifies the Heed/Garfield
LUT provenance and reports REAL wall seconds. The synthetic regression
runs before the real 200-hit fold. No new Geant4/Garfield computation.
After a PASS, measure one COMPLETE Jpsi MS0 50k-parent chunk:

    python3 test/mwpc-stage-c/pi_2m_native_six_geometry_fold.py fold \
      --channel Jpsi --station 0 --chunks 1

That output is resumable; recorded wall time gives a measured folding
throughput BEFORE committing to all 36 real native station shards.
If memory/performance is acceptable, run:

    bash test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh fold
    bash test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh merge

New HTML portal, distinct from the old 3k-parent dashboard:

    test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/
      pages/pi_pitch_fieldwire_timing_dashboard_v3.html

No actual 100k-source Heed/LUT fold has yet been observed on WSL,
so these instructions ARE NOT a claim that high-statistics time PDFs
or event-level errors have been completed. The PDF and mean/RMS
currently omit parent-event bootstrap errors; the native Heed/Garfield
LUT trials also carry separate finite-sampling uncertainty.
Do not run FULL 2M/channel Geant4 extension until representative
Heed/LUT folding wall time, memory, and intermediate file sizes
have been measured and reviewed.


### 2026-10-23: matched physical geometry/binomial vs 6 Garfield LUT timing study

The user has now measured a COMPLETE first source/station 50k-parent
Heed+six-frozen-Garfield fold:
  Jpsi, chunk0001, MS0: 92,017 actual native muon gas hits;
  18,738 grouped response histograms; 4.5 MiB compressed;
  wall time 158.0 seconds, checkpoint PASS.
This is distinct from the upstream 50k-parent Geant4 elapsed time
(~4 minutes for this 50k Jpsi chunk) and from ORIGINAL microscopic
Garfield map construction.

The PI also requested a first-order geometric/control model and
F-A-I-R wall-clock comparison with the microscopic Garfield approach.
Existing original October 2026 prototype:
  c_pi_geant_wire_split.py
  analyze_c_pi_geometric_split.py
  fold_c_pi_geometric_clusters_with_transport.py
This prototype operated on the previous geometry and should not be
mistaken for the current symmetric 3+3 mm six-field comparison.

NEW on the CURRENT actual native Geant4 benchmark data:
  pi_geometric_vs_lut_benchmark.py
  test_pi_geometric_vs_lut_benchmark.py
  run_pi_geometric_vs_lut_benchmark.sh

It samples genuine station muon gas segments uniformly by deterministic
reservoir, preserving source/chunk/local_event/track identity. Each hit
uses the ORIGINAL seed and RNG sequence from the 100k/source six-field
Heed fold (proven byte-identical for primary positions and quantiles by
test), draws a shared Ncl~Poisson(lambda*L), and independently marked
cluster multiplicities from empirical Heed. These exact cluster births
are shared by BOTH geometry and LUT methods.

GEOMETRY: assign each COMPLETE primary-electron CLUSTER to the closest
finite physical anode. Exact track arc-length fractions through
nearest-anode Voronoi cells determine binomial/multinomial cluster
counts. For exactly two cells with p=L1/L:
  N1|Ncl ~ Binomial(Ncl,p);
  P(both cells occupied | Ncl)=1-p**Ncl-(1-p)**Ncl.
With Poisson marginalization:
  P(both)= (1-exp(-lambda*L1))*(1-exp(-lambda*L2)).
No geometric loss, avalanches, electronic thresholds or intrinsic
drift timing are asserted. Unlike LUT, the two field-wire variants at
the SAME pitch have identical geometric nearest-wire predictions.

LUT: for the SAME individual Heed primary electrons and same q, call
the original six independently computed Garfield lookup-table kernels.
Real finite-wire collection, migration and noncollection are allowed.
Compare for each configuration:
  P(Ndistinct_anodes>=2) with thresholds >=1,2,3,5 collected primary
  electrons PER anode, LUT collection and migration fractions, and
  conditional on native gas hits crossing exactly two geometric cells
  analytic binomial and Poisson expectations vs sampled Geom and LUT.

Measured wall time is split rather than conflated:
  - source native CSV disk read and reproducible selection,
  - LUT load from disk,
  - SHARED Poisson/Heed cluster generation (identical both models),
  - GEOMETRY path split and per-anode electron count, PER config,
  - GARFIELD LUT lookup and per-anode electron count, PER config,
  - original microscopic Garfield++ LUT creation, IF verified per-column
    wall-time sidecars survived from old pilot,
  - original upstream 50k Geant4 elapsed from completed.json.
This never claims a new direct microscopic Garfield++ simulation was
performed, or calls the 6 LUT fold itself an independent full
detector physics simulator.

FIRST confirm synthetic physics tests plus 200 real crossing benchmark:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_geometric_vs_lut_benchmark.sh smoke

Then with 2,000 uniformly sampled real Jpsi MS0 crossing segments:

    bash test/mwpc-stage-c/run_pi_geometric_vs_lut_benchmark.sh benchmark

Native result JSON (English keys) appears at:
  test_runs/mwpc_stage_c/pi_geometric_vs_lut_benchmark_v1/
    jpsi_chunk0001_MS0_sample2000.json

Separately COMPLETE six-field Heed+LUT stage already has a validated
Jpsi MS0 chunk0001 checkpoint. Resume other 35 source-station-chunk
physics shards without recomputing that checkpoint:

    bash test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh fold
    bash test/mwpc-stage-c/run_pi_2m_heed_six_geometry_fold.sh merge

The current final merge still requires all 36 native shards before
any 100k GENERATED parents/source high-stat V3 pages can be labelled
complete. English localisation of V3 and PI slide deck are deferred
until full statistical comparisons and physics claims are verified.


### 2026-10-23: user verified all six-source/geometry native 100k pages

User reports the 36 independent Heed+six Garfield response shards AND
the strict two-chunk/source MERGE all completed and the high-stat
Russian V3 pages opened and worked correctly in Windows browser.
This is genuine 100,000 GENERATED parent events for EACH of the 3
sources, not 2,000,000 yet. Event-level parent bootstrap errors still
have NOT been calculated. No new 2M/channel transport should be
initiated until further scaling and scientific review.

English-localisation stage (read-only on physics):
  test/mwpc-stage-c/build_pi_100k_english_dashboard.py
  test/mwpc-stage-c/test_pi_100k_english_dashboard.py
  test/mwpc-stage-c/run_pi_2m_english_dashboard.sh

Run from the existing checkout (DO NOT rerun Geant4, Heed, Garfield,
ROOT or the expensive native fold, and do not need to re-merge):

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_2m_english_dashboard.sh check
    bash test/mwpc-stage-c/run_pi_2m_english_dashboard.sh build

The check mode first runs a synthetic fixture for English portal,
interactive V3, and detailed per-anode dashboard with Python
no-JavaScript static fallback, then validates the 36 already merged
Russian HTMLs and their true 100k parent exposure without writing
translated files. Build mode writes English copies next to all 36
original HTMLs without overwriting any physics JSON or original
Russian pages. It localises ALL dynamic JavaScript selection text,
labels, legends, plot axes, warnings, event-exposure notices,
detailed anode spectra, clickable station navigation, and the
static Python fallback.
The embedded <script type=application/json id=payload/data> physics
payload must match EXACTLY byte-for-byte between original and English.
Every output must contain no Cyrillic outside physics JSON and must
retain six configs, logarithmic timing mode, and angular observables.
An English page refusing an untranslated Russian term is a deliberate
QA FAIL: adjust localisation, DO NOT change the numeric physics data.

The English entrypoint is:

    test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/pages/
      pi_pitch_fieldwire_timing_dashboard_v3_en.html

Open in Windows from WSL:

    FILE="$HOME/na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/pages/pi_pitch_fieldwire_timing_dashboard_v3_en.html"
    explorer.exe "$(wslpath -w "$FILE")"

A detailed physics slide deck should wait for parent-level uncertainty
analysis and controlled normalisation/efficiency comparisons.


### 2026-10-23: scientific interpretation of the merged 100k anode timing

The user confirmed that both native six-config 100k/parent/source merge
and English V3 navigation work in their WSL/Windows browser. They now
request a scientific analysis of WHICH ANODES collect at WHICH electron
drift times, pitch and grounded-field-wire effects, 10 ns distributions,
and possible finite-edge effects.

Real merged station JSON files are user-local and NOT committed to
GitHub, so numerical station-level effects cannot be asserted from
source-code inspection alone. A read-only scientific summary/interactive
English dashboard has been added:

    test/mwpc-stage-c/analyze_pi_100k_anode_timing.py
    test/mwpc-stage-c/test_pi_100k_anode_timing.py
    test/mwpc-stage-c/run_pi_100k_anode_timing_analysis.sh

Run in WSL:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_100k_anode_timing_analysis.sh

This does NOT run Geant4/Heed/Garfield++/ROOT, nor re-fold six LUTs.
It streams all 18 existing verified source/station merged JSONs with
bounded memory and checks:
  - source, station, six pinned configs, 100k parents/source
  - real gas-crossing exposure, counts and all 31 time bins
  - station=chambers=anodes hist, known+unknown electron conservation
  - unrounded sum(t)/sum(t²) moment consistency
  - center vs outer anode time differences as electron-weighted means
  - per-anode conditional time histogram profiles normalised to each
    pitch's native physical anode lattice

Computes per source/station/config:
  - total gas hits (NOT parent count)
  - primary Heed, anode collected and known drift time electron counts
  - anode collection probability / ALL born primary electrons
  - P(e collected with t<100 ns) / ALL born primary electrons
  - P(t<100 ns|collected, known drift time)
  - mean and centered RMS from raw exact moments
  - >=300 ns known-electron tail fraction
  - mean timing and F100 vs normalized physical anode location
  - station- and source-pooled field-wire minus no-wire differences.

Physical caution: the per-anode groups have NO denominator of born
primary electrons assigned to their physical anode, so anode-by-anode
collection efficiency is undefined from these data. Do NOT equate
timing differences (conditioned on collection) to true electronics
timing resolution. Do NOT call a multiwire fraction a result of the
100k timing merge; multiwire benchmark was a separate 2000
gas-hit comparison. Do NOT compare same integer wire IDs at
different anode pitches (their physical positions change). No
parent-event confidence intervals for the current 100k analysis yet;
finite microscopic Garfield LUT trials are an additional uncertainty.

Read-only new files:
    test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/
      anode_timing_analysis_en/
        summary.txt                  # plain-text science summary
        station_summary.csv          # all 18×6 station/config metrics
        station_anode_timing.json    # compact 18×6×20 anode-position bins
        anode_timing_report_en.html  # six-curve time PDF and spatial profiles

Open independent English analysis:

    explorer.exe "$(wslpath -w "$HOME/na6p/src/NA6PRoot-mwpc/test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/anode_timing_report_en.html")"

User can share summary.txt or station_summary.csv to enable a
rigorous numerical interpretation and discussion of physical tradeoffs
based on actual measured results. A separate parent-event bootstrap
and independent finite-LUT Monte Carlo systematic study should follow.


### 2026-10-23: mandatory per-muon truth retention BEFORE 6M launch

The user wants to avoid a future full Geant4 rerun if prototype readout,
anode layout, Garfield drift, event reconstruction or electronics change.
The canonical detailed retention specification now lives at:
  test/mwpc-stage-c/PI_2M_MUON_DATA_RETENTION.md

DO NOT press start before:
  bash test/mwpc-stage-c/run_pi_weekend_2m.sh check

CHECK now calls:
  run_pi_2m_truth_preservation_check.sh
  audit_pi_muon_truth_bank.C
  freeze_pi_2m_transport_inputs.py

These scripts verify six original 50k-parent ROOT benchmarks with
event-indexed TParticle/hits, inspect 200 spread events per chunk
for track associations, require original raw worker ROOT, positive
Edep gas steps, virtual plane trees, source+chunk ROOT-to-CSV SHA
consistency; they also archive actual original field/cuts/source and
simulator-binary file bytes with SHA and immutable snapshot manifest.
This is only a SAMPLE per-object audit, not a complete simulation proof.
Every future 50k chunk still passes the original complete ROOT tree
and SHA validator. Do NOT change frozen Geant4 input source/binaries
to add diagnostics now: that would conflict with successful 100k chunks.

Stored input is 3 forced dimuon types, NOT realistic minimum-bias
heavy-ion backgrounds or beam-time occupancy/pileup. MCKine remapped
track IDs and raw worker original step IDs must not be joined blindly.
Reference planes are virtual interpolations, not exact sensitive hits.
The 20-normalized-position-bin station_anode_timing.json uploaded by
user is a summary, NOT per-muon event data. Original ROOT and native
CSV are the reusable track bank.

After check PASS:
  bash test/mwpc-stage-c/run_pi_weekend_2m.sh start
  bash test/mwpc-stage-c/run_pi_weekend_2m.sh status

Three production stages checkpoint independently; however the current
40-chunk full merge and parent-event bootstraps are still downstream
development. Always distinguish user-local simulation clock from
geometry vs precomputed six-Garfield-LUT postprocessing timing. Separate
the precomputed LUT creation clock, Heed, Geant4, native ROOT export and
folding costs in eventual study.


### 2026-10-24: 6M campaign resumes + connected Google Drive backup

User reports 120 verified native Geant4 CSV exports COMPLETE, all 120
Geant4 50k-parent chunks complete, and a later Heed/LUT fold failure:
Jpsi chunk0009 MS3, event 4040, hit index 3, track 2, detector 110.
The independent frozen CSV audit finds exactly one sub-eV zero-length
chord with identical gas entry/exit points, Edep 0.0465400322447 eV,
out of 78,928 muon crossings; 78,927 were valid.

Do NOT modify the original Heed sampler or fold driver: their SHAs
are frozen in all earlier manifests. The new compatibility wrapper
pi_2m_zero_path_compatible_resume.py keeps such an exact zero chord
in gas-hit exposure but samples Poisson(lambda*0)=0 Heed electrons,
and writes a separate per-shard exception ledger. Any other invalid
path is rejected by the original code. A separate regression test and
one-pass audit_pi_2m_all_native_paths.py on all 120 verified exports
are now mandatory in run_pi_weekend_2m.sh check/start, so all future
rare native path anomalies are identified before resuming long folds.
The main work mode calls this wrapper, and phase 4 runs paired
Poisson/binomial geometry-vs-six-frozen-Garfield-LUT timing benchmarks
on matching Heed births for each source/station.

A connected Google Drive project folder was created with ID
1Svbp3E230YmVk53D0M1AHwTGo6cBxy38 and title NA60plus_MWPC_Research.
URL: https://drive.google.com/drive/folders/1Svbp3E230YmVk53D0M1AHwTGo6cBxy38
The user must authorize WSL rclone separately (remote na60drive,
advanced root_folder_id set to that ID) to upload local WSL results.
Script run_pi_drive_archive.sh supports check/copy/verify; copy is
idempotent and non-deleting; verify checks local vs Drive checksums.
Cloud is an offsite archive and a connector-accessible reporting surface,
NOT the live native ROOT/Geant4 working filesystem. Preserve all
original MCKine/Hit/ReferencePlane ROOT and raw worker data locally.

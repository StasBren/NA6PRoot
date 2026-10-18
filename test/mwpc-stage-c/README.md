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

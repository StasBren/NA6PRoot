# MWPC Stage-C local data inventory

This file is the **canonical inventory of large/local simulation inputs and
runtime outputs** used by the current MWPC Stage-C / PI study.

The files listed here are intentionally **not committed to Git**.  The purpose
of this document is to make their locations and provenance easy to recover
from GitHub.

> Local checkout assumed below:
>
> ```text
> $HOME/na6p/src/NA6PRoot-mwpc
> ```
>
> If the checkout is moved, paths below `$REPO/test_runs/...` move with it.

---

## Path aliases

For the commands and paths below:

```bash
REPO="$HOME/na6p/src/NA6PRoot-mwpc"

BORE_BASE="$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/mnp33-bore-runs/20260912-202112-0wrxhdg_"

STAGEC="$REPO/test_runs/mwpc_stage_c"
```

---

## 1. Large historical Geant4 production: three physics sources

Root directory:

```text
$BORE_BASE/
```

Absolute path:

```text
$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/
  mnp33-bore-runs/20260912-202112-0wrxhdg_/
```

This is the large MNP33-bore production made in the separate
`na60-dice-mwpc-simulation` study.

Sources:

```text
J/psi
Phi
Omega
```

Statistics:

```text
200,000 parent dimuon events per source
= two independent 100,000-event Geant4 chunks per source
```

So the six production chunks are:

```text
$BORE_BASE/jpsi40-bore-chunk001/
$BORE_BASE/jpsi40-bore-chunk002/

$BORE_BASE/phi40-bore-chunk001/
$BORE_BASE/phi40-bore-chunk002/

$BORE_BASE/omega40-bore-chunk001/
$BORE_BASE/omega40-bore-chunk002/
```

Important persistent ROOT files in every chunk:

```text
LocalTrajectoryAudit.root
MCKine.root
```

For example:

```text
$BORE_BASE/jpsi40-bore-chunk001/LocalTrajectoryAudit.root
$BORE_BASE/jpsi40-bore-chunk001/MCKine.root
```

### LocalTrajectoryAudit content relevant to Stage C

`LocalTrajectoryAudit.root` contains a `Crossings` tree with transported
track states at fixed reference planes.  It stores quantities such as:

```text
event, track, originalTrack, parent, pdg, plane,
x_cm, y_cm, z_cm,
px_GeV, py_GeV, pz_GeV, p_GeV
```

It does **not** store the actual gas-volume entry/exit pair or gas
energy deposition for the historical run.

Observed reference planes:

```text
plane 0 : z = 300 cm   MS0
plane 1 : z = 340 cm   MS1
plane 2 : z = 365 cm   MNP33 reference
plane 3 : z = 430 cm   MNP33 reference
plane 4 : z = 495 cm   MNP33 reference
plane 5 : z = 550 cm   MS2
plane 6 : z = 590 cm   MS3
plane 7 : z = 810 cm   MS4
plane 8 : z = 850 cm   MS5
```

For the six MWPC stations the historical Stage-C sample uses planes:

```text
0, 1, 5, 6, 7, 8
```

---

## 2. Historical crossing CSVs derived from the big production

Already-produced per-source crossing tables:

```text
$BORE_BASE/analysis/jpsi40-bore/bore_muon_crossings.csv
$BORE_BASE/analysis/phi40-bore/bore_muon_crossings.csv
$BORE_BASE/analysis/omega40-bore/bore_muon_crossings.csv
```

These contain:

```text
job,event,track,pdg,station,
x_cm,y_cm,z_cm,
px_GeV,py_GeV,pz_GeV,p_GeV,
theta_x_deg,theta_y_deg
```

Stage-C mapping used for the PI study:

```text
global Y -> local u  (across wires)
global Z -> local v  (chamber normal)
global X -> local w  (along wires, with installed-sign conventions handled
                       where needed)

theta_u = atan(py/pz)
theta_w = atan(px/pz)
```

---

## 3. Compact 36k realistic GEANT track sample

Prepared from the three sources above with 2000 tracks per
(source, station) combination:

```text
3 sources x 6 stations x 2000 = 36,000 tracks
```

Files:

```text
$STAGEC/c_pi_tracks/c_pi_track_sample.csv
$STAGEC/c_pi_tracks/c_pi_track_distribution_summary.csv
```

The sample is deliberately balanced by source and station.  It is therefore a
response-study sample, not a naturally weighted physics-rate sample.

---

## 4. Heed ionisation reference library for the PI 2+4 mm geometry

Current PI geometry used in this branch:

```text
gas gap          : v = -2 ... +4 mm
total gas gap    : 6 mm
wire pitch       : 4 mm
wire diameter    : 30 um
HV               : 1800 V
gas              : Ar/CO2 70:30
B                : 0 T baseline
```

Heed library:

```text
$STAGEC/c_heed_library_2p4/heed_2p4_tracks.csv
$STAGEC/c_heed_library_2p4/heed_2p4_clusters.csv
```

Reference statistics from the 20k-track Heed sample:

```text
tracks                         : 20,000
clusters                       : 438,898
mean clusters / sampled track  : 21.9449
mean primary e- / track        : 54.424
mean primary e- / cluster      : 2.4800
cluster density used by PI     : 3.719475 / mm
```

The empirical Heed cluster-electron multiplicity distribution is used rather
than a Gaussian approximation.

---

## 5. Baseline single-electron Garfield transport box map

One 4-mm wire cell is divided into:

```text
u = -2 ... +2 mm
v = -2 ... +4 mm
du = dv = 0.2 mm

20 x 30 = 600 boxes
```

Baseline statistics:

```text
50 microscopic electrons / box
600 x 50 = 30,000 Garfield electron drifts
```

Raw map:

```text
$STAGEC/c_transport_box_map_gm2p0_gp4p0/c_transport_box_electrons.csv
```

Box summary:

```text
$STAGEC/c_transport_box_map_gm2p0_gp4p0/c_transport_box_summary.csv
```

Column checkpoints and run log:

```text
$STAGEC/c_transport_box_map_gm2p0_gp4p0/columns/
$STAGEC/c_transport_box_map_gm2p0_gp4p0/scan_progress.log
```

Each raw electron trial stores the starting box and microscopic response,
including destination wire, collection/loss, drift time, `phi_hit`,
hit side, last microscopic direction, and incidence angle.

---

## 6. Earlier targeted high-stat timing study

Selected u-columns:

```text
u = -1.9, +0.1, +1.5, +1.9 mm
```

Current existing timing run used:

```text
500 electrons / box
4 u-columns x 30 v-bins x 500 = 60,000 microscopic drifts
```

Location:

```text
$STAGEC/c_transport_timing_highstat_2p4/
```

Important files:

```text
$STAGEC/c_transport_timing_highstat_2p4/c_transport_timing_highstat_electrons.csv
$STAGEC/c_transport_timing_highstat_2p4/c_transport_timing_highstat_boxes.csv
$STAGEC/c_transport_timing_highstat_2p4/c_transport_timing_highstat_box_wires.csv
$STAGEC/c_transport_timing_highstat_2p4/c_transport_timing_highstat_representative.csv
```

---

## 7. Current targeted high-stat boundary run

Purpose: improve the Garfield statistics exactly where the midpoint
approximation and 100-ns timing are most sensitive.

Selected u-columns:

```text
u = -1.9, -1.7, +1.7, +1.9 mm
```

Default statistics:

```text
2000 electrons / box
4 u-columns x 30 v-bins x 2000
= 240,000 microscopic Garfield drifts
```

Run directory:

```text
$STAGEC/c_transport_boundary_highstat_2p4/
```

Column checkpoints:

```text
$STAGEC/c_transport_boundary_highstat_2p4/columns/
```

Combined targeted map, after the run completes:

```text
$STAGEC/c_transport_boundary_highstat_2p4/c_transport_boundary_highstat_electrons.csv
```

Timing summaries:

```text
$STAGEC/c_transport_boundary_highstat_2p4/c_transport_boundary_highstat_boxes.csv
$STAGEC/c_transport_boundary_highstat_2p4/c_transport_boundary_highstat_box_wires.csv
$STAGEC/c_transport_boundary_highstat_2p4/c_transport_boundary_highstat_representative.csv
```

### Hybrid full transport map

After the targeted run completes, its 120 high-stat boundary boxes replace the
same boxes in the old full 50-electron map.

Hybrid map:

```text
$STAGEC/c_transport_box_map_hybrid_boundary_2p4/c_transport_box_electrons.csv
```

Expected statistics in that hybrid map:

```text
480 ordinary boxes  : 50 trials / box
120 boundary boxes  : 2000 trials / box
```

---

## 8. PI nominal gas segments reconstructed from historical GEANT tracks

The old production does not contain true gas entry/exit hits, so for the
historical response study a nominal 6-mm gas segment is extrapolated around the
saved station reference state.

File:

```text
$STAGEC/c_pi_geometric_split/c_pi_nominal_segments.csv
```

Important limitation:

```text
This is a nominal intrinsic 6-mm segment reconstruction.
It is not the exact staggered current chamber crossing and it has no historical
event-by-event gas Edep.
```

For a new Geant4 production with modern `NA6PMuonSpecHit`, use the real
entry/exit and accumulated gas energy loss instead.

---

## 9. PI geometric Poisson/binomial + Heed outputs

Current geometric PI model:

```text
realistic GEANT track
  -> split path length between 4-mm wire cells
  -> Ncl ~ Poisson(lambda * L)
  -> multinomial/binomial cluster assignment by path fraction
  -> empirical Heed primary-electron multiplicity per cluster
```

Outputs:

```text
$STAGEC/c_pi_geometric_split/c_pi_wire_split_crossings.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_split_wires.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_split_clusters.csv
```

The cluster file contains the generated `(u,v,w)` position of every cluster
and its primary-electron multiplicity.  This is the source that is folded with
the Garfield box map.

Geometric-analysis summaries:

```text
$STAGEC/c_pi_geometric_split/c_pi_geometric_split_summary.txt
$STAGEC/c_pi_geometric_split/c_pi_geometric_split_binned.csv
```

---

## 10. PI folding through the baseline Garfield map

Baseline folding output:

```text
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_crossings.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_clusters.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_wires.csv
```

This combines the synthetic ionisation clusters with the finite-statistics
single-electron Garfield response map.

Reference result from the baseline 50-electron/box map:

```text
36,000 crossings
2,015,952 primary electrons
1,628,087 collected electrons
collection fraction = 0.807602
migration fraction  = 0.008992
```

Here "migration" means that a collected electron ends on a wire different from
the hard geometric midpoint assignment.

---

## 11. PI folding through the high-stat hybrid Garfield map

After the current boundary run is complete, use:

```text
$STAGEC/c_transport_box_map_hybrid_boundary_2p4/c_transport_box_electrons.csv
```

and write a separate result set:

```text
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_highstat_crossings.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_highstat_clusters.csv
$STAGEC/c_pi_geometric_split/c_pi_wire_transport_highstat_wires.csv
```

Do not overwrite the baseline folding; keeping both makes the effect of the
high-stat boundary boxes directly checkable.

---

## 12. PI-facing plots

Baseline figures:

```text
$STAGEC/c_pi_geometric_split/pi_plots/
```

High-stat figures, after the hybrid refolding:

```text
$STAGEC/c_pi_geometric_split/pi_plots_highstat/
```

Current PI-facing plot set includes:

```text
pi_01_theta_u_distribution.png
pi_02_theta_w_distribution.png
pi_03_path_length_zoom.png
pi_04_geometric_vs_garfield_sharing.png
pi_05_two_wire_probability.png
pi_06_garfield_migration_vs_u.png
pi_07_timing_100_200_vs_abs_u.png
pi_08_event_wire_multiplicity.png
pi_09_multiwire_threshold_scan.png
pi_plot_summary.txt
```

---

## 13. Modern Geant4 crossing exporter

For new/current NA6PRoot Geant4 runs with modern `NA6PMuonSpecHit`, use:

```text
test/mwpc-stage-c/export_c_geant_muon_crossings.C
```

Unlike the old `LocalTrajectoryAudit` production, the modern hit object can
provide the actual gas crossing information needed by the PI workflow,
including entry/exit positions and accumulated energy deposition.

No single canonical local path for a new full production is recorded here yet.
When a new production is chosen as the reference, add its absolute root
directory to this inventory.

---

## 14. Quick existence check

From the NA6PRoot checkout:

```bash
cd "$HOME/na6p/src/NA6PRoot-mwpc"

BORE_BASE="$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/mnp33-bore-runs/20260912-202112-0wrxhdg_"
STAGEC="$PWD/test_runs/mwpc_stage_c"

ls -lh \
  "$BORE_BASE/analysis/jpsi40-bore/bore_muon_crossings.csv" \
  "$BORE_BASE/analysis/phi40-bore/bore_muon_crossings.csv" \
  "$BORE_BASE/analysis/omega40-bore/bore_muon_crossings.csv" \
  "$STAGEC/c_pi_tracks/c_pi_track_sample.csv" \
  "$STAGEC/c_heed_library_2p4/heed_2p4_clusters.csv" \
  "$STAGEC/c_transport_box_map_gm2p0_gp4p0/c_transport_box_electrons.csv" \
  "$STAGEC/c_pi_geometric_split/c_pi_wire_split_clusters.csv"
```

---

## Maintenance rule

Whenever a simulation becomes a **reference input** for subsequent Stage-C
work, record here:

1. its absolute or `$HOME`-relative local directory;
2. what was simulated;
3. event/trial statistics;
4. which files are the persistent inputs;
5. which downstream analysis depends on them.

This file should be updated before an old local run is renamed, moved, or
deleted.

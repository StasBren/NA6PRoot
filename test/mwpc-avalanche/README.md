# MWPC microscopic response sandbox — Stages A and B

This directory contains the **validated microscopic-response work that Stage C
reuses**.  It is intentionally standalone from the normal NA6PRoot build;
Garfield++ is still not a production dependency of the main detector code.

New event-level integration belongs in `../mwpc-stage-c/`.
Do not create a second coordinate transform, Heed chain, avalanche model, or
strip-response model there unless an existing implementation is first factored
into a shared helper.

## Coordinate convention

The detector convention is

- detector/global **X**: horizontal, along the anode wires;
- detector/global **Y**: vertical, across the wires;
- detector/global **Z**: along the beam / nominal chamber normal.

The chamber-local convention is

- `u`: across wires;
- `v`: chamber normal;
- `w`: along wires.

Garfield uses

```text
(x_G, y_G, z_G) = (u, v, w)
```

For the nominal NA60+/DiCE orientation this is

```text
u = global Y
v = global Z
w = global X
```

The reusable rigid-frame implementation is `MWPCCoordinateFrame.h`.
`phaseA_coordinate_test.cxx` validates point/vector round trips and basis
directions.

## Reference microscopic cell

The current reference operating point is Prototype-3-like rather than a claim
of frozen final chamber geometry:

- Ar/CO2 = 70:30;
- anode-wire diameter = 30 um;
- wire pitch = 4 mm;
- asymmetric cathode gaps = 2 mm + 4 mm;
- anode voltage = +1.8 kV;
- reference strip pitch/width = 1.7 mm / 1.7 mm in the later Stage-B studies.

Most executables expose the relevant geometry/field values as runtime options.

## Stage A — gas, transport and avalanche

Key executables/source files:

- `phaseA_coordinate_test.cxx` — global/local frame validation;
- `phaseA_collection_scan.cxx` — wire collection / separatrix behaviour;
- `phaseA_heed_muon.cxx` — Heed ionisation clusters and conduction electrons;
- `phaseA_heed_avalanche.cxx` — full Heed -> microscopic avalanche chain;
- `phaseA_heed_statistics.cxx` — ionisation/energy-loss statistics;
- `phaseA_gain_fluctuations.cxx` — stochastic single-electron gain;
- `phaseA_magnetic_drift.cxx` — B-field drift sensitivity;
- `phaseA_avalanche_visualization.cxx`, `heed_muon_visualization.cxx` —
  presentation/diagnostic visualisations.

Established qualitative/quantitative behaviour includes sensible symmetric
wire collection at a cell boundary, Poisson-like primary-cluster multiplicity,
long-tailed energy loss, broad avalanche-gain fluctuations, and mainly
along-wire Lorentz displacement for the nominal field orientation.

Absolute gas gain is **not** prototype-calibrated.

## Stage B — induced strip response and reconstruction

Key response pieces:

- `phaseB_strip_weighting.cxx` — strip weighting potential/field;
- `phaseB_single_ion_signal.cxx` — single-ion Shockley-Ramo test;
- `phaseB_electron_ion_pair.cxx` — electron/ion signal decomposition;
- `phaseB_ion_cloud_signal.cxx` and
  `phaseB_electron_ion_cloud_signal.cxx` — compact cloud studies;
- `phaseB_real_avalanche_signal.cxx` — real microscopic avalanche -> strip signal;
- `phaseB_resolution_ensemble.cxx` — repeated-response ensemble;
- `phaseB_strip_geometry_scan.cxx` — pitch/width trade-off;
- `phaseB_tilt_scan.cxx` — tilted strip families;
- `phaseB_opposite_cathode_stereo.cxx` — two-family/opposite-cathode study;
- `phaseB_avalanche_angle_scan.cxx` — avalanche azimuth around the wire.

Associated `analyze_*.py`, `plot_*.py` and `run_*.sh` files are analysis and
presentation utilities for these studies.

Stage B established that a real avalanche gives a prompt electron component plus
a long ion tail, while normalized spatial sharing/CoG stabilises much earlier
than the absolute integrated charge.  Five-strip CoG was already close to the
full exported cluster in the ideal weighting model.

Stage B still excludes final electronics/noise/threshold/shaping and therefore
does not claim final chamber resolution.

## Build

Load the NA6P/ROOT and Garfield environments, then build this standalone
directory:

```bash
source ~/na6p/env.sh
source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh

cd ~/na6p/src/NA6PRoot-mwpc/test/mwpc-avalanche

cmake -S . -B build \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"

cmake --build build -j4
```

The exact target names are listed in `CMakeLists.txt`.

## Stage C handoff

Stage C starts from a full muon crossing and combines the existing pieces:

```text
Geant4 crossing
  -> ChamberFrame global-to-local transform
  -> Heed ionisation
  -> electron drift + many avalanches
  -> Stage-B strip signals
  -> event-level strip cluster / CoG / timing
```

The canonical Stage-C workspace and current C0 acceptance analysis are in:

```text
test/mwpc-stage-c/
```

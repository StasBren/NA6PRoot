# Stage A/B reuse inventory for Stage C

This note records the pieces already established before Stage C so that they
are reused rather than reimplemented.

## Coordinate transform

Canonical implementation:

```text
test/mwpc-avalanche/MWPCCoordinateFrame.h
```

Definitions:

```text
u = across wires
v = chamber normal
w = along wires
Garfield (x,y,z) = (u,v,w)
```

Nominal station mapping:

```text
global Y -> u
global Z -> v
global X -> w
```

The transform supports points and vectors separately and is round-trip tested by
`phaseA_coordinate_test.cxx`.

## Stage A reusable physics

Use the existing code paths for:

- Heed primary ionisation: `phaseA_heed_muon.cxx`;
- Heed + full microscopic avalanche: `phaseA_heed_avalanche.cxx`;
- collection/separatrix checks: `phaseA_collection_scan.cxx`;
- gain fluctuations: `phaseA_gain_fluctuations.cxx`;
- magnetic drift: `phaseA_magnetic_drift.cxx`;
- ionisation statistics: `phaseA_heed_statistics.cxx`.

Reference results from the completed studies:

- 10 GeV muons through the 6 mm reference gas gap:
  `<N_cluster> ~= 21.77`, variance `~= 20.11`;
- energy loss is strongly right-tailed:
  median `~= 1.16 keV`, mean `~= 2.06 keV`, with rare large delta-electron events;
- large-statistics wire-cell boundary collection is symmetric
  (`~50.23%` vs `~49.77%` for the two neighbouring wires);
- the current reference gain model is broad and uncalibrated:
  positive-gain mean `~= 2780 e-`, median `~= 2126 e-`,
  10--90% `~= 491--5887 e-`, and `~15.2%` zero-collected trials;
- for the nominal local frame, the B-field produces mainly along-wire
  displacement, approximately `<Delta w> ~= -0.93 B[T] mm` in the studied setup.

These values are diagnostics of the current sandbox, not final detector
calibration constants.

## Stage B reusable physics

Use the existing code paths for:

- weighting field/potential: `phaseB_strip_weighting.cxx`;
- Shockley-Ramo single-charge checks: `phaseB_single_ion_signal.cxx`;
- prompt electron + ion tail: `phaseB_electron_ion_pair.cxx`;
- compact ion/electron-ion clouds: `phaseB_ion_cloud_signal.cxx`,
  `phaseB_electron_ion_cloud_signal.cxx`;
- real microscopic avalanche -> strip waveform:
  `phaseB_real_avalanche_signal.cxx`;
- response ensemble: `phaseB_resolution_ensemble.cxx`;
- strip geometry: `phaseB_strip_geometry_scan.cxx`;
- tilted/two-family readout: `phaseB_tilt_scan.cxx`,
  `phaseB_opposite_cathode_stereo.cxx`;
- around-wire avalanche azimuth:
  `phaseB_avalanche_angle_scan.cxx`.

Reference Stage-B observations:

- one real avalanche has a ns prompt electron component and a long ion tail;
- in the representative event used for the Stage-B story, the ion fraction of
  integrated signal was about 78% at 25 ns, 94% at 100 ns and 96% at 200 ns;
- normalized three-strip sharing was essentially unchanged between 25 and
  200 ns;
- CoG at 25 ns differed from the 200 ns reference by only a few microns in the
  time-stability test;
- three strips contained about 95.38% of the ideal-model signal and gave
  `~116.2 um` projected spread;
- five strips contained about 99.67% and gave `~72.1 um`, close to the full
  exported cluster at `~71.2 um`.

Again, these are Stage-B model diagnostics before electronics/noise/thresholds,
not a final chamber-resolution claim.

## Stage-C rule

When C1 needs one of the operations above, factor it out of the existing
Stage-A/B executable into a shared helper if necessary.  Do not copy the
physics implementation into a second independent version.

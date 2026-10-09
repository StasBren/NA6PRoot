# C2c across-wire cell scan — first response curve

Reference scan:

- momentum: 5 GeV/c
- theta_u = theta_w = 0
- w0 = 0
- u0 = -1.5, -0.75, 0, +0.75, +1.5 mm
- 5 events per point
- same seed block 130001..130005 reused at every u0
- wire pitch = 4 mm
- numerical Elmer strip weighting
- strip bank k = -4..+4

## Numerical stability

All 25 events pass:

```text
readout_valid               = 1 for 25/25
ion_transport_failures      = 0 for 25/25
nonfinite_weighting_queries = 0 for 25/25
wire_cog_u                  = 0 mm for 25/25
```

The complete scan therefore stays on the central anode wire while the true track
crossing moves through most of the +/-2-mm wire cell.

## Mean raw stereo response

```text
u0 [mm]    <u_reco> [mm]    sample std [mm]
-1.50       -0.14187          0.00829
-0.75       -0.10169          0.00844
 0.00       +0.04330          0.00746
+0.75       +0.15961          0.01569
+1.50       +0.21262          0.00761
```

The response is strongly compressed but it is not flat.  A first linear
description over this restricted range gives

```text
u_reco ~= 0.12937 * u0 + 0.03439 mm
R^2 ~= 0.965
```

so the raw stereo coordinate carries sub-wire information with a response slope
of only about 0.13.  The inverse slope is about 7.7, meaning that a calibration
would amplify both statistical fluctuations and any FEM/reconstruction
systematics by the same factor.

The mean deviations of the five response points from this straight line are at
the 0.01--0.04 mm raw-coordinate level, so a production calibration should not
assume exact linearity.

## Gain-weighted seed truth

At normal incidence, the gain-weighted primary-ionisation coordinate follows
the requested crossing extremely closely:

```text
u0 [mm]    <u_gain> [mm]
-1.50       -1.50036
-0.75       -0.75023
 0.00       -0.00025
+0.75       +0.74952
+1.50       +1.49978
```

Thus the compressed response is not caused by stochastic displacement of the
primary ionisation centroid.  It arises downstream in electrostatic focusing /
avalanche formation / strip induction.

## Interpretation

This scan rejects the simplest "all sub-wire information is lost once every
electron reaches the same wire" picture.  All avalanches are collected on the
central wire, but the two stereo cathode charge patterns still retain a
monotonic dependence on the original crossing position.

At the same time, the dependence is weak: a 3-mm displacement of the true
crossing changes the mean raw stereo coordinate by only about 0.35 mm.

Likely information carriers include:

- the side from which primary electrons approach the anode;
- the microscopic avalanche distribution around the wire;
- the electron-drift induced component;
- the ion cloud's azimuthal distribution around the wire.

A 3-D event visualization is therefore especially useful now: it can show
whether the avalanche birth cloud visibly migrates around the wire as u0 is
changed.

## Next physics scan

Before assigning a calibrated SHR to this coordinate:

1. increase u0 sampling density across nearly the full +/-2-mm wire cell;
2. include points approaching the wire-assignment boundary;
3. increase statistics per point;
4. fit a monotonic calibration curve rather than forcing a straight line;
5. repeat selected angles only after the normal-incidence transfer curve is
   established.

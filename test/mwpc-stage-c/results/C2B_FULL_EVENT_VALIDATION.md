# C2b full-event numerical weighting — first validation event

Reference event:

- random seed: `120001`
- particle: `mu-`
- momentum: 5 GeV/c
- truth at wire plane: `u0 = 1.0 mm`, `w0 = 0.0 mm`
- theta_u = 10 deg, theta_w = 25 deg
- symmetric gap: 2.5 + 2.5 mm
- wire pitch / diameter: 4.0 mm / 30 um
- strip pitch / width: 1.7 / 1.7 mm
- tan(alpha) = 0.1
- numerical readout: translated central-strip Elmer map
- strip bank: k = -4..+4
- observation time: 100 ns
- all Heed electrons processed

## Microscopic chain

```text
clusters                    = 29
primary electrons           = 47
processed seeds             = 47
zero avalanches             = 5
avalanche electrons         = 177048
avalanche ions              = 177022
recorded ion births         = 177022
collected electrons         = 177107
active wires                = 1
ion transport failures      = 0
late electron segments      = 0
non-finite weighting queries= 0
readout valid               = 1
```

The microscopic event is the same reference event used in the C1 hybrid
comparison.

Gain-weighted Heed truth:

```text
<u_seed>_gain = +1.00638 mm
<v_seed>_gain = +0.0369995 mm
<w_seed>_gain = +0.0179259 mm
```

All collected avalanche charge lands on the central anode wire
(`wire_u = 0`).

## Direct numerical strip charges

No hybrid rescaling is applied:

```text
minus scale = 1
plus scale  = 1

Q_cathode_minus = 3.95895 fC
Q_cathode_plus  = 3.88069 fC

sum numerical minus strips = 3.83009 fC
sum numerical plus strips  = 3.76763 fC

minus coverage = 0.967451
plus coverage  = 0.970866
```

The 9-strip numerical bank therefore contains about 97% of the corresponding
complete-cathode induced charge.

Every exported strip charge is positive for this event.  The pathological
bipolar strip pattern and O(100) normalization factors of the full-gap hybrid
model are gone.

## Reconstructed coordinates

Track truth at v=0:

```text
u0 = 1.00000 mm
w0 = 0.00000 mm
```

Numerical stereo CoG:

```text
x_plus  = -0.0670667 mm
x_minus = +0.0521916 mm

u_stereo = +0.599266 mm
w_stereo = -0.00747466 mm
```

Residual to the requested track intercept:

```text
delta u = -0.400734 mm
delta w = -0.00747466 mm
```

Residual to the avalanche-gain-weighted Heed seed position:

```text
delta u_gain = -0.407117 mm
delta w_gain = -0.0254006 mm
```

The important detector-axis interpretation is:

- local `w` = global X = along the horizontal anode wires = across the
  vertical MNP33 magnetic field;
- local `u` = global Y = across the wires = along the magnetic field.

Thus `w` is the high-precision bending coordinate, while `u` is the looser
second coordinate.

The current design requirements are approximately 0.1--0.2 mm in the
high-precision coordinate and 0.5--1 mm in the other coordinate.  A single
event is not a resolution measurement, but this event is on the expected scale:
the numerical model gives a few-10-um w residual and a 0.4-mm u residual.

## Strip-bank convergence from the same event

Recomputing the CoG from nested subsets of the exported 9-strip bank gives:

```text
half-width      u_reco [mm]    w_reco [mm]
+/-1 strips       0.40484       -0.00124
+/-2 strips       0.55394       -0.00474
+/-3 strips       0.58954       -0.00700
+/-4 strips       0.59927       -0.00747
```

The change from 7 to 9 strips is only about 0.010 mm in u and 0.0005 mm in w.
Therefore the remaining ~0.4-mm u offset is not primarily caused by truncating
the strip bank at k=+/-4.

## Interpretation

This is the first full C-stage event in which the strip response is physically
well behaved without an empirical cathode normalization.  It strongly supports
the numerical weighting-field approach.

Do not interpret the single-event residuals as detector resolution.  The next
physics validation should be an event ensemble and a within-wire position scan
to measure bias, stochastic width and possible non-linearity of the raw stereo
CoG estimator.

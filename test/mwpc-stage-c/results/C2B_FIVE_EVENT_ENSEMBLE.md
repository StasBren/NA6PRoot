# C2b five-event numerical-weighting ensemble

Reference point:

- particle: `mu-`
- momentum: 5 GeV/c
- `u0 = 1 mm`, `w0 = 0 mm` at the wire plane
- `theta_u = 10 deg`, `theta_w = 25 deg`
- symmetric 2.5 + 2.5 mm gap
- wire pitch: 4 mm
- strip pitch/width: 1.7 / 1.7 mm
- `tan(alpha)=0.1`
- numerical Elmer weighting, translated central-strip map
- strip bank `k=-4..+4`
- five events, seeds `120001..120005`

## Numerical stability

All five events are valid:

```text
readout_valid                 = 1 for 5/5
ion_transport_failures        = 0 for 5/5
late_electron_segments        = 0 for 5/5
nonfinite_weighting_queries   = 0 for 5/5
active wires                  = 1 for 5/5
collected wire                = central wire u = 0 mm for 5/5
```

No strip charge is negative in any of the 90 exported strip entries.

The nine-strip bank captures a stable fraction of the complete cathode signal:

```text
minus coverage mean = 0.967641, sample std = 0.001074
plus  coverage mean = 0.971880, sample std = 0.007267
```

There is no empirical cathode normalization: `minus_scale=plus_scale=1`.

## Event-by-event response

```text
seed    u_reco     w_reco      u_gain      w_gain      du_gain     dw_gain
        [mm]       [mm]        [mm]        [mm]        [mm]        [mm]
120001  +0.59927   -0.00747    +1.00638    +0.01793    -0.40712    -0.02540
120002  +0.56341   -0.11496    +0.95765    -0.11264    -0.39424    -0.00232
120003  +0.54375   +0.06273    +1.01619    +0.04237    -0.47244    +0.02037
120004  +0.46685   -0.08033    +0.97669    -0.05826    -0.50984    -0.02207
120005  +0.33460   +0.20798    +1.08327    +0.21082    -0.74867    -0.00285
```

Relative to the fixed wire-plane track intercept:

```text
<u_reco-u0> = -0.49843 mm
sample std  =  0.10516 mm

<w_reco-w0> = +0.01359 mm
sample std  =  0.12844 mm
```

The latter spread is not mainly a weighting-field error.  At theta_w = 25 deg,
the stochastic avalanche charge centroid moves along the finite gas thickness.
The gain-weighted Heed position itself has

```text
<w_gain>    = +0.02004 mm
sample std  =  0.12312 mm
```

and the numerical readout follows it closely:

```text
<w_reco-w_gain> = -0.00646 mm
sample std      =  0.01839 mm
correlation(w_reco,w_gain) = 0.990
```

A five-point linear fit gives approximately

```text
w_reco = 1.033 * w_gain - 0.0071 mm
```

so the along-wire stereo response is behaving consistently with the microscopic
charge centroid in this small sample.

## Across-wire coordinate: do not label the 0.5-mm shift as a final detector bias yet

All five events focus onto the same central anode wire at `u_wire=0`, while
the requested track intercept is `u0=1 mm`.  Therefore three different
coordinates must be kept distinct:

1. the geometric track crossing at the wire plane;
2. the pre-avalanche/gain-weighted ionisation coordinate;
3. the avalanche/wire coordinate after electrostatic focusing.

The current chamber documentation explicitly discusses chamber resolution using
wire positions reconstructed from strip measurements.  Therefore the correct
transfer function for the coarse/across-wire coordinate has to be established
with a controlled scan through one wire cell before deciding which residual is
the physically relevant single-hit residual.

For the present five events,

```text
<u_reco> = +0.50157 mm
sample std = 0.10516 mm
```

but this should be treated as one point on the detector response curve, not yet
as a calibrated coordinate measurement.

## Strip-bank convergence

Recomputing the same five events from nested strip subsets gives:

```text
bank        mean u_reco   mean w_reco
            [mm]          [mm]
+/-1        0.33215       0.00508
+/-2        0.45870       0.00990
+/-3        0.49177       0.01255
+/-4        0.50158       0.01359
```

The mean change from 7 strips to 9 strips is only about 0.010 mm in u and
0.001 mm in w.  The remaining response structure is therefore not dominated by
the finite k=-4..+4 bank.

## Decision

The numerical weighting model is stable enough to move from integration
debugging to detector-response characterization.

Next measurement:

- scan `u0` across a single 4-mm wire cell;
- start at normal incidence to isolate the electrostatic/strip transfer
  function from finite-gap angular effects;
- measure wire assignment, raw stereo coordinate, bias/non-linearity and
  stochastic spread;
- only after that repeat representative nonzero incidence angles.

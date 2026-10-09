# C1 full-event diagnostic — full-cathode-gap sharing

Reference event:

- seed: 120001
- mu-, 5 GeV/c
- u0 = 1 mm, w0 = 0 mm at v = 0
- theta_u = 10 deg, theta_w = 25 deg
- symmetric 2.5 + 2.5 mm gap
- tan(alpha) = 0.1
- full-cathode-gap segmented weighting
- observation time = 100 ns
- all Heed electrons processed

## Microscopic chain

```text
clusters                 = 29
primary electrons        = 47
processed seeds          = 47
zero avalanches          = 5
energy loss              = 1274.62 eV
avalanche ions           = 177022
recorded ion births      = 177022
collected electrons      = 177107
active wires             = 1
ion transport failures   = 0
late electron segments   = 0
readout valid            = 1
```

The transport/integration chain is numerically stable.

## Event truth

From the per-seed CSV, avalanche-gain-weighted Heed seed coordinates are

```text
<u_seed>_gain = +1.00638 mm
<v_seed>_gain = +0.03700 mm
<w_seed>_gain = +0.01793 mm
```

which are consistent with the intended symmetric track intercept
`(u0,w0) = (1,0) mm`.

The microscopic avalanche ion-birth mean is instead

```text
<u_birth> = +0.01347 mm
<w_birth> = +0.00241 mm
```

because electrons focus onto the central anode wire in u.  Therefore
`<u_birth>` is **not** an appropriate truth coordinate for the sub-wire track
position; it is an avalanche-location diagnostic.  For w it remains useful
because the anode wire extends along w.

## Current readout result

```text
x_plus  CoG = -0.467765 mm
x_minus CoG = -0.606861 mm
u_stereo    = -0.698951 mm
w_stereo    = -0.539993 mm
wire CoG u  =  0.000000 mm
```

Compared with the physical track truth, the full-gap stereo estimator is poor:

```text
u_stereo - u0 = -1.69895 mm
w_stereo - w0 = -0.53999 mm
```

The wire CoG being zero is expected because every avalanche is collected by the
central wire; wire identity alone cannot recover the sub-wire track coordinate
u0 = 1 mm.

## Hybrid-weighting warning

```text
Q_cathode_minus = +3.95895 fC
Q_cathode_plus  = +3.88069 fC

Q_local_minus   = +0.0391265 fC
Q_local_plus    = -0.0391265 fC

minus scale     = +101.183
plus scale      =  -99.1833
```

The raw segmented sums nearly cancel while the wire-aware complete-cathode
signals remain O(4 fC).  The resulting ~100 rescaling is a warning that the
full-gap ideal planar segmented model is not a reliable approximation for the
complete full-track signal.

The exported strip pattern is also bipolar.  Using absolute net strip charge for
the CoG gives a strong common displacement of both stereo families.

Capture around the maximum absolute strip:

```text
              3 strips   5 strips   7 strips   9 strips
minus side     67.63%     88.68%     96.90%     99.06%
plus side      67.24%     89.56%     97.32%     99.21%
```

This is much broader than the symmetric truncated smoke event.

## Likely modelling limitation

The full-cathode-gap planar weighting construction removes the anode-wire
boundary from the segmented weighting problem.  In the real MWPC, cathode-strip
weighting potential must satisfy the conductor boundary condition on the anode
wires.  The historical half-gap Stage-B model did enforce a zero-weighting
boundary at the wire plane, but could not represent carriers from the opposite
half-gap.

The next controlled comparison is therefore **not** an angle scan.  It is the
same complete event with the historical half-gap/legacy sharing model, now that
both gas halves are populated.  If that reconstructs the track substantially
better, it isolates the full-gap planar extension as the problem.

The final production solution may require a weighting-field treatment that
includes the discrete anode wires explicitly rather than either planar
approximation.

# C1 smoke validation — symmetric two-sided baseline

Reference run:

- branch: `feature/mwpc-stage-c`
- random seed: 120001
- particle: mu-
- momentum: 5 GeV/c
- local track intercept: u0 = 1 mm, w0 = 0 mm at v = 0
- theta_u = 10 deg
- theta_w = 25 deg
- gap = 2.5 + 2.5 mm
- strip pitch/width = 1.7 / 1.7 mm
- tan(alpha) = 0.1
- B = 0 T
- observation time = 100 ns
- first 8 Heed electrons only (intentional smoke truncation)

## Result

The Stage-C1 chain is numerically valid end-to-end for the symmetric two-sided
geometry:

```text
clusters                 = 4
primary electrons        = 12
processed seeds          = 8 (truncated)
zero avalanches          = 1
avalanche ions           = 30006
collected electrons      = 30017
active wires             = 1
ion transport failures   = 0
readout valid            = 1
```

Microscopic avalanche-charge-weighted truth from ion birth positions:

```text
<u>_ion-birth = +0.00830 mm
<w>_ion-birth = -1.04469 mm
```

Readout estimators:

```text
u_stereo = -0.06358 mm
w_stereo = -1.16609 mm
u_wire   =  0.00000 mm
```

Differences to microscopic avalanche truth:

```text
u_stereo - <u>_ion-birth = -0.07188 mm
w_stereo - <w>_ion-birth = -0.12140 mm
u_wire   - <u>_ion-birth = -0.00830 mm
```

This is a strong integration sanity check.  In particular, the large ~1.3 mm
stereo-u bias seen with the asymmetric 2+4 mm geometry collapses to ~72 um in
the symmetric two-sided configuration.

## Cathode / strip diagnostics

```text
Q_cathode_minus = 0.694709 fC
Q_cathode_plus  = 0.604104 fC
Q_local_minus   = +0.045303 fC
Q_local_plus    = -0.045303 fC
minus scale     = +15.3349
plus scale      = -13.3349
```

The equal-and-opposite raw local sums are consistent with the symmetric
full-gap planar weighting construction.  The remaining side-dependent scale
factors reflect the hybrid nature of the model: wire-aware full-cathode
coupling combined with ideal planar segmented sharing.

Absolute-strip capture around the maximum strip:

```text
              3 strips   5 strips   7 strips
minus side     91.39%     99.12%     99.77%
plus side      89.43%     98.55%     99.77%
```

For this symmetric smoke event, five strips again capture essentially the full
ideal-model cluster, close to the behaviour established in Stage B.

## Interpretation boundary

This is not a detector-resolution measurement:

- only the first 8 Heed electrons were processed;
- the event is intentionally truncated in v;
- the absolute gain is not calibrated for the symmetric 2.5+2.5 mm geometry;
- electronics/noise/threshold/shaping remain absent;
- segmented weighting is still the hybrid ideal-planar approximation.

Because the event is truncated, the meaningful truth comparison is to the
avalanche ion-birth mean, not directly to the requested track intercept
(u0,w0).

## Next check

Run the same Geant4/Heed phase-space point with the **same random seed 120001**
and `--max-seeds 0`.  This gives an exact paired truncated/full-event
comparison before starting an angle or momentum ensemble.

# Stage C0 findings — finite working-area angular tails

Source: completed high-statistics MNP33-bore analysis, 200k parent dimuons/source.
This note records only the `useful_bore_working_area` selection.

Coordinate convention:

```text
theta_u = atan2(p_u, p_v) = atan2(pY_global, pZ_global)   # across wires
theta_w = atan2(p_w, p_v) = atan2(pX_global, pZ_global)   # along wires
```

## Main result

The finite station working-area cut removes only a small fraction of the
useful-bore population at MS2/MS3, but it does reduce the most extreme
downstream `theta_w` tails.

For the low-mass channels:

| source | station | retained vs useful bore | q95 |theta_u| | q99 |theta_u| | q95 |theta_w| | q99 |theta_w| |
|---|---:|---:|---:|---:|---:|---:|
| Phi | MS2 | 99.57% | 10.68 deg | 12.87 deg | 21.45 deg | 39.30 deg |
| Phi | MS3 | 98.90% | 10.59 deg | 12.77 deg | 20.72 deg | 36.05 deg |
| Omega | MS2 | 99.37% | 10.51 deg | 12.93 deg | 24.93 deg | 44.30 deg |
| Omega | MS3 | 98.52% | 10.40 deg | 12.81 deg | 23.95 deg | 40.90 deg |

So the large magnetic-bending tail is real for muons that still lie inside the
finite station benchmark.  It is not merely a population far outside the
station footprint.

## Tail fractions inside the finite working area

| source | station | |theta_w| > 30 deg | > 40 deg | > 50 deg |
|---|---:|---:|---:|---:|
| Phi | MS2 | 2.328% (5010) | 0.942% (2028) | 0.381% (821) |
| Phi | MS3 | 1.899% (4058) | 0.671% (1434) | 0.243% (519) |
| Omega | MS2 | 3.262% (6636) | 1.463% (2977) | 0.585% (1191) |
| Omega | MS3 | 2.771% (5587) | 1.106% (2230) | 0.361% (727) |

J/psi has essentially no corresponding large tail after the same cuts.

## Stage-C implication

For the first detailed C1 response scan:

- treat `theta_u` as a compact variable; about +/-15 deg already covers the
  99% scale for the accepted population;
- use a dense/core `theta_w` coverage through roughly +/-30 deg;
- include explicit tail validation points around 40 and 50 deg;
- do not build the full 2-D scan as a uniform +/-50 deg square, because the
  large tail is strongly anisotropic and mainly along `w`.

The >50 deg population is sub-percent and should be treated as a tail/stress
sample rather than the core grid.  The 30--40 deg population is still large
enough to be physically relevant.

## Important boundary

The working-area condition is applied at the nominal station plane.  It is not
yet a full finite-module-edge / stagger-layer intersection test.  Stage C1 can
therefore use these distributions to define the local response phase space,
while a later integration check should verify any very oblique tracks close to
actual chamber edges.

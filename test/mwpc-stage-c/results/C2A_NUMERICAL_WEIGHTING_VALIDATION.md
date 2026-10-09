# C2a numerical weighting validation — first successful map

Reference run:

- branch: `feature/mwpc-stage-c`
- symmetric gas gap: 2.5 + 2.5 mm
- wire pitch: 4.0 mm
- wire diameter: 30 um
- strip pitch/width: 1.7 / 1.7 mm
- tan(alpha): 0.1
- explicit wires: 5
- local patch: u = +/-12 mm, w = +/-8.5 mm
- Gmsh quadratic tetrahedral mesh
- Elmer electrostatic solve
- Garfield++ `ComponentElmer` import

## Boundary probe

```text
minus target cathode phi      = 0.98384148
minus opposite cathode phi    = 0.00068772
plus target cathode phi       = 0.98363079
plus opposite cathode phi     = 0.00068571
minus phi near central wire   = 0.03154289
plus phi near central wire    = 0.04010706
boundary sanity               = PASS
```

This is qualitatively the expected Shockley-Ramo weighting behaviour: the map
approaches one near the selected strip, zero near the opposite cathode and is
strongly suppressed near the grounded anode wire.

## Mesh inventory

The first mesh contains:

```text
nodes                         = 456378
quadratic tetrahedra (type11) = 279923
quadratic triangles (type9)   = 84978
```

The physical surface groups are present for both cathodes, the side walls and
the explicit anode wires.

## Central-strip profile checks

All exported profile values are finite.  Over the sampled profile:

```text
phi range                         = 0.00460591 .. 0.332590
max minus/plus reflected mismatch
  near selected cathode           = 0.005519
  far from selected cathode       = 0.0003114
```

The largest reflected-map mismatch is below 0.6 percentage points in absolute
weighting potential and is compatible with the present FEM discretisation.

At the central profile point, for example:

```text
minus phi(v=-1.25 mm) = 0.332590
minus phi(v=+0.05 mm) = 0.0315429
minus phi(v=+1.25 mm) = 0.0435719

plus  phi(v=-1.25 mm) = 0.0433000
plus  phi(v=+0.05 mm) = 0.0401071
plus  phi(v=+1.25 mm) = 0.331701
```

The expected cathode-exchange symmetry is therefore visible.

## Remaining interpolation warning

Garfield++ emitted two isolated `ComponentElmer::Coordinates13` convergence
warnings at

```text
(u,v,w) = (0, 0.05, 1.7) mm
```

during the profile scan.  The corresponding exported values remained finite,
so the first C2a map passes the coarse boundary test, but this warning should
not be accepted for production event integration: full avalanche transport will
query the map many times very close to the wire.

The next revision therefore keeps quadratic FEM basis functions but uses a
linear isoparametric geometry mapping for the quadratic mesh
(`Mesh.SecondOrderLinear = 1`).  This removes unnecessary curved-element
inversion around the 30-um wire while preserving the explicitly parameterised
wire boundary.

The revised probe also reports profile bounds and reflected-map symmetry
directly in the terminal.

## Decision

C2a has validated the numerical modelling route.  Do not connect this first map
to the full event yet.  First rerun the stability revision and confirm that the
`Coordinates13` warnings disappear; after that perform mesh/patch convergence
and generate a neighboring-strip bank.


## Stability rerun after linear isoparametric geometry mapping

The revised mesh with `Mesh.SecondOrderLinear = 1` was rerun successfully.

```text
minus target cathode phi      = 0.98384560
minus opposite cathode phi    = 0.00069177
plus target cathode phi       = 0.98363487
plus opposite cathode phi     = 0.00068976
minus phi near central wire   = 0.03300738
plus phi near central wire    = 0.04174073

profile phi range             = 0.00483802 .. 0.33288576
max reflected diff (near)     = 0.00552075
max reflected diff (far)      = 0.00031177

boundary sanity               = PASS
profile sanity                = PASS
```

No `ComponentElmer::Coordinates13` convergence warnings were observed.
C2a is therefore considered stable enough to proceed to event integration.

## C2b strategy

Neighbouring parallel strips do not require independent Elmer solves.  In the
ideal local chamber the electrodes and wires are translationally invariant
along local `w`.  For a strip family

```text
x = w cos(alpha) +/- u sin(alpha)
```

strip `k` is obtained from the central-strip weighting map by

```text
w -> w - k * strip_pitch / cos(alpha).
```

C2b therefore loads only the central minus/plus FEM maps and evaluates a
parameterised translated strip bank.  This preserves the explicit anode-wire
boundary while avoiding 2*N+1 separate Elmer solves per cathode.

The finite `w` patch must be wider than the translated strip bank plus the
event footprint.  The first C2b reference config therefore enlarges the patch
from +/-8.5 mm to +/-12 mm and uses strips k=-4..+4.

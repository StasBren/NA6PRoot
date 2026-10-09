# Numerical MWPC weighting fields (Stage C2)

This directory contains the first **numerical** replacement for the analytic
hybrid strip-sharing model used in Stage B/C1.

The goal is deliberately narrow:

```text
keep the validated analytic/Heed/Garfield microscopic transport
replace only the cathode-strip weighting potential
```

The numerical model is generated with **Gmsh + Elmer** and imported back into
Garfield++ with `ComponentElmer`.

## Why this is needed

The C1 full-event comparison showed that neither analytic approximation is
stable for reconstructing the weak stereo difference:

- the Stage-C full-gap planar model omits the anode-wire conductor boundary;
- the historical half-gap model enforces the wire-plane boundary but cannot
  describe carriers in the opposite half-gap.

For an exact Shockley-Ramo weighting solve, the selected strip is held at
1 V and every other conductor, including every anode wire, is held at 0 V.

## Parameterisation

All geometry/mesh inputs live in:

```text
weighting_config.json
```

The current reference values are:

```text
gap minus / plus     2.5 / 2.5 mm
wire pitch           4.0 mm
wire diameter        30 um
strip pitch / width  1.7 / 1.7 mm
tan(alpha)           0.1
explicit wires       5  (half_wires = 2)
local patch          u = +/-12 mm, w = +/-8.5 mm
```

Nothing in the generator hard-codes these values. In particular the following
can be varied independently:

- `gap_minus_mm`, `gap_plus_mm`;
- `wire_pitch_mm`, `wire_diameter_um`, number of explicit wires;
- strip pitch, strip width and stereo angle;
- local patch size;
- anode voltage used for the reference physical-field solution;
- mesh sizes/refinement distances;
- strip indices and cathode sides for which weighting maps are solved.

To scan another geometry, copy the JSON file, edit the parameters and pass the
new file to the runner.

## Geometry and coordinates

The numerical coordinates are exactly the Stage-C local coordinates:

```text
Elmer/Gmsh x = u  across wires
Elmer/Gmsh y = v  chamber normal
Elmer/Gmsh z = w  along wires
```

The gas volume is a finite local box. Cylindrical anode wires are explicitly
subtracted from this volume, so their surfaces are real conductor boundaries.

The artificial side walls have natural zero-normal-flux boundary conditions.
They are **not** detector electrodes. Patch-size convergence therefore has to
be checked before production use.

The cathode-strip pattern is imposed as a coordinate-dependent Dirichlet
condition on the cathode plane. This avoids an unnecessary geometric
fragmentation of the plane. For the central minus strip, for example,

```text
| w cos(alpha) - u sin(alpha) | <= strip_width / 2
```

is set to 1 V; the rest of the cathode, the opposite cathode and every anode
wire are at 0 V.

For `strip_width < strip_pitch`, the present prototype treats the uncovered
cathode area as grounded. A later PCB model can replace this with explicit
dielectric gaps/FR4 if needed.

## C2a smoke test

The first smoke test solves:

1. the physical reference field (anode at the configured HV, cathodes at 0 V);
2. the central minus-strip weighting potential;
3. the central plus-strip weighting potential;
4. imports all three maps into Garfield++;
5. checks the weighting potential near the target cathode, opposite cathode
   and central wire;
6. exports a central-strip weighting profile.

Run:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

source ~/na6p/env.sh
source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh

bash test/mwpc-stage-c/run_c2_weighting_smoke.sh
```

Or use another parameter file:

```bash
bash test/mwpc-stage-c/run_c2_weighting_smoke.sh \
  /path/to/my_weighting_config.json
```

Runtime maps go to:

```text
test_runs/mwpc_stage_c/c2_weighting_smoke/
```

and are intentionally not committed because second-order FEM meshes/results can
be large.

## Dependencies

This step adds two external command-line tools beyond Garfield++:

- Gmsh;
- Elmer (`ElmerGrid` and `ElmerSolver`).

The runner performs a preflight check and exits before doing any work if one of
them is unavailable.

## What comes after C2a

Do **not** immediately connect the numerical map to the full muon event.

First validate the map itself:

- target strip approaches weighting potential 1 near its cathode;
- opposite cathode approaches 0;
- weighting potential approaches 0 on the anode wire;
- map is finite and smooth in the gas;
- repeat with larger local patch / finer mesh to establish convergence.

After this passes, C2b will solve a small bank of neighboring strip maps
(`strip_indices` is already parameterised) and replace `MWPCStereoReadout`
inside the existing C1 event accumulator. The Heed, avalanche, ion transport
and event-level logic remain unchanged.

# DiCE MWPC implementation handoff

This branch (`feature/mwpc-dice`) replaces the placeholder silicon muon-spectrometer planes with a detailed MWPC chamber geometry and a configurable six-station staggered layout. This document summarizes the current working point, implementation structure and validation status for review.

## Coordinate convention

Detector coordinates are:

- X: horizontal, positive toward Jura when looking downstream;
- Y: vertical, positive upward;
- Z: along the beam.

The MNP33 field is vertical. The MWPC wires therefore run horizontally, along detector X.

The detailed chamber builder preserves the prototype-local mechanical convention (`local x = short side`, `local y = long side`) and the complete chamber is installed with a +90 degree rotation around Z. This avoids swapping only the body dimensions while leaving the asymmetric electronics/divider pieces in the wrong orientation.

## Current working geometry

Sensitive gas dimensions in installed detector coordinates:

```text
standard chamber : 66.56 x 51.20 cm  (X x Y)
MS0 chamber      : 66.56 x 30.72 cm  (X x Y)
```

Current grid:

```text
station   Nx x Ny   chambers   active span X x Y [cm]
MS0        4 x 8       32       257.24 x 224.76
MS1        4 x 5       20       257.24 x 244.00
MS2        6 x 7       42       384.36 x 340.40
MS3        6 x 7       42       384.36 x 340.40
MS4        8 x 9       72       511.48 x 436.80
MS5        9 x 10      90       575.04 x 485.00
                              ----------------
                              298 chambers total
```

Common layout parameters:

```text
active overlap Ox = Oy = 3 cm
adjacent stagger step Delta z = 4 cm
A/B/C/D gas-centre offsets = -6, -2, +2, +6 cm
station Z = 300, 340, 530, 590, 810, 850 cm
```

Checkerboard assignment is

```text
q = 2*(row % 2) + (col % 2)
q = 0,1,2,3 -> A,B,C,D
```

## Main implementation files

`sim/include/NA6PMWPCParam.h`
: User-facing acceptance/layout controls only; detailed MWPC construction constants are kept in `NA6PMWPCChamber`.

: Configurable MWPC geometry/layout parameters: detailed chamber dimensions, rotation, gas mixture, station working areas, `Nx,Ny`, active overlaps, stagger spacing and MS0 narrow dimension.

`sim/include/NA6PMWPCChamber.h`, `sim/src/NA6PMWPCChamber.cxx`
: Detailed reusable chamber builder, including the layered mechanical/material model and Ar/CO2 sensitive gas.

`sim/src/NA6PMuonSpec.cxx`
: Builds the six stations, places chambers on the regular X/Y grid, assigns A/B/C/D z offsets, rotates each detailed chamber as a rigid object, and keeps sequential sensitive-volume IDs.

`base/include/NA6PLayoutParam.h`
: Current MS station Z positions. MS1 is at 340 cm in this branch.

`NA6PSim.cxx`

## Validation

The current 298-chamber working point has been checked with the branch build/install and a fresh full NA6PRoot geometry.

The following tests pass:

```text
- full branch compilation/install
- real Geant4 transport smoke test
- checkMWPCStaggeredLayout.C
- checkMWPCWorkingAreaCoverage.C
- ROOT CheckOverlaps(1e-4): 0 illegal overlaps/extrusions
- 20,000 J/psi events through the dimuon pipeline
- 20,000 omega events through the dimuon pipeline
- 20,000 phi events through the dimuon pipeline
```

The geometry checkers do not hard-code the present `Nx,Ny`. They infer the grid from the generated geometry, so later changes to chamber counts can be checked without rewriting the tests.

For the full generated geometry:

```bash
root -l -b -q 'test/checkMWPCStaggeredLayout.C("geometry.root")'
root -l -b -q 'test/checkMWPCWorkingAreaCoverage.C("geometry.root")'
```

A ROOT overlap check should return zero illegal overlaps/extrusions.

## Dimuon transport validation

The intended validation pipeline is documented in `test/mwpc-dimuon-validation/README.md`.

Typical run:

```bash
bash test/mwpc-dimuon-validation/run_validation.sh 20000 4
```

This runs J/psi, omega and phi samples through the actual NA6PRoot/Geant4 geometry. The validation hook forces the signal parents to the dimuon channel for this test only.

The direct-hit plots use the actual sensitive-gas entrance coordinates from `HitsMuonSpec.root`; no track extrapolation is used. Maps are approximately 1 cm x 1 cm and show the complete configured chamber grid with an outer margin.

## Isolated chamber test

`test/mwpc-particle-gun/` provides a small standalone Geant4VMC transport smoke test of one detailed chamber. It uses the same `NA6PMWPCChamber` builder as the production geometry.

## Present scope

The branch models geometry, materials and Geant4 transport. It does not yet model MWPC signal formation: no electron/ion drift, avalanche gain, induced strip/wire signals, electronics threshold, detector inefficiency or reconstruction response is included.

# MWPC dimuon validation

This directory contains the full-NA6PRoot transport validation used for the current DiCE MWPC implementation on `feature/mwpc-dice`.

The test uses the existing `test/genDimuonBgEvent.C` generator for three signal channels:

- `Jpsi` (PDG 443)
- `Omega` (PDG 223)
- `Phi` (PDG 333)

`mwpcDimuonHooks.C` forces each signal parent to decay to `mu+ mu-` with 100% branching fraction for this validation only and keeps both daughter muons in `MCKine.root`, including daughters that do not hit an MWPC sensitive volume. Production decay tables are not changed in the source code.

## Current MWPC working point

Detector coordinates are X horizontal / along the wires, Y vertical / across the wires, Z along the beam. The detailed chamber model keeps the prototype-local mechanical axes and is installed with a +90 degree rotation around Z.

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

Common geometry parameters:

```text
standard active gas : 66.56 x 51.20 cm in detector X,Y
MS0 active gas      : 66.56 x 30.72 cm in detector X,Y
active overlap      : Ox = Oy = 3 cm
stagger step        : Delta z = 4 cm
A/B/C/D gas z       : -6, -2, +2, +6 cm relative to station centre
MS z positions      : 300, 340, 530, 590, 810, 850 cm
```

The grid is a current working point, not a frozen final optimization. `NA6PMWPCParam` keeps the acceptance/layout knobs (`Nx`, `Ny`, overlap, chamber size and stagger spacing) configurable; detailed chamber construction stays fixed in `NA6PMWPCChamber`.

## Run

The default validation uses 40 GeV/nucleon beam energy, one signal parent per event, `0 < pT < 5 GeV/c`, `0 < y < 6`, no hadronic background and 20,000 events per channel.

After building/installing this branch and selecting the branch installation in the environment:

```bash
bash test/mwpc-dimuon-validation/run_validation.sh 20000 4
```

The second argument is the number of independent `na6psim` worker processes. `NA6PSim` uses Geant4 in single-application mode; `na6psim_parallel` supplies process-level parallelism. The script first performs a one-event compile/smoke run and then runs J/psi, omega and phi sequentially.

Outputs are written to:

```text
test_runs/mwpc_dimuon/Jpsi/
test_runs/mwpc_dimuon/Omega/
test_runs/mwpc_dimuon/Phi/
```

## Outputs

`analyzeMWPCDimuon.C` checks the generated decay and writes:

```text
plots/kinematics_<channel>.png
plots/kinematics_<channel>.root
```

The analysis contains parent |p|, parent pT, parent rapidity and daughter-muon |p|/pT distributions, and reports four-momentum closure of the forced two-body decay.

`plotMWPCHitDensity.C` reads the actual `HitsMuonSpec.root` Geant4 sensitive-gas hits. It does not extrapolate tracks to nominal planes. Every stored chamber crossing is filled at `hit.getXIn(), hit.getYIn()` after subtracting the station centre.

For each station it writes:

```text
plots_direct_hits/direct_hits_MS<n>_<channel>.png
plots_direct_hits/hit_multiplicity_MS<n>_<channel>.png
plots_direct_hits/direct_hit_validation_<channel>.root
```

The direct-hit maps use approximately 1 cm x 1 cm bins and are scaled to hit density in counts/cm^2. Each station has its own color scale. The axes are derived from the complete configured active grid plus an outer margin, so the full outermost chamber row and column are visible rather than clipped.

The multiplicity histogram is per generated direct daughter muon. A value of zero means that daughter produced no sensitive-gas hit in that station; values above one occur when the trajectory crosses more than one staggered chamber because of overlap/stagger geometry.

`plotMWPCHitMultiplicityMap.C` is an optional diagnostic for the local mean number of chamber hits per contributing muon. It is not part of the default validation script.

## Geometry regression checks

The full generated geometry should also pass:

```bash
root -l -b -q 'test/checkMWPCStaggeredLayout.C("geometry.root")'
root -l -b -q 'test/checkMWPCWorkingAreaCoverage.C("geometry.root")'
```

`checkMWPCStaggeredLayout.C` infers the configured `Nx,Ny` from the generated geometry and checks the installed active dimensions, 3 cm overlap, checkerboard parity and A/B/C/D z offsets. It therefore remains valid when the grid counts are changed.

`checkMWPCWorkingAreaCoverage.C` verifies that the active footprint covers the configured station working rectangle. Minimality is reported only as information; deliberately adding rows or columns is allowed.

A full ROOT `CheckOverlaps(1e-4)` should report zero illegal overlaps/extrusions before a geometry is treated as validated.

## Scope

This is a transport-level geometry validation. It includes the full detailed chamber material geometry and Geant4 transport, but no MWPC avalanche, electron/ion drift, wire response, threshold, electronics response, intrinsic detector inefficiency or reconstruction model.

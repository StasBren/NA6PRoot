# MWPC particle-gun smoke test

This is a deliberately small transport test for the detailed `NA6PMWPCChamber` geometry. It is the NA6PRoot-native analogue of the earlier isolated-chamber muon tests in `na60-dice-mwpc-simulation`.

The test creates one complete MWPC chamber at the origin, shoots one controlled forward muon per event through the centre, and lets Geant4 transport the particle through the ROOT/TGeo geometry via Geant4VMC. It is intentionally independent of the full station layout so that the chamber implementation can be checked in isolation.

## Current production integration

The same chamber builder is integrated into `NA6PMuonSpec`, which constructs the current provisional six-station staggered MWPC layout:

```text
MS0   4 x 8   = 32
MS1   4 x 5   = 20
MS2   6 x 7   = 42
MS3   6 x 7   = 42
MS4   8 x 9   = 72
MS5   9 x 10  = 90
                 ---
                 298
```

Detector X is horizontal / along the wires and detector Y is vertical / across the wires. The detailed prototype chamber keeps its local mechanical axes and is installed with a +90 degree rotation around Z.

The installed sensitive-gas sizes are:

```text
standard chamber : 66.56 x 51.20 cm in detector X,Y
MS0 chamber      : 66.56 x 30.72 cm in detector X,Y
```

Active gas overlaps are `Ox = Oy = 3 cm`. Chambers are assigned to four checkerboard parity classes A/B/C/D and their gas centres are staggered at `-6, -2, +2, +6 cm` relative to the station centre, corresponding to an adjacent layer spacing of `Delta z = 4 cm`.

The current station Z positions are `300, 340, 530, 590, 810, 850 cm` for MS0..MS5.

`test/checkMWPCStaggeredLayout.C` checks the generated full-station geometry without hard-coding `Nx,Ny`: it infers the grid from chamber centres and verifies installed gas dimensions, overlap, parity and four-level stagger. `test/checkMWPCWorkingAreaCoverage.C` verifies coverage of the configured working rectangles.

The grid above is a working point for the present acceptance study, not a frozen final optimization.

## What the isolated particle-gun test checks

- the installed branch can link and instantiate `NA6PMWPCChamber`;
- the complete detailed chamber is passed to ROOT/TGeo and Geant4VMC;
- ROOT reports no overlaps in the isolated geometry;
- the chamber gas is recognized as a sensitive NA6PRoot volume;
- the primary muon enters the gas and produces Geant4 transport steps there;
- gas energy deposition is recorded as a simple transport-level diagnostic.

This is **not** an MWPC signal-response simulation. There is no avalanche gain, drift, wire response, threshold, electronics, intrinsic inefficiency or reconstruction.

## Prerequisite

Build and install the current `feature/mwpc-dice` branch first. The test links to the installed `simLib`, `baseLib` and `utilsLib`, so this also checks that the new chamber code is part of a real NA6PRoot build rather than only a standalone stub.

With the usual O2/NA6PRoot environment loaded:

```bash
cd NA6PRoot
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/NA6PRoot/install"
cmake --build build --parallel 5 --target install
```

Select that installation in the environment before running the test.

## Run

From the repository root:

```bash
bash test/mwpc-particle-gun/run.sh
```

Defaults are 100 mu- events at 20 GeV/c with seed 20260914.

Arguments are:

```text
run.sh [events] [momentum_GeV] [PDG] [seed]
```

Examples:

```bash
# 100 central 20 GeV/c mu-
bash test/mwpc-particle-gun/run.sh 100 20 13

# 50 central 4 GeV/c mu+
bash test/mwpc-particle-gun/run.sh 50 4 -13
```

Each invocation creates a new time-stamped directory under `runs/`.

## Output

`summary.txt` contains the overall PASS/FAIL result. The run passes only when the primary muon enters the sensitive gas and has gas transport steps in every event.

`events.csv` contains one row per event:

```text
event,gas_steps,gas_entries,gas_exits,primary_steps,primary_entries,energy_deposit_keV
```

`geometry.root` is the exact isolated TGeo geometry given to Geant4VMC. `geometry_summary.txt` records the gas dimensions/centre and full chamber envelope. `geometry_checks.txt` records ROOT's overlap result. `MCKine.root` is the usual NA6PRoot Monte Carlo stack output. Empty hit files from other registered NA6PRoot modules may also appear because the standard `NA6PMC::init()` lifecycle is reused; they are not part of this test.

## Code flow

```text
mwpc_gun.C
    -> one known forward muon/event

MWPCTestModule
    -> calls the real NA6PMWPCChamber builder
    -> gas gets a normal NA6PRoot sensitive-volume ID

MWPCParticleGunApplication
    -> creates an air world + one chamber only
    -> exports geometry.root
    -> hands the TGeo geometry to Geant4VMC with geomRoot

Geant4
    -> transports the muon
    -> NA6PMC::Stepping() routes sensitive gas steps to MWPCTestModule::stepManager()
    -> events.csv / summary.txt
```

The local `MWPCTestModule` remains test scaffolding only. Production geometry uses the same `NA6PMWPCChamber` builder through `NA6PMuonSpec`.

# MWPC Phase A: microscopic wire / avalanche sandbox

This directory is intentionally standalone from the normal NA6PRoot build.
The validated `feature/mwpc-dice` geometry is therefore not affected by the
microscopic-response work.

## Goal

Phase A establishes the microscopic MWPC model before adding cathode strips
or reconstruction:

1. analytic wire/cathode electric field;
2. Ar/CO2 electron transport with Magboltz;
3. microscopic electron drift and avalanche;
4. configurable magnetic field;
5. basic wire-collection and avalanche-size checks.

No strip weighting fields, Shockley-Ramo signal calculation, electronics, or
NA6PRoot digitisation are included yet.

## Coordinate mapping

`ComponentAnalyticField` uses wires parallel to its local z axis. We map this
to the detector convention as

- Garfield x = detector/global Y (vertical, across the wires),
- Garfield y = detector/global Z (through the gas gap),
- Garfield z = detector/global X (horizontal, along the wires).

Thus a physical vertical MNP33 field (+global Y) is supplied as +Garfield x.

## Initial benchmark parameters

The first benchmark deliberately reproduces the documented Prototype-3
microscopic parameters rather than claiming to be the final chamber design:

- Ar/CO2 = 70:30,
- wire diameter = 30 um,
- wire pitch = 4 mm,
- asymmetric cathode gaps = 2 + 4 mm,
- anode voltage = 1.8 kV.

After this benchmark works, the next scan should compare the symmetric
3+3 mm geometry used by the present 6 mm gas model and the proposed
2.5+2.5 mm option.

## Build

First source the Garfield++ environment, then build this directory by itself:

```bash
source /path/to/Garfield/install/share/Garfield/setupGarfield.sh

cd ~/NA6PRoot/test/mwpc-avalanche
cmake -S . -B build
cmake --build build -j4
```

Run at B=0 and HV=1800 V:

```bash
./build/mwpc_phase_a 0.0 1800
```

The program writes `phase_a_results.csv`.

The first argument is the physical vertical magnetic field in tesla; the
second is the anode voltage in volts. For example:

```bash
./build/mwpc_phase_a 0.5 1800
```

## First acceptance checks

For B=0:

- electrons launched within one 4 mm pitch should be collected by the
  geometrically expected neighbouring wire;
- the response should be symmetric about x=0;
- the avalanche should remain finite and reproducible statistically.

Then turn on B and compare the endpoint/wire assignment and avalanche
statistics. Plotting and drift-line visualisation are the next step once the
basic executable is confirmed to compile and run.

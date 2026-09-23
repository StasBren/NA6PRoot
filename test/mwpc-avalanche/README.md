# MWPC microscopic response sandbox — Phase A

This directory is intentionally **standalone** from the main NA6PRoot build.
Garfield++ is not added as a production dependency of NA6PRoot at this stage.

The goal of Phase A is to understand the microscopic gas response of a small
MWPC cell before adding cathode-strip readout or integrating a parameterized
response into the full detector simulation.

## Coordinate convention

The detector convention is

- detector **X**: horizontal, along the anode wires;
- detector **Y**: vertical;
- detector **Z**: along the beam.

`ComponentAnalyticField` uses wires parallel to its local z axis. Therefore
the sandbox uses

- Garfield x = detector Y,
- Garfield y = detector Z,
- Garfield z = detector X.

The MNP33 magnetic field +Y therefore maps to Garfield +x.

## Geometry is intentionally configurable

The chamber/readout design is not frozen, so the microscopic model must not
assume a single final cathode spacing or strip topology.

The current default operating point is only a Prototype-3 reference:

- Ar/CO2 = 70:30;
- anode-wire diameter = 30 um;
- wire pitch = 4 mm;
- wire plane between cathodes with 2 mm and 4 mm gaps;
- anode voltage = +1.8 kV;
- magnetic field = 0 T for the first smoke test.

The two wire-to-cathode gaps are independent runtime parameters. Examples:

```bash
# Prototype-3-like asymmetric configuration
./mwpc_phase_a --gap-minus-mm 2 --gap-plus-mm 4

# symmetric 3+3 mm
./mwpc_phase_a --gap-minus-mm 3 --gap-plus-mm 3

# symmetric 2.5+2.5 mm
./mwpc_phase_a --gap-minus-mm 2.5 --gap-plus-mm 2.5
```

Wire pitch and diameter are configurable as well:

```bash
./mwpc_phase_a --pitch-mm 4 --wire-diam-um 30
```

This means Phase A can compare field, drift and avalanche behaviour for
different gap choices without changing the source code.

## Readout topology to be tested in Phase B

The strip/readout design is also deliberately not frozen. Phase B will support
at least two distinct configurations rather than baking one into the geometry:

1. **two-sided readout** — signal pickup on both cathode sides, with the two
   strip coordinates distributed between the two cathodes;
2. **single-sided readout** — both strip-coordinate patterns are placed on one
   readout cathode while the opposite cathode is not used for strip pickup.

These will share the same Phase-A gas/wire transport model. The difference
enters when strip electrodes, weighting fields, Shockley-Ramo induced signals,
charge sharing and reconstruction are added.

Keeping the gas geometry and readout topology separate is intentional: we want
to be able to scan, for example, symmetric versus asymmetric wire placement
independently of one-sided versus two-sided strip readout.

## Prerequisite

Garfield++ must be installed and its environment loaded. For an installed
Garfield++ tree this is normally

```bash
source /path/to/garfield/install/share/Garfield/setupGarfield.sh
```

The official Garfield++ build requires ROOT 6, GSL, CMake, a compatible C++
compiler and a Fortran compiler.

## Build

From this directory:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j4
```

If CMake cannot find Garfield++, prepend its install prefix to
`CMAKE_PREFIX_PATH`, for example:

```bash
cmake .. -DCMAKE_PREFIX_PATH="$GARFIELD_HOME/install"
```

## First test: one electron, B = 0

```bash
./mwpc_phase_a
```

The program constructs nine parallel wires between two cathode planes, launches
one low-energy electron, runs microscopic transport and avalanche multiplication,
and prints the avalanche size and electron endpoints.

Useful options:

```bash
./mwpc_phase_a --hv 1800 --b 0
./mwpc_phase_a --gap-minus-mm 3 --gap-plus-mm 3
./mwpc_phase_a --x0 0.10 --y0 0.20
```

## What this test proves — and what it does not

A PASS means only that

1. the Ar/CO2 Magboltz medium is usable,
2. the analytic wire/cathode field is valid,
3. an electron is transported,
4. microscopic multiplication can be generated.

It does **not** yet validate the physical gas gain. Penning transfer,
gas conditions, exact gap geometry and high-voltage operating point will need
to be tied to prototype measurements before interpreting the avalanche charge
quantitatively.

## Next Phase-A checks

After the smoke test works:

1. scan starting position across one wire pitch at B = 0;
2. scan HV and record the gain distribution;
3. scan the two cathode gaps, including symmetric and asymmetric cases;
4. switch on a constant vertical MNP33 field and measure the arrival shift;
5. replace a single starting electron by ionisation clusters from a muon track.

Cathode-strip weighting fields and Shockley-Ramo signals belong to Phase B.

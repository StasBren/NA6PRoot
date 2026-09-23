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

## Phase-A starting geometry

The initial model follows the Prototype-3 operating point as a controlled
starting configuration:

- Ar/CO2 = 70:30;
- anode-wire diameter = 30 um;
- wire pitch = 4 mm;
- wire plane between grounded cathodes with 2 mm and 4 mm gaps;
- anode voltage = +1.8 kV;
- magnetic field = 0 T for the first smoke test.

All of these are model parameters, not assumptions about the final chamber.
In particular, later scans should include symmetric gaps and the actual MNP33
field.

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
./mwpc_phase_a --x0 0.10 --y0 0.30
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

1. scan starting position across one 4 mm wire pitch at B = 0;
2. scan HV and record the gain distribution;
3. switch on a constant vertical MNP33 field and measure the arrival shift;
4. replace a single starting electron by ionisation clusters from a muon track.

Cathode-strip weighting fields and Shockley-Ramo signals belong to Phase B.

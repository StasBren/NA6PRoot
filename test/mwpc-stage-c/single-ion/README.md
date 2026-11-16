# NA60+ MWPC: single-ion current on a stereo cathode strip

**Status as of 2026-11-16: GitHub source committed; Garfield++/Elmer numerical execution is NOT verified.**

This is an isolated first response benchmark. It is NOT a validated
muon waveform, a production simulation, or a selected operating HV.

## Physics and scope

The earlier `phaseB_strip_weighting.cxx` models strips running along u,
segmented in w: that is not the planned stereo strip. The earlier
`phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx` calculates whole-plane
cathode currents, not individual readout strips.

The new benchmark reuses the existing 3D Gmsh/Elmer generator under
`test/mwpc-stage-c/numerical-weighting/`. Explicit wire surfaces are held
at zero weighting potential; the **selected physical oblique cathode
strip** is set to 1 V, and the other conductors to zero. A second map
selects the central strip on the opposite cathode. Garfield++
`ComponentElmer` provides weighting potentials, while
`ComponentAnalyticField` provides the *real drift field* and
`DriftLineRKF` propagates one positive ion.

The separate JSON here sets P3 noFW: anode pitch 3 mm, 30 um diameter,
cathodes at ±3 mm, Ar/CO2 70:30, B=0 and **1800 V TEST ONLY**;
1.7 mm strip width/pitch and tan(alpha)=0.1. Actual strip dimensions
and stereo sign remain provisional.

## First WSL execution

Use a separate checkout; do not switch the laptop HV-scan checkout.

```bash
git clone --branch feature/mwpc-single-ion-stereo-20261116 \
  https://github.com/StasBren/NA6PRoot.git \
  ~/na6p/mwpc-cathode-single-ion

bash ~/na6p/mwpc-cathode-single-ion/test/mwpc-stage-c/single-ion/run_single_ion.sh
```

The runner requires Gmsh, ElmerGrid, ElmerSolver, CMake, Garfield++ and ROOT.
It automatically looks for the provisional Garfield data file
`IonMobility_Ar+_Ar.txt`. If unavailable, set
`MWPC_ION_MOBILITY=/full/path/to/mobility.txt`.

The runner never alters historical references, the ongoing HV scan or
other build folders. It **refuses to overwrite an existing run directory**.
For another run choose a fresh path:

```bash
SINGLE_ION_OUT="$HOME/na6p/new_single_ion_test" \
  bash ~/na6p/mwpc-cathode-single-ion/test/mwpc-stage-c/single-ion/run_single_ion.sh
```

To study another HV, copy and edit the JSON configuration, then set
`SINGLE_ION_CONFIG=/absolute/path/new.json` AND a new
`SINGLE_ION_OUT`. No future operating voltage has yet been selected.

## Acceptance gates and outputs

One +e ion starts at local (u,v,w)=(0,0.04,0.2) mm.
Ion diffusion is disabled solely for a deterministic numerical test.
The C++ executable checks selected/opposite cathode weighting potentials,
physical cathode arrival, nontruncated time windows, and the signed
Shockley–Ramo charge identity:

```text
Q_sensor = sum_b I_b * dt
Q_Ramo   = (+e) * [phi_w(endpoint) - phi_w(birth)]
```

It repeats for dt=10,5,2 ns and maximum drift steps 0.04,0.02 mm,
checking both strips. The code also compares cumulative charge at fixed
physical times (0.1, 1, 10, 50 and 100 microseconds) across time-bin
choices using a 5%-of-Ramo-charge tolerance (floor 0.01 elementary charges).
The endpoint charge tolerance is
max(3% of |Q_Ramo|, 0.005 elementary charges), with additional
endpoint/step stability gates. Failures return nonzero and must not be
described as validation.

Outputs include `single_ion_summary.csv`, high-resolution ion-current
and cumulative-charge CSVs, Gmsh/Elmer logs, map manifests and SHA256 sums.
ROOT produces ATLAS-like `single_ion_minus/plus.png` and `.pdf` plots
labelled **NA60+ Simulation**.

## Important unresolved limitations

- Finite FEM domain and Neumann artificial side walls; mesh/domain-size
  convergence is still required. A passing Ramo integral does not prove
  that the weighting field is a correct representation of hardware.
- Inter-strip gaps are idealized grounded boundaries, not full PCB/FR4
  material, traces and guard structures.
- The Ar+ mobility in pure Ar is only a provisional approximation for
  Ar/CO2 70:30. Ion drift times require validated species/mobility.
- There is no microscopic avalanche, prompt electron current,
  electronics, diffusion, space charge or muon signal in this benchmark.
- Different start positions and mirror-side tests are still required.

## Run/dataset status update from the user, 2026-11-16

The group is reconsidering the anode voltage **upward for higher gain**.
The higher-HV operating point is undecided pending gain, fluctuation,
avalanche-cap and space-charge-model checks.

The user reports a separate HPC calculation of **approximately 8,000
muons at 1800 V** is nearly finished; this has not been independently
verified. Preserve its manifests/raw files and use it, after audit,
as a historical independent test of Fast4 and potentially a donor/holdout
resource for Fast5 development. Do not confuse it with the archived
30k dataset, completed 1000 per-wire sample or current HV scans.
Verify whether it stores ion-birth positions and times; summary moments
alone cannot provide physical strip waveforms.

Next gate: actual WSL compile and Garfield/ROOT run, then mesh
refinement, opposite hemisphere and independent ion-start tests.

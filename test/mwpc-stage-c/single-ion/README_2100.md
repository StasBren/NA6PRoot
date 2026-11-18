# 2100 V P3 two-sided *physical stereo strip* one-ion benchmark (WSL)

This is the existing `stereo_single_ion.cxx` and `run_single_ion.sh`, not a new physical model. Source branch `feature/mwpc-single-ion-stereo-20261116` supplies C++ Garfield++ ion RKF/3D Gmsh/Elmer strip weighting and ROOT plots; new wrapper `run_single_ion_2100.sh` supplies the chosen 2100 V P3 geometry and the known local Ar+ mobility path. **Status: static source checked, but numerical runtime/Gmsh-Elmer/ROOT convergence has not yet run or passed on the laptop.**

## Reproduce without changing frozen laptop and cluster archives

```bash
git clone --depth 1 --branch feature/p3-2100v-stereo-single-ion-20261118 \
  https://github.com/StasBren/NA6PRoot.git \
  ~/na6p/cathode-stereo-2100
cd ~/na6p/cathode-stereo-2100
bash test/mwpc-stage-c/single-ion/run_single_ion_2100.sh
```

The wrapper requires `gmsh`, `ElmerGrid`, `ElmerSolver`, `cmake`, ROOT and Garfield++. Its physics calculations and statistical checks are in C++/ROOT; Python is used **only to create the FEM geometry/solver inputs**. It explicitly exports `LD_LIBRARY_PATH` for the known WSL Garfield++ install, solving the prior `libHeed.so` runtime issue. The provisional ion-mobility table defaults to `~/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt`, and may be overridden with `MWPC_ION_MOBILITY`.

## Input geometry and resource implications

- Nine anodes (local indices -4..+4), diameter 30 um, pitch 3 mm, at 2100 V
- Both cathodes at v=+/-3 mm, readout central positive/negative *stereo strips*
- Initial ion position (u,v,w)=(0,0.04,0.2) mm; radius from wire center 40 um (>15 um wire radius)
- Provisional strip pitch/width 1.7 mm, tan(alpha)=0.1; not a final detector specification
- Explicit conducting wire surfaces in FEM weighting solves; 3D mesh u extent +/-15.5 mm, w extent +/-11 mm; artificial side walls
- Ion transport physical field: Garfield `ComponentAnalyticField`; signal weighting: Garfield `ComponentElmer` imported from Gmsh/Elmer
- Six drift runs: time bins 10/5/2 ns for RKF max step 0.04/0.02 mm, 250 us observation window

The Gmsh second-order FEM mesh can take substantial RAM/disk space. The original C2a baseline at 4 mm pitch and +/-2.5 mm gaps is **not** a substitute for this new 3 mm pitch, +/-3 mm gap geometry. The first run is a numerical pilot, not a validated full-chamber/PCB map.

## Acceptance gates

Look for `PASS single-ion current/charge and convergence gates` in `results/physics.log`, `GATES=PASS` in `results/single_ion_metadata.txt`, sensible selected/opposite strip potentials, cathode arrival, complete recording window, and charge differences between sensor integral and Ramo endpoint. Additionally inspect early charge at 0.1, 1, 10, 50 and 100 us, and convergence across dt and RKF steps. Verify `sha256sum -c SHA256SUMS` from the run root. **Important limitation:** `ion.UseWeightingPotential(true)` uses the weighting potential to calculate current, so endpoint/integral agreement is not an independent test of the FEM field being correct.

Successful output path: `test_runs/mwpc_stage_c/p3_2100_stereo_9wire_v1/`. Contains maps, model config, `results/single_ion_summary.csv`, waveforms and ROOT-generated PNG/PDF/ROOT plots, plus manifest and checksums. The old microscopic data and earlier Stage C directory are not accessed or overwritten.

## Storage cleanup AFTER the numerical result passes

First inspect size:

```bash
RUN=test_runs/mwpc_stage_c/p3_2100_stereo_9wire_v1
du -sh "$RUN" "$RUN/map" "$RUN/results"
(cd "$RUN" && sha256sum -c SHA256SUMS)
```

The Gmsh transfer file `map/chamber.msh` is an intermediate duplicated by Elmer's imported `map/chamber/mesh.*`; it is not part of the current SHA manifest. Once Elmer produces all three `.result` files and the ROOT checks pass, it is safe to **preview** and optionally delete this one bulky intermediate:

```bash
ls -lh "$RUN/map/chamber.msh"
rm -i -- "$RUN/map/chamber.msh"
```

Keep `map/chamber/mesh.*`, `map/chamber/drift.result`, `map/chamber/weight_minus_0.result`, `map/chamber/weight_plus_0.result` until subsequent strip signal work. Do NOT delete current `~/na6p/cathode-ion-2100-wsl` weighted archives or muon 0 full ion archive. If plot/waveform CSV files need compression after inspection, use lossless gzip and record this separately; those are not needed by the later C++ ion-birth source itself.

## Still not physics-calibrated

The ionic transport uses provisional Ar+ in pure Ar mobility (not measured Ar/CO2 ion species/mobility). Mesh/patch convergence, voltage-independent weighting-field check, alternative ion start positions, prompt electron signal, neighbouring strips, PCB/FR4/gaps, and actual electronics response remain separate acceptance tasks.


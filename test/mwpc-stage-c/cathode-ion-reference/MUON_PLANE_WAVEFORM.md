# Stage C3: complete-muon ion-only response of both WHOLE cathode planes

**Status:** code committed, unbuilt in user WSL as of 2026-11-18. Numerical outcome unknown until user runs it.

## Why this is the next experiment

Previously we found 1.34-1.36% mean ion-level response interpolation L2 in 48 held-out positive-side ion starts, with large tail cases. Further tuning a coordinate metric is not the immediate detector goal. We now want an actual complete Geant4 muon crossing, all its preserved microscopic ion-birth points (weighted per Heed primary) and both time-dependent cathode currents.

Inputs are **already existing**, no microscopic Garfield avalanches should be rerun:

- The `weighted_k512/muon_0000/ion_births.root` sealed muon archive, retaining 3072 weighted macroions representing 554943 microscopic ion births, with actual per-ion creation time, physical start coordinates and per-primary inclusion weight.
- The `polar_kernel_v3/train256_test48_basis.root` from the previous K256 training run (the code reads the stored 256 RAW ion response templates; it does NOT rely on the polar or PCA fit). Those templates have 100-ns bins, 100-us duration, positive ions and two whole-plane currents.
- The `muon0_weighted512_direct_waveforms.csv` produced with direct `DriftLineRKF` for exactly the SAME K512 muon archive.

Method: load training bank; per archived weighted macroion, find four nearest same-side training responses with the previously frozen Cartesian/radial IDW metric; shift each response by the actual `t_birth_ns` with charge-conserving fractional-bin deposition; add both cathodes. Compare full time series with the directly drifted K512 response. Store absolute charge errors, separate plus/minus and combined L2, early cumulative charge errors and peak time. Plot ROOT 2x2 current/cumulative response, overlaying direct RKF in black and kernel in blue. No diffusion, electronics, or prompt electron current.

**Error attribution:** this comparison isolates kernel approximation plus discrete time shifting on the SAME fixed weighted macroion population. It does NOT test sampling versus ALL individual ion births. K256 training and K512 testing use the same original underlying Geant4 muon, and could share callback births, so this is an event-level smoke test, NOT held-out-muon generalization.

## WSL commands: reuse current compact-kernel checkout, NOT frozen runner

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
git status --short
git fetch origin feature/p3-2100v-muon-plane-waveform-20261118:refs/remotes/origin/feature/p3-2100v-muon-plane-waveform-20261118
git switch --track origin/feature/p3-2100v-muon-plane-waveform-20261118
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
cmake --build test/mwpc-avalanche/build-kernel --target mwpc_plane_muon_kernel_compare -j2

BIN=./test/mwpc-avalanche/build-kernel/mwpc_plane_muon_kernel_compare
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
TRAIN="$PWD/test_runs/mwpc_stage_c/polar_kernel_v3/train256_test48_basis.root"
REF="$HOME/na6p/cathode-plane-response-2100/test_runs/mwpc_stage_c/plane_response_pilot_v1/muon0_weighted512_direct_waveforms.csv"
OUT="$PWD/test_runs/mwpc_stage_c/muon_plane_waveform_v1"
mkdir -p "$OUT"
test -s "$TRAIN" && test -s "$REF" && test -s "$DATA/weighted_k512/muon_0000/ion_births.root"

"$BIN" \
  --ion-file "$DATA/weighted_k512/muon_0000/ion_births.root" \
  --training-root "$TRAIN" \
  --reference-waveforms "$REF" \
  --out-prefix "$OUT/jpsi_muon0_K512" \
  --neighbors 4 --dt-ns 100 --window-us 100

cat "$OUT/jpsi_muon0_K512_summary.csv"
ls -lh "$OUT/jpsi_muon0_K512_muon_overlay.png"
~~~

Output: `<prefix>_summary.csv`, `<prefix>_waveforms.csv`, `<prefix>_response.root` and `<prefix>_muon_overlay.png`.

### Interpretation

- All currents are **ion induced**, not measured voltages and not full electron+ion currents.
- The reference is `direct RKF + statistical K512 sample`, NOT full-ion truth.
- The ion mobility is `IonMobility_Ar+_Ar.txt`, provisionally used with Ar/CO2 70:30; absolute drift times are not calibrated.
- The geometry is the analytic nine-wire 2100V P3 with two *whole* cathode planes, NOT real stereo strips.
- Check any kernel charge truncated outside the 100-us window. A nonzero `truncated_kernel...` field must be investigated.
- Compare both cathodes independently, especially the bipolar lower-cathode current, and absolute early cumulative charge (not only signed final charge).

## Project presentation uncertainty table

The living reference table is `docs/mwpc-cathode-response/APPROXIMATION_ERROR_REGISTER.md`. Update that file when actual muon-level comparisons become available; preserve a clearly labeled `not quantified` status for all untested physical approximations.

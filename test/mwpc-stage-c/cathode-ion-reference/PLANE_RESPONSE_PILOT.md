# First 2100-V two-plane response pilot: archived ion births → RKF current

**Status:** C++/ROOT/Garfield code committed but not yet compiled or run in the user's WSL. This is a pilot, not the final validated kernel or electronic response.

## Input and method

Uses existing sealed one-crossing ion birth ROOT archives. The real field is the P3 noFW nine-wire Garfield analytic field (3 mm pitch, 30 um anode, ±3 mm cathodes, 2100 V, zero B). Ion mobility file is PROVISIONAL Ar+ in Ar, not calibrated for Ar/CO2 70:30. Upper/lower weighting electrodes are entire cathode planes, not U/V strips.

Program has two intentionally different modes:

- `--mode direct`: accepts only a `weighted` archive and integrates a full `DriftLineRKF` positive-ion trajectory / time-dependent whole-plane induced current for **every retained macroion**, with its exact inverse-inclusion weight and actual birth time. Checks cathode arrival and complete time window; refuses to silently truncate on the number of macroions. Saves two current histograms in ROOT plus CSV waveforms.
- `--mode charge`: accepts the full archive and streams **every callback birth** but calculates only an integrated-charge *endpoint-side estimate*: assumes an ion ultimately reaches the cathode on its initial sign(v) side; `Q_k=e*[1_on_ending_cathode - phi_w,k(birth)]`. It does NOT drift the ion, does NOT calculate a complete waveform and must be validated against `direct` endpoints in the same geometry.

The fixed-field model is linear. For a true full-ion waveform, and for a precomputed time-kernel, additional work is needed. The kernel is not claimed implemented by this module. The expensive direct current is measured for sampled macroions; check against direct microscopic ion trajectories on a feasible subset and validate a future lookup grid before treating the full-plane muon waveform as a verified high-precision reference. Prompt electron currents are not yet summed; full-source ROOT archives hold a direct whole-plane electron current for first Heed primary only.

## Isolated WSL checkout

```bash
cd ~/na6p
git clone --depth 1 --branch feature/p3-2100v-plane-response-pilot-20261118 \
  https://github.com/StasBren/NA6PRoot.git cathode-plane-response-2100
cd cathode-plane-response-2100
if [ -f ~/na6p/env.sh ]; then source ~/na6p/env.sh; fi
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-plane-response \
  -DCMAKE_PREFIX_PATH="$HOME/na6p/install/garfieldpp"
cmake --build test/mwpc-avalanche/build-plane-response \
  --target mwpc_plane_response_from_ions_root -j2

BIN=./test/mwpc-avalanche/build-plane-response/mwpc_plane_response_from_ions_root
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
OUT="$PWD/test_runs/mwpc_stage_c/plane_response_pilot_v1"
mkdir -p "$OUT"

# Direct weighted macroion waveform for the same muon 0:
"$BIN" --ion-file "$DATA/weighted_k256/muon_0000/ion_births.root" \
  --out-prefix "$OUT/muon0_weighted256_direct" --mobility "$MOBILITY" \
  --mode direct --dt-ns 100 --window-us 100 --max-step-mm 0.05

# Stream ALL 554,943 true ion births for a fast endpoint-side charge estimate.
# This is NOT an all-ion direct waveform.
"$BIN" --ion-file "$DATA/full_k0/muon_0000/ion_births.root" \
  --out-prefix "$OUT/muon0_full_charge" --mobility "$MOBILITY" \
  --mode charge --dt-ns 100 --window-us 100

cat "$OUT/muon0_weighted256_direct_summary.csv"
cat "$OUT/muon0_full_charge_summary.csv"
```

Outputs: `<prefix>_summary.csv`, `<prefix>_response.root` and `<prefix>_waveforms.csv` (latter only for `direct`). Source ROOT archives and existing physics locks are read-only. Refuses to overwrite earlier result prefixes. Resource safeguards: default `--max-direct-ions 20000`, `--window-us` limited to 500 and number of time bins <=10000; these are analysis limits, NOT Garfield avalanche gain caps.

## Next scientific gates

1. Confirm direct RKF cathode endpoints agree with the initial sign(v) hypothesis for *all retained macroions*; otherwise the full charge estimate is not valid.
2. Compare full-archive endpoint-charge sum to sample endpoint-charge sum for K256/512 across multiple muons, including multi-wire/large-angle stress tests.
3. Check dt, RKF max-step and ion mobility sensitivity. Both direct and charge paths currently use the same weighting model; only the direct path actually transports ions.
4. Build a single-ion `K_{plane}(u_0,v_0,t-t_birth)` waveform LUT over physically populated start coordinates, with tested interpolation and numerical errors.
5. With verified kernels construct *both* full-ion and macroion-sampled entire muon waveforms, then replace whole-plane weighting with physical U/V 3D stereo strip weighting.

## Historical 371 MB C2 FEM directory

Old numerical weighting maps are **not used here**. Their README default describes a different 4 mm-pitch/±2.5 mm geometry, but the actual saved `manifest.json` must be inspected before concluding what those 371 MB contain. Do not delete them as a prerequisite for this pilot; retain for stereo-stage provenance or regenerate after recording their hashes/configuration.

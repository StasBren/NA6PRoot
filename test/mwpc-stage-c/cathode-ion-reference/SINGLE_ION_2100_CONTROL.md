# Stage M2a: one positive ion at 2100 V — whole-plane control

**Status: code committed and source-reviewed. It has NOT been compiled/executed with the user's Garfield++ installation.** The outputs must not be interpreted as a validated individual stereo-strip current.

This is an independent 2100-V fork of `phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx`; no historical or cluster simulation is changed. Only the physical anode voltage is made configurable (`--hv`, default 2100 V). The ion drift uses Garfield++ `DriftLineRKF` with the P3 no-field-wire 9-anode analytic field (3 mm pitch, 30 um diameter, cathodes at v = +/-3 mm). The readouts here are exactly the two complete cathode planes, with real wire-aware weighting fields.

## Use a SEPARATE checkout

Do NOT `git pull` in your frozen `~/na6p/cathode-ion-2100-wsl` clone: its previous microscopic results are tied to its Git SHA/config.lock. Use this new worktree:

```bash
git clone --depth 1 --branch feature/p3-2100v-single-ion-control-20261118 \
  https://github.com/StasBren/NA6PRoot.git ~/na6p/cathode-single-ion-2100
cd ~/na6p/cathode-single-ion-2100
[[ ! -f ~/na6p/env.sh ]] || source ~/na6p/env.sh
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-single-ion-2100
cmake --build test/mwpc-avalanche/build-single-ion-2100 \
  --target mwpc_phase_c_p3_noFW_2100_single_ion_control_v1 -j2
BIN=./test/mwpc-avalanche/build-single-ion-2100/mwpc_phase_c_p3_noFW_2100_single_ion_control_v1
"$BIN" --selftest
```

If CMake cannot find Garfield++, supply `-DGarfield_DIR` using the actual location from your prior successful WSL build (`cmake -LA -N ~/na6p/cathode-ion-2100-wsl/test/mwpc-avalanche/build-cathode-2100-wsl | grep Garfield_DIR`).

## Input: exactly one physical positive ion

```bash
OUT="$PWD/test_runs/mwpc_stage_c/single_ion_2100_v1"
mkdir -p "$OUT"
printf 'u_mm,v_mm,w_mm,t_ns,weight\n0,0.030,0,0,1\n' > "$OUT/one_ion_plus.csv"
printf 'u_mm,v_mm,w_mm,t_ns,weight\n0,-0.030,0,0,1\n' > "$OUT/one_ion_minus.csv"

# Find the mobility table actually installed with Garfield++
find "$HOME/na6p" -type f -name 'IonMobility_Ar+_Ar.txt' -print

# Replace this path with the file found by the preceding command
export MWPC_ION_MOBILITY_FILE='/absolute/path/to/IonMobility_Ar+_Ar.txt'

"$BIN" --mobility "$MWPC_ION_MOBILITY_FILE" \
  --input "$OUT/one_ion_plus.csv" --hv 2100 \
  --dt-ns 100 --tmax-us 500 --out-prefix "$OUT/ion_plus"
"$BIN" --mobility "$MWPC_ION_MOBILITY_FILE" \
  --input "$OUT/one_ion_minus.csv" --hv 2100 \
  --dt-ns 100 --tmax-us 500 --out-prefix "$OUT/ion_minus"
head -n 2 "$OUT/ion_plus_ions.csv"
head -n 2 "$OUT/ion_minus_ions.csv"
```

The two CSVs report both `plus_ramo_delta_Q_fC`/`minus_ramo_delta_Q_fC` (weighting-potential endpoints) and `plus_window_Q_fC`/`minus_window_Q_fC` (time-integrated induced currents); they must agree within numerical and finite-window errors after sign conventions are verified. The `_waveforms.csv` files contain the induced-current histories.

**Ion mobility limitation:** `IonMobility_Ar+_Ar.txt` is a provisional Ar+ in pure Ar transport model for a 70:30 Ar/CO2 gas. Its drift-time predictions are not validated mixture predictions. No ions species inference or measured electronics is claimed.

**No physical stereo strips yet:** the old `phaseB_strip_weighting.cxx` is a deliberately simplified plane-strip model and does not include the nine-anode electrostatics and real U/V stereo strip configuration. After the whole-plane Ramo charge control passes, implement and validate the correct two-sided physical strip weighting fields before summing muon clouds.

**Storage:** single-ion outputs are tiny. Keep the previously archived `full_k0/muon_0000/ion_births.root` and `weighted_k256/muon_0000/ion_births.root` for the forthcoming full-vs-weighted avalanche signal comparison.


# Whole-plane complete muon e + ion currents at 2100 V

Status as of 2026-11-18: C++/ROOT/Garfield++ source committed but NOT YET COMPILED OR EXECUTED IN USER WSL. This is an instrumented one-muon pilot, not a final simulation validated against detector electronics.

## What the detector sources actually say

- NA60+/DiCE Proposal `SPSC-P-373.pdf`, Section 2.4.6, specifies VMM/VMM3a readout for the planned muon-spectrometer strips, with on-chamber electronics and FPGA/ROC aggregation. The proposal describes a **25 ns VMM timestamp grouping** downstream; that is a DAQ time grouping, **not** an integration/peaking time or physical detector-current filter.
- `New_chamber_design(1).pdf`, Table 1 and Section 2, lists 2023-26 chamber prototypes, gas Ar/CO2 70:30 and experimental SRS/STAR DAQ technologies. A prototype's DAQ chain need not match the final VMM setup.
- VMM3a User Guide (Iakovidis, CERN ATL-MUON-PUB-2022-002; 2022), Section 1: a charge-sensitive amplifier, third-order DDF shaper with selectable **peaking times 25/50/100/200 ns**, followed by discriminator, peak detector and timing circuit; gain options 0.5/1/3/4.5/6/9/12/16 mV/fC. The 100 ns setting is possible but not confirmed for NA60+/DiCE.
- Thus **100 ns peaking != a 100 ns boxcar integral**. A physically meaningful VMM prediction requires frontend polarity, gain, peaking selection, pulse transfer function, threshold, etc., and strip capacitance/load. The code introduced in this branch does **not** invent those settings.

## Why a microscopic electron replay is necessary

`CathodeIonArchive.hh` saved the full muon ion births, but `electron_endpoints` and direct `electron_control` cover **only the FIRST Heed primary, and only in FULL archive mode**. Weighted K512 mode has no full-muon direct electron current. Reconstructing all prompt electrons from ion births and endpoint *totals* is not physically justified.

The independent `mwpc_p3_2100_full_muon_electron_planes` target copies the 2100V uncapped microscopic Geant4/Heed/Garfield logic into an isolated new source file (not a modification of the historical frozen runner). It enables `AvalancheMicroscopic::EnableSignalCalculation(true)` for EVERY Heed primary, with the same real 2100V nine-wire analytic field and both complete cathode electrodes.

- Input is the exact immutable one-muon Geant4 crossing `weighted_k512/muon_0000/input.csv` and sealed reference `weighted_k512/muon_0000/ion_births.root`.
- Per primary, it checks **original Heed start coordinates/time, Garfield seed, exact ion birth callback count, electron collection count, and success flag**. Any mismatch aborts. A successful replay calculates the full electron current, **without changing the ion archive or recomputing macroion weights**.
- Each primary's 0.05ns-binned 200ns electron current is summed into a new, small ROOT file, with reference identity and physics audit metadata. The 200ns window is checked against electron endpoint times; out-of-window endpoints fail rather than being ignored.
- No positive-ion current is generated in this electron ROOT; existing K512 direct weighted-RKF current remains the separate ion reference. Electron + ion superposition is a linear fixed-field approximation. At high gain, space charge remains unmodeled.

The companion `mwpc_assemble_electron_ion_planes` reads the audited fine electron ROOT, existing **same-muon** K512 direct/fast-ion waveform CSV, and its identity summary. It bins electron current conservatively to the existing **100ns ion grid**, forms direct electron+ion and kernel electron+ion currents and outputs CSV, ROOT and a 2x2 PNG with an early-2us zoom. First-100ns electron and direct-ion cumulative charges are reported separately.

### Important timing limitation

Although microscopic electrons are retained at **0.05ns** in the independent ROOT archive, the combined plotting grid uses the existing **100ns ion bins**. Those coarse ion bins already show a strong first-bin kernel discrepancy; they cannot resolve a 25/50ns VMM shaper or reliably predict a 100ns peaking detector. After verifying the complete raw current, the next timing step is **recomputing the direct K512 ion current on a finer grid (for example 10ns) and using a finer-grained ion kernel**, before any electronics shaping study.

## Reproducible laptop instructions

Use the existing non-frozen compact kernel checkout, fetch the new child branch, and build both independent targets:

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
git status --short
git fetch origin feature/p3-2100v-complete-electron-plane-20261118:refs/remotes/origin/feature/p3-2100v-complete-electron-plane-20261118
git switch --track origin/feature/p3-2100v-complete-electron-plane-20261118
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-kernel -DCMAKE_PREFIX_PATH="$HOME/na6p/install/garfieldpp"
cmake --build test/mwpc-avalanche/build-kernel --target mwpc_p3_2100_full_muon_electron_planes mwpc_assemble_electron_ion_planes -j2
~~~

Now replay the full microscopic muon for its prompt electron current only. The source ion ROOT file is read-only:

~~~bash
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
OUT="$PWD/test_runs/mwpc_stage_c/whole_muon_electron_plus_ion_v1"
mkdir -p "$OUT"
ELEC=./test/mwpc-avalanche/build-kernel/mwpc_p3_2100_full_muon_electron_planes
"$ELEC" --whole-muon-electron-signal \
  "$DATA/weighted_k512/muon_0000/input.csv" \
  "$OUT/electron_replay_muon0" \
  "$DATA/weighted_k512/muon_0000/ion_births.root" \
  "$OUT/muon0_direct_electron.root"
~~~

Check for `PASS audited COMPLETE-MUON direct ELECTRON-only plane currents`; the printed per-primary audit must cover all Heed primaries. If microscopic counts differ, stop and inspect the failure. The ROOT file is not created until all microscopic primary audits pass.

Then assemble the raw currents and plots:

~~~bash
MERGE=./test/mwpc-avalanche/build-kernel/mwpc_assemble_electron_ion_planes
ION="$PWD/test_runs/mwpc_stage_c/muon_plane_waveform_v1/jpsi_muon0_K512"
"$MERGE" \
  --electron-root "$OUT/muon0_direct_electron.root" \
  --ion-waveforms "${ION}_waveforms.csv" \
  --ion-summary "${ION}_summary.csv" \
  --out-prefix "$OUT/muon0_complete_planes_coarse100ns"
cat "$OUT/muon0_complete_planes_coarse100ns_summary.csv"
ls -lh "$OUT/muon0_complete_planes_coarse100ns_electron_ion_planes.png"
~~~

OUTPUT: `muon0_direct_electron.root` fine 0.05ns electron and metadata; coarse combined `<prefix>_waveforms.csv`, `<prefix>_summary.csv`, `<prefix>_response.root`, `<prefix>_electron_ion_planes.png`. Parent data and already validated ion waveforms are not overwritten.

## How to report uncertainty

Update `docs/mwpc-cathode-response/APPROXIMATION_ERROR_REGISTER.md` after observing actual measured numbers. Distinguish:

- **Muon source/census fidelity**: exact Geant4 crossing and all Heed primary Garfield counts; fail on mismatch.
- **Microscopic electron current**: direct Garfield model, not sampled or inferred. It is still only the analytic nine-wire / whole-plane geometry.
- **Positive-ion sampling**: K512 weighted macroions versus all 554943 births remains a statistical approximation.
- **Kernel interpolation**: direct versus fast-ion same-reservoir disagreement, concentrated in the first 100ns for the previous run.
- **Signal combination**: fixed electric field + linear Shockley-Ramo superposition, no ion-cloud space-charge correction.
- **Electronics**: **not modeled**; 100ns is one possible VMM3a peaking setting, not a readout box/window.

## Reference-vs-hardware geometry caveat

The inspected `New_chamber_design(1).pdf` Table 1 lists the real **2026 Prototype-3** at **4 mm pitch, 2+4 mm gaps, 1.8 kV, one-side SRS readout**. The present *frozen* internal `p3_noFW` simulation is **3 mm pitch, 3+3 mm gaps, 2.1 kV, two ideal whole-plane cathodes**. These are NOT the same detector. We preserve internal apples-to-apples comparisons now, and will align the simulation to target as-built/final geometry before detector-performance claims. This is an important unquantified physical systematic independent of kernel or macroion sampling error.

## Runtime environment fix — TrackHeed database (2026-11-18)

A WSL run stopped at `TrackHeed::Initialise: Cannot retrieve database path`, before avalanche replay. This is a missing environment setup, not a physical failure. Garfield++ explicitly documents `HEED_DATABASE` as the installed `share/Heed/database` directory. Preflight before launching the expensive uncapped microscopic replay:

~~~bash
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export HEED_DATABASE="$GARFIELD_INSTALL/share/Heed/database"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
if [ ! -f "$HEED_DATABASE/cbdel.dat" ]; then
  echo "Heed database missing at $HEED_DATABASE"
  find "$HOME/na6p/install/garfieldpp" "$HOME/na6p/src/garfieldpp" \
    -type f -name cbdel.dat -print 2>/dev/null
  echo "Do not launch the replay until HEED_DATABASE is fixed."
fi
echo "HEED_DATABASE=$HEED_DATABASE"
~~~

If the `find` lists a different directory, point `HEED_DATABASE` to the **parent directory of cbdel.dat**. Because the failed run exited before replay, the audited electron ROOT should not exist; if it does, inspect provenance before running again. The replay source invokes `AvalancheMicroscopic::EnableSignalCalculation(true)` for *every* Heed primary and reads `Sensor::GetElectronSignal`; this comprises both the initial electron drift AND all multiplication (secondary) electrons in each microscopic avalanche. The separate positive-ion birth archive contains the ionization-callback ions; it does not include the original positive ions left behind by Heed primary ionization, an additional small omitted component to register.

## 2026-11-18 uploaded full-muon result and ROOT plotting repair

The uploaded same-muon electron+ion summary shows all 7 Heed primaries, 554943 original avalanche ion births, 555139 collected electrons. Both whole-cathode induced-current CSVs have finite nonzero values. Full-window combined electron + K512 direct-RKF ion charges are +88.700555 fC and +0.271157 fC. Electron endpoints finish by 40.6543ns and all their current is in the first 100ns coarse ion bin. The total first-100ns direct ion+electron charges are +10.909134/+9.453484 fC; the simple fast ion template yields 18.21/19.56% less during that bin. The underprediction matters for later shaping and cannot be declared harmless.

The first `*_electron_ion_planes.png` was BLANK except legends because `MakePlot` constructed TH1D histograms on the stack inside a loop and destroyed them before TCanvas::SaveAs. Fixed by keeping `std::unique_ptr<TH1D>` histograms alive and calling canvas.Update before SaveAs; all existing physics outputs remain valid.

Only rerun the **assembler** with a FRESH output prefix. The microscopic electron replay and weighted-ion RKF simulations do not need to be rerun:

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
git fetch origin feature/p3-2100v-complete-electron-plane-20261118
git merge --ff-only FETCH_HEAD
cmake --build test/mwpc-avalanche/build-kernel --target mwpc_assemble_electron_ion_planes -j2
OUT="$PWD/test_runs/mwpc_stage_c/whole_muon_electron_plus_ion_v1"
ION="$PWD/test_runs/mwpc_stage_c/muon_plane_waveform_v1/jpsi_muon0_K512"
MERGE=./test/mwpc-avalanche/build-kernel/mwpc_assemble_electron_ion_planes
"$MERGE" --electron-root "$OUT/muon0_direct_electron_retry1.root" \
  --ion-waveforms "${ION}_waveforms.csv" --ion-summary "${ION}_summary.csv" \
  --out-prefix "$OUT/muon0_complete_planes_plotfix_v2"
ls -lh "$OUT/muon0_complete_planes_plotfix_v2_electron_ion_planes.png"
~~~

Next physics timing gate: preserve the 0.05ns direct electron ROOT and regenerate direct K512 ion current at a finer (e.g. 10ns) time grid, then compare against similarly fine kernel/fast-ion response. **Do not invent a VMM3a shaper response from 100ns box bins.**


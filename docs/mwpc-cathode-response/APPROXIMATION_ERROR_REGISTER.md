# NA60+ MWPC cathode-response approximation and uncertainty register

**Reference for final presentation.** Version 1 — 2026-11-18. The table is the authoritative project-wide record of what the reduced detector-response model omits, how the approximation may bias observables, and how well that bias has actually been quantified. Update the evidence column only after a corresponding reproducible numerical or measured comparison. In the final presentation, include an up-to-date form of this table; do not present pilot numerical errors as detector-systematic uncertainties.

| Simplification | Consequence | What we actually know | Status / next evidence |
|---|---|---|---|
| Weighted macroions rather than individually transporting all ions | Sampling noise in waveform shapes and rare ion-arrival tails; signed bipolar currents may cancel and amplify the relative uncertainty of small net charges | For muon 0 at 2100 V, K256 vs K512 direct weighted RKF: upper plane 87.598630 vs 87.608213 fC, lower plane -0.832667 vs -0.743271 fC. These are two sampled versions of the SAME microscopic avalanche, not two independent muons. | Event-level waveform comparisons; compare against complete-ion reference or converged charge-basis, and test repeated reservoirs |
| Interpolated single-ion response kernel instead of direct RKF | Numerical interpolation error in the time-dependent induced current, including initial charge rise and bipolar tail | On 48 HELD-OUT positive-side ion positions from one muon, 256 training templates yielded Cartesian raw-IDW mean L2 1.364%, p90 2.837%, worst 8.314%; polar raw-IDW mean 1.339%, p90 2.526%, worst 8.889%. These are **single-ion** errors only; event-level error unknown. | Compare whole-muon kernel waveform with direct weighted RKF for the EXACT SAME reservoir |
| Nine-wire analytic electric field and full-plane weighting | Excludes finite PCB, guard traces, edge fields, mechanical support and actual segmented readout. Whole-plane response cannot predict stereo strip sharing. | No physical-geometry systematic quantified. Historical FEM `c2_weighting_smoke` may have a DIFFERENT wire pitch and cathode gaps; its files are not a 2100V P3 validation by themselves. | Obtain physical strip layout and validate 3D FEM, mesh and patch-size convergence |
| Pure-argon Ar+ mobility file used with Ar/CO2 70:30 | Ion drift speed and signal duration may be biased; the relevant ion species/mobility in the mixture is not established by this input file | Uncalibrated ion-transport times: sample latest arrivals 61.33 us (K256) and 82.41 us (K512) but these differ also because of sampling. Cannot interpret differences as gas-model error. | Find/measure mixture-appropriate positive-ion transport; sensitivity scan of waveform duration |
| Fixed external electric field, no avalanche space charge | At large gain, ion cloud can modify drift/avalanche field and signal | Unquantified at 2100 V; successful uncapped microscopic avalanche transport does not by itself prove space charge negligible. | Quantify ion density/gain and assess threshold for field perturbation |
| Positive-ion signal only | Omits prompt electron-induced current, especially the leading edge, and possibly early readout timing | Existing full ROOT reference contains direct electron current only for FIRST Heed primary. The new **independent all-primary full-muon electron replay** and electron+ion assembly are implemented but NOT YET executed; no numerical electron fraction can be claimed. | Run the exact Garfield seed/callback/electron census replay against the frozen K512 muon; then compare 100ns-bin electron, ion and total contributions |
| No electronics response | Raw induced current is NOT shaped voltage, threshold or measured hit time | Electronics shaping/noise/threshold/digitization parameters are not yet validated; measured-response error unknown. | Obtain front-end transfer function, measured VMM/electronics settings; then convolve and digitize |

## Readout and error bookkeeping

1. Always label two **whole cathode planes** separately from individual stereo U/V strips.
2. Report `kernel error | same exact archived macroion reservoir` separately from `reservoir sampling error | same microscopic avalanche` and `physical model error | different field/mobility/electronics`. Do not add unrelated percentages as if they measured the same effect.
3. For comparisons, retain waveform L2 on each cathode, peak magnitude/time, early cumulative charge, bipolar sign crossing, total integrated charge, and recording-window completeness. Use absolute errors when a signed quantity crosses zero.
4. Pilot ion-response templates are trained on a single positive-side muon. Models applied to a separate reservoir of the SAME muon are useful event-level smoke tests but not independent scientific generalization.
5. Current reference: 2100 V, P3 noFW, 3 mm anode pitch, 30 um wire diameter, cathode gaps ±3 mm, Ar/CO2 70:30 760 Torr and 293.15 K, B=0, provisional Ar+ in pure Ar mobility.

## Presentation instruction

Place the updated table on a dedicated slide titled **Model simplifications and known uncertainties**, using at most a short evidence line per row. If space is insufficient, split into (A) numerical approximations and (B) physical/detector omissions. Retain an explicit `not yet quantified` status; do not replace it with guessed percentages. Link slide notes to the corresponding run manifests and analysis comparisons.

### Symmetry-based response extension (2026-11-18)

The K256 single-ion kernel bank contained only positive-v starts; K512 has at least one ion with an opposite-side start. For the existing analytically symmetric P3 noFW nine-wire / ±3 mm whole-cathode model and B=0, reflection `v -> -v` and interchange of the two cathode readout channels is an exact **model symmetry** (no extra approximation within this idealized geometry). We apply it to the single-hemisphere template bank; it is not valid for asymmetric cathodes, magnetic-field effects or individual U/V stereo strips. The event-level K512 comparison against direct Garfield RKF remains the numerical validation gate.

### Electronics source check (2026-11-18)

- The attached NA60+/DiCE proposal `SPSC-P-373.pdf`, Section 2.4.6, specifies VMM-based strip electronics. Its 25-ns time-stamp grouping describes the DAQ packet clock, not analog integration.
- The attached `New_chamber_design(1).pdf`, Section 2 / Table 1, gives SRS/STAR for 2023-26 tested prototypes. That does not establish a final VMM shaping configuration.
- The VMM3a User Guide by G. Iakovidis (2022, Section 1) states third-order DDF shaping with selectable 25/50/100/200-ns peaking times, eight gain settings and peak/timing detectors. **100-ns peaking time is not a 100-ns boxcar integration**. No final NA60+/DiCE VMM peaking/gain/polarity/threshold setting was found in the inspected project documents.
- The existing 100-ns-ion-bin waveform shows first-bin interpolation discrepancy. Before convolving with a 100-ns-peaking circuit, we need a finer ion-time grid plus the direct prompt-electron waveform; avoid converting raw currents directly into peak/threshold predictions.

### Explicit distinction between numerical P3-noFW reference and actual Prototype-3

The attached 2026 chamber-design note, Table 1, lists its tested **Prototype-3 (2026)** as **30 um wires, 4 mm wire pitch, asymmetric 2+4 mm gas gap, 1.8 kV, single-sided SRS readout** (1.2 mm strip dimension). Our frozen numerical reference is instead **30 um wires, 3 mm pitch, 3+3 mm symmetric gaps, 2.1 kV and two ideal whole-cathode readouts**. The source file naming uses `p3_noFW`, but this is a deliberate simplified reference configuration **NOT the as-built Prototype-3 geometry**. This mismatch is a separate unquantified geometry/operation systematic, not an interpolation error. Keep the frozen 2100V configuration for internal model comparisons; agree the eventual target geometry with the detector team before claiming as-built electronics or strip response.

### Small primary-ion component not yet included

The archived microscopic Garfield++ ion callback counts ion-electron pairs created **during avalanches**. TrackHeed's original muon ionizations also leave one positive ion per original primary electron in the gas. Those **Heed primary positive ions are not contained in the avalanche-ion archive and are currently absent from the cathode ion current**. Their count is likely much smaller than the amplified avalanche ion population at 2100 V, but this contribution has not been calculated; do not call the electron+avalanche-ion sum fully charge-complete until it has been assessed.

### Full-muon electron+avalanche-ion plane response — first actual result (2026-11-18)

Uploaded `muon0_complete_planes_coarse100ns_summary.csv` and waveform CSV correspond to the same crossing `jpsi40-bore|1|8953|3|1|1|39`, 7 Heed primaries, 554943 avalanche-ion births, 555139 collected avalanche electrons. Replayed microscopic electron endpoints end by 40.6543 ns and 0.05ns electron samples are conservatively rebinned into existing 100ns ion samples.

- Integrated microscopic electron signals: +1.092342 fC (+ plane), +1.014428 fC (- plane).
- Direct weighted RKF ion responses: +87.608213 fC (+), -0.743271 fC (-).
- Combined electron + direct weighted RKF ions: +88.700555 fC (+), +0.271157 fC (-). These are two ideal complete planes with provisional ion mobility, not shaped voltage or real stereo strips.
- In the first 100 ns: electron +1.092342/+1.014428 fC, ion +9.816792/+8.439056 fC, total +10.909134/+9.453484 fC (plus/minus). Electrons are 10.0% and 10.7% of the first-bin total charges. Current fast-ion template underestimates first-bin totals by 18.21% (+) and 19.56% (-). This early-time kernel discrepancy is NOT negligible for a VMM3a shaping study and must be addressed with finer ion bins.
- Total-ion fast-kernel vs direct-RKF integrated-charge error remains +0.024465/+0.023445 fC; small final-charge error does NOT prove first-bin waveform accuracy.
- The initially exported ROOT PNG was blank because histogram objects were stack-local and destroyed before rendering. The ROOT plot lifetime bug is now corrected in C++ source; this affected visualization, not the archived CSV currents.

**Physical scope:** only avalanche ions (not original Heed positive ions), unsegmented two-plane analytic field, no VMM3a frontend, and only one muon event. These numbers must NOT be presented as detector-wide systematic uncertainties.


**Attribution warning:** The fast event response starts from 100ns-bin *integrated* single-ion templates and uses a fractional-bin time shift for each true ion birth. Therefore the measured 18–20% first-bin discrepancy combines **spatial template interpolation and temporal discretization/shift error**; it is NOT justified to attribute all of it to spatial kernel physics. In fact, the uploaded waveform shows the upper-plane first-bin residual -1.98664 fC, second-bin residual +1.77446 fC (net -0.21217 fC after 200 ns). The strong compensation suggests resampling/shift sensitivity. Distinguish these errors with fine-bin single-ion impulses and a separate like-for-like rebinned RKF test before changing the kernel model or predicting VMM peaking.

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

### Controlled ion-current time grid test (run on WSL 2026-11-20; first muon K512)

Previous first-100ns kernel-minus-direct K512 full-muon ion response was approximately -1.987 fC on the upper plane, while the ~100us integrated ion charge differed by only +0.0245 fC. A 100ns-bin impulse library shifted by fractional ion birth times can distort the early transient even if total charge is conserved; spatial interpolant error is not isolated by this observation.

Stage C 10ns test in `test/mwpc-stage-c/cathode-ion-reference/PLANE_TIMING_10NS_20261120.md` holds fixed (i) microscopic muon crossing, (ii) original K512 sampled ion births and statistical weights, (iii) the same original 256 RKF template start positions, (iv) nine-wire analytic field/ion mobility/geometry. It recomputes direct and fast responses at 10ns and uses charge-conserving rebinning to compare against the old 100ns currents. Separately records old-vs-fine direct reference numerical error and old-vs-fine fast interpolation error for early 100/200/400ns windows, and optionally adds the already-audited full-muon electron current for context. The WSL run **completed successfully** for the original K512 benchmark; quantitative results follow. No VMM shaper, noise, realistic strips or hardware calibration are included.

**Observed results, exact same ion reservoir and training sites:**

| Quantity | Plus plane | Minus plane |
|---|---:|---:|
| Direct RKF first 100ns (both grids, fC) | +9.81679168838 | +8.43905632153 |
| Old fast first 100ns (fC) | +7.83015489848 | +6.58952528012 |
| 10ns fast integrated first 100ns (fC) | +9.83697554425 | +8.45397206411 |
| Old fast minus direct first 100ns (fC) | -1.98663678990 | -1.84953104140 |
| 10ns fast minus direct first 100ns (fC) | +0.02018385587 | +0.01491574258 |
| Old relative early integrated error | -20.24% | -21.92% |
| New relative early integrated error | +0.206% | +0.177% |
| Old current-waveform L2 relative to direct, native 100ns bins | 22.88% | 27.37% |
| New current-waveform L2 relative to direct, native 10ns bins | 2.943% | 3.447% |
| New response *rebinned to 100ns* L2 relative to direct rebinned to 100ns | 0.226% | 0.179% |
| Old direct vs 10ns direct rebinned to 100ns, waveform L2 | ~1.3e-15 | ~1.1e-15 |
| Final 100us fast-vs-direct induced-charge bias (fC; unchanged) | +0.02446462114 | +0.02344520447 |

Including the **same** previously calculated microscopic full-muon electron current, the first-100ns combined direct electron+ion charges are +10.909133967 fC / +9.453483978 fC, and fast10 electron+ion totals are +10.929317823 fC / +9.468399721 fC (+0.185% / +0.158% errors). Electrons are common to both, so this comparison tests ion-kernel error only.

**Interpretation:** 100ns direct RKF exactly matches charge-conserving rebin of 10ns direct RKF to numerical precision; the 100ns reference numerical grid was not the cause. Time-shifting *already averaged* 100ns single-ion impulse templates produced the dominant spurious first-bin error. Using 10ns impulse sampling and shifts nearly eliminates the early **integrated charge** bias for this same K512 muon while keeping all spatial interpolation/sampling unchanged. This is an empirical attribution to coarse temporal response handling, not a general proof across muons or real detector geometry. The 2.94–3.45% **native 10ns pointwise waveform L2 error** is considerably larger than the ~0.18–0.23% rebinned 100ns L2: bin averaging hides high-frequency mismatches. For ~25ns-peaking VMM3a, do not equate integrated-charge agreement with shaped peak amplitude/threshold timing precision.

**Still unknown:** independent-muon response-kernel generalization; actual strip weighting-field response, ion mobility in the operating gas, model-vs-hardware systematics, stochastic ion-sampling variance, and actual VMM3a shaping/noise. No need to increase the template count or add a neural network on the basis of the old first-100ns error. Preserve these distinctions in group-meeting slides.

### Independent-muon validation and candidate two-sided stereo strip geometry (initiated 2026-11-20; subsequently completed below)

Three distinct ROOT-audited K256 muons from the first-571-event snapshot were chosen for held-out event-level kernel validation: phi source index 301 (station 0; 4096 saved ions), jpsi source index 408 (station 3; 7424 ions), omega source index 446 (station 2; 5376 ions). Their frozen archived K256 birth positions and per-primary weights will each be drifted directly using 10ns RKF and also evaluated by the **same original 256 site** four-neighbor spatial interpolant. For each event, report native 10ns waveform L2, same data rebinned to 100ns L2 and early 50/100/200/400ns induced-charge errors on both ideal planes. Do not confuse these with ensemble-wide systematics, physical electrode behavior, or sampling-vs-full truth. **Historical pre-run note; the test was subsequently executed successfully. See measured independent-muon results at the end of this register.** C++/ROOT and no-overwrite execution script in test/mwpc-stage-c/cathode-ion-reference/INDEPENDENT_3MUON_STEREO_READOUT_PILOT.md.

The preferred future detector geometry supplied for study and described in New_chamber_design(1).pdf uses simpler two-sided readout strips of width ~1.7 mm **measured along the wires**, crossing at pi/2 +/- alpha with tan(alpha)=0.1. Symmetric 3+3mm is a **design candidate, not finalized**; 2.5+2.5mm is an alternate in the document. Real copper gaps/pitch are unspecified. A new analytic wire-free two-plate stereo-strip weighting-potential geometry pilot assumes pitch=width=1.7mm (ideal contiguous strips), with correct opposed inclinations, v-mirror/u-mirror symmetry, conductor boundary values and analytic Laplace tests. **It omits all anode wires**, PCB dielectric, strip copper gaps/traces and finite edges, therefore it is NOT the physical 3D MWPC strip weighting field and does not provide valid real strip currents or charge sharing. The required next stage is a 3D wire-aware BEM/FEM/neBEM weighting solution, with all other conductors including anodes grounded, mesh convergence, boundary and sum-rule checks. No electronics calculations are authorized or included.

### Independent muon 0301 RKF guard stop (2026-11-20; historical diagnostic, resolved below)

Separate-worktree WSL run of first selected held-out muon muon_0301 passed immutable input integrity preflight with indices 0301, 0408, 0446 fully sealed and audited. During the weighted direct ion DriftLineRKF pass for 0301, executable aborted with **ion drift failure or time window truncated**, before producing valid direct/fast result CSVs. **This was the initial failure; the controlled finite-window treatment was later implemented and all three tests completed successfully.** The guard combines transport boolean, endpoint time outside the 100000 ns acquisition window, and possible invalid/no-progress endpoint. **The exact cause and original sampled ion row were not printed by old code.** Do not infer whether ion drift physically failed, whether it exceeded the capture time, or whether an alternate status (e.g. electrode/gas boundary) occurred. No changing acceptance, dropping failed ions, changing mobility or waveforms, or accepting partial results without separate physics justification. Diagnostic additions on the current branch report event, K256 ion row, Heed primary, original birth coordinates/time and weight, Garfield drift return, endpoint code/coords/time, and recording-window bound. Rerun first event in *diagnostic mode only*, then resolve policy on evidence. **The other two held-out muons were subsequently tested successfully.** Independently, the candidate **wire-free** stereo geometry self-test and PNG/ROOT maps passed; no claim of physical wire-aware strip currents.

### Root cause confirmed: first independent-muon ion arrives AFTER the 100us recording window

Garfield diagnostic for archived phi muon 0301 weighted ion row 257 (Heed primary 2, sample weight 225.9296875) reports successful DriftIon `ok=1`, boundary status -1, valid +3 mm cathode endpoint at time 131747.107ns, but fixed Sensor time window ends at 100000ns. Thus the original guard failed due to **late physical cathode arrival**, not failed transport. This is possible for near-anode-symmetry-plane ions with long curved drift trajectories. The original first benchmark happened not to expose this condition.

Correction is **opt-in** via --record-late-ion-endpoints for the specific 100us windowed independent-muon study: do not exclude late ions or artificially change their time, instead integrate their valid induced current *inside* 0-100us and separately preserve exact full-electrode endpoint charge and the unrecorded charge tail. Log late ion rows, weighted late-ion population, largest late arrival, and both cathodes' missing tails. All ions must still transport successfully and end on the physically correct cathode; complete-in-window ions still require original per-ion integrated-vs-Ramo closure. Only the over-restrictive assumption that endpoint arrival precedes the observation cutoff has changed. Whole-event 100us charge, in presence of a late tail, is not equal to final infinite-time charge. **Historical statement superseded: all three independent events were completed on WSL on 2026-11-20; see the measured table below.**

### Measured held-out 3-muon ion-kernel generalization (2026-11-20, WSL success)

From user-supplied, successfully completed `independent3_windowed.log`, `independent3_events.csv` and `independent3_windows.csv`: all three SHA256/ROOT-audited K256 microscopic muon archives passed, direct ion RKF ran to physical cathode endpoints, and four-nearest-neighbor interpolation of the exact same sampled ion reservoir using the original frozen 256-site 10 ns kernel produced waveforms for both ideal complete cathode planes. Combined three-muon two-plane **native 10 ns waveform relative L2 = 0.04810785 (4.81%)**. These results are within-model, NOT comparisons to hardware or independent avalanche sampling.

| Muon (source/station) | Upper 100ns charge bias | Lower 100ns charge bias | Upper 10ns L2 | Lower 10ns L2 | Upper final 100us charge bias | Lower final 100us charge bias | Late ion rows / K256 rows | Max late arrival |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0301 φ / 0 | +1.9822% | +0.8664% | 6.9928% | 5.2827% | +0.9630% | +2.1883% | 328/4096 | 191.078 μs |
| 0408 J/ψ / 3 | +0.5243% | +1.8501% | 2.2522% | 4.4945% | +1.3531% | +0.4477% | 248/7424 | 188.779 μs |
| 0446 ω / 2 | +2.3709% | +2.6553% | 4.1850% | 5.1790% | +1.3945% | +1.3992% | 359/5376 | 193.223 μs |

Relative L2 after rebinning the fine 10 ns currents into 100 ns bins, by plus/minus side: 0301 5.4331% / 2.2285%; 0408 1.3538% / 4.5605%; 0446 3.8200% / 5.0681%. The first benchmark's ~0.2% first-100ns charge errors and ~3% 10ns L2 were therefore **not universal**. No independent test was used to retrain/tune the 256 template sites.

**Late trajectories, finite 100us acquisition window:** weighted ion population arriving after 100us was 67204.171875 / 886727 = 7.579% for 0301; 32720.953125 / 1013413 = 3.229% for 0408; 57195.058594 / 871236 = 6.565% for 0446. The positive cathode unrecorded post-window induced-charge tails were +0.7129673, +0.2518651, +0.6228628 fC, respectively; negative cathode tails were +0.6315717, +0.3642773, +0.5573880 fC. With these valid late-arriving ions, direct-RKF full-endpoint charge and 100us-window recorded charge are **not equal**; the 100us fast/direct comparisons are windowed induced-charge observables, not complete charge endpoints.

**Scientific decisions:** this is evidence that a fixed 256-site plane response kernel has a useful degree of event-level generalization to three other muons (three sources, stations 0/3/2), but errors of 2–7% in the 10ns current waveforms and up to 2.66% in 100ns induced charge remain. It is insufficient to declare final 256-template adequacy across all regions/strips, to claim ±0.2% global response accuracy, or to establish event-level uncertainty from K256 reservoir sampling. The proposed stereo readout still requires a full three-dimensional wire-aware weighting field including anodes, electrode segmentation and actual strip copper gaps. No electronics calculations are authorized in this phase.

### Stage C candidate 3D U/V weighting-field solver with grounded anodes (implemented 2026-11-20, WSL tests pending)

Physics motivation: the earlier ideal stereo strip Laplace boundary test is **wire-free** and thus cannot describe MWPC strip induced signals near anode avalanches. New `WireAwareStereoWeighting.hh` implements a 3D electrostatic **thin-wire approximation**: add to the exact two-infinite-plate strip Dirichlet solution a sum of line-element sources along each of nine conducting anodes. The sources use the Green function for the region between both cathode planes with homogeneous Dirichlet boundaries; reusable dense LU coefficients are chosen to ground the **four-angle mean** of potential on each 30μm wire cylinder segment. Exact conductor-surface potential is approximated, not assumed; the additional 12-angle residual is measured. This uses the contemplated ±tan(α)=0.1, width 1.7mm measured along anode wires, symmetric ±3mm as a *candidate*, plus provisional pitch=width contiguous strips. No PCB dielectric, copper gaps, trace patterns, guard structures, finite readout edges or self-consistent microscopic fields. No induced-current calculations or electronics.

New independent diagnostics at `wire_aware_stereo_strip_pilot.cc` compute both upper adjacent strips, one lower strip, the screened whole-plane weighting potentials, selected strip grounded-wire residual, stereo mirror check, sum-of-strips vs screened whole plane and comparison with **Garfield::ComponentAnalyticField::WeightingPotential** for the same nine-anode analytic reference. With `--convergence`, independently halves longitudinal charge-segment length and expands correction support by ±4mm, reporting absolute potential deviations. Need test image-sum truncation and full cylinder surface BEM, and impose quantitative tolerance to Garfield. New C++ source **not yet compiled or executed on user WSL**. A separate Python numerical derivation sanity check with the same 9×20 source elements produced max midpatch wire-cylinder potential residual ~0.0038 (dimensionless) and sum-of-23-strip error vs screened whole plane ~4.6×10^-6 at one point; these are not execution results for the committed C++ solver or general performance claims.

Design source: New_chamber_design(1).pdf sections 4.2, 5 and 6 states 1.7mm and tan(α)=0.1, symmetric 3+3mm contemplated and 2.5+2.5mm option, and that strips on opposed cathodes are preferable to a complex interlaced one-sided layout. See `test/mwpc-stage-c/cathode-ion-reference/NA60_STEREO_3D_WIRE_FIELD_STAGE1.md`. Do not claim a validated physical stereo-strip field until the independent whole-plane Garfield control, wire-surface residuals and longitudinal/image convergence are observed and appropriately accepted on the WSL installation.

### 3D wire-aware thin-wire stereo-strip field: first WSL numerical result (2026-11-20)

User compiled `mwpc_wire_aware_stereo_strip_pilot` successfully and executed the candidate **3+3mm**, **1.7mm strip width and provisional strip pitch**, **tan(α)=0.1**, nine **30μm diameter** anodes of 3mm pitch. Initial M=12 Dirichlet image cutoff, ±10mm longitudinal 180-unknown line-element solve, and optional segment/extent/image convergence all completed with program `PASS`.

| Diagnostic | Measured result |
|---|---:|
| Max additional 12-angle wire-surface potential residual on U / V / whole plate | 0.003788416 / 0.003789194 / 0.002538913 |
| U↔V mirrored stereo potential discrepancy at selected locations | 6.533331e-5 |
| 23 contiguous upper-strip potentials minus screened whole cathode | 5.026515e-6 |
| Interior finite-difference 3D Laplace residual | 3.258118e-5 mm^-2 |
| Max abs 3D screened whole-cathode potential difference vs Garfield analytic 9-wire whole-plane at six selected interior points | **8.872376e-4** |
| Max change in monitored potentials: w segmentation 1→0.5mm | 1.592416e-4 |
| Max change in monitored potentials: w source half-span 10→14mm | 5.785199e-4 |
| Max change in monitored potentials: Green image cutoff M=12→20 | 6.713097e-4 |

Interpretation: this is promising preliminary interior *weighting-potential* agreement and consistent U/V geometry. Do NOT present it as a <0.1% calibrated strip-current error: Garfield comparison sampled only six bulk points, not a near-anode sweep or full gradient; the anode surface residual of ~0.0038 dimensionless can be comparable to the entire local strip weighting potential near a wire. Gradient/detailed induced currents are much more sensitive than potential agreement at bulk sample points. The strip partition is primarily an INTERNAL consistency check; it is not an independent physical reference. Absolute differences also do not imply small relative errors for channels with weak induced response. Actual strip pitches/copper gaps, PCB dielectrics, edge guards, module finite extent and ion mobility remain unmodeled.

Follow-up code now offers `--near-wire`: scans weighting potentials and finite-difference gradients in the region 1.1–20 anode radii from wire centers, across multiple wires/angles/w locations, comparing the **screened whole-plane fields** against independent Garfield analytic weighting potentials. It writes `_near_wire.csv`, max/RMS potential differences and a global vector-gradient L2; additionally outputs the candidate selected-strip potential at the same near-anode positions. This is diagnostic, not an automatic accuracy certification, and has NOT been run on user WSL yet. Do not compute real strip currents from this thin-wire solver before near-anode field checks and refinement criteria are assessed.

### Measured near-anode screened whole-cathode field versus Garfield reference (CSV 2026-11-20)

User supplied `wire_aware_nearwire_v1_near_wire.csv`, 240 rows and 15 columns, no NaN entries. This is an *independent whole-cathode* Garfield comparison (not a reference for the stereo-strip field): three anodes (-1,0,+1), radial shells r/R = 1.1,1.5,2,5,20 (R=0.015mm), four polar angles, w=0/1.1mm and upper/lower cathodes.

Combined vector weighting-gradient relative L2 = **0.001109 (~0.1109%)**; maximum pointwise gradient-relative error = **0.0054795 (~0.548%)**, at r=20R, on the lower cathode where the reference gradient norm is comparatively small. Median pointwise relative gradient error = 0.001133 (~0.113%). RMS weighting potential difference = 0.000151, maximum absolute = 0.000385, max pointwise relative potential error = 0.00905 (~0.905%) in a tiny-potential region near r=1.1R on the lower cathode.

| Radial shell | Max absolute potential difference | Shell gradient relative L2 | Max pointwise relative gradient error |
|---|---:|---:|---:|
| 1.1R (16.5μm center radius) | 4.24e-5 | 0.1105% | 0.1359% |
| 1.5R (22.5μm) | 6.62e-5 | 0.1107% | 0.1401% |
| 2R (30μm) | 8.89e-5 | 0.1111% | 0.1456% |
| 5R (75μm) | 1.832e-4 | 0.1147% | 0.1825% |
| 20R (300μm) | 3.849e-4 | 0.1458% | 0.5480% |

These validate that *screened whole-plane* gradients match Garfield's independent analytic weighting solution well at these sampled locations. No divergent near-anode relative-gradient error is apparent in these particular shells. **CRITICAL LIMITATION:** the uploaded `selected_strip0_potential` column has NO independent Garfield (or other) *selected individual strip* reference or gradient. At 1.1R, local selected-strip potential reaches values as low as ~0.00375, comparable in scale to the prior ~0.00379 max cylinder-surface boundary residual (different sampled positions); this makes extrapolating whole-plane accuracy to individual stereo-strip signals unjustified. No PCB copper gaps/dielectrics, physical edges, or field-level strip convergence has been validated. Next scientific control should be near-wire **individual-strip** vector gradient sensitivity to angular collocation, longitudinal segmentation, finite w correction span and Green-image cutoff, followed by explicit physical induced-strip charge closure. No electronics.

### Candidate U/V near-wire *strip-gradient* convergence target (committed, not yet executed)

The new 3D individual-strip convergence program evaluates both U and V selected-strip three-component gradients at 1.1R, 2R, 5R near three wire centers, four angles and two longitudinal w locations. It compares the baseline 4-angle wire-cylinder mean, 1mm longitudinal segments, ±10mm longitudinal wire-source support and 12 Dirichlet image pairs separately to 8-angle cylinder mean, 0.5mm segments, ±14mm extent, and 20 image pairs. It outputs separate global relative vector-gradient L2, max local relative gradient deviation, max absolute gradient difference, plus potential sensitivities for each variation/plane. A full-plane Garfield agreement of 0.11% in vector gradient does not establish individual-strip-field accuracy. The interstrip gap=0, 3+3mm and other candidate geometry limitations persist. No new avalanches, no RKF drift, no electronics. WSL validation is pending.

### Measured 3D individual U/V near-wire gradient parameter sensitivity (2026-11-20 WSL)

User successfully ran `mwpc_stereo_strip_gradient_convergence` at candidate gap ±3mm, strip width and nominal zero-gap pitch 1.7mm, tan(α)=0.1, 9 anodes of radius 0.015mm and 3mm pitch. Baseline line-source wire segmentation is **1mm** along w, finite source half-length 10mm, 12 Green image pairs, and four-angle **mean** cylinder collocation. 72 sample points per side around three wires at r/R=1.1, 2 and 5. The executable completed, but this does **NOT** establish that the individual strip weighting gradients are converged.

| Independently changed parameter | U strip 3D gradient relative L2 change | V strip relative L2 change | Max local relative gradient change (U/V) | Max abs potential change |
|---|---:|---:|---:|---:|
| Wire-cylinder angle-average 4→8 angles | 6.63e-11 (~0) | 6.63e-11 (~0) | ~1.46e-10 | ~8.26e-12 |
| **Longitudinal segment 1→0.5mm** | **6.03347%** | **6.03344%** | **11.7396% / 11.7408%** | **0.00410861** |
| Source half-span 10→14mm | 0.0002018% | 0.0002018% | 0.0002996% / 0.0003336% | ~7.05e-8 |
| Image cutoff 12→20 pairs | 0.036223% | 0.036223% | 0.0486% / 0.0654% | ~1.3122e-5 |

**Critical interpretation:** the first 1mm longitudinal segmentation is *not* sufficiently stable for individual strip gradients near the anode. The 6.03% is a difference between numerical discretizations, not the true unknown strip-field error; full absolute physical accuracy may be worse. Previous 0.1109% Garfield reference agreement applies only to WHOLE cathode, not individual strips. The 4→8 angle-average virtually unchanged result is not a successful check of wire-cylinder pointwise boundary conditions: angularly averaged zero-potential constraints with only axial monopole line sources cannot reproduce arbitrary angular variation at the 30μm physical cylinder surface. We still need a finer segmentation check and potentially angular surface modes/full conductor BEM. The earlier maximum U/V 12-angle grounded-wire potential residual ~0.00379 is comparable to small near-wire strip weighting potential ~0.00375, so physical strip charge sharing is not certified. Do not proceed to ion-induced U/V strip currents or optimize electronics on these fields.

Follow-up branch code adds `--only-variation segment_half` to the convergence executable so the next test can isolate **0.5→0.25mm** longitudinal segments, using a new output prefix and without redundantly rebuilding the insensitive angular, image and span variants. It prints `NOT_CONVERGED` if the requested gradient changes exceed a 2% internal-screening threshold. That threshold is a diagnostic, **not** a physical acceptance criterion; gradients stable under just one segment change do not prove strip field correctness. This additional test is committed, not yet compiled or run in user WSL.

### Measured second longitudinal mesh refinement and adaptive follow-up (2026-11-20)

The user's WSL 0.5→0.25mm uniform wire line-segment convergence run completed with the original P3 noFW ideal nine-wire reference and candidate 1.7mm U/V stereo strip boundaries. Per cathode 72 near-wire points r/R=1.1,2,5. Baseline 360 dense unknowns (0.5mm) and refined 720 unknowns (0.25mm). This still **did not meet** the program's 2% internal screening criterion; no physical strip current should be calculated from the existing coarse mesh.

| Uniform longitudinal segment refinement | U gradient vector L2 difference | V gradient vector L2 difference | Max local relative gradient change (U/V) | Max abs potential change |
|---|---:|---:|---:|---:|
| 1→0.5mm (previous) | 6.03347% | 6.03344% | 11.7396% / 11.7408% | 0.00410861 |
| **0.5→0.25mm (new)** | **2.571503%** | **2.571493%** | **5.09373% / 5.09397%** | **0.00166804** |

The refinement sensitivity fell to ~42.62% of its earlier value, suggesting convergence but not proving the solution is accurate or that the physical cylinder boundary is enforced. Current uniform solver has a conservative 800-unknown guard; halving again to 0.125mm across all 9×20mm wires would require 1440 line-charge elements, which is not an appropriate next action without a computational strategy.

New **adaptive central w-mesh** capability is now committed in `WireAwareStereoWeighting.hh` and `mwpc_stereo_strip_gradient_convergence`: optional `localHalfSpanWmm`/`localStepMm` refine the central wire region while keeping the outer segmentation coarse, with exact junction alignment, symmetric contiguous source elements and the same maximum 800 dense unknowns. First candidate: global source half-length 10mm, **1mm outer segments**, refined central region |w|≤2mm with **0.125mm** elements (432 unknowns) compared to **0.0625mm** central elements (720 unknowns) using `--only-variation local_segment_half`. The actual separate U and V vector-field gradients are evaluated at the same near-wire locations as before and the independent grounded-wire 12-angle **surface residual** is logged for both models. Candidate controls additionally check halving outer mesh 1→0.5mm and expanding central fine region |w|≤2→3mm, since a numerically stable central mesh is not sufficient if coarse outer elements introduce bias. The Garfield whole-cathode pilot now also accepts the new adaptive settings.

**This adaptive code has NOT yet compiled or run in user's WSL.** The chosen 2mm local refinement region (and 1mm coarse exterior) is an efficiency hypothesis, not a demonstrated accurate mesh. Original uniform-grid and K256 microscopic reference work remains unchanged. If strip gradients remain sensitive or physical 12-angle wire-cylinder boundary residuals stay large, future work should consider pointwise surface boundary enforcement with angular multipoles or full 3D BEM/FEM, rather than repeatedly increasing 1D axial line-source densities. No electronics, no strip induced-current physics claims.

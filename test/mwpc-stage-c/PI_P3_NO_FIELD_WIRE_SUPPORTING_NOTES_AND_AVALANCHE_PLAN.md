# NA60+ MWPC — provisional 3 mm anode-only baseline and avalanche follow-up
**Supporting notes for discussion with A. Milov | 28 October 2026 | NA60+ Simulation — Work in progress**

**Status:** internal working note; NOT a validated detector-design recommendation or a calibrated gas-gain prediction.

## Scope lock for the NEW avalanche campaign (28 Oct 2026)

**Do not transfer previous 2+4 mm Stage A/B numerical results into the new 3+3 mm geometry.** The user's local inventory confirmed legacy files including `gain_fluctuations.csv`, `gain_fluctuations_5000.csv`, `heed_avalanche_seed_summary.csv`, `heed_avalanche_track_summary.csv`, `stageB_real_avalanche_event_summary.csv`, `stageB_real_avalanche_analysis_summary.csv`, and `stageB_avalanche_angle_scan_summary_by_angle.csv`. These validate *workflow and output-column interpretation* but are **not measurements of the new four-geometry campaign**. For example, the old Stage B integrated-charge fraction changes markedly with observation window, and its asymmetric gap necessarily changes cathode weighting fields; these values must not be quoted for ±3 mm.

All NEW avalanche gain/ion birth/induced charge runs shall fix **cathode planes at v=−3 and +3 mm**, and compare the following exactly four configurations under otherwise matched gas/HV/wire diameter/B/seed conditions:

| Configuration | Anode pitch | Grounded midpoint field wire | Purpose |
|:--|--:|:--:|:--|
| **P3_NOFW** | **3.0 mm** | no | **primary reference baseline** |
| P3_FW75 | 3.0 mm | 75 μm | assess FW effect at baseline pitch |
| P4_NOFW | 4.0 mm | no | pitch effect without field wire |
| P4_FW75 | 4.0 mm | 75 μm | pitch–field wire interaction |

**Code audit:** `phaseA_gain_fluctuations.cxx`, `phaseA_heed_avalanche.cxx` and `phaseB_avalanche_angle_scan.cxx` already expose pitch and both cathode-gap options, but their current `ComponentAnalyticField` construction adds **anode wires only** and does not include any grounded midpoint field-wire option. Therefore setting `--gap-minus-mm 3 --gap-plus-mm 3` is sufficient to start a **no-FW 3+3** methodological smoke; it is NOT sufficient to simulate P3_FW75/P4_FW75. Extend the field builder with a consistent grounded wire diameter/position, assert actual geometry and wire labels, and add no-FW/FW75 regression checks **before** doing the four-way quantitative gain/charge comparison. Changes must preserve old 2+4 standalone defaults and legacy result filenames.

**Definitions fixed for subsequent runs:** record **both** the full hemisphere and cathode-facing ±45° sectors for each of (i) primary-electron anode arrival position, (ii) avalanche ionization-birth angles, (iii) avalanche electron endpoint position. Separately compute `A_Q(T)` from upper/lower cathode weighting fields and positive-ion drift. Geometrical sector A is **not** a measured cathode-charge asymmetry. Report zero-response/censored avalanches and finite observation windows.

**Selected working configuration:** anode pitch **3.0 mm**, **no field wires**, anode diameter **30 μm**, cathode planes at **v = ±3 mm**, anode potential **+1800 V**, **Ar/CO₂ (70:30)**, T=293.15 K, P=760 Torr, **B=0 T**. Compare with a midpoint **grounded 75 μm field wire**. The current working selection is conditional on the **full-upper/lower-hemisphere electron endpoint asymmetry**, which is *not* an induced-cathode-charge observable.

## 1. Research objective and definitions

For a physical muon crossing of one chamber gas volume, Heed samples primary ionisation clusters and their electron multiplicities along the Geant4 muon trajectory. For each primary electron the frozen Garfield electron-transport LUT maps the initial (u,v) to a wire endpoint, drift time and anode-surface impact-position azimuth \`φ_hit\`, or an electron loss. This calculation is **before avalanche multiplication**.

Define **hemisphere-count** asymmetry per physical gas crossing:
\[
A_{\mathrm{hemi},i}=\frac{N_{i,\phi>0}-N_{i,\phi<0}}{N_{i,\phi>0}+N_{i,\phi<0}}.
\]
Exact horizontal angular boundary hits are excluded. An event with zero denominator has **undefined** asymmetry and is not assigned A=0. Each valid crossing gets equal weight in P(A), regardless of electron multiplicity.

**Do not mix statistical units:** the earlier p=3 mm LUT500 MS0 screenshot had **99,968 parents with crossing**, and measures one A **after pooling every gas crossing of that J/ψ parent within a station**; its RMS 0.3106 (no FW) vs 0.3500 (FW75) is a **parent-event** result. The 100k-parent LUT20 \`full34\` study below instead reports *one A per gas crossing* and thus different RMS values. These are independently useful, but not interchangeable.

A distinct **cathode-facing sector** asymmetry uses only ±45° sectors about +v and −v:
\[
A_{\mathrm{TB},i}=(N_{\mathrm{top},i}-N_{\mathrm{bottom},i})/
                  (N_{\mathrm{top},i}+N_{\mathrm{bottom},i}).
\]
The two observables can have different trends, and neither is yet the measured cathode-charge asymmetry.

## 2. Evidence supporting the present no-FW selection

### 2.1 Primary electron transport: independent 126k Garfield full-cell pilot (21 Oct)

Same pitch and full-cell source, **18,000 Garfield microscopic electron starts/geometry** at p=3 mm. Original output documented in \`test/mwpc-stage-c/README.md\` under “2026-10-21: full 126k Garfield pitch pilot completed” and saved under \`test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/pilot/\`.

| p=3 mm microscopic full-cell pilot | No field wire | 75 μm grounded FW | FW − noFW |
|:--|--:|--:|--:|
| Anode-collected electrons | 14,875 / 18,000 | 14,258 / 18,000 | −617 |
| Anode collection fraction | 82.64% | 79.21% | **−3.43 pp** |
| Mean drift time **conditional on collection** | 27.152 ns | 26.614 ns | **−0.538 ns** |
| Centered drift-time RMS **conditional on collection** | 14.949 ns | 14.727 ns | −0.222 ns |

Interpretation: adding field wires produces a small shift toward faster *surviving-electron* drift (~0.54 ns or 2.0% of the noFW mean), **not** a demonstrated improvement in signal timing resolution. It also reduces the anode collection fraction by ~3.43 pp in this run; this is a measured model trade-off, not an automatically “negligible” loss. Field-wire geometries can change which slow electrons survive; the conditional time difference is **not** necessarily a paired physical drift-time acceleration.

The paired exact-birth Garfield comparison exists in \`pi_pitch_fieldwire_paired_analysis.py\` and its read-only runner; preserve the special exclusion for source electrons artificially relocated out of metal volume.

### 2.2 Complete 100k-parent J/ψ Geant4 → Heed → six LUT20 response (uploaded \`full34.zip\`)

989,695 physical **muon gas crossings** across MS0–MS5 from 100,000 generated J/ψ parents; for each crossing, use matched Heed primary electrons and the frozen p=3 mm noFW/FW75 LUT20. Both variants have 55,784,197 **primary electron instances** in the folded gas crossings (not independent microscopic Garfield histories).

| p=3 mm, 100k-parent full34, LUT20 | No field wire | FW75 | FW − noFW |
|:--|--:|--:|--:|
| Primary electron instances | 55,784,197 | 55,784,197 | same |
| Collected at a physical anode | 46,066,785 | 44,168,063 | −1,898,722 |
| Collection / primary electron | **82.580%** | **79.177%** | **−3.404 pp** |
| RMS of crossing \`A_hemi\` | **0.36490** | **0.40066** | **+0.03576** |
| Valid crossing \`A_hemi\` | 988,899 | 988,830 | −69 |
| Electron fraction landing in ±45° cathode-facing sectors | 41.563% | **64.408%** | **+22.845 pp** |
| RMS of crossing \`A_TB\` among defined crossings | 0.61871 | **0.45675** | **−0.16196** |

This demonstrates that **no FW has a smaller RMS for the specified full-hemisphere crossing-count asymmetry**, and **higher primary anode collection**. It also demonstrates the **opposite trend** for the restricted cathode-facing-sector A: with FW, fewer extreme differences appear between those sectors because many more electrons reach them. Therefore the data **do not justify the general claim** “no FW gives more balanced upper/lower cathode signals”. A full cathode signal requires avalanche ion-birth positions, ion transport, Shockley–Ramo weighting fields, time integration and strip response.

The 100k-parent **LUT20** has finite-20-histories-per-cell re-use and nontrivial parent-event correlations. No independent finite-LUT error budget or complete matched-parent confidence interval is claimed for these numbers.

### 2.3 Separate parent-level p=3 mm LUT500 comparison

In the earlier LUT500 dashboard, MS0 had 99,968 J/ψ **parents with a station crossing**:
- **No FW:** mean A=+0.0008, RMS(A)=0.3106.
- **FW75:** mean A=−0.0014, RMS(A)=0.3500.
- These are *parent-aggregated* A values and **must not** be overlaid with the full34 crossing A values as though they were the same measurement.

The original p=3 mm independent microscopic LUT500 maps remain frozen and are the natural check of finite-LUT uncertainty in the new crossing-level study. The ongoing/finished p=4 mm LUT500 is separate, and **not needed** to document the present p=3 mm decision.

## 3. Decision for the PI and open conditions

**Working decision:** continue with **p=3 mm, no field wire** for the **first avalanche and gain simulations**. This keeps the current geometry simple, is consistent with the lower observed *full-hemisphere* event-by-event A RMS and higher primary-electron anode collection, and avoids sacrificing ~3.4 pp of collection for a shift of only ~0.54 ns in the uniform Garfield pilot's *conditional mean drift time*.

**Important qualifications:**
1. \`0.54 ns\` is a small *magnitude*, not yet a formal statistically non-significant difference; obtain uncertainties from paired Garfield and native 100k timing reports.
2. This is a **research baseline, not final hardware selection**. Field wires increase the electron population in cathode-facing sectors and reduce restricted-sector A RMS; the downstream cathode signal could prefer a different design.
3. No induced cathode charge, gas gain, front-end shaping, thresholds, noise, occupancy or double-hit separation is included.
4. Baseline gain has not been calibrated against prototype data; do not present absolute gain as confirmed.

**Before freezing text for Overleaf:** obtain the existing **100k-parent native six-geometry drift-time report** \`test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/station_summary.csv\` (or \`summary.txt\`), and the p=3 mm LUT500 crossing-level summary if complete. Verify consistency of the actual J/ψ timing selection; no new Geant4 or Garfield production is required to access these files.

## 4. Avalanche work already completed in the repository

The existing reference research code is in \`test/mwpc-avalanche/\`, documented in its \`README.md\`. These are previous validated *methodological prototypes*, largely using **p=4 mm and asymmetric 2+4 mm gaps** by default; they **are not yet production measurements of the present p=3 mm, ±3 mm, noFW baseline**.

- \`phaseA_gain_fluctuations.cxx\` / target \`mwpc_phase_a_gain\`: independent **single-electron AvalancheMicroscopic** ensembles. Outputs \`avalanche_electrons\`, \`avalanche_ions\`, **\`collected_on_wire\`**, attached and other endpoints, eventwise collection times; prints unconditional effective gain mean/RMS/quantiles and zero-collected fraction. Has explicit avalanche-size cap and detects capped events.
- \`phaseA_heed_avalanche.cxx\`: **Heed muon clusters -> all primary electrons -> microscopic avalanches -> collected charge per anode**. Already captures seedwise and trackwise electron/ion counts and wire multiplicity, but its default reference geometry is p=4 mm, gap 2+4 mm, no FW.
- \`phaseB_real_avalanche_signal.cxx\`: actual single-electron microscopic avalanche with **ionisation-birth callback**, real positive-ion drift and induced electron/ion strip currents. A one-seed demonstrator, **not** a full muon hit.
- \`phaseB_avalanche_angle_scan.cxx\` / \`analyze_stageB_avalanche_angle_scan.py\`: scans avalanche-seed azimuth, records ionisation birth angles/centroids, ion multiplicity, and integrated **upper/lower complete-cathode induced charges** versus seed azimuth. The source defaults to **p=4 mm, gap 2+4 mm**. It skips zero/failed avalanches in the output; **do not** derive unbiased gain from its accepted-only subset.
- \`MWPCSignalUtils.h\` and \`MWPCStereoReadout.h\` implement reusable signal/strip infrastructure.
- \`test/mwpc-stage-c/PI_FAST_LUT_TRANSPORT_ONLY_README.md\` documents prior fast **AvalancheMC** transport-only (attachment OFF, 5 μm step) comparisons; **these are not avalanche-gain simulations**. Do not infer avalanche multiplicative gain or signal charges from this transport benchmark.

The existence and prior validation of these workflows is supported by repository documentation, but the original numerical CSV/run manifests for Phase-A gain or Stage-B ion signal **have not been rechecked against the local WSL outputs in this review**. No numeric gain is claimed here.

## 5. What does a Geant4 muon hit tell us about avalanche gain?

Geant4's muon gas hit gives the **entry/exit trajectory, gas path and deposited energy**, but no avalanche multiplication. Heed turns a crossing into **ionisation clusters and N_primary primary conduction electrons**, with start space/time coordinates. Garfield **AvalancheMicroscopic** must then track each primary seed in the chamber's high field, where it can create additional e⁻–ion⁺ pairs. All offspring form one microscopic avalanche for that seed.

For one seed k:
- \`N_pairs,k\`: count of avalanche ionisation electron–ion **pairs generated**, best monitored using the ionisation callback; **not** initial Heed primary pairs.
- \`N_e,aval,k\`, \`N_i,aval,k\`: bookkeeping numbers returned by \`GetAvalancheSize(ne, ni)\`, whose detailed counting conventions/attachment effects must be inspected.
- \`N_e,anode,k\`: **number of avalanche electron endpoints actually collected at an anode**; a different number, and the observable used for effective charge gain.
- \`G_eff,k=N_e,anode,k/1\` for one initial primary electron, **including zero-response seeds**.

For a muon gas crossing with \`N_primary\` Heed electrons:
\[
G_{\mathrm{eff,hit}} =
\frac{\sum_{k=1}^{N_{\rm primary}}N_{e,\mathrm{anode},k}}
     {N_{\rm primary}}
\quad(N_{\rm primary}>0).
\]
Average gain across crossings **must specify the weighting**: pooled collected electrons / pooled primary seeds vs equal-hit average of \`G_eff,hit\`. Also report conditional gain for primaries surviving to the multiplication region separately. Do **not** use \`N_e+N_i\` (a double counting of ionisation products) as a gain denominator. Do not treat a Garfield upper avalanche size cap as a physical gain value.

## 6. Proposed avalanche asymmetries — three successive physics levels

**(A) Angular location of ionisation births**, for each primary-electron avalanche (and ultimately each full muon crossing). Classify every microscopic ionisation-birth point in the near-anode u–v plane:
\[
A_{\mathrm{ion,hit}} =
\frac{N_{\rm ion,birth}(v>0)-N_{\rm ion,birth}(v<0)}
     {N_{\rm ion,birth}(v>0)+N_{\rm ion,birth}(v<0)}.
\]
Also compute restricted ±45° sectors and the circular mean of the ion-birth angle. This accounts for the *gain-weighted* avalanche-ion population, not just one \`φ_hit\` per initial electron. Preserve zero-ion/undefined events.

**(B) Amplified electron endpoint asymmetry**, counting ALL final avalanche-electron anode endpoints, not one count per primary seed:
\[
A_{e,\mathrm{aval,hit}} =
\frac{N_{e,\mathrm{anode,up}}-N_{e,\mathrm{anode,down}}}
     {N_{e,\mathrm{anode,up}}+N_{e,\mathrm{anode,down}}}.
\]
This differs from primary-electron A because every seed has a stochastic gain, and attached/failed avalanches contribute zero.

**(C) Cathode integrated-charge asymmetry at observation time T**:
\[
A_Q(T)=\frac{|Q_{+}(T)|-|Q_{-}(T)|}
             {|Q_{+}(T)|+|Q_{-}(T)|},
\]
using positive-ion and electron charge transport and each cathode's weighting field. This is the first stage that addresses the **actual upper/lower readout**, and differs fundamentally from A_ion and A_e. For real U/V strips, compute Q per strip, then define plane-level totals consistently.

## 7. Minimal safe execution plan

1. **Inventory (no physics run):** \`find test_runs/mwpc_stage_c -iname '*gain*' -o -iname '*avalanche*'\`, and inspect any existing Phase A / Stage B output directories, manifests, timestamps and number of events. Verify gas composition, wire diameter/pitch, HV, gap, B, ion mobility, Penning settings, and any size-limit censoring. Never overwrite old output.
2. **Rebuild pre-existing Phase A/B targets** (standalone \`test/mwpc-avalanche/CMakeLists.txt\`). Check \`mwpc_phase_a_gain\`, \`mwpc_phase_b_avalanche_angle_scan\`, \`mwpc_phase_b_real_avalanche\`. Reuse Stage-A gas and Stage-B ionization callback; do not copy those models into independent implementations.
3. **Small p=3 mm, noFW, ±3 mm gain smoke**: run an unconditional, location-resolved single-electron ensemble (top, bottom, side, several seed distances) with Garfield \`AvalancheMicroscopic\` and explicit cap flags. Existing \`phaseA_gain_fluctuations.cxx\` already accepts \`--pitch-mm 3 --gap-minus-mm 3 --gap-plus-mm 3\`, etc. Track 100%-population mean, zeros, RMS, distribution tails and attachment. A one-position smoke is **not** representative gain over a Heed muon hit.
4. **Same-geometry angle/ion-birth smoke:** use \`mwpc_phase_b_avalanche_angle_scan\` with \`--wire-pitch-mm 3 --gap-minus-mm 3 --gap-plus-mm 3\`. It already records seed azimuth, \`birth_mean_theta_deg\`, ion counts and upper/lower cathode Q. Its acceptance filtering makes it an **angle-response diagnostic, not an unbiased gain ensemble**.
5. **New full-hit integration on a bounded subset** of *original* Geant4 muon gas crossings: deterministic shared Heed primary births; independent microscopic avalanche per seed; accumulate raw \`N_ionbirth, N_e,anode\` and 3 asymmetry definitions for ONE physical crossing. Run a small provenance-checked smoke first, then size-tested representative station/DET samples. Record computational time and size-limit fraction before scale-up. Do **not** attempt all 989,695 crossings with full microscopic avalanche before runtime review.
6. **Stage-B signal:** for selected representative hits, drift the positive ions (or use a verified response kernel) and calculate \`Q_±(T)\` for multiple time windows; later the real cathode strip Q and ASIC shaping. Study correlations between \`A_primary\`, \`A_ion\`, \`A_e,aval\` and \`A_Q(T)\`. Bootstrap by *independent J/ψ parent*, and separately assess finite microscopic sampling.

## 8. Planned figures for eventual Overleaf supporting note

1. **Figure 1 — Geometry and definitions.** Schematic of p=3 mm, gap ±3 mm, noFW vs grounded FW75; anode hemispheres versus cathode-facing ±45° sectors (different definitions).
2. **Figure 2 — Primary collection and time.** Anode-collected fractions and drift-time PDFs (paired 18k-cell pilot; in final version prefer 100k native time report). Show observation population and uncertainty.
3. **Figure 3 — Native J/ψ event-by-event primary asymmetry.** p=3 mm \`P(A_hemi)\` on one **muon gas crossing** (LUT20 full34, blue solid noFW/red dashed FW75; RMS/valid counts). Separate parent-level LUT500 only if labelled distinctly.
4. **Figure 4 — Counterpoint relevant to future cathodes.** p=3 mm \`P(A_TB)\` and cathode-facing electron fraction; explicitly show opposite RMS trend.
5. **Figure 5 — Future avalanche physics.** Unconditional \`P(G_eff)\`, zero-gain fraction, ion-birth angular histograms, then \`A_ion\`, \`A_e,aval\`, and \`A_Q(T)\` vs each other from same primary-crossing ensemble.

**Presentation:** ATLAS-inspired professional plot style (white canvas, black axes, complete units/normalization; blue solid noFW, red dashed FW75; sample and error definitions); label **“NA60+ Simulation — Work in progress”**, not “ATLAS”, and avoid calling numerical gain prototype-calibrated.

## 9. Repository evidence / source files

- \`test/mwpc-stage-c/README.md\`: 21 Oct 126k six-geometry Garfield pilot and conditional drift times.
- \`test/mwpc-stage-c/pi_pitch_fieldwire_paired_analysis.py\`: matched-birth collection/time audit.
- \`test/mwpc-stage-c/pi_lut500_muon_crossing_asymmetry.py\`: exact crossing-level A, raw IDs.
- \`test/mwpc-stage-c/PI_FAST_LUT_TRANSPORT_ONLY_README.md\`: AvalancheMC is a **transport-only** candidate, not multiplication.
- \`test/mwpc-avalanche/README.md\`: existing Phase-A and Phase-B scope; older p=4 mm and 2+4 mm reference geometry.
- \`test/mwpc-avalanche/phaseA_gain_fluctuations.cxx\`, \`phaseA_heed_avalanche.cxx\`, \`phaseB_real_avalanche_signal.cxx\`, \`phaseB_avalanche_angle_scan.cxx\`: existing tested methods.
- Garfield++ User Guide and ROOT/Garfield maintainers' clarification of \`GetAvalancheSize\` vs electron endpoint-count and the role of attachment: https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf and https://root-forum.cern.ch/t/difference-between-getavalanchesize-and-getnumberofelectronendpoints/50834

All tables here use the existing full-cell Garfield pilot numbers from repository documentation or the exact full34 ZIP provided earlier. **No new Garfield avalanche physics run was carried out in producing these notes.**

## 10. First symmetric 3+3 mm avalanche smoke5 (28 October 2026)

Google Drive frozen run:
https://drive.google.com/drive/folders/1-UCTer9rh2tuiZBkL_OWQxldFeUBPNHe

All four source geometry logs confirm **cathode gaps 3+3 mm**, anode 30 μm, +1800 V, Ar/CO₂ (70:30), B=0 T, and grounded midpoint FW diameter 0 or 75 μm. The smoke sampled just **five single-electron avalanches per geometry**, starting at (u,v)=(0,+1 mm), with base seed 280103. All four logs report **zero size-limit-like events** (limit 200,000).

| 3+3 mm single-seed smoke | 3 mm noFW | 3 mm FW75 | 4 mm noFW | 4 mm FW75 |
| --- | ---: | ---: | ---: | ---: |
| Effective collected gain mean | 600.2 | 8,958.4 | 1,429.6 | 8,098.6 |
| Population RMS | 495.7 | 10,201.6 | 505.3 | 8,107.7 |
| Zero collected avalanches | 1/5 | 1/5 | 0/5 | 1/5 |
| Size-limit-like | 0/5 | 0/5 | 0/5 | 0/5 |

These results establish that our **new geometries run** and suggest a large *local* FW-induced change in gain near the chosen anode. However, five highly fluctuating avalanches do **not** estimate the representative gas gain or precision FW effect. In particular, the p3 FW75 sample contains avalanches with 23,724 and 18,877 collected electrons, while p3 noFW peaked at 1,371. Do not interpret mean ratios as precise physics and do not change the baseline selection based on this smoke test.

**Metadata warning**: the historical immutable smoke5 `run_config.json` incorrectly retains status=`incomplete`, even though `summary.json`/`manifest.json` say complete, all four result CSVs contain five entries, and rclone reported 12 matching cloud files. The new writer has been fixed so later runs finalize `run_config.json` and SHA-pin it; *do not mutate the historical run*.

Next: run the 100-seed fixed-position pilot in the same four geometries, then extend to multiple physically relevant primary-electron locations. Only after this should the work proceed to avalanche ion-birth `A_hemi`/`A_TB` and the actual cathode induced-charge asymmetry `A_Q(T)`.

## 11. Completed 100-seed symmetric 3+3 gain pilot (28 October 2026)

Original immutable cloud run (all four CSV, logs, run_config, SHA-validated manifest):
https://drive.google.com/drive/folders/1TQvDjcowzlSdN1GBw50-sLWIziGsvVny

Same four geometries and common parameters as smoke5: Ar/CO₂ 70:30,
+1800 V, 30 μm anodes, symmetric ±3 mm cathodes, B=0 T,
one seed per avalanche, source at (u,v)=(0,+1) mm, 100
independent-event RNG seeds per geometry. `run_config.status=complete`;
all four logs report **zero avalanche-size-limit-like events** (cap=200,000).
Means include failed/zero-collection seeds.

| Single-seed observable (100 per geometry) | p3 noFW | p3 FW75 | p4 noFW | p4 FW75 |
| --- | ---: | ---: | ---: | ---: |
| Mean collected gain | **554.45** | 6,418.18 | 1,648.90 | 5,039.70 |
| RMS (population) | 585.72 | 6,513.07 | 1,740.67 | 5,884.51 |
| Median collected gain | 411.5 | 4,499.5 | 1,238.0 | 3,666.0 |
| Zero collected seeds | 18/100 | 15/100 | 10/100 | 18/100 |
| Maximum of sample | 2,962 | 28,716 | 12,912 | 35,796 |
| Size-limit-like events | 0/100 | 0/100 | 0/100 | 0/100 |

**Observation, not extrapolation:** At this exact seed point, the FW75/noFW
mean-gain ratio is 11.58 for p3 and 3.06 for p4. Broad distributions with
RMS comparable to means require larger ensembles and spatial sampling
before characterizing the chamber. Matching seed integers across
geometries are *not* a guarantee of physically identical microscopic
avalanche paths (RNG consumption and fields differ).

**Limits:** These are per-single-seed electron collections, not Heed primary
cluster-averaged or J/ψ hit-averaged gas gains; there are no ion-birth
angular asymmetries or induced cathode charges in these outputs. The
provisional **p3 noFW** reference remains unchanged pending downstream
avalanche asymmetry and signal evaluation. The earlier 2+4 mm
results must never be substituted for this ±3 mm campaign.

## 12. Literature interpretation of FW gain and Stage B1 protocol

**Question raised 28 October 2026:** At a fixed +1800 V, the p3 noFW to
p3 FW75 gain changed from 554.45 to 6418.18 (11.58×), and at p4 from
1648.9 to 5039.7 (3.06×), using 100 microscopic single-electron
avalanches per geometry with one seed at (u,v)=(0,+1 mm).

**External support, not a numerical calibration:**
- H. J. Hilke, “Time projection chambers,” *Rep. Prog. Phys.* **73**
  (2010) 116201, section 3.2.1.1, reports that chamber operation
  **without field wires requires higher anode voltage for the same gain**.
  It also notes a substantial dependence of induced pad charge on
  field wires (ALEPH comparison: ~20% with field wires vs ~35% without,
  under different geometry/gas); see
  https://lhcb.web.cern.ch/lhcb_page/lhcb_documents/presentations/TalksNotDirectlyRelatedLHCb/Postcript/Hilke.pdf
- “High-voltage switching multiwire proportional chamber,”
  *Prog. Theor. Exp. Phys.* (2017) 023C01, uses alternating anode and
  potential wires and Garfield++, with a ±3 mm cathode gap. Gain can
  be strongly controlled by potential-wire bias; this paper studies
  **potential-wire voltage variation**, NOT a quantitative p3
  noFW-versus-FW75 factor:
  https://academic.oup.com/ptep/article/2017/2/023C01/2996387
- Gerald Eigen, “Gaseous Tracking Detectors” (2025), section 4.1,
  explains the near-anode field, how anode capacitance/charge
  controls gain, and how changes of geometry can amplify gain changes:
  https://link.springer.com/chapter/10.1007/978-3-031-67336-8_4

**Interpretation:** a large increase in gain at fixed HV is physically
plausible; **11.6× is not a universal or literature-certified FW
multiplier**. At p3 the *conditional positive-response* mean changed
from 676.16 (noFW) to 7550.8 (FW75), ~11.17×, so the 11.6×
difference is not simply the fraction of zero-response seeds
(18% vs 15%). Absolute gain depends on field geometry, gas
transport, Penning model and gain calculation conventions.

**Required electric-field audit before treating this as an intrinsic
detector result:** for all four geometries with identical +1800 V,
sample Garfield |E| and radial component at 20, 30, 50, 100,
and 200 µm from the central physical anode along top, bottom and
side directions; verify cathodes at v=±3 mm, FW at grounded midpoint
and near-wire field convergence on increasing the number of analytic
wires. Add a small HV scan to separate *fixed-voltage effects*
from geometry dependence *at equal effective gain*. Neither
gain nor cathode induced charge should be judged only at +1800 V
if the chosen configurations have dissimilar multiplication.

### Stage B1: primary single-seed ion-birth angular asymmetry

First implement **only the ion-birth and avalanche-electron endpoint
observables** by deriving a separate, non-overwriting program from
`test/mwpc-avalanche/phaseB_avalanche_angle_scan.cxx`.
Keep the historical 2+4 mm code/outputs immutable; switch the
four new configurations to ±3 mm and grounded midpoint FW75
where appropriate.

- For each independent initial electron seed, call
  `AvalancheMicroscopic::SetUserHandleIonisation` with the
  existing `RecordIonisation(x,y,z,t,type,level,Medium*)` callback,
  then `AvalancheElectron`. The callback records each
  simulated ionising collision (including ionisations initiated by
  Penning transfer if enabled). Count these as new avalanche
  positive-ion **births**, not as drifted ions or a cathode current.
  Compare callback counts with Garfield bookkeeping and flag
  discrepancies rather than silently equating conventions.
- For each birth, find the associated physical anode center
  `u_a` and calculate
  `phi=atan2(v-v_a,u-u_a)`; histogram full [-180,180) azimuth
  AND radius `r` from the anode. Track births far from the
  multiplication region separately rather than silently assuming
  all are close to the anode.
- Record `N_ion_up` at 0<phi<180,
  `N_ion_down` at −180<phi<0, and exact boundaries;
  separately `N_ion_top` for 45<=phi<135 and
  `N_ion_bottom` for −135<=phi<−45, plus side counts.
  Compute **two distinct ion-birth asymmetries** per seed:
  `A_ion_hemi=(N_up-N_down)/(N_up+N_down)` and
  `A_ion_TB=(N_top-N_bottom)/(N_top+N_bottom)`.
  Zero denominators are **undefined/null**, not assigned 0.
- Retain ALL requested seeds, including zero-gain, attachment,
  no-anode endpoint, time-window/cap failure, with explicit
  statuses and the original seeded geometry/initial position.
  Distinguish equal-seed averaging from pooled-ion counting
  (different weighting by stochastic gain).
- Begin with the already validated source (u,v)=(0,+1mm),
  then mirror (0,−1mm) to test top/bottom symmetry
  and add ±u off-axis starts and a near-anode radial/angle
  scan. The full six-station Geant4+Heed muon-crossing
  weighting is a **later step**; a single fixed top-side
  birth point does not give an unbiased detector A.
- Use a no-ion-drift B1 mode first for efficient counting.
  The **B3** induced cathode charge requires ion transport
  and the cathode weighting potentials; do NOT infer
  Q_upper/Q_lower directly from ion-birth sector counts.
  Later calculate `A_Q(T)` for several observation windows
  and, for controlled comparisons, both fixed 1800 V and
  gain-matched voltages.

**Status:** only the protocol and original legacy callback were
reviewed here. This is NOT a claim that the new four-geometry
B1 implementation has compiled or produced data.

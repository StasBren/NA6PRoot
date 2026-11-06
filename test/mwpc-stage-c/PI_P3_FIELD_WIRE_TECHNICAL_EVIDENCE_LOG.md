# NA60+ MWPC — technical evidence and reproducibility log
**Internal simulation record | Version 1.0, 6 November 2026 | Living document**

**Status:** a record of existing simulation outputs and their provenance. Values copied from existing notes or dashboard screenshots are explicitly identified. This document is *not* a final chamber-design specification, detector-calibration report, or independent validation of the numerical results.

**Companion documents (three active documents in total):**
- [Supporting note: provisional no-FW decision](PI_P3_NO_FIELD_WIRE_SUPPORTING_NOTES_AND_AVALANCHE_PLAN.md) — concise physics argument and principal whole-MS0 result; send this one for a decision discussion.
- **This document** — methods, statistics, numerical detail, historical studies, traceability and validation issues.
- [Follow-up research and signal-response plan](PI_P3_FIELD_WIRE_SIGNAL_RESPONSE_FOLLOWUP_PLAN.md) — future measurements, decision gates and unresolved questions.
- [Immutable editorial archive of former combined document](PI_P3_NOFW_COMBINED_WORKING_NOTE_ARCHIVE_2026-11-06.md) — preserves the previous text verbatim, including superseded plans. The archive is **not** an extra active note.

## 1. Fixed scope and definitions

The *current* no-field-wire comparison is at **3.0 mm anode pitch** with a symmetric **3+3 mm** cathode gap, **30 μm anode diameter**, **+1800 V anodes**, **Ar/CO₂ 70:30**, **293.15 K**, **760 Torr**, **B=0 T**, compared with a **grounded 75 μm midpoint field wire**. These are the simulation conditions, not an experimental hardware recommendation. The comparison of **3 versus 4 mm anode pitch is intentionally deferred until signal response is studied**; the two pitches are not ranked here.

Three distinct measurements must never be conflated:

| Level | Observation unit | Observable | What it excludes |
|:--|:--|:--|:--|
| MS0 station parent | One distinct generated J/ψ parent with ≥1 MS0 muon gas crossing | A after **summing classified up/down electron endpoints over every MS0 crossing and module of that parent** | Avalanche amplification, cathode Q, actual spatial resolution |
| Physical detector parent | One distinct parent with ≥1 crossing in a selected physical DET | A after counting endpoints across this parent's crossings in that DET | Whole-station parent weighting |
| Gas crossing | One physical muon traversal of a chamber gas volume | A from endpoints associated with **that crossing alone** | Parent-level pooled statistics |

For any of these, define `A_hemi=(N_up-N_down)/(N_up+N_down)` by **anode-surface endpoint position azimuth** `φ_hit`, not arrival velocity, electron avalanche gain or induced charge. Up: `0<φ_hit<180°`; down: `-180°<φ_hit<0`; exact boundaries excluded. If the denominator is zero, A is **undefined**, never A=0. Statistics are *equal-weight observations* at the designated level (population RMS). A separate cathode-facing ±45° sector asymmetry `A_TB` is a different observable. Signal-charge asymmetry `A_Q(T)` requires weighting fields, avalanche and ion transport.

### Input provenance and model-specific qualifications

- Native Geant4 crossings provide the physical muon tracks/chords; the displayed p3/LUT500 and full34 studies fold a **Heed-calibrated marked-Poisson primary-ionisation sampler** along those chords, not a direct per-parent `TrackHeed` run.
- The frozen p3 **LUT500** maps contain **500 microscopic Garfield histories per spatial cell** (900 cells/configuration in the dashboard builder). Folding tens of millions of primaries through them **does not** create that many independent Garfield trials.
- The separate `full34` analysis used **LUT20**, sampling only 20 microscopic histories per spatial cell. Its crossing-level RMS values must not be treated as direct uncertainty checks on the LUT500 station-parent distribution.
- Reused lookup maps, shared parent input tracks and correlated event populations require dedicated uncertainty treatment. None of the displayed RMS differences has a validated finite-LUT-plus-parent confidence interval.
- The source files in `test_runs/` are local/cloud outputs. They are not automatically available merely because the plotting scripts are present in GitHub.

## 2. Primary result audit: MS0 entire station, p3/LUT500 (6 Nov)

Original output path: `test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/summary.json` (and `event_counts.csv.gz` with SHA-pinned `manifest.json`). Dashboard generator: [`pi_lut500_event_asymmetry_root_dashboard.py`](pi_lut500_event_asymmetry_root_dashboard.py). See [LUT500 dashboard documentation](PI_LUT500_EVENT_ASYMMETRY_ROOT_DASHBOARD.md).

| Whole-MS0 parent-level observation | No field wire | Grounded FW75 |
|:--|--:|--:|
| Input generated parents | 100,000 | 100,000 |
| Parents with at least one MS0 crossing | 99,968 | 99,968 |
| Defined A | 99,968 | 99,968 |
| Undefined A | 0 | 0 |
| Mean A | +0.0008 | −0.0014 |
| Population RMS(A) | **0.3106** | **0.3500** |
| N classified electron endpoints | 8,524,668 | 8,242,284 |
| N up | 4,274,033 | 4,125,035 |
| N down | 4,250,635 | 4,117,249 |

**Observed width difference:** `0.3500−0.3106=0.0394`, corresponding to **11.26% smaller RMS in noFW relative to FW75**. The paired mean `A_FW75−A_noFW ≃−0.0022` is a mean shift, not the width difference. The populations have different numbers of collected endpoints (~3.43% more in noFW), which can change counting noise. A simple independent-binomial estimate would suggest a much smaller ~1.7% change in statistical RMS from the aggregate electron-count difference alone, but this is *not* a corrected detector result; correlations and parent multiplicity must be retained.

**Provenance status:** these rounded table values are transcribed from the user's **6 November 2026 complete MS0 dashboard screenshot**. The saved histogram `summary.json` was **not independently accessed** during this note update. The traced illustration in the supporting note is therefore **not** a numerical source. Render the real histogram bins with:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
python test/mwpc-stage-c/plot_pi_ms0_station_parent_asymmetry.py \
  --summary test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/summary.json \
  --output test/mwpc-stage-c/figures/ms0_station_parent_A_p3_LUT500.png
```

The script selects **`stations["0"]["station_parent"]`**, verifies complete source exposure and bin-count conservation, and runs **no Geant4, Heed or Garfield**. The resulting exact-bin PNG should replace the temporary screenshot-traced SVG when available; record its checksum and dataset manifest.

## 3. Independent electron-transport and crossing-level evidence

### 3.1 Primary electron transport: independent 126k Garfield full-cell pilot (21 Oct)

Same pitch and full-cell source, **18,000 Garfield microscopic electron starts/geometry** at p=3 mm. Original output documented in `test/mwpc-stage-c/README.md` under “2026-10-21: full 126k Garfield pitch pilot completed” and saved under `test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/pilot/`.

| p=3 mm microscopic full-cell pilot | No field wire | 75 μm grounded FW | FW − noFW |
|:--|--:|--:|--:|
| Anode-collected electrons | 14,875 / 18,000 | 14,258 / 18,000 | −617 |
| Anode collection fraction | 82.64% | 79.21% | **−3.43 pp** |
| Mean drift time **conditional on collection** | 27.152 ns | 26.614 ns | **−0.538 ns** |
| Centered drift-time RMS **conditional on collection** | 14.949 ns | 14.727 ns | −0.222 ns |

Interpretation: adding field wires produces a small shift toward faster *surviving-electron* drift (~0.54 ns or 2.0% of the noFW mean), **not** a demonstrated improvement in signal timing resolution. It also reduces the anode collection fraction by ~3.43 pp in this run; this is a measured model trade-off, not an automatically “negligible” loss. Field-wire geometries can change which slow electrons survive; the conditional time difference is **not** necessarily a paired physical drift-time acceleration.

The paired exact-birth Garfield comparison exists in `pi_pitch_fieldwire_paired_analysis.py` and its read-only runner; preserve the special exclusion for source electrons artificially relocated out of metal volume.


### 3.2 Complete 100k-parent J/ψ Geant4 → Heed → six LUT20 response (uploaded `full34.zip`)

989,695 physical **muon gas crossings** across MS0–MS5 from 100,000 generated J/ψ parents; for each crossing, use matched Heed primary electrons and the frozen p=3 mm noFW/FW75 LUT20. Both variants have 55,784,197 **primary electron instances** in the folded gas crossings (not independent microscopic Garfield histories).

| p=3 mm, 100k-parent full34, LUT20 | No field wire | FW75 | FW − noFW |
|:--|--:|--:|--:|
| Primary electron instances | 55,784,197 | 55,784,197 | same |
| Collected at a physical anode | 46,066,785 | 44,168,063 | −1,898,722 |
| Collection / primary electron | **82.580%** | **79.177%** | **−3.404 pp** |
| RMS of crossing `A_hemi` | **0.36490** | **0.40066** | **+0.03576** |
| Valid crossing `A_hemi` | 988,899 | 988,830 | −69 |
| Electron fraction landing in ±45° cathode-facing sectors | 41.563% | **64.408%** | **+22.845 pp** |
| RMS of crossing `A_TB` among defined crossings | 0.61871 | **0.45675** | **−0.16196** |

This demonstrates that **no FW has a smaller RMS for the specified full-hemisphere crossing-count asymmetry**, and **higher primary anode collection**. It also demonstrates the **opposite trend** for the restricted cathode-facing-sector A: with FW, fewer extreme differences appear between those sectors because many more electrons reach them. Therefore the data **do not justify the general claim** “no FW gives more balanced upper/lower cathode signals”. A full cathode signal requires avalanche ion-birth positions, ion transport, Shockley–Ramo weighting fields, time integration and strip response.

The 100k-parent **LUT20** has finite-20-histories-per-cell re-use and nontrivial parent-event correlations. No independent finite-LUT error budget or complete matched-parent confidence interval is claimed for these numbers.


## 4. Module-level cross-check: MS0/DET21 (not representative of the entire station)

### 4.1 Local MS0 / physical DET21 parent-level LUT500 cross-check (6 Nov 2026)

**Purpose:** strengthen the **provisional choice of p=3 mm with no field wire on the basis of the smaller RMS of full-hemisphere electron-count asymmetry**. This cross-check concerns a **single physical chamber (MS0 / DET21)**, not the full station, and it is distinct from the full34 *crossing-level* values in §3.2 and the *station-parent* values in §2.

The complete selected source consists of **100,000 generated J/ψ parent slots in two native Geant4 input chunks**. Only **10,803 parents actually crossed this DET**. For each such parent, aggregate classified upper/lower primary-electron endpoints across all of its gas crossings in DET21, then compute **one** valid `A_hemi = (N_up-N_down)/(N_up+N_down)` for that parent. A zero denominator is undefined, not A=0. Geometry uses p=3 mm and the matched noFW/FW75 **LUT500** maps; no avalanche multiplication or cathode signal is included.

| MS0 / DET21, parent-level P(A), p3 LUT500 | No field wire | Grounded FW75 |
|:--|--:|--:|
| Parents with gas crossing | 10,803 | 10,803 |
| Valid parent A | **10,767** | **10,761** |
| Undefined parent A | 36 | 42 |
| Mean A | −0.0017 | −0.0030 |
| **Population RMS(A)** | **0.3616** | 0.3972 |

**Numerical provenance:** the DET21 numbers below were transcribed from the 6 November 2026 dashboard screenshot. The original frozen `summary.json` and bin-level histograms were **not independently re-read in preparing this update**; the plotting command below is provided to regenerate and check the exact distribution locally.

**RMS contrast:** `RMS_FW75 − RMS_noFW = 0.0356`; noFW gives a **~8.96% lower RMS relative to FW75** (equivalently FW75 is ~9.85% higher than noFW). The two means are both near zero, so the observed difference is primarily in **width**, not a sizeable bias in the mean. The dashboard reports **10,761 paired parents with both A defined** and mean paired `ΔA = A_FW75 − A_noFW ≈ −0.0011`; that paired *mean shift* is not the RMS effect. The unequal valid populations do not invalidate the descriptive comparison, but a matched-parent uncertainty analysis should precede a significance claim.

**Decision relevance and limits:** this local P(A) corroborates the selected no-FW *electron-transport* baseline because a narrower full-hemisphere A distribution is desirable for balanced upper/lower primary-electron counting. It **does not** establish smaller *reconstructed spatial-resolution RMS*, smaller measured cathode-charge asymmetry, or better readout performance. Do not extend the DET21 result to every MS0 module or all stations without checking their own distributions. The restricted cathode-facing A_TB still exhibits the opposite trend in §3.2.

**Reproducible Figure 3b — MS0 / DET21, parent-level P(A).** Generate the *actual 41-bin histograms* directly from the frozen `summary.json` using the new source-only plotter; do not digitize this screenshot or approximate either distribution by a Gaussian:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/plot_pi_ms0_det21_parent_asymmetry.py \
  --summary test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/summary.json \
  --station 0 --det 21 \
  --output ms0_det21_parent_A_p3_LUT500.png
```

The script reads existing JSON only and refuses partial-exposure data; **no Geant4, Heed, Garfield or LUT production is rerun**. Once its PNG is reviewed, add it to the notes with a caption that states *parent-level MS0/DET21*, **10,803 exposed parents**, **p3 LUT500**, and the two valid samples. The image itself is not yet part of this commit; plotting code and exact summary inputs are supplied rather than a synthetic/guessed histogram.


## 5. Avalanche model reference code and bookkeeping

### 5.1 Existing reference workflows already completed in the repository

The existing reference research code is in `test/mwpc-avalanche/`, documented in its `README.md`. These are previous validated *methodological prototypes*, largely using **p=4 mm and asymmetric 2+4 mm gaps** by default; they **are not yet production measurements of the present p=3 mm, ±3 mm, noFW baseline**.

- `phaseA_gain_fluctuations.cxx` / target `mwpc_phase_a_gain`: independent **single-electron AvalancheMicroscopic** ensembles. Outputs `avalanche_electrons`, `avalanche_ions`, **`collected_on_wire`**, attached and other endpoints, eventwise collection times; prints unconditional effective gain mean/RMS/quantiles and zero-collected fraction. Has explicit avalanche-size cap and detects capped events.
- `phaseA_heed_avalanche.cxx`: **Heed muon clusters -> all primary electrons -> microscopic avalanches -> collected charge per anode**. Already captures seedwise and trackwise electron/ion counts and wire multiplicity, but its default reference geometry is p=4 mm, gap 2+4 mm, no FW.
- `phaseB_real_avalanche_signal.cxx`: actual single-electron microscopic avalanche with **ionisation-birth callback**, real positive-ion drift and induced electron/ion strip currents. A one-seed demonstrator, **not** a full muon hit.
- `phaseB_avalanche_angle_scan.cxx` / `analyze_stageB_avalanche_angle_scan.py`: scans avalanche-seed azimuth, records ionisation birth angles/centroids, ion multiplicity, and integrated **upper/lower complete-cathode induced charges** versus seed azimuth. The source defaults to **p=4 mm, gap 2+4 mm**. It skips zero/failed avalanches in the output; **do not** derive unbiased gain from its accepted-only subset.
- `MWPCSignalUtils.h` and `MWPCStereoReadout.h` implement reusable signal/strip infrastructure.
- `test/mwpc-stage-c/PI_FAST_LUT_TRANSPORT_ONLY_README.md` documents prior fast **AvalancheMC** transport-only (attachment OFF, 5 μm step) comparisons; **these are not avalanche-gain simulations**. Do not infer avalanche multiplicative gain or signal charges from this transport benchmark.

The existence and prior validation of these workflows is supported by repository documentation, but the original numerical CSV/run manifests for Phase-A gain or Stage-B ion signal **have not been rechecked against the local WSL outputs in this review**. No numeric gain is claimed here.


### 5.2 What does a Geant4 muon hit tell us about avalanche gain?

Geant4's muon gas hit gives the **entry/exit trajectory, gas path and deposited energy**, but no avalanche multiplication. Heed turns a crossing into **ionisation clusters and N_primary primary conduction electrons**, with start space/time coordinates. Garfield **AvalancheMicroscopic** must then track each primary seed in the chamber's high field, where it can create additional e⁻–ion⁺ pairs. All offspring form one microscopic avalanche for that seed.

For one seed k:
- `N_pairs,k`: count of avalanche ionisation electron–ion **pairs generated**, best monitored using the ionisation callback; **not** initial Heed primary pairs.
- `N_e,aval,k`, `N_i,aval,k`: bookkeeping numbers returned by `GetAvalancheSize(ne, ni)`, whose detailed counting conventions/attachment effects must be inspected.
- `N_e,anode,k`: **number of avalanche electron endpoints actually collected at an anode**; a different number, and the observable used for effective charge gain.
- `G_eff,k=N_e,anode,k/1` for one initial primary electron, **including zero-response seeds**.

For a muon gas crossing with `N_primary` Heed electrons:
\[
G_{\mathrm{eff,hit}} =
\frac{\sum_{k=1}^{N_{\rm primary}}N_{e,\mathrm{anode},k}}
     {N_{\rm primary}}
\quad(N_{\rm primary}>0).
\]
Average gain across crossings **must specify the weighting**: pooled collected electrons / pooled primary seeds vs equal-hit average of `G_eff,hit`. Also report conditional gain for primaries surviving to the multiplication region separately. Do **not** use `N_e+N_i` (a double counting of ionisation products) as a gain denominator. Do not treat a Garfield upper avalanche size cap as a physical gain value.


## 6. Historical symmetric-gap avalanche pilots (NOT a chamber-design verdict)

These are very small *single-seed, fixed-position* tests at ±3 mm cathodes, **not** a native muon-crossing response nor a gain-matched geometry comparison. Earlier 2+4 mm code/output belongs to a different geometry and must not be used for the new ±3 mm simulations.

### 6.1 First symmetric 3+3 mm avalanche smoke5 (28 October 2026)

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


### 6.2 Completed 100-seed symmetric 3+3 gain pilot (28 October 2026)

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


## 7. Literature and interpretation of the field-wire gain shift

### 7.1 Literature interpretation of FW gain (single-position pilot only)

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


## 8. Validation debt and reproducibility checklist

Before any quantitative recommendation goes beyond the provisional statement in the supporting note:

- [ ] Re-read frozen `summary.json` and `manifest.json` and compare all MS0 values, raw up/down totals and 41 bins against the screenshot.
- [ ] Regenerate exact whole-MS0 distribution; preserve checksum and source identifiers. Keep the traced SVG marked illustrative.
- [ ] Obtain parent-level bootstrap confidence intervals for RMS(FW75)−RMS(noFW) under **matched J/ψ parent resampling**, accounting separately for finite LUT histories.
- [ ] Compare asymmetry distributions conditional on similar electron multiplicity; investigate parent-position/angle dependence and potential selection effects.
- [ ] Verify provenance and timing definition of the paired 18k-cell Garfield study and the completed native timing report before quoting a timing-precision improvement.
- [ ] Maintain distinct labels for LUT500 parent, LUT20 crossing, and single-electron avalanche observables.
- [ ] Verify actual field geometry and near-anode field strength on the new ±3 mm four-geometry gas-gain studies.
- [ ] Do not extrapolate the small 100-seed gain tests to signal charge or final effective gain without spatial/gain-matched studies.

## 9. Source registry

### Original source references

- `test/mwpc-stage-c/README.md`: 21 Oct 126k six-geometry Garfield pilot and conditional drift times.
- `test/mwpc-stage-c/pi_pitch_fieldwire_paired_analysis.py`: matched-birth collection/time audit.
- `test/mwpc-stage-c/pi_lut500_muon_crossing_asymmetry.py`: exact crossing-level A, raw IDs.
- `test/mwpc-stage-c/PI_FAST_LUT_TRANSPORT_ONLY_README.md`: AvalancheMC is a **transport-only** candidate, not multiplication.
- `test/mwpc-avalanche/README.md`: existing Phase-A and Phase-B scope; older p=4 mm and 2+4 mm reference geometry.
- `test/mwpc-avalanche/phaseA_gain_fluctuations.cxx`, `phaseA_heed_avalanche.cxx`, `phaseB_real_avalanche_signal.cxx`, `phaseB_avalanche_angle_scan.cxx`: existing tested methods.
- Garfield++ User Guide and ROOT/Garfield maintainers' clarification of `GetAvalancheSize` vs electron endpoint-count and the role of attachment: https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf and https://root-forum.cern.ch/t/difference-between-getavalanchesize-and-getnumberofelectronendpoints/50834

All tables here use the existing full-cell Garfield pilot numbers from repository documentation or the exact full34 ZIP provided earlier. **No new Garfield avalanche physics run was carried out in producing these notes.**


Other existing source documents that are **not merged or overwritten**:
- [Parent-event geometry selector](PI_PARENT_ASYMMETRY_GEOMETRY_SELECTOR.md).
- [Muon-driven B1 detailed protocol](PI_MUON_DRIVEN_AVALANCHE_B1_PLAN.md).
- [Provisional original combined note, preserved as archive](PI_P3_NOFW_COMBINED_WORKING_NOTE_ARCHIVE_2026-11-06.md).

## 10. Append-only results register for future updates

Maintain one row per **completed, auditable dataset** and link to its local/cloud manifest; do not replace historical observations with later plots.

| Date | Geometry / statistical unit | Observable | Dataset and SHA / manifest | Validation status | Relation to decision |
|:--|:--|:--|:--|:--|:--|
| 2026-10-21 | p3 noFW / FW75; 18k microscopic starts/config | collection and conditional drift time | `pi_pitch_fieldwire_pilot_v1/pilot/` | reported in README; original full dataset not independently rechecked here | supports transport-level baseline |
| 2026-10-28 | p3/p4 noFW/FW75; fixed one-seed source | avalanche gain, smoke5 | cloud folder in §6.1 | geometry smoke; n=5/config | no design inference |
| 2026-10-28 | p3/p4 noFW/FW75; fixed one-seed source | avalanche gain, n=100 | cloud folder in §6.2 | very small pilot; not gain matched | no design inference |
| 2026-11-06 | p3 noFW/FW75; complete MS0 parent | P(A), 41 bins, RMS | `pi_lut500_parent_asymmetry_root_v1/` | dashboard screenshot inspected, exact bins pending verification | primary provisional evidence |
| TBD | — | — | — | — | — |

**Update rule:** append completed run metadata, exact observable definition, observed numbers, uncertainties, caveats and the decision impact; leave historical entries intact. Analysis tasks and future decision gates belong in the [follow-up plan](PI_P3_FIELD_WIRE_SIGNAL_RESPONSE_FOLLOWUP_PLAN.md), not in this log.

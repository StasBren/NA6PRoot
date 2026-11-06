# NA60+ MWPC — field-wire follow-up and signal-response research plan
**Internal living work programme | Version 1.0, 6 November 2026 | NA60+ Simulation — Work in progress**

**Status:** planned checks and acceptance criteria; NOT evidence that a task has been performed or a detector design has been approved.

**Three-document set:**
1. [Supporting note: provisional no-FW selection](PI_P3_NO_FIELD_WIRE_SUPPORTING_NOTES_AND_AVALANCHE_PLAN.md) — the scientific question, whole-MS0 result and current decision.
2. [Technical evidence / reproducibility log](PI_P3_FIELD_WIRE_TECHNICAL_EVIDENCE_LOG.md) — data provenance, all numeric checks and prior work.
3. **This document** — what remains to be established and how future results change the decision.

[Historical combined working note (archive)](PI_P3_NOFW_COMBINED_WORKING_NOTE_ARCHIVE_2026-11-06.md). Further pre-existing task-specific guides, especially [the true muon-driven B1 protocol](PI_MUON_DRIVEN_AVALANCHE_B1_PLAN.md), remain as implementation references and are not superseded.

## 1. Scientific decision tree and intentionally unresolved questions

**Current limited decision:** use **p=3.0 mm with no field wires** as the *working reference geometry* for the next microscopic avalanche and induced-signal calculations. The main evidence is the smaller **parent-event RMS of primary-electron hemisphere-count asymmetry** in the complete MS0/LUT500 result, supported by higher anode collection in the transport studies. These observations alone do not establish better real cathode signal, spatial resolution, double-hit resolution, triggering, or efficiency.

**Question 1 — field wire at the fixed p=3 mm baseline:** does adding a grounded 75 μm midpoint field wire improve or degrade *effective readout performance* after avalanche multiplication, positive-ion drift, cathode weighting fields, strip geometry and front-end response? The full-hemisphere primary count RMS favors noFW, but restricted cathode-facing-sector A_TB gives the opposite trend. A realistic charge study must resolve this.

**Question 2 — 3 mm versus 4 mm anode pitch:** **OPEN / DEFERRED UNTIL SIGNAL RESPONSE.** The working p3 baseline is *not* a validated preference over p4. No pitch ranking will be made from the old isolated fixed-position gain ratios, electron-count asymmetry alone, or a single fixed-voltage test. Return to p3 versus p4 when the signal response (including comparable effective gain) can be computed for both and when a relevant reconstruction metric is available.

**Question 3 — sufficient detector performance:** establish whether a candidate meets the experiment's required *position resolution, efficiency and occupancy/multiple-hit behaviour*, not merely a qualitative symmetry criterion. Treat the detector requirements in [A. Milov, *Readout for the MS modules of NA60+/DiCE* (30 Aug 2026)](../../../README.md) and the collaboration experiment proposal as separate design requirements to be verified against the appropriate reconstructed observables. No current A RMS is a single-hit position resolution.

### Geometry and physics comparability lock

All follow-up ±3 mm cathode studies should identify **gas, P, T, B, anode diameter, HV, pitch, field-wire diameter/bias, readout geometry, ion mobility, Penning setting, time window and avalanche cap**. Initial matrix:

| ID | Anode pitch | Grounded midpoint FW | Role |
|:--|--:|:--|:--|
| P3_NOFW | 3.0 mm | absent | working primary baseline |
| P3_FW75 | 3.0 mm | 75 μm | fixed-pitch field-wire comparison |
| P4_NOFW | 4.0 mm | absent | **reserve for signal-response pitch study** |
| P4_FW75 | 4.0 mm | 75 μm | pitch–field-wire interaction after response readiness |

Baseline simulation settings inherited from the previous study: 30 μm anode, symmetric v=±3 mm cathodes, +1800 V, Ar/CO₂ 70:30, 293.15 K, 760 Torr, B=0 T. The **+1800 V case compares at equal voltage, not equal gain**. For meaningful signal comparison include a *gain-matched HV scan* in addition to the fixed-voltage benchmark, and do not confuse preliminary gain with a prototype-calibrated measurement.

## 2. Observables, populations and non-negotiable distinctions

All asymmetries with zero denominator are *undefined/null*, not zero. Preserve parent IDs, gas crossing IDs, initial track angle, path length and the exact population used in each comparison.

| Stage | Observation | Definition / core outputs | Why it matters |
|:--|:--|:--|:--|
| Primary transport | `A_primary,hemi`; `A_primary,TB`; collection probability; drift time | Count **pre-avalanche** anode-arrival surface angles, classified by hemisphere versus ±45° cathode-facing sectors | existing screening criterion; cannot predict measured cathode charge |
| Avalanche ion creation | `A_ion,hemi`; `A_ion,TB`; `N_ionbirth` | Use microscopic ionisation-birth coordinates, not final ion-drift positions | gain-weighted spatial distribution of positive-ion origins |
| Avalanche electron endpoints | `A_e,aval`; `G_eff`; zeros and tail | Count **all** collected avalanche offspring on anodes, not one count per initial Heed seed | incorporates stochastic multiplication and attachment |
| Induced cathode charge | `Q_+(T)`, `Q_-(T)`, `A_Q(T)` | `A_Q(T)=(|Q_+|-|Q_-|)/(|Q_+|+|Q_-|)` from Shockley–Ramo weighting fields and electron/positive-ion trajectories | first direct counterpart to upper/lower cathode signal balance |
| Strips + electronics | `Q_U(t)`, `Q_V(t)`, cluster amplitude/timing, reconstructed hit | Real U/V strip weighting, shaping, threshold, noise and time window | first stage where reconstructed resolution and detection criteria can be assessed |

Do not call ion-birth angles cathode induced signals, or the RMS of the primary A distribution the MWPC spatial resolution.

## 3. Work packages and gates

### WP0 — freeze and validate present transport evidence (read-only)

**Input:** completed original 100k-parent p3/LUT500 `summary.json`, `event_counts.csv.gz`, `manifest.json`; original native crossing sources and frozen Garfield LUTs.

**Tasks:** Verify SHA and full exposure; regenerate whole-MS0 41-bin P(A) from **`station_parent`** (not a single DET), and check every printed RMS/mean, classified electron count, valid/undefined count and paired population. Keep p3/LUT20 crossing distributions separate from p3/LUT500 station-parent. Review a multiplicity-matched control for the observed ~3.43% difference in endpoint counts.

**Acceptance gate:** provenance and bin checksum pass, no partial events mislabelled as complete, and uncertainty method documented. **Output:** exact-data figure replacing the provisional screenshot trace, plus a validation entry in the technical log. No physics run is required to start.

### WP1 — controlled electric-field and gain audit

**Inputs:** verified symmetric ±3 mm analytic field cells for the two fixed-pitch p3 geometries, eventually the p4 cases. The earlier **2+4 mm** demonstration is *methodological only*; its numerical gain is not transferrable.

**Checks:** field-wire bias, wire/cathode positions and labels, symmetry, near-anode `|E|` at radial distances 20/30/50/100/200 μm, convergence with number of analytic wires; gas/transport setup and avalanche cap status.

**Comparison:** fixed +1800 V first, then scan voltage to compare at similar collected effective gain. Keep zero-collected avalanches in unconditional gain and attach/error fractions; do not calculate mean gain from success-only subset. For source positions sample top, bottom and off-axis starts, and vary distance to avoid single-point overinterpretation.

**Acceptance gate:** all four geometries physically audited, source coordinate conventions verified and effective-gain distributions (with tail/cap flags) reproducible. The existing 5- and 100-seed pilots **already happened** as feasibility checks; they do not pass this gate.

### WP2 — microscopic avalanche, muon-driven B1/B2 (scientific priority)

**Actual unit of analysis:** one **real native Geant4 muon gas crossing** with true track position and direction, whose primary ionisation is generated with **TrackHeed** (or clearly labelled validated equivalent, not silently substituted). For each primary electron, use Garfield **AvalancheMicroscopic** and record all offspring.

**First measurable level (B1):** callback-recorded avalanche **ionisation births**, `(u,v,z,t)`, radius/angle relative to the correct physical anode and counts in the upper/lower hemispheres and cathode-facing ±45° sectors. Compute `A_ion,hemi` and `A_ion,TB` per crossing, not only per seed. Include empty/failed/attachment/capped events.

**Second level (B2):** count avalanche **electron anode endpoints** and derive `A_e,aval`, `G_eff,hit`, collected charge and per-wire multiplicity. Report both pooled-electron and equal-crossing-weighted averages with their denominators. `GetAvalancheSize` bookkeeping is not automatically the same as collected anode electron count.

**Control sample:** begin with bounded original muon-crossing subsets spanning track angle, entry point, station and position, keeping the exact same input tracks and primary births for competing geometries. Reuse existing `phaseA_heed_avalanche.cxx`, `phaseB_avalanche_angle_scan.cxx` and callbacks as methodological references without overwriting historical files. Before a full native campaign, inspect run time, timeouts, avalanche-size caps and ion bookkeeping. Follow the detailed [true-muon B1 protocol](PI_MUON_DRIVEN_AVALANCHE_B1_PLAN.md).

**Acceptance gate:** eventwise `A_ion` and `A_e,aval` outputs are complete for a bounded, provenance-checked muon set, with zeros, cap/selection rules and independent physical observations identified. Fixed single-seed angle scans are diagnostics **only**.

### WP3 — induced cathode signal B3

For the **same muon-crossing ensemble** as WP2, drift positive ions (plus relevant electron contributions), compute weighting potentials/fields, and integrate induced **upper/lower cathode charge** `Q_±(T)` for several explicitly specified observation times. Check numerical convergence, charge conservation conventions, ion mobility and sensitivity to truncating ion tails.

**Acceptance gate:** `A_Q(T)` is calculated for valid signal populations; zero/undefined cases, gain matching, time windows and corrections documented. Compare correlations with `A_primary`, `A_ion`, `A_e,aval`. A favorable primary count asymmetry need not imply favorable `A_Q`.

### WP4 — realistic strip/ASIC signal response and geometry decision

Only once cathode plane signals are credible, proceed to **actual U/V strip weighting and response**, charge-sharing across strips, shaping, digitization, thresholds, noise and representative occupancy. Quantify (as available) hit efficiency, cluster-size distribution, timing, hit-position bias and RMS, charge-signal asymmetry and false associations/double-hit separation.

**This is where the p=3 mm versus p=4 mm decision is reopened.** Compare P3_NOFW, P3_FW75, P4_NOFW and P4_FW75 using common muon input and, where feasible, **gain-matched** operating points. Document any cost/complexity differences separately from the physical response. A stable decision should rely on a metric tied to tracking/readout requirements rather than minimizing an auxiliary A by itself.

**Acceptance gate:** a side-by-side four-configuration signal-response study with explicit physical metrics, uncertainties and source provenance. Until then, **3 mm remains a reference pitch, not a selected optimum**.

## 4. Analysis controls and uncertainty strategy

- **Matched parent/bootstrap:** for station-level P(A), resample independent *J/ψ parents* with both FW/noFW configurations evaluated in each bootstrap replicate; report CI for `RMS_FW−RMS_noFW` and how invalid A is handled. For crossing-level studies, consider clustering crossings within parent rather than pretending independence.
- **Finite LUT error:** separate repeated 500-history LUT sampling uncertainty from genuinely new microscopic histories and parent-event variation; test selected source cells and cross-validations rather than naive `1/sqrt(N_parent)` errors.
- **Electron multiplicity:** compare fixed/stratified `N_primary` or `N_collected` and native track-angle/position subsets before attributing an RMS difference purely to field topology.
- **Conditional selection:** distinguish all produced primary electrons, collected primaries, positive-gain avalanches, complete ion-drift events and threshold-passing reconstructed hits.
- **Physics checks:** top-bottom symmetry, electron endpoint boundaries, gain/charge conservation conventions, zero-response handling, time-window dependence, wire identities, field convergence and angle mirroring.
- **No apples-to-oranges:** station-parent p3/LUT500 versus crossing p3/LUT20, 3+3 versus 2+4 gaps, fixed-voltage versus gain-matched operating point, and ion-birth counts versus induced charge must always be separately labelled.

## 5. Evidence update contract

Each future completed result should be appended to the [technical log](PI_P3_FIELD_WIRE_TECHNICAL_EVIDENCE_LOG.md) with **date, dataset SHA/manifest, software revision, scope, geometry, sample denominator, full selection definition, observable and uncertainties**. Then revise only the short supporting-note conclusion if the interpretation genuinely changes; retain superseded findings with dates.

For each work package use: **planned → running → data available → QA passed → interpretation reviewed**. Only *QA passed* evidence should update the main decision. Do not silently promote feasibility smokes or a pending output to complete science.

| Work package | Present status (as of 6 Nov 2026) | Required next result |
|:--|:--|:--|
| WP0 | Complete dashboard screenshot available; exact-bin source verification pending | independently checked original MS0 41-bin plot and data provenance |
| WP1 | fixed-position 5- and 100-seed feasibility studies documented; not a spatial/gain-matched performance study | spatial/gain comparison and electric-field QA |
| WP2 | methodology and other campaign notes exist; this three-document review did **not** verify a complete corresponding comparison | verified true-muon crossing-level ion-birth and avalanche-endpoint distributions |
| WP3 | earlier one-seed methodology exists only | new symmetric-gap upper/lower cathode `A_Q(T)` on muon ensemble |
| WP4 | not established in current source set | first matched signal-response result, including **3 versus 4 mm** |

## 6. Figures to add when supported

1. **Already available** — whole-MS0 parent-event `P(A_primary,hemi)`, p3/LUT500 noFW/FW75, with original-bin plot replacing the interim screenshot trace.
2. **Add after WP0/WP1** — anode collection and conditional drift-time comparison with uncertainties and source population.
3. **Add after WP2** — gain distributions, zero collection, `A_ion`, `A_e,aval` on actual matched muon crossings.
4. **Add after WP3** — `Q_+(t)`, `Q_-(t)` and `A_Q(T)` time-window dependence.
5. **Add after WP4** — p3 vs p4 and noFW vs FW75 **signal response** and hit-reconstruction performance.

Follow established plot conventions (white canvas, labelled axes/units/normalization, blue solid noFW and red dashed FW75, **NA60+ Simulation — Work in progress**). Do not use ATLAS experiment branding, even if plots are ATLAS-inspired.

## 7. Decisions and open issues register

| Date | Decision / open question | Basis | Revisit when |
|:--|:--|:--|:--|
| 2026-11-06 | **Provisionally use p3 noFW for subsequent simulations** | whole-MS0 parent-level electron-count RMS 0.3106 vs 0.3500, plus primary collection comparison | matched avalanche and `A_Q(T)` results |
| 2026-11-06 | **Do not choose between p3 and p4** | present P(A) and single-point gain do not establish the relevant strip/ASIC performance | WP4 signal-response comparison, preferably at gain-matched operating points |
| 2026-11-06 | **Do not claim final cathode-charge symmetry** | full-hemisphere and restricted-sector primary asymmetries move in opposite directions | WP3 induced-charge transport |
| TBD | — | — | — |

This is a *living plan*. Extend the register and checkpoints as results arrive; do not rewrite historical observations as if they had been known on 6 November 2026.

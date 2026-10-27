# PI follow-up — p=4 mm LUT500 and per-muon gas-crossing asymmetry

Branch: `feature/mwpc-stage-c`. These are **independent, non-destructive
sidecar calculations**. Never replace the previous p=3 mm LUT500, the
six-geometry LUT20, parent-event dashboard or original 2M source results.

## Part A — Garfield p=4 mm, no field wire and field wire 75 μm

New isolated directory:
`test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p4_v1/pilot/`.

Frozen physics: Ar/CO2 70:30, 293.15 K, 760 Torr; B=0 T;
+1800 V 30 μm anodes; 3+3 mm cathode gaps; pitch 4.0 mm;
grounded 75 μm field wires between anodes or none.
Grid: 40 u × 30 v boxes; 500 independent microscopic Garfield
drifts per box per configuration (600,000 trials/config; 1.2 million total).
Field-wire diameter is a **provisional study setting**, not final hardware.

The existing `pi_pitch_fieldwire_pilot.py` is invoked with
`--only-pitch 4`, an isolated output and a new fixed seed.
Completed 15,000-trial u-columns are checked and reused; interrupted
columns are safely recomputed. Per-column SHA and timing sidecars are
retained. The old p=3 mm 500-trial LUT is not read or overwritten.

Run:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500_p4.sh check
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500_p4.sh start
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500_p4.sh status
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500_p4.sh verify
```

**Start triggers a genuine Garfield computation** on the user's WSL
machine (low nice priority). Do not assume completion from GitHub alone.
Avoid competing for CPU with the existing 6M production; `status`
is observational, whereas `verify` inspects every completed u-column
and every assembled electron row. The original environment activation,
Garfield linker check and build logic are inherited from
`run_pi_pitch_fieldwire_pilot.sh`.

## Part B — true per-muon gas-crossing asymmetry

New analysis:
`pi_lut500_muon_crossing_asymmetry.py`, with synthetic fixture:
`test_pi_lut500_muon_crossing_asymmetry.py`.

The existing 100k-parent station+DET J/ψ dashboard is
**parent-event** statistics. For this follow-up, the statistical unit is
instead **one physical muon crossing of one chamber gas volume**.
Original `source,chunk,event,track_id,hit_index,station,det_id` and
signed `theta_u,theta_w` remain attached to each crossing.

Read-only inputs:

* Native original Geant4 muon gas crossing CSVs (SHA provenance);
* unchanged empirical Heed cluster library and stable source+chunk+event+hit+
  track primary birth sampler;
* p=3 mm noFW/FW75 frozen LUT500 (already completed);
* optional new p=4 mm noFW/FW75 frozen LUT500 (only after Part A verifies).

Both field-wire geometries get the same birth positions and LUT lookup
quantiles for each crossing. Their microscopic random trajectories
are separate precomputed LUT trials, reused across events. No Garfield,
Geant4 or fresh Heed library jobs are launched by Part B.

### Distinct observables

`phi_hit` is **anode endpoint position azimuth**; NOT electron
approach-velocity direction and NOT a simulated cathode signal.

1. Original half-circumference asymmetry:
   `A_hemisphere=(N_up-N_down)/(N_up+N_down)`, with
   `up = 0<phi<180`, `down = -180<phi<0`. Exact horizontal
   boundary hits are excluded.
2. Illustrated cathode-facing sectors:
   `top=[45,135)`, `bottom=[-135,-45)` degrees,
   remaining surface = side.
   `A_top_bottom=(N_top-N_bottom)/(N_top+N_bottom)`.
   `f_top_bottom=(N_top+N_bottom)/N_anode_collected`.

**A_top_bottom is a geometric endpoint proxy, not an induced-charge
ratio.** This ±45° definition is a study choice and needs sector
sensitivity tests before making hardware claims. The two asymmetries
MUST NOT be silently interchanged. Undefined denominators are
stored as blank in raw CSV and null in summaries, NEVER zero.

### Commands

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only

# Always start with tests and a bounded p3-only pilot.
bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh check
bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh pilot3
bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh full3

# Only once independent p4 LUT500 'verify' PASSES:
bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh pilot34
bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh full34

bash test/mwpc-stage-c/run_pi_lut500_muon_crossing_asymmetry.sh status
```

The pilots take the first 200 accepted muon gas crossings **per
station** from one 50k-parent input chunk, and are explicitly partial.
The full modes read 2 × 50k generated J/ψ parents, every gas crossing
at all six stations. A 100k-generated-parent exposure is NOT the same
as 100k independent chamber crossings.

Separate output folders under:
`test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/`
(`pilot3`, `full3`, `pilot34`, `full34`).
An existing dataset or interrupted `.incomplete` is **never
overwritten**. Inspect/backup before a restart.

### Deliverables

* `crossings.csv.gz`: one crossing per geometry, preserving parent and
  muon-track IDs, local entrance angles, hit counts in both hemispheres,
  top/bottom/side, A values and their denominators.
* `summary.json`: station and physical DET distributions; valid/
  undefined populations, means, RMS; field-wire matched
  `A_FW-A_noFW`; original angle-bin conditional distributions;
  input Garfield LUT SHA and native CSV SHA.
* `hemisphere_A_p3.png`, `top_bottom_A_p3.png` and corresponding
  p4 figures when p4 is present: valid-crossing fraction per ΔA bin.
* `up_down_joint_p3_noFW.png` and `...FW75.png`: per-crossing
  upper versus lower electron-count joint distributions.
* `paired_delta_A_p3.png`: within-crossing FW-minus-noFW hemisphere A.
* `RMS_A_vs_theta_u_p3.png`, `RMS_A_vs_theta_w_p3.png`:
  RMS versus signed entrance-angle bin, with sparse bins suppressed.
* Same figures for p4 if requested.
* `manifest.json`: files SHA256 and input provenance. Descriptive
  means and RMS are NOT model or sampling error estimates.

### Presentation (ATLAS plotting CONVENTIONS, never ATLAS branding)

White plot canvas; black axes and high-contrast labels; axis quantities
and units, top/right ticks; blue solid for no field wire and red dashed
for grounded FW75; shared scales for six station panels; legends and
actual denominator/bin width; no ROOT automatic plot-statistics boxes.
Label results `NA60+ Simulation — Work in progress`.
The plots are NOT official ATLAS figures.

### Physics boundaries

The simulation is primary-electron *anode-collection position*
statistics, **not** avalanche gain, induced charge or real cathode
strip amplitudes. In particular, an electron counted on the "top" of
the anode is not equivalent to an ion-induced current confined to the
top cathode. True cathode asymmetry requires Shockley–Ramo weighting
fields, microscopic avalanche ion-birth locations, ion drift,
time integration and readout segmentation.

Do not interpret the large raw sample of collected electrons as
independent Garfield microscopic trials. Finite LUT500 uncertainty,
Geant4/Heed parent-event correlations, field-wire geometry uncertainty
and sector-boundary dependence still require separate validation.

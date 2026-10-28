# NA60+ muon-crossing asymmetry: rebuild, not parent pooling

**Scientific definition:** One native Geant4 physical muon crossing of ONE
chamber gas volume = ONE observable A. The 100k generated J/ψ parents
are **input exposure**, not histogram entries. The old
`pi_lut500_event_asymmetry_root_dashboard.py` intentionally pools
crossings per J/ψ parent and station/DET, so its P(A) is a **different**
observable and must not be relabelled as muon-crossing P(A).

The new physical fold was already implemented in
`pi_lut500_muon_crossing_asymmetry.py`. For each crossing:
Geant4 entry/exit & truth IDs -> stable Heed clustered ionization ->
one frozen Garfield LUT20/LUT500 response per primary electron ->
count upper/lower endpoint `phi_hit` -> compute A. Repeat over
all valid crossings, equally weighted, with exact original
source+chunk+event+track+hit+DET identity retained in
`crossings.csv.gz`.

- `A_hemisphere=(N_up-N_down)/(N_up+N_down)`
  for `0<phi<180` vs `-180<phi<0`.
- `A_top_bottom=(N_top-N_bottom)/(N_top+N_bottom)`
  for geometric cathode-facing sectors `[45,135)` and `[-135,-45)`
  degrees. Electrons arriving on side sectors are not in this
  second denominator.
- Undefined denominators are excluded from P(A) and counted explicitly.
  The equal-weight crossing histogram uses 41 bins over [-1,+1];
  each configuration is normalized by its own number of valid crossings.
- These count primary-electron endpoint **positions on anodes**,
  NOT avalanche gain, ion-induced charge, cathode currents, or strip signals.
- Correlations across crossings of the same parent, and reuse of the
  same finite microscopic LUT, must be accounted for before claiming
  independent statistical errors.

## New interactive dashboard

`pi_crossing_asymmetry_dashboard.py` renders
`crossing_dashboard.html` **from saved** crossing-level `summary.json`,
`crossings.csv.gz`, and the original SHA-pinned `manifest.json`.
It does not run Geant4, Heed or Garfield and never changes these
existing physics files. The dashboard includes:
MS0–MS5, physical DET map (clickable), anode pitch, hemisphere vs
cathode-facing sector A, probability or counts, and chamber table.
It explicitly displays "Muon gas crossings", not parents, for
histogram denominators. It follows ATLAS plotting conventions
(white graph, black axes and tick marks, blue solid no FW, red dashed
FW75), but is clearly branded `NA60+ Simulation — Work in progress`,
never ATLAS.

### A. Use the already completed LUT20 p3+p4 full34

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

bash test/mwpc-stage-c/run_pi_crossing_asymmetry_dashboard.sh check
bash test/mwpc-stage-c/run_pi_crossing_asymmetry_dashboard.sh lut20-full34
```

Open (Windows Explorer via `wslpath -w` printed by runner):

`test_runs/mwpc_stage_c/pi_lut20_p4_muon_crossing_asymmetry_v1/full34/crossing_dashboard.html`

This command only renders a dashboard from **existing 989,695
muon gas crossings** (2×50k generated J/ψ input). It will not
repeat the nearly million crossing folds.

### Compatibility with earlier LUT500 crossing summaries

Some pre-LUT20-versioned `full3/summary.json` files have the correct
crossing-level `LUT500` schema, but no explicit
`lut_trials_per_cell` field. The dashboard renderer now infers
500 **only** from the exact versioned schema, requires the matching
`LUT500` provenance and validates the original manifest SHA256 of
both the summary and raw CSV. It does not edit the old physics files.

To retry a previously failed dashboard publication, simply update the
branch and rerun the same `lut500-p3` command. If another metadata
mismatch remains, the validator reports exactly which field differs;
do not delete or recompute the physical fold merely to satisfy a
dashboard metadata check.

### B. Recompute p=3 mm LUT500 with the correct statistical unit

```bash
bash test/mwpc-stage-c/run_pi_crossing_asymmetry_dashboard.sh lut500-p3
```

If the correct crossing-level LUT500 output already exists, it is
reused and not overwritten; otherwise the existing pilot3 and full3
scripts fold the original native Geant4+Heed input with the
already completed p3 LUT500 map. This is an actual crossing-level
recalculation but NOT a new microscopic Garfield LUT job.

Output:

`test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full3/crossing_dashboard.html`

### C. When independent p4 LUT500 is complete & verified

```bash
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500_p4.sh verify
bash test/mwpc-stage-c/run_pi_crossing_asymmetry_dashboard.sh lut500-p34
```

This produces p3+p4 at equal LUT500 resolution under:

`test_runs/mwpc_stage_c/pi_lut500_muon_crossing_asymmetry_v1/full34/crossing_dashboard.html`

The runner never overwrites existing physics results or old dashboards.
An incomplete crossing fold is intentionally NOT silently resumed;
inspect its saved incomplete artifacts first.

### D. Dashboard-only HTML refresh

For a completed source directory `$INPUT`:

```bash
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_crossing_asymmetry_dashboard.py --input-dir "$INPUT"
# To replace ONLY the new HTML when styling changes, never physics:
"$PY" test/mwpc-stage-c/pi_crossing_asymmetry_dashboard.py \
  --input-dir "$INPUT" --replace-dashboard
```

The renderer verifies the physics manifest and SHA256 for both the
summary and raw CSV (streaming read, no materialized rows) and checks
station/physical-DET event exposure and hist bin conservation.
It generates a sidecar `crossing_dashboard.provenance.json`.

## How to interpret the old dashboard screenshot

The old LUT500 p=3 mm MS0 chart with `Parents with crossing: 99,968`
and RMS `0.3106 / 0.3500` is a **parent-event-level station** distribution.
Those RMS values CANNOT be treated as per-crossing asymmetry.
It is scientifically legitimate in its own right, but answers a
different question from the PI's requested per-muon crossing statistics.
The original dashboard remains available for a provenance-labelled
comparison; the corrected new dashboard uses crossing-level numbers only.

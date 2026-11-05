# Full-Garfield J/psi 1000 muons: azimuth and event-level asymmetry

**Geometry:** `p3_noFW` ONLY: 3-mm anode pitch, no field wire,
cathodes at ±3 mm, Ar/CO₂ (70:30), +1800 V, B=0.

**Frozen physics reference:** `test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/`.
The user completed and fully SHA-audited **1000** Geant4 J/psi
parent-independent selected gas crossings, **54,912** actual microscopic
Garfield++ avalanches; **0 size caps**, **2 unsuccessful calls** (both
0 ion births and 1 electron collected, starting microns from anode),
in 50 sealed native chunks. An unsuccessful `AvalancheElectron` boolean
is not silently interpreted as a cap or used as a reconstructed charge.
No new physics is run for these plots.

## Run after pulling branch

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi1000_ion_angles_v1.py

"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
"$PY" -u "$SCRIPT" plot | tee ~/na6p/p3_noFW_jpsi1000_ion_angles.log
~~~

The plotter refuses to run until all 50 native chunks are fully rehashed
and the expected 1000 independent muon IDs are present. It also verifies
`Nbirth/point` reservoir weights and conservation against each event and
the original SHA checkpoint census.

Plots are **new, non-overwriting** output in
`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/plots_ion_birth_azimuth_muon1000_v1/`:

1. `01_ion_birth_azimuth.png`: angle around nearest 3mm anode,
   ion-weighted and equally weighted per nonzero muon, 72×5° bins.
   **Approximate** because up to 128 ion positions saved per avalanche.
2. `02_exact_event_asymmetry_distributions.png`: full-muon
   **exact** full-ion sector asymmetry `A_hemi` and `A_TB` histograms.
   Each valid muon has equal weight.
3. `03_asymmetry_each_of_1000_muons.png`: separate event-indexed
   asymmetries for every original J/psi crossing; undefined omitted.
4. `04_each_muon_azimuth_heatmap.png`: reservoir-estimated azimuth
   histogram for each of 1000 muons (equal event normalization).
5. `05_asymmetry_vs_true_incidence.png`: exact muon asymmetry against
   native Geant4 signed `theta_u,theta_w`.
6. `06_single_avalanche_vs_whole_muon_asymmetry.png`: diagnostic contrast
   of per-primary and whole-muon asymmetry, each separately normalized.
   These are **different statistical units**, never equated.

CSV:
- `muon_exact_asymmetry_1000.csv`: per-parent event identifiers,
  station/angles, Heed multiplicity, full ion counts and both exact A.
- `ion_birth_azimuth_5deg.csv`: probability per angle bin, estimated
  ion-weighted and per-muon weighted.
- `manifest.json`: SHA of all 150 source physics CSVs,
  analysis script, plot+CSV artifacts, source counts and exact A definitions.

Definitions match the **native C++ BirthAccumulator**:
`phi = atan2(v, u - 3mm * nearest_anode_index)`, with C++ `std::round`
ties-away-from-zero behavior:
- `A_hemi = (N_up - N_down)/(N_up+N_down)`, strict upper/lower
  half planes: `phi>0` and `phi<180`, or `phi<0` and
  `phi>-180`; exact boundary excluded.
- `A_TB = (N_top - N_bottom)/(N_top+N_bottom)`, cathode-facing
  angular sectors `[45,135)` vs `[-135,-45)`; side excluded.
- Zero denominator -> **undefined**, not 0.
- `phi` histogram uses per-avalanche uniformly sampled ion births
  weighted by `N_ion_births / min(128, N_ion_births)` to undo
  reservoir truncation. It is not an exact per-bin ion census;
  only sector asymmetries/counts are exact.

**Not calculated:** positive-ion motion, cathode-induced charge,
U/V strip charge sharing, electronics, position resolution.
The original 1000 are one selected crossing per generated J/psi parent
from Geant4 chunk 3, **not the true unstratified all-crossings flux**.
The plots describe ionization **birth positions**, not electron impact
angles and not a Shockley–Ramo cathode signal.

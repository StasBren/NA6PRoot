# NA60+ MWPC — PI-targeted **MUON-DRIVEN** avalanche B1 (3+3 mm)

**Status 28 October 2026:** design specification, not yet executed. This plan **supersedes the single-fixed-electron Stage B1 as the primary physics objective**. Single-seed simulations (including smoke5 / pilot100) are engineering and microscopic validation only.

## 1. Actual scientific question

For muons arriving at the NA60+ MWPC with the **physical Geant4 ensemble of impact positions, signed entrance angles, trajectories, charges and momenta**, what primary-electron avalanches result from their Heed ionization clusters, and how do the avalanche ion-birth/electron-endpoint angular distributions differ across **four MWPC geometries**, all with symmetric **v=±3 mm cathodes**?

| ID | Anode pitch | Midpoint field wire | Priority |
|---|---:|---:|---|
| p3_noFW | 3 mm | none | **baseline** |
| p3_FW75 | 3 mm | grounded 75 µm | control |
| p4_noFW | 4 mm | none | control |
| p4_FW75 | 4 mm | grounded 75 µm | control |

Fixed preliminary gas/HV/B for first comparison: Ar/CO₂ 70:30, 293.15 K, 760 Torr, 30 µm anodes, +1800 V, B=0 T. The local baseline choice is **provisional**; do not choose a detector design from fixed-point gain alone.

## 2. Existing validated parts and important distinction

**Native Geant4 parent tracks:** original frozen J/ψ muon gas crossings from the 100,000-parent native simulation, by `crossing_provenance("Jpsi", chunk)` in `test/mwpc-stage-c/pi_2m_native_six_geometry_fold.py`. A crossing row includes:
`source,chunk,event,hit_index,track_id,station,det_id,pdg,path_length_mm,u_in_mm,v_in_mm,w_in_mm,u_out_mm,v_out_mm,w_out_mm,theta_u_in_deg,theta_w_in_deg,edep_eV` (some fields supplied outside a single CSV). Preserve signed angle and parent identifiers. **Do not reconstruct independent tracks by random angles**.

**The earlier `full34` source is NOT fully microscopic Heed per crossing.** In `pi_six_geometry_native_timing.py::birth_sample`, for each native Geant4 chord, primary clusters are sampled as `N_clusters ~ Poisson(lambda * length)`, multiplicities are drawn from empirical Heed library marks (which originated in `c_heed_library_2p4`), and electron positions are placed along the actual 3D Geant4 chord. This is a reproducible **Heed-calibrated statistical fold**, preserving true muon directions and lengths, but it does not replay the exact `TrackHeed::GetCluster/GetElectron` process for each native muon. Preserve earlier frozen transport outputs without mislabeling them.

**The existing true microscopic Heed→avalanche chain exists already** in `test/mwpc-avalanche/phaseA_heed_avalanche.cxx`: `TrackHeed::NewTrack/GetCluster/GetElectron` followed by per-electron `AvalancheMicroscopic` and anode collected charge. But its current `NewTrack(x0,yStart,0,t, 0,-1,0)` fixes **normal incidence** and its default gaps are 2+4 mm with no midpoint field wire. It needs real-muon import and geometry support; do not use its legacy CSV as a final NA60+ angular ensemble.

**The existing Stage B** `phaseB_avalanche_angle_scan.cxx` already has `RecordIonisation` callback, ion drift and cathode weighting potentials. These are reusable components, not a representative muon-hit event model. The existing Stage A single-seed `phaseA_gain_fluctuations.cxx` was adapted to p3/p4 ±3 mm with FW support. Its fixed electron start `u=0,v=+1mm` is an instrument/debug point only.

## 3. Main implementation: ONE PHYSICAL GEANT4 MUON GAS CROSSING = ONE OBSERVATION

1. Read SHA-validated, immutable native crossing rows; select original gas hits across MS0–MS5, preserve `(source,chunk,event,hit_index,track_id,station,det_id)`. Read actual entrance and exit 3D positions; normalized `d=(r_out-r_in)/|r_out-r_in|` in local `u,v,w`, and `theta_u_in_deg`, `theta_w_in_deg`. Verify local coordinate conventions and entry/exit orientation (not arbitrary midpoint assumptions).
2. Obtain **true muon momentum and charge for Heed**, not one invented shared momentum. **The existing immutable native Geant4 CSV already contains** `pu_in_GeV,pv_in_GeV,pw_in_GeV` (and `px_in_GeV,py_in_GeV,pz_in_GeV`) produced by `export_c_geant_muon_crossings.C`; the earlier Python `REQUIRED_HIT` was only a minimal subset and should not be mistaken for the full schema. Compute `|p_in|` from the existing local components, compare incoming momentum to native chamber chord, and fail validation rather than silently substituting 10 GeV/c. The existing PDG+MCKine joining remains in force. **No extra ROOT re-export or Geant4 production is needed for the small B1 smoke.**
3. Create `TrackHeed` for each **selected native crossing** with `SetParticle("mu-"/"mu+")` and `SetMomentum(p_eV)`; `NewTrack` with that **real** local entrance and unit direction. Generate actual Heed clusters and primary electrons via `GetCluster/GetElectron`. Enforce/verify gas path clipping to the original Geant4 **crossing segment** and physical chamber boundaries, including edge/zero-length hits. Match path length, number of clusters/electrons, and any delta-electron transport assumptions. Retain Heed `(u,v,w,t,E,direction)` electron records.
4. Freeze this **one shared microscopic primary-ionization realization per physical muon crossing**, with a deterministic Heed seed independent of geometry. Retain identical Geant4/Heed primaries across the four pitch/FW candidates. Use separate reproducible microscopic avalanche RNG streams for each `(crossing, geometry, seed)` and preserve caps, failed/zero outcomes. Do not rerun Heed independently per geometry and then mistake random ionization variation for geometry effects. Compare on physically identical selected crossing IDs.
5. For each Heed primary electron, call `AvalancheMicroscopic` in the selected geometry and attach the existing `SetUserHandleIonisation(RecordIonisation)` callback. Verify how the Heed-returned initial electron energy/direction and delta electron transport are passed to avalanche simulation (the legacy `phaseA_heed_avalanche.cxx` hard-codes 0.1 eV at avalanche seed); do not assume this approximation is automatically valid. Record every ionization-birth position and every terminal avalanche electron endpoint. Gas gain is `N_e,collected / N_primary,Heed` and must include all zero-response primary seeds.
6. For ONE PHYSICAL MUON CROSSING, sum upper/down **raw counts across every primary-electron avalanche** first, then form the ratio. A large avalanche contributes proportionally more ion births than a small one, as required for ion-population asymmetry.

```
crossing_id := (source, chunk, event, hit_index, track_id, station, det_id)
Nion_upper(crossing) := SUM_{Heed primary seeds k} Nion_upper(k)
Nion_lower(crossing) := SUM_{Heed primary seeds k} Nion_lower(k)
A_ion,hemi(crossing) := (Nion_upper-Nion_lower)/(Nion_upper+Nion_lower)
A_ion,TB(crossing) := (Nion_top45-Nion_bottom45)/(Nion_top45+Nion_bottom45)
```

In cathode facing-sector counts use exact `phi` angular bins `[45°,135°)` and `[-135°,-45°)` in the local anode-centered `u,v` plane. Use the **nearest physically struck anode**, not a hypothetical field wire center. Test boundary conventions. A denominator of 0 is **undefined**, not A=0.

Retain also per-primary-electron birth count and individual `A_k` for diagnosis, but **do not average A_k to obtain per-muon asymmetry**. Also record amplified electron endpoint `A_e,hemi`/`A_e,TB` and `G_eff,hit`, and keep downstream cathode induced-charge `A_Q(T)` explicitly separate.

## 4. Muon-angle dependence is central, not optional

For each physical chamber crossing, keep `theta_u_in_deg` and `theta_w_in_deg` plus local 3-vector `d` and physical gas path length. Report:

- `P(A_ion,hemi | theta_u, theta_w)`, `P(A_ion,TB | theta_u, theta_w)`, `P(G_eff,hit | angles)` by station and (sampled) DET.
- Mean and RMS versus **signed** angle, with **both** positive and negative sides, and coarse 2D signed-angle bins only when sufficient occupancy.
- Equal-weight individual physical gas crossings as the primary statistical unit. Crossings within a single J/ψ event and multiple stations are correlated; parent-clustered bootstrap is required for global intervals.
- Angle, position, path length and momentum are correlated in the Geant4 ensemble. Preserve the joint distribution; **do not simply draw uniform or independent angle grids** and call it realistic.

Also test negative/positive angle symmetry only for geometrically appropriate mirrored samples, with station/geometry acceptance accounted for. Separate genuine MWPC effects from J/ψ phase-space/acceptance effects.

## 5. Bounded execution sequence (no million-crossing microscopic run)

**Phase A: minimal unit tests.** Fixed single electron, mirrored seeds and field probe for validating geometry, callbacks, counting and sign. This is never a PI physics plot.

**Phase B: real-muon chain smoke.** Run 5–10 SHA-selected real crossing rows with different signed entrance angles and path lengths, each for all four geometries, generating Heed clusters **from the actual track** and all their microscopic avalanches. Compare underlying primary-ionization identity across geometries. Keep complete seed and crossing audit. Measure runtime, cap rates and ion-birth output sizes.

**Phase C: representative real-muon angular ensemble.** Stratified selection from true native J/ψ crossings by station and signed angles (plus positions and momentum, with sampling weights), ideally starting with tens of native crossings per station/angle bin to assess runtime and variance, then scaling on demand. Keep explicit sampling probability/weight; do not mistake stratified equal-size bins for the true angle PDF. Report per-crossing A and dependence on `theta_u,theta_w` plus physical detector acceptance.

**Phase D: Cathode signal.** For a bounded, representative selected set of muon crossings, drift all positive avalanche ions or use a validated response kernel, compute cathode/strip induced charge via Shockley–Ramo for several integration windows. Study correlations between primary transport `A_hemi/A_TB`, ion-birth `A_ion` and true `A_Q(T)`. Match or scan gain separately from fixed +1800V geometry comparisons.

**Performance:** start with a capped small representative sample and estimate CPU/memory before scaling. The prior 100k-parent `full34` involved roughly 989,695 muon crossings and ~55.8 million Heed-calibrated primary-electron draws, but re-running *all* of these with microscopic avalanche multiplication is not an appropriate initial job.

## 6. Deliverables for P.I. / Overleaf

1. True muon angular/path-length distributions in MS0–MS5, with selection definition and parent/chamber exposure.
2. One illustrated native muon gas chord with genuine Heed cluster positions, primary-electron seeds and separate secondary avalanche ion-birth clouds (gas local u-v-w, with cathode ±3 mm).
3. Four-geometry histograms of true crossing-level `N_primary`, `G_eff`, `N_ionBirth`, `A_ion,hemi`, `A_ion,TB`, *including undefined/zero* counts.
4. Mean and RMS `A_ion` versus signed muon incidence angles and position. Show statistical and finite-MC limitations.
5. The measured cathode-readout follow-up is `A_Q(T)`, separately computed in Stage B3 and never asserted from the angular sector proxy.

**Physics and provenance checkpoint:** before running the muon-driven production study, explicitly validate the *already present* Geant4 momentum columns, TrackHeed segment boundaries, real primary Heed generation versus the old marked-Poisson surrogate, and correct coordinate transformations for the 3+3 mm chamber. The 2+4 mm Stage A/B legacy numbers remain *methodological only*.

## 7. Implementation checkpoint: true-muon B1 smoke (28 Oct 2026)

New non-overwriting code in `feature/mwpc-stage-c`:

- `pi_prepare_muon_avalanche_b1.py` reads the frozen, SHA-validated `Jpsi/chunk0001/geant_muon_crossings.csv` and selects 2 or 8 **genuine** gas crossing rows across physical stations and signs of `theta_u`. It uses existing `pu_in_GeV,pv_in_GeV,pw_in_GeV` and rejects edge/grazing/misaligned rows; selection is **stratified for validation and is NOT the true angle distribution**.
- `phaseB_native_muon_3p3_avalanche.cxx` generates true `TrackHeed` clusters once per native crossing using its incoming momentum direction. It reuses that same electron list across all four p3/p4 noFW/FW75 variants, then drives independent `AvalancheMicroscopic` seeds with `SetUserHandleIonisation(OnIonisation)`. It collects ions by full half-plane and ±45° cathode-facing sector *and* amplified electron anode endpoints per full muon crossing. An observation is one physical gas crossing, never one seed electron.
- `pi_run_native_muon_b1.py` audits all 4×N crossings, shared Heed primaries, angular count conservation, IDs, per-crossing true momentum, undefined denominators, zero/capped/failed avalanches, source hashes and immutable output manifests.
- `run_pi_native_muon_b1.sh` builds just the new CMake target, selects original gas crossings and runs small non-overwriting smoke tests.
- `test_pi_native_muon_b1.py` contains synthetic no-physics tests of momentum/chord direction and crossing-level count conservation.

**Current status:** code committed but Garfield compilation, execution on the actual WSL machine and physics validity **remain unverified**. Do not cite B1 numerical results before the CMake build and real smoke pass; a separate, independent-RNG reservoir holds at most 100 unbiased ion-birth coordinates per seed, with explicit inverse-inclusion weights for spatial visualizations; all full-hit A values use the complete ion-birth census. The physics fields still use the earlier bounded 9-wire local analytic proxy and exclude physical u-edge crossings. Extend to full station-angle acceptance only after validating the small chain.


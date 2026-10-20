# PI: centred-anode angular and time-response specification

Scope: NA60+/DiCE MWPC; branch feature/mwpc-stage-c, 2026-10-19.
This is an IMPLEMENTATION SPEC, not a claim that all analyses are done.

## Baseline and physically distinct quantities

The already finished native Geant4 transport has 200k generated events for each of Jpsi, phi and omega (67 chunks each). All native MWPC hit positions, entry/exit momenta, deposited energies and reference planes are frozen. Never rerun or modify the frozen transport for this task.

The detector/Garfield local axes are:
- u: across wires, nominally global Y
- v: chamber normal / gas thickness, nominally global Z
- w: along wires, nominally global X
- Garfield (x,y,z)=(u,v,w); actual transformed native coordinates must use the saved TGeo frame from export_c_geant_muon_crossings.C.
- Signed entrance angles are theta_u=atan2(pu,pv), theta_w=atan2(pw,pv). Analyse backward-going tracks separately, avoid silently taking absolute angle.

The chamber has 6mm of gas, with cathodes at v=-3mm and +3mm and anode wire plane at v=0. Wire pitch is 4mm. MS0 has 77 finite physical anodes; MS1..MS5 each have 128. Provisional field parameters from current Stage-C reference are Ar/CO2 70:30, 30um anode diameter, +1800V. Confirm the voltage and actual chamber geometry with PI before claiming calibrated predictions.

Existing 2+4mm Garfield B=0 response tables cannot be re-used by merely removing the historical +1mm local-v translation. A FRESH Garfield electric field and drift-table computation is necessary, with a NEW geometry identity. Old 2+4 results must retain a visible LEGACY label and remain unmodified.

## Output A: centered drift response map

Source: test/mwpc-stage-c/c_transport_box_map.cxx through run_c_transport_box_map.sh, setting C_TRANSPORT_GAP_MINUS_MM=3.0 and C_TRANSPORT_GAP_PLUS_MM=3.0. New output dir c_transport_box_map_gm3p0_gp3p0.

Parameters: u=-2..+2mm, v=-3..+3mm, du=dv=0.2mm, B=0 reference, 20x30=600 bins, >=50 single-electron trials/bin for pilot. Record collection status, relative physical wire index and drift_time_ns, including uncollected and timing-unknown trials. The compatible response map cell centres are:
 u_center=-1.9+0.2*iu, iu=0..19
 v_center=-2.9+0.2*iv, iv=0..29.
NO +1 mm offset in the new centered sampler. Preserve a separately named legacy 2+4 adapter for historical reproduction.

New strict metadata file must record geometry, B, gas, HV, wire diameter, pitch, transport seed and all input-code checksums. Fail closed when a file without a validated 3+3 geometry ID is used with centered mode. Validate full 600-box coverage, number of drift trials per box and local u/v map placement; include v and u mirror-symmetry *statistical* control comparisons at B=0, not exact identity.

## Output B: native muon hits

File: muon_gas_hits.csv, one row per TRUE Geant4 gas-hit SEGMENT:
source,chunk,event,hit_index,track_id,pdg,mother_track_id,
station,det_id,layer,row,col,chamber_center_x_cm,chamber_center_y_cm,
x_in_cm,y_in_cm,z_in_cm,x_out_cm,y_out_cm,z_out_cm,
u_in_mm,v_in_mm,w_in_mm,u_out_mm,v_out_mm,w_out_mm,
pu_in_GeV,pv_in_GeV,pw_in_GeV,pu_out_GeV,pv_out_GeV,pw_out_GeV,
theta_u_in_deg,theta_w_in_deg,theta_u_out_deg,theta_w_out_deg,
gas_path_mm,edep_eV,geant_hit_time_ns,valid_segment.

Unique segment ID: (source,chunk,event,hit_index), not merely event ID.
Muon track ID for same-chunk MCKine remapping is (source,chunk,event,track_id).
Different overlapping chambers and repeated segments of the same track in one chamber MUST be identified separately. Gas hits are not incoming-station-muon denominators; calculate incident counts from MuonReferencePlanes.root with verified track-ID reconciliation.

## Output C: pilot electron endpoints, preserving origin coordinates

File: primary_electrons.csv.zst (pilot only; stream aggregates in full production):
source,chunk,event,hit_index,track_id,station,det_id,cluster_id,electron_id,
u_birth_mm,v_birth_mm,w_birth_mm,t_birth_rel_muon_ns,
theta_u_in_deg,theta_w_in_deg,
collected,collection_status,anode_index,t_drift_ns,
t_arrival_rel_muon_ns,edge_loss,garfield_sample_id.

Sample Heed primary clusters along ACTUAL TGeo-local three-dimensional native gas segments. Do not convert Geant4 Edep to a fixed number of electrons: Geant4 energy loss and independently sampled Heed ionisation are two different models. Drift time is relative to electron creation; Geant4 time_ns is particle flight time; electronics discriminator crossing time is a third distinct physical time. If ionisation-creation offset is unavailable, leave t_birth_rel_muon_ns and t_arrival_rel_muon_ns unknown, not zero. Include attached/uncollected electrons and unknown drift time (not prompt).

## Output D: position evolving with time (requires NEW trajectory data)

File: electron_drift_snapshots.csv.zst
source,chunk,event,hit_index,electron_id,station,det_id,
theta_u_bin,theta_w_bin,t_rel_birth_ns,u_mm,v_mm,w_mm,
alive,attached,collected.

Save genuine Garfield intermediate drift positions at time snapshots (pilot suggestion: 0,5,10,20,50,100,200,400 ns). The current 600-box map contains only endpoint/time data and CANNOT generate genuine P(u(t)), P(v(t)), P(w(t)). Never use simple linear interpolation and label it Garfield drift. Distributions must identify survival/collection as a function of t. For the current translationally invariant 2-D ideal wire field, a physical w dependence from the electrostatic field itself is absent; sample correlation in w can arise from track geometry/acceptance or B-field. Include B=0 reference and only run B-on if actual field is known.

## Output E: anode multiplicity time response

File: anode_multiplicity_by_hit.csv:
source,chunk,event,hit_index,track_id,station,det_id,
theta_u_bin,theta_w_bin,gas_path_mm,edep_eV,
n_anodes_t50,n_anodes_t100,n_anodes_t200,n_anodes_t400,
n_anodes_t800,n_anodes_all,n_collected_electrons,
n_unknown_drift_times,n_edge_lost_electrons.

Each n_anodes_tT counts DISTINCT finite physical anodes with >=1 collected primary electron and t_drift<=T, not electrons and not hardware ASIC channels. N=0 gas hits remain in denominators. For every angular, chamber, station, source and time bin, emit N0,N1,N2,N3,N4,N5plus, denominator and P>=2 with uncertainty; assert sum(N0..N5plus)==gas_hit_segments for each selection. Check monotonic mean N(T) as T grows, thresholds <= all-collected, anode IDs within 0..76 or 0..127, no more distinct anodes than collected primary electrons.

## PI plot set and exact definitions

1. Position snapshots: 1-D normalised P(u;t), P(v;t), P(w;t) at the same set of times, separate by entrance theta_u and theta_w bins. Show survival populations/collected fraction at each time and non-selected versus selected sampling. Distinguish moving electrons from initial birth distributions.
2. Electron arrival: 2-D density (u_birth,t_drift), (v_birth,t_drift), (w_birth,t_drift) and corresponding 1-D conditional P(t_drift | birth-coordinate bin,muon-angle bin), normalised over collected electrons with known times.
3. Angle acceptance vs module position: XY map for each station and source, with individual chamber PDFs of signed theta_u_in and theta_w_in, and optional outgoing/delta angles. Map P(theta_u>0), P(theta_w>0), P(abs(theta_u)>=10deg), P(abs(theta_w)>=10deg), median and 16/84 quantiles. Existing full QA angle hist bins can support probability fractions at exact bin boundaries. Exact moment and quantile calculations require streaming native hits; histogram midpoint estimates MUST be visibly labelled approximations.
4. Angle comparison: centre versus upper/lower edge, left/right edge and corner chamber; label each by DET, layer, row/col, centre global X/Y and sample size.
5. Anode timing-gate study: P(N_active_anodes=0,1,2,3,4,5plus | chamber,theta bins,T) with T=50,100,200,400,800 ns and all, plus P(N>=2) and a sign-sensitive angle split. Normalise per REAL native gas-hit segment; provide a separate analysis if segment coalescing by track×det is needed.
6. Each anode occupancy: count activations by source×station×det×physical_anode; separate per gas segment, per generated parent event and per unique incident-station muon once virtual-plane mapping is verified.
7. Native material truth Edep and path-length: retain existing full-campaign spectra and physics flags as diagnostic overlays; no hidden forced link from Edep to sampled Heed cluster counts.
8. Plot error bars via event-level bootstrap for correlated hits/electrons and meaningful low-n markings; no spurious per-electron independence.

Angle bins: retain signed angles, starting with [-180,-30,-10,-5,-2,0,2,5,10,30,180] degrees; also offer occupancy-balanced exploratory quantile bins if rare high-angle strata too sparse. Always display exposure and avoid comparing different selections as if equal acceptance.

## Implementation steps (separate scripts, do not modify frozen campaign)

1. run_pi_centered_transport_box_map.sh: launch new 3+3 Garfield response into new output directory. Validate and version source; separate old map.
2. pi_validate_centered_box_map.py and test_pi_centered_geometry.py: strict map identity, bin grid, symmetries, no +1 mapping and deterministic seeds.
3. pi_centered_electron_sampler.py: native Heed ionisation with actual (u,v,w), zero offset, new map, raw pilot and aggregated full mode; never mix geometry epochs.
4. pi_native_angle_position.py: immediately build coarse module angular maps from full QA summary; in a second pass compute exact means/quantiles from native hits; do not rerun transport.
5. pi_native_centered_timing_pilot.py: full first chunk per source, gate multiplicity, electron drift-time distributions conditioned on angle and birth location; ensure N0, missing-time and edge losses correct.
6. pi_centered_drift_snapshot_pilot.cxx: reuse Garfield analytic field and step callback to store full u(t),v(t),w(t) for sampled Heed electron origins across angle and chamber-position strata.
7. pi_native_centered_full_fold.py: resumable checkpoints for all 201 chunks, one ROOT export at a time; per-module/anode/angle/time counters. Freeze manifest and verify the source model hashes before reuse.
8. build_pi_centered_pi_dashboard.py: one offline HTML with XY module-angle map, PDFs, time-sliced coordinate distributions, angular/time-gated anode multiplicity, legends explaining actual Geant4 vs Heed/Garfield model.
9. tests: true native hit conservation against 201-chunk QA, symmetric 3+3 field smoke and v+-mirror, near-wire and far-wire edge cases, all N0..N5plus totals and timing monotonicity.

Acceptance gates:
A. 3+3 centered Garfield map validated and new map ID pinned (requires new microscopic Garfield computation).
B. Complete native angular map (894 module/source records), source/station hit counts exactly match saved Geant4 QA.
C. Pilot timing counts accurate with all zero-hit and missing-timing cases preserved.
D. Real Garfield trajectory time slices available before claiming evolving-coordinate PDFs.
E. Full 201-chunk primary-electron response complete with correct incoming-muon denominators and no cathode strips.
F. Old 2+4 results remain reproducible and visibly identified as legacy.

PI to confirm before calibrated detector-firing predictions: final HV, gas mixture, 3+3 anode-centering across all modules, actual magnetic field at each station, and electronics threshold/shaping. Anode "collected >=1 electron" is never automatically a fired discriminator output.


## PI additional requirement: anode-wire pitch scan 3–4 mm

Requested 2026-10-19: study drift-time and angular response as the **distance between neighbouring anode wires** changes from 3 to 4 mm. This is wire PITCH in the transverse u direction, not anode–cathode spacing. Cathodes stay at local v=-3mm/+3mm and the anode plane at v=0; Ar/CO2, wire diameter, HV, B and physical gas aperture are held fixed for a controlled comparison.

Suggested FIRST scan grid p={3.0,3.5,4.0} mm (expand to 0.25mm increments after initial validation). All pitch points require fresh Garfield ComponentAnalyticField drift/endpoint tables. The current centered manifest, sampler and response lookup are HARD-CODED to p=4mm, u in [-2,+2]mm and 20 u bins; never reuse the 4mm response by relabelling p or simply rescaling endpoints. Set u in [-p/2,+p/2] and use an exact tiling, e.g. du=0.1mm for p=3,3.5,4, and retain v in [-3,+3]mm with dv=0.2mm. Refactor both generator and loading/sampling/validation to dynamic pitch and dynamic number of u bins. Version each map by physical pitch, geometry, B, HV and sampled-response provenance; refuse cross-pitch mixing.

IMPORTANT controlled comparison: the PHYSICAL gas active width in u stays unchanged. Changing p changes the number of *finite* anodes and end positions. Do NOT retain MS0=77 and MS1..5=128 at p=3mm (that artificially shrinks the instrumented aperture). Derive integer physical wire count and centered grid phase for each pitch and original module active width; explicitly document how partial edge cells are treated. Show both central unbounded-cell response and full finite-module occupancy to separate field/pitch physics from edge coverage.

Two DIFFERENT possible meanings of "angles" must be kept separate:
1. **Muon entrance angles** theta_u and theta_w are already stored in native Geant4 and, for unchanged gas geometry/selection, their incoming distributions DO NOT change with anode pitch. Plot time response and P(N>=2) *conditioned on the same native muon entrance-angle bins* for each pitch; this measures the geometry-induced response to oblique tracks.
2. **Electron arrival/incidence angle at anode wire**, e.g. saved Garfield phi_hit_deg / incidence_normal_deg, DOES depend on the electric field, which depends on p. Compare signed/absolute electron incidence angle and t_drift across pitch, and joint p(t_drift,theta_e_hit | pitch, initial u/v), not mix it with muon entrance angles.

Deliverables for p=3.0,3.5,4.0mm:
- Overlaid electron drift-time PDFs and CDFs, e.g. median, 90th/95th and fraction arriving <=50/100/200ns, at controlled same starting positions (and global native-angle-weighted Heed ensemble).
- Garfield electron incidence-angle (signed phi_hit_deg and/or well-documented incidence_normal_deg) distributions and conditional t_drift vs electron incidence angle, plus theta_e_hit vs pitch.
- Native Geant4/Heed full gas-hit sample reused with IDENTICAL deterministic cluster starting locations and correlated random trial quantiles across pitch where possible to lower comparison variance; do not regenerate a different muon sample.
- P(N=0,1,2,3,4,5+) and P(N>=2) per native muon gas-hit segment, vs muon theta_u/theta_w and drift-time cutoff T; propagate finite-grid edge losses and actual wire counts for each pitch.
- Explicit separation of field response vs geometric sharing effects; report per-pitch map-statistical uncertainty and gas-hit event correlation. Re-check B=0 first and only later physical B if known.

Implementation order: (i) keep current tested centered p=4mm baseline; (ii) create pitch-parametrised standalone Garfield runner and map manifest validator; (iii) dynamic-grid native Heed sampler using same tracks; (iv) one PI comparative dashboard with pitch selector/overlays; (v) pilot first, then full 201-chunk re-fold if comparison proves stable. No Geant4 transport rerun and no cathode-strip processing.


## PI follow-up, 2026-10-20: statistics, quadrant generation, field wires, benchmark

### Statistical reporting (supersedes earlier median-based proposals)

For drift/arrival times, electron arrival-angle distributions, and timing-gated
responses, report **arithmetic mean, centered RMS and statistical uncertainties**.
Define RMS = sqrt(sum((x_i - mean)^2)/N); clearly label this as a width, not
sqrt(mean(x_i^2)) and not the uncertainty on the mean. Give population N,
collection/unknown-time fractions, and error estimates. For full physical
response, resample source Geant4 **parent events** (clusters/electrons remain
together) so uncertainty is not biased by falsely independent electrons.
For single-electron LUT calibration, resample independent Garfield trials
inside each starting-position cell instead. Sparse cells must be marked.
Report neither medians nor percentile curves as primary PI deliverables.

### Exact versus approximate quadrant symmetries

The installed NA6PMuonSpec module placement is not invariant under general
(X,Y)->(abs(X),abs(Y)): ABCD layers have different z offsets (-6,-2,+2,+6 cm),
set by row/column parity. At least one reflection changes the parity-layer
for each of MS0..MS5. Dipole field, material, muon charge, vertex source,
acceptance and correlated dimuon daughters introduce further constraints.
DO NOT restrict Geant4 production to one quadrant or quadruple the event
count/statistical independence without independent validation of the whole
transport problem. Restricting BOTH muons to Q1 also changes the dimuon
correlation/acceptance. Safer alternative: preserve the full GEANT4 sample,
use the local periodic MWPC unit-cell symmetry for B=0 LUT reduction only
after comparing direct Garfield trajectories, and/or use symmetry for
descriptive acceptance diagnostics without pretending new independent events.

Read-only baseline before enabling any generation change:

    bash test/mwpc-stage-c/run_pi_stagec_preflight.sh

The script checks exact reflected layer/z geometry, compares existing full
201-chunk gas-hit segment count per chamber and reflected signed theta_u/theta_w
marginals by quadrant, and leaves a JSON + CSV. Counts and total-variation
distances are descriptive; full acceptance and event-level bootstrap require
per-event virtual-plane muon counts (not merely correlated gas-hit rows).
Never infer quadrant-transport equivalence solely from one hit-count ratio.

### Six geometries in the physics pitch scan

For each anode-to-anode pitch p={3,3.5,4} mm, compare without field wires
and with a grounded field wire at each anode midpoint. V_field=0 V
is **potential**, not zero electric field. Freeze cathode spacing 3+3 mm,
gas, nominal anode voltage, diameter, B and gas aperture for controlled
comparisons. Parameterize field-wire diameter and material, and record
collected-on-field-wire, cathode loss, anode collection and unknown-time
statuses separately. Rebuild the electric field and Garfield lookup for
EVERY geometry; no post-hoc relabeling of the p=4 mm LUT. Keep the physical
module aperture fixed; change finite anode count and edge treatment as p
changes. With alternating anode-field wire pattern, validate spatial cell
periodicity afresh. First study B=0; magnetic field requires independent
physical station conditions.

Principal outputs per p/field-wire configuration and signed incident muon-angle
bin: mean drift time +/- error, RMS drift time +/- error, collection fraction,
mean/RMS electron arrival-angle +/- errors (with angle wrapping/convention
explicit), angular/time correlations, and time-gated active-anode multiplicity.
For angle-of-impact, distinguish electron hit position phi_hit_deg from
final-velocity incidence_normal_deg: the latter is based on last microscopic
step and is sensitive to collision sampling. Retain the same Geant4/Heed
origins for all six fields, and keep muon-angle sign.

### Runtime and lookup accuracy benchmark

The frozen old Geant4 200k/source completed.json checkpoints did not record
independent wall-time durations. Do NOT estimate these from file timestamps.
New future Geant4 chunks now stamp monotonic wall seconds for the simulation,
ROOT validator, and entire chunk. The preflight reads old checkpoints
honestly as unmeasured. Timing reports give per-event mean/RMS and errors
bootstrap-resampled across independent completed chunks.

For the existing centered p=4mm, B=0 600-cell response map, preflight
also reports independent-trial mean and RMS of collected electron drift
time with bootstrap errors for EACH grid cell, how many cells have sparse
known times, and an in-memory NumPy lookup microbenchmark. Geometry-only
nearest-anode assignment is a SPEED comparator, **not a physical Garfield
replacement**. The full comparison of direct Garfield vs LUT accuracy,
field-wire scan, realistic memory/CPU profiling, and ultimately the new
2,000,000 generated parent events/source must follow this preflight.
Do not extend the original 200k/source suite: use a new frozen suite and
a fresh benchmark/smoke before long production.

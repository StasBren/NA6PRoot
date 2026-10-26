# NA60+ PI — ROOT-style parent-event A dashboard (p=3 mm, LUT500)

**Branch:** feature/mwpc-stage-c.
**New independent output:** test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/
**Do NOT mix with the existing complete J/psi 2-million-parent
six-geometry LUT20 dashboard.**

## What was inspected and reused

The latest existing dashboard is the offline Full J/psi 2M six-station
browser at:
  test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1/pages/index.html

Its interactive XY physical DET selector, native timing, phi-hit polar
atlas, 3/3.5/4-mm geometry switch and six stations are sourced from
**20 microscopic Garfield drifts per spatial cell** and are not event
asymmetries. Its merged counters do not retain source parent IDs,
so true per-parent A cannot be reconstructed by summing the existing
histograms. We therefore re-fold ORIGINAL sha-validated native gas
crossing rows with the already-frozen **p=3mm 500-trial** microscopic
Garfield response maps (noFW / FW 75um only), and the already-specified
stable-seeded Poisson-Heed primary birth sampler. No new transport or
Geant4 simulation. Every parent event has its original (chunk,event)
identifier, per station and physical DET.

For a parent with one or more classified anode electrons:
   A = (N_up - N_down)/(N_up + N_down).
Classify *anode surface impact POSITION* phi_hit sign, strict
0<phi<180 vs -180<phi<0; exact 0,+/-180 is unclassified.
Do NOT use velocity incidence azimuth or electron charge weight.
If N_up+N_down=0, A is UNDEFINED and only counted in the exposure.
Each parent's two J/psi muons may traverse multiple chambers/stations.

**Station:** sum electron up/down counts from every crossing in that
station for THIS parent first, then calculate one A.

**Chamber:** sum all the parent's crossing electrons in a single
physical DET first, then calculate one A per parent per DET. Do NOT
compute a station distribution by adding individual chamber A
histograms, or by averaging crossing-level A.

Histograms are P(A) over valid parent events only (41 native bins);
same initial Heed positions and same LUT lookup quantiles are used
for the two field-wire options; n_valid, n_undefined, mean(A),
population RMS(A), paired mean FW-minus-noFW A are all shown.
Raw event-level (up,down) counters are also saved in gzip CSV.

## 1. Update and perform ZERO-physics regression

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python

    "$PY" -m py_compile \
      test/mwpc-stage-c/pi_lut500_event_asymmetry_root_dashboard.py \
      test/mwpc-stage-c/test_pi_lut500_event_asymmetry_root_dashboard.py

    "$PY" test/mwpc-stage-c/test_pi_lut500_event_asymmetry_root_dashboard.py

## 2. Small, explicitly partial all-station smoke pilot

    BASE=test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1

    "$PY" test/mwpc-stage-c/pi_lut500_event_asymmetry_root_dashboard.py \
      --chunks 1 --stations 0 1 2 3 4 5 \
      --max-crossings 200 \
      --out "$BASE/pilot200"

The pilot samples ONLY the first 200 accepted muon gas crossings
PER STATION from one 50k-parent chunk. Its outputs MUST be labelled
PARTIAL, not scientific full exposure. It uses no Garfield process:
only existing frozen LUT500 lookups.

Open:
   test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/pilot200/dashboard.html

Check station selectors MS0..MS5, physical DET IDs, canvas and 41-bin
ROOT-like P(A) step histograms, entries/mean/RMS, paired geometry stats,
undefined A and the interactive physical XY map.

## 3. Full initial matched LUT500 J/psi response, 2 x 50k chunks

Only after the pilot checks succeed:

    "$PY" test/mwpc-stage-c/pi_lut500_event_asymmetry_root_dashboard.py \
      --chunks 2 --stations 0 1 2 3 4 5 \
      --max-crossings 0 \
      --out "$BASE"

The complete dashboard above uses up to 100,000 generated J/psi
parents in two input chunks. This is NOT the separate existing
2,000,000-parent LUT20 dashboard, nor a 2M LUT500 claim.

Outputs:
 - dashboard.html : ROOT-inspired light histogram canvas with stepped
   overlaid P(A) for p3/noFW and p3/FW75, interactive DET XY selector,
   clickable per-DET table, station selection, visual RMS(A) maps.
 - summary.json : station + physical chamber complete eventwise P(A)
   counts and moments, valid/undefined counts, same-parent geometry
   delta moments, native geometry layout.
 - event_counts.csv.gz : independently auditable true raw per-parent
   up/down counts at each station and physical chamber. Columns:
   source,chunk,event,station,det_id,gas_crossings,
   up_noFW,down_noFW,up_FW75,down_FW75.
   det_id=-1 = entire station, >=0 = physical DET.
 - manifest.json : SHA256 of generated files and referenced original
   native crossing SHA256s.

The builder refuses to overwrite existing outputs.
Every station's chamber electron totals must reproduce the station's
original summed electron counts for BOTH geometries. The nonlinear
eventwise A histogram is deliberately NOT summed across chambers.

## 4. Optional: apply ROOT-like theme to the latest 2M dashboard

The browser CSS (both XY chamber/time/angular page, per-anode detail
and polar atlas) is restyled to ROOT-like white histogram canvases
with clearer dark axes and two-tone series. This is VISUAL ONLY.
The full 2M backend still correctly remains six-geometry LUT20.
Re-rendering existing SHA-verified merged station JSON now rebuilds
all four station HTML variants if the new UI theme marker is absent,
without rerunning 240 native Geant4/Heed/20-trial Garfield shards.

    bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh check
    bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh start

The runner runs the read-only dashboard builder from already
checked local files. Monitor its STATUS with the existing script.
The portal will show a verified link to the companion LUT500
event-asymmetry dashboard once the latter's manifest and HTML SHA
exist. DO NOT modify the original publication manifest by hand.

## Physics interpretation limits

- This is electron endpoint **count** asymmetry, not induced U/V
  charge, gain, threshold, drift time or reconstructed muon A.
- The 500 stored microscopic trials are independently simulated per
  spatial cell and then REUSED. Parent-level Heed random seeds are
  stable, but reuse of LUT trajectories limits systematic precision.
- Station/chamber parent samples overlap; error bars need a parent-ID
  aware resampling method. RMS(A) is an event dispersion, NOT standard
  error of mean(A).
- This sidecar does not change full 2M LUT20 histograms or imply any
  3.5/4mm LUT500 exists.

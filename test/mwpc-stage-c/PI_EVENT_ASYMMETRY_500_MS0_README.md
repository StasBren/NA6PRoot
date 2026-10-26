# PI question 2: eventwise MWPC up/down collection asymmetry

This is a new isolated analysis for real native Geant4 J/psi muon gas
crossings in MS0, pitch = 3 mm, B = 0, without and with 75 um field wire.
Both geometries use the completed 900-cell microscopic Garfield LUT with
500 independent drift trials per cell.

Each gas crossing reuses the EXACT stable-seeded original Poisson-Heed
primary-electron birth sampler. The Garfield response is selected using
the original shared per-electron lookup quantile. Both configurations
receive identical initial electron positions and random quantiles.
No new Garfield/Geant4 physics and no changes to frozen 720 old folds.

## Observable

Up and down are classified by the sign of the Garfield anode endpoint
'phi_hit_deg = atan2(v_end, u_end-u_anode)':

- upper: 0 < phi_hit < 180 degrees
- lower: -180 < phi_hit < 0 degrees
- exact phi = 0 or +/-180: boundary, excluded from the asymmetry

Also exclude electrons missing an actual finite physical anode. No
avalanche multiplication, induced charge, threshold or electronics.
For N_up+N_down>0 define A = (N_up-N_down)/(N_up+N_down).
When N_up+N_down=0, A is undefined (do not incorrectly set A=0).

Two disjoint analysis units:

1. PARENT event: sum all selected MS0 muon gas crossings belonging to
   the same (source, chunk, local Geant4 event index). Report P(A),
   mean and population RMS across parent events. A J/psi parent has
   two muons, therefore it has no unique muon incidence angle.
2. GAS crossing: one native TGeo chamber crossing and one signed
   muon incidence theta_u_in_deg or theta_w_in_deg. Report the mean
   and population RMS of crossing A by entrance-angle bin, also
   showing the number of valid and empty crossings. The plotted RMS
   omits angle bins with fewer than 30 valid crossings, but the JSON
   retains every bin.

RMS(A) = sqrt(<A^2> - <A>^2): physical dispersion ACROSS sampled
parent events/crossings, NOT the uncertainty on its mean. Garfield LUT
finite statistics and cluster correlations remain possible systematics.
The earlier analytic expected Heed occupancy cannot be used to compute
event-level fluctuations; here we use the seeded individual Heed births.

## Run the small independent pilot

Run inside the existing project tree using the Python environment already
configured for PNG/PPTX slides:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only

    test_runs/mwpc_stage_c/.venv_pi_slides/bin/python -m py_compile \
      test/mwpc-stage-c/pi_event_asymmetry_500_ms0.py \
      test/mwpc-stage-c/test_pi_event_asymmetry_500_ms0.py

    test_runs/mwpc_stage_c/.venv_pi_slides/bin/python \
      test/mwpc-stage-c/test_pi_event_asymmetry_500_ms0.py

    test_runs/mwpc_stage_c/.venv_pi_slides/bin/python \
      test/mwpc-stage-c/pi_event_asymmetry_500_ms0.py \
      --chunks 1 --max-hits 5000

Outputs live independently in:
test_runs/mwpc_stage_c/pi_event_asymmetry_500_ms0_v1/

- asymmetry_summary.json: immutable LUT and native CSV hashes, event/crossing
  mean and RMS, signed entrance angle response and electron conservation.
- parent_event_A_distribution.png: normalized P(A) for parent events.
- crossing_RMS_vs_muon_angles.png: RMS versus individual muon crossing
  theta_u and theta_w, two Garfield geometry configurations.

## Full first-two-chunk sample

Only after the small pilot passes, use both already existing 50k parent
Geant4 chunks, without rerunning the detector simulation:

    test_runs/mwpc_stage_c/.venv_pi_slides/bin/python \
      test/mwpc-stage-c/pi_event_asymmetry_500_ms0.py \
      --chunks 2 --max-hits 0

This is heavier than the mean-Heed LUT convergence plot because true
per-event Poisson-Heed and sampled Garfield response fluctuations matter.
A partial pilot with --max-hits is NOT a complete generated-parent
exposure and should not be cited as such.

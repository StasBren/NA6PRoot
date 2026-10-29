# B2: full muon response versus LUT500 on shared real Heed primaries

Goal: take 20 PHYSICALLY DISTINCT real Geant4 J/psi muon tracks from original frozen NA60+ Muon Spectrometer gas-crossing export, generate their TrackHeed clusters and primary electrons ONCE, and compare the four p3/p4, field-wire/no-field-wire geometries with (a) full microscopic Garfield++ and (b) empirical 20-avalanche/cell overnight LUT500.

**This is a new stage beyond the 84,000 single-electron avalanche nights.**

## Input and geometry

- Source: original frozen Geant4 J/psi Muon Spec crossing, chunk 1; no new Geant4 simulation, no invented incidence angles.
- Deterministic hash-ranked sampling across six MS stations and both signs of theta_u, with 16mm clearance from actual u edge and chord/3D-momentum consistency.
- 20 unique PHYSICAL muon tracks, keyed by (source,chunk,event,track_id); repeated station hits of one muon forbidden. Other archived 2/8 B1 input and physics files are never changed.
- TrackHeed is called ONCE for each muon. Store accepted clusters (location, time, original/retained electron counts) and every retained electron with exact coordinate/time, cluster ID, primary index.
- Same electron list is then used for each full Garfield field. Ar/CO2 70:30, 293.15K, 760Torr; 1800V anode, ±3mm cathodes, B=0; each primary thermalized at 0.1eV as in original B1. All ion births/angular sectors preserved; ion POSITION reservoir limited to 5 per primary to save file size.
- Same immutable cached Heed electrons are replayed 200× per real muon/geometry with a new independent random draw from the 20 original complete LUT500 avalanches in each nearby 0.1mm-u × 0.2mm-v spatial bin. Anode phase uses actual station-specific wire-grid parity. Old drift 'collected' is NOT used as a gate; never double-count original attachment.
- Compare full Garfield's SINGLE independent realization per muon and field against LUT conditional distribution of 200 independently sampled gains and ion angular sectors. Report 5/50/95 percentiles, fraction of real-muon full reference falling in conditional 90% LUT interval, mean ratios, zero-birth counts and ion-weighted A_hemi/A_TB. Do not claim the two stochastic realizations agree seed by seed.

**Important boundaries:** 20 muons stratified by station/sign of angle is NOT the true incidence-angle PDF; each LUT cell has only 20 sampled avalanches; conditioning on exact Heed ionization does not quantify Heed-source fluctuations; this comparison is not yet a validated continuous/interpolated LUT, ion 3D spatial density or birth-time replay. Ion births are NOT induced cathode strip signals and do not include Shockley-Ramo current or electronics.

## Exact WSL commands

From existing branch, initial compile and source/Heed physics preflights:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
    python3 test/mwpc-stage-c/pi_run_muon_lut_validation.py prepare20
    python3 test/mwpc-stage-c/pi_run_muon_lut_validation.py check20

Only if these PASS, create full microscopic Garfield reference (unattended runtime not guaranteed; full physics benchmark must complete successfully):

    nice -n 15 python3 -u test/mwpc-stage-c/pi_run_muon_lut_validation.py reference20 | tee ~/na6p/muon_reference20_20261029.log

Next calculate 200 conditional empirical LUT responses (no Garfield or new TrackHeed):

    python3 -u test/mwpc-stage-c/pi_run_muon_lut_validation.py compare20

Finally independently publish full real-muon reference AND LUT comparisons under two separate immutable and SHA-verified folders on Google Drive:

    python3 -u test/mwpc-stage-c/pi_run_muon_lut_validation.py publish20 | tee ~/na6p/muon_lut_publication_20261029.log

Existing project Drive archive: https://drive.google.com/drive/folders/1Xbh45V86BhVu-WLYSM5qtosqPEG363M7.

Local immutable outputs:
- test_runs/mwpc_stage_c/pi_muon_lut_validation_v1/input_20_native_geant4_v1
- test_runs/mwpc_stage_c/pi_muon_lut_validation_v1/reference_20_muon_full_garfield_v1
- test_runs/mwpc_stage_c/pi_muon_lut_validation_v1/compare_20_muon_lut500_v1

Inside the full reference: original native Geant4 selection and SHA input manifest, real Heed cluster and electron CSVs, one Garfield full avalanche per cached electron/geometry, real muon-level angular sum, ion birth samples, trace and SHA manifests.

Inside the LUT comparison: muon_comparison.csv (80 original muon/geometry rows), muon_lut_replicas.csv (16000 resampled muon/geometry replicas), source SHA and method+limitations, summary and SHA manifest.

All old source LUT500/Geant4 and completed B1 runs stay immutable, and incomplete new run directories are retained rather than overwritten. Do not claim physical reference success or LUT agreement until the actual user WSL calculation and Google Drive readback have passed.

Check state with:

    python3 test/mwpc-stage-c/pi_run_muon_lut_validation.py status

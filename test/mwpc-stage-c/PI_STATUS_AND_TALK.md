# PI MWPC Stage C — status and review notes (October 2026)

Scope: NA60+/DiCE 298-chamber MWPC; six stations, J/psi, phi and omega at 40 AGeV. Baseline MS2 at 530 cm.

## Critical interpretation

The large old Geant4 sample stores muon kinematics at station **reference planes**, not actual current-MWPC gas entry/exit hits, Edep or exit momentum. We extrapolate these saved states along straight rays to the chamber gas slabs. Primary Heed ionization and sampled Garfield drift are then applied to these nominal segments. This is a **consistency check**, not full Geant4-to-electronics efficiency.

## PI requested items

| Request | Status | Needed |
|---|---|---|
| Activity of every wire in every module relative to muon run | Complete per-wire counts and fractions conditional on chamber gas crossing; module and station images | Additional normalization by **unique incident muons at station** and optionally total generated muons; clarify which is intended |
| Entry and exit particle angles per module | Reference-plane theta_u and theta_w available; entry/exit positions geometrically extrapolated | New Geant4 gas hits recording input AND exit directions; old straight segment cannot describe scattering |
| Energy deposited in different modules | Heed primary cluster count from modeled gas path, no real event Edep | New actual Geant4 gas Edep per hit/module, with consistency checks |
| One/two/three/four wire probability per single module hit | Older per-module random-wire-phase multiplicity available; modern fixed-grid multiplicity per **station** | New fixed-grid per-module N0/N1/N2/N3/N4/N5plus summary to align with individual wire maps |
| Sanity of number of wires | Geometry vs Garfield overlay, finite-wire checks, phase scans, 1e/5e/100ns thresholds | Tail and edge loss analysis against real new Geant4 hits and electronics |
| Physical module muon distribution | Full 298 chamber reference-to-gas projection with exact independent count checks | Cross-check with new actual Geant4 gas-hit topology and losses |
| Are angle and time distributions useful? | theta_u and theta_w distribution on accepted stations, theta_u vs multiplicity, 20x30 Garfield timing map, 100ns cutoff | Determine actual front-end integration/shaping threshold before making efficiency claims |

## Distinct denominators: do not conflate

For station s, module m and fixed wire j:
- N_station_muons(s): unique incoming reference-plane transported muons.
- N_chamber(s,m): gas slab crossings. A muon may cross multiple overlapping modules.
- N_on(s,m,j): crossings where that anode collects >=1 primary electron.

The chamber profile plots show N_on/N_chamber. They are NOT the
rate per incoming station muon, N_on/N_station_muons, requested by one reading of PI wording. V2 station colors are 10000*N_on/(sum of module crossings in that station): still not per unique muon. Add explicit dual normalization.

## Physical wire numbering

Four-mm nominal pitch. MS0 has **77 anodes, j=0..76** (30.72 cm gas height);
MS1..5 have **128 anodes, j=0..127** (51.20 cm gas height).
Index increases with global Y across wires; wires extend along global X.
Coordinates are *nominal centred wire-grid assumptions*, not a verified
mechanical first-wire drawing.

Example: MS0, DET17, layer B, centre X=-31.78 cm, Y=13.86 cm.
For J/psi: 13,086 chamber intersections; wire j=0 is at global Y=-1.34 cm,
j=76 at +29.06 cm. Wire j=36 has about 1.58% probability of >=1 collected
primary electron per DET17 gas crossing.

One crossing may activate several distinct wires, so sum over j of
P(wire j active per module crossing) can exceed unity.

## Example 0/1/2/>=3 multiplicities on MS0 (modern fixed grid, >=1 primary e-)

- J/psi: N=155,998 crossings, N0=53, N1=124,422, N2=31,523, N>=3=0.
  P(N>=2)=20.207%.
- phi: N=152,359, N0=22, N1=128,829, N2=23,508, N>=3=0.
  P(N>=2)=15.429%.
- omega: N=128,342, N0=22, N1=108,901, N2=19,419, N>=3=0.
  P(N>=2)=15.131%.

These are conditional accepted gas-crossing response results, NOT absolute
production rates or measured electronic discriminator multiplicities.

## Reproducible PI detector diagram

From the already computed fixed wire table, without any Geant4 or Garfield run:

    cd ~/na6p/src/NA6PRoot-mwpc
    git switch feature/mwpc-stage-c
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_station_numbered_anodes.sh

Outputs: station geometry with all 32 MS0 module IDs and one zoom of
DET17; second very detailed file with all 77 anodes literally numbered
0..76 with three probabilities at each anode. Optional different detector
and channel: pass an integer det_id and then source label.

## Prioritized work

1. Clarify with PI *which incoming muon count* should normalize raw
   individual-wire firings. Report both per-chamber-crossing and per-station
   incident muon rates to avoid ambiguity.
2. Add per-module **fixed-wire** N0..N4+ multiplicity probabilities and
   zero-hit fractions, preserving single module exposure, phase and edge data.
3. Use a modern full-chamber Geant4 sample. Export actual gas entry/exit
   coordinates and momenta, particle identity, path, accumulated Edep and time.
   Existing exporter: test/mwpc-stage-c/export_c_geant_muon_crossings.C.
   It currently exports incoming momentum, not explicitly outgoing momentum:
   check NA6PMuonSpecHit and extend it if necessary.
4. Cross-check real hits against the reference-plane geometric projection,
   and derive Edep distributions per station, module and incidence angle.
5. Only after a true gain and frontend integration/threshold model, decide
   whether wire arrival time or long-wire incidence angle should affect
   a **detectable** hit definition.

## Compact PI explanation

We mapped accepted Geant4 reference-plane muons to the six station
checkerboard MWPC layout and checked all per-module crossing counts. We then
folded Heed primary ionization and a 600-box Garfield single-electron
transport map onto those gas intersections. In the current nominal fixed
4-mm anode grid, the MS0 probability to collect primary electrons on
at least two anodes per gas crossing is about 20.2% for J/psi, 15.4% for
phi and 15.1% for omega. This is consistent with broader accepted J/psi
angular kinematics. The actual XY track footprint is compared with
the geometric anode cells and the Garfield-collected wire profile.
But these results still use an extrapolated gas segment and collected
primary electrons: we have not yet demonstrated actual 298-module Geant4
gas Edep, exit-angle scattering, gas avalanche gain or readout threshold.
These are the next tests, not results to claim as completed.

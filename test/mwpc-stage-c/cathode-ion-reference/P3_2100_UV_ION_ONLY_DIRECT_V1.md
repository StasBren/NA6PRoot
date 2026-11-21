# P3 2100V: direct ion-only U0/V0 preparation (2026-11-21)

**STATUS: Code committed, NOT YET compiled or tested on user's WSL, no physical U/V waveform established.** Safe light tests only while the original 120-event microscopic recovery is still running. Do not rerun original Geant4, Heed or Garfield AvalancheMicroscopic; the campaign is isolated in `~/na6p/cathode-overnight-2100`.

## Implemented (separate development branch)

1. `TimeBinnedRamoCharge.hh`: ROOT/Garfield-free time-window accumulator. For a simulated ion between consecutive RKF samples `(r_j,t_j)` and `(r_{j+1},t_{j+1})`, it deposits `+e * sampling_weight * [phi_w(r_{j+1})-phi_w(r_j)]` (Garfield Sensor convention) into 10ns time bins, assuming linear `phi_w(t)` inside each RKF segment. An O(1)-per-segment difference array handles long segments; charges outside [0,100us] are recorded as unobserved tail. Per-channel arithmetic telescopes to the full electrode endpoint Ramo charge. This proves bin arithmetic, not FEM accuracy or trajectory time-resolution convergence.
2. `stereo_ion_response_from_ions_root.cc`: imports `ComponentElmer` plus and minus selected central-strip maps; constructs exact analytic nine-wire 2100V ArCO2 physical drift geometry and RKF positive-ion paths with the historical provisional mobility table. The physical drift is independent of weighting strip segmentation. Reads the existing trusted `Plane2100Archive.hh` K256 per-primary weighted ion births; default processes every saved ion and audits weighted census. Explicit `--ion-limit 1` pilot processes only a single ion and MUST NOT be called a complete muon. Reads RKF vertices via `GetNumberOfDriftLinePoints` / `GetDriftLinePoint`, uses `phi_w` differences, analytic 0/1 endpoint conductor values, and enforces per-strip charge-bin closure. Outputs `*_summary.csv`, `*_channels.csv`, `*_waveforms.csv` and either `*_pilot.txt` or `*_complete.txt`, never overwriting existing output.
3. Noncentral translated strips are DISABLED by default; opt-in `--half-strips N --allow-translated-strips` is explicitly UNVALIDATED on the finite-Neumann-domain map. The first result only U_plus_0 and V_minus_0; the full U/V cluster, independent adjacent-strip FEM validation, electron-induced current, detector electronics and fast response library remain future work.
4. `run_p3_2100_stereo_ion_response.sh` checks source `complete.txt`, original ROOT physics audit, uncapped microscopic physics evidence, SHA256 event integrity, FEM map SHA256 and candidate geometry. Does not change frozen original.
5. The existing faceted FEM baseline is **PROVISIONAL**: old quadratic mid-edge wire surface nodes deviate inward as much as 2.222μm compared to only 1.5μm near-wire evaluation clearance. Do not present computed U/V signals as physically validated until curved-P2 mesh, Garfield curved-element interpolation, mesh/patch sensitivity and adjacent-strip checks pass.

## Low-memory preparation while muon recovery runs

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git merge --ff-only FETCH_HEAD
g++ -std=c++17 -O2 -Wall -Wextra -pedantic \
  test/mwpc-stage-c/cathode-ion-reference/test_time_binned_ramo_charge.cc \
  -o /tmp/p3_time_binned_ramo_test
/tmp/p3_time_binned_ramo_test
bash -n test/mwpc-stage-c/cathode-ion-reference/run_p3_2100_stereo_ion_response.sh
bash test/mwpc-stage-c/cathode-ion-reference/run_p3_2100_stereo_ion_response.sh preflight
~~~

Optional **compile only** (does NOT load the 1.85M-node FEM mesh; no ion trajectories computed):

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build-fem-strip-probe \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL" \
  -DROOT_DIR="$HOME/na6p/install/root/cmake"
cmake --build test/mwpc-stage-c/build-fem-strip-probe \
  --target mwpc_stage_c_stereo_ion_response_2100 -j1
~~~

## Only after original microscopic recovery finishes and RAM is available

Existing unvalidated faceted map can be used for a **provisional algorithm test**; it must not establish physical strip accuracy. The map can later be replaced via an independent completed curved-FEM run directory `MWPC_STEREO_MAP_RUN` without rewriting any muon archive:

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
RUNNER=test/mwpc-stage-c/cathode-ion-reference/run_p3_2100_stereo_ion_response.sh
MWPC_STEREO_EVENT=0301 bash "$RUNNER" pilot
# If the one-ion pilot succeeds:
MWPC_STEREO_EVENT=0301 bash "$RUNNER" muon
~~~

If the sample w coordinate falls outside the finite patch, or the FEM returns nonfinite/out-of-range values, the run explicitly fails with an ion-row diagnostic. It must not substitute zero for an invalid map point. Garfield ComponentElmer may return zero for domain misses and a dedicated location diagnostic will be needed before production.

## Required next validations

- Compare charge and early bins under halved ion step limit; don't mistake a telescoping `Delta phi_w` identity for independent accuracy.
- Verify the curved-wire P2 map actual radius and Garfield interpolation before using real U/V fields.
- Generate explicitly displaced strip FEM maps k=±1 and compare to translated central maps: finite Neumann boundaries break perfect symmetry.
- Choose cluster width from actual signal tails; reuse the exact same source reservoir for direct-vs-fast comparisons.
- Add prompt microscopic electron signal separately; never call current release full e+ion.
- No electronics/noise/VMM.

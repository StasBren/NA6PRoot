# Joint development decision: physically validated U/V weighting and low-error deterministic response operator
*2026-11-21. Integrates independent conversation "Gemini Claude Model Prompt".*

## Frozen input status (user WSL output)
Overnight recovery completed at 2026-11-21T16:12:27+02:00:
`STATUS VERIFIED=997/1000 FAILED=1 UNSEALED=2 ABSENT=0`.
All 120 previously loader-failed late slots were recovered; `newly_failed=0`.
Read-only SHA and event-physics verification remain mandatory. Do not modify frozen
worktree, rescued failure evidence or 997 archives. Exclude the three unverified
events from physical response statistics unless separately recovered and audited.

## Conclusions from the independent whole-plane operator study
The other worktree/branch is
`~/na6p/response-rank-pilot-20261121`,
`feature/mwpc-avalanche-response-rank-pilot-20261121`.
These are TWO WHOLE, UNSEGMENTED cathode planes, not U/V strip results.

- Previous 256-start IDW4 waveforms on the development muon0301 gave joint
  native-10ns waveform L2 about 6.2304% against direct ion RKF.
- The original bank had a severe `(ln radius, polar angle)` support mismatch:
  71.001% weighted source ions were outside its intrinsic combined marginal
  coverage. Global dipole/quadrupole moments are NOT an adequate physically
  validated substitute for a spatial response operator.
- A deterministic uncompressed 12x24 polar bank gave 4.21196% L2.
- Paired EXACT-SPATIAL control demonstrated that even ion starts at their exact
  true positions incur 3.99763% whole-plane waveform L2 under the old
  10ns fractional current-histogram shift (4.29363% in 0-100ns).
- Discretizing the birth-time phase into four 2.5ns phases reduced exact-spatial
  shift L2 to 0.924936% total, 0.993818% in 0-100ns.
- Latest results reported in that independent chat for the 12x24x4 polar
  phase-resolved bank: 1152 trajectory/phase responses, single development
  muon0301 joint whole-plane L2 about **1.3545%**; exact-spatial phase4
  about **0.9249%**, reported spatial component about **0.8516%**.
  These values have not been independently rerun in THIS U/V worktree and are
  not results on independent muons or real strip weighting maps. The
  0.8516 and 0.9249 errors cannot be combined in quadrature without checking
  covariance/cross term.

**Architectural conclusion:** retain the direct single-ion RKF trajectory
and cumulative induced charge `Q_s(t)=q[phi_w,s(r(t))-phi_w,s(r_birth)]` as
the reference physics object. Time-bin currents are derived from increments in
Q, not by shifting coarse current histogram samples. A deterministic
geometry-aware polar trajectory/operator bank is the current leading fast
candidate. SVD is a STORAGE compression layer to test AFTER uncompressed
spatial/temporal interpolation works, NOT a replacement for the correct
spatial coordinates or physical FEM strip weighting. No fixed moment-regression
or unverified nine-moment exactness claims.

## Crucial distinction: DIRECT strip pilot vs FAST library

The existing `mwpc_stage_c_stereo_ion_response_2100` executable in this branch
reads one audited K256 real-ion birth archive, solves the actual physical
Garfield RKF drift for each saved ion, then takes `ComponentElmer` U0/V0
strip-weighting potential differences on the path. It does NOT use IDW,
the 256-template bank, the polar bank, fractional template shifts, or SVD.
Therefore the measured 6.23%, 4.21%, 3.998% fast-operator errors DO NOT
invalidate the purpose of this **direct reference/pipeline smoke test**.

However, this direct pilot itself is NOT yet a numerical ground truth because:
(a) the original FEM 15um anode wire has mid-edge radial errors up to 2.222um
and a new curved P2 mesh has not completed;
(b) the Garfield RKF path-vertex time spacing and piecewise-linear
phi(t) integration still require a step-size/early-bin convergence check;
(c) the initially solved U0/V0 maps are finite-box central strips; a
multi-strip translated cluster requires explicit k=+/-1 map tests,
physical channel origin and finite boundary patch validation;
(d) only avalanche positive ions are represented; microscopic prompt
electrons and original Heed positive ions are not yet included.
The pilot's telescoping Q identity is an INTERNAL arithmetic check,
not a field-accuracy certificate.

## Execution decision after 997-muon completion

1. Keep two isolated development worktrees and existing SHA-sealed datasets;
   **do not merge the two branches wholesale or overwrite any old output**.
2. Now safe to perform the read-only one-real-ion **pilot** against the saved
   faceted baseline ONLY as a software/integrity smoke test, after
   `test_time_binned_ramo_charge` passes and compiling the latest code.
   It may find unsupported strip phase/w coordinate or Garfield map lookup.
   Never infer strip physical accuracy from its amplitude.
3. Run new, no-overwrite physical **curved 12um P2** FEM map in
   `p3_2100_fem_curved012_v2` and audit actual mid-edge radius, Gmsh
   high-order mesh Jacobians, Elmer 510 import, Garfield point location,
   U/V symmetry and potential changes against baseline. Observe laptop
   7.4GiB RAM; do not start the large FEM solve concurrently with an
   ongoing ion pilot or rank-bank generation.
4. At identical one-ion paths compare baseline versus curved map induced
   `Delta Q`, and independently halve RKF max step (e.g. .05 to .025mm)
   to quantify **trajectory time discretization** separately from
   **FEM electrostatic geometry**. An exact coaxial Ramo test is an
   independent validation target for the full chain.
5. With field gate met, process the full weighted ion reservoir for one
   muon, **still only U0/V0 ion-only**. Add independently computed
   neighboring k strips and verify finite box translation hypotheses;
   do not call two channels a complete physical U/V cluster.
6. For fast strips, compute deterministic RKF trajectories at a polar
   grid in transverse (u,v), with w coordinate and stereo strip phase
   handled in the **weighting-potential evaluation** along each
   trajectory. Precompute a time- and strip-index-resolved cumulative
   induced-charge operator `Q_{s,k}(t; ln(r),theta,phi_strip)`, not an
   interpolation/shift of coarse current bins. Preserve actual ion
   birth time by querying/interpolating cumulative Q at shifted bin
   *boundaries*, test 4x2.5ns phase version as an independent simpler
   control. Distinguish finite-bank edge-wire and local coordinate
   classes. Reuse same K256 ions for direct-vs-fast comparisons.
7. Freeze model parameters on the development sample, then test true
   independent muons 0408 and 0446, and later source/station-balanced
   cases. Report native 10ns waveform L2, 0-100ns charge, 100us charge,
   ion tails, per-strip signed charge/centroid differences, incomplete
   strip coverage and exact 0% out-of-grid before applying fast operator
   to **997 verified** microscopic archives.
8. Compression: only after uncompressed U/V operator proves accurate,
   evaluate SVD basis rank and reconstructive error versus speed/RAM.
   Do not use Gemini's unvalidated 'nine global moments' as the fast model.

## Priority status
The independent fast whole-plane improvement is impressive but has not
completed true U/V strip accuracy or independent-muon validation.
The original single-ion U0/V0 pilot is still scientifically useful as a
**diagnostic**, whereas the curved FEM geometry audit is the critical
physical-accuracy gate. No claim is made of validated strip
signals, physically calibrated ion mobility, electronics or full 1000
verified events.

# Requirements, provenance and current status

## A. Explicit group guidance (confirmed in earlier correspondence)

- In the 2026-09-30 MWPC simulation/reconstruction correspondence, Sasha Milov said the **detector-response Monte Carlo should remain within the group's implementation/maintenance responsibility** because mechanical/readout parameters are changing and require optimization.
- Provisional strip width **about 1.7–2 mm** and strip-to-wire stereo inclination **tan(alpha)=0.1** were stated. These are NOT a final geometry freeze.
- A later November 2026 group meeting concerned alternative downstream MS4/MS5 readouts and possible mechanical/electronics changes. Do not hardcode the current strip segmentation for all stations.

**Privacy/provenance:** those are paraphrases of group correspondence, not pasted private mail. Do not publish participant email addresses or private meeting details.

## B. User-relayed Sasha feedback (reported 2026-11-15; not independently verified as his exact words)

1. Investigate **gain vs anode voltage** and **gain fluctuations vs voltage**.
2. Characterize gain fluctuations with **RMS/mean**, not Gaussian-fit sigma/mean.
3. Show complete gain **distributions**.
4. Understand channel differences in ion azimuth (J/psi, phi, omega), especially phi; check whether input geometry/station mix explains them.
5. Analyze observables **per station**, with per-channel station occupancies and detector IDs.
6. Revisit slide-9 avalanche vs muon asymmetry at altered voltage; distinguish meaningful RMS(A) from RMS(A)/mean(A) (signed mean near zero).
7. **ROOT/C++ only for physics/statistical calculations, and ROOT with ATLAS-like plot styling**; do not use Python for physics analysis or plot generation. Label plots honestly as NA60+ simulation.
8. Document exact **selection of 30k muon crossings**, including station/det_id, incident angle, momentum and path distributions; don't conflate equal per-channel with equal per-station sampling.

These are research requirements; do not claim any voltage trend, channel effect or station bias until computed with verified code.

## C. Confirmed current simulation facts

### Original Garfield++ 30k ensemble

- 10,000 genuine Geant4 crossings from each of J/psi, phi, omega; exact original selection algorithm and station composition still require audit.
- TrackHeed -> AvalancheMicroscopic in P3 noFW at nominal 1800 V.
- 1,554,258 primary-electron avalanches and ~5.073e9 ion births; 14.69% zero-collected primary avalanches.
- Source library built on 10k J/psi; phi/omega reserved as physics-channel validation samples.
- Original 30k export includes seed/muon summaries and ion birth reservoirs, **but not true physical collecting-wire IDs**.
- Original fast simulation Fast4 uses 900 LUT500 initial-position cells x 500 drift trials/cell; samples collecting wire / arrival azimuth / delta-w / drift time jointly; J/psi donors supply gain, ion birth count and angular/longitudinal distributions.
- Fast4 matches global gain and global signed A fairly well, but local ion-cloud correlations need improvement.

### True physical-wire reference (completed)

- 1000 *Geant4 crossings* pseudo-randomly selected from original metadata with 334 J/psi, 333 phi, 333 omega; all six stations; **fresh Heed RNG streams** (not exact replay of old Heed primaries).
- 52,997 actual Heed primaries; 172,548,342 ion births.
- 246/1000 muons >1 actual collecting wire with nonzero collected charge.
- 192/1000 have a SECOND wire above 5% of event collected charge.
- In those 192, the two leading wires have opposite signed ion-birth hemisphere asymmetry in 143/192 (74.5%).
- \`A_mu = sum_k I_k A_k / sum_k I_k\` using ion births associated with nearest physical wire. Local |A_k| can be large even when |A_mu| is small. Never treat this asymmetry as an electronics charge measurement.
- Most recent preliminary Fast4 verification used same **new** Heed inputs; do not claim calibration on them. Its apparent good global multi-wire fraction does not validate cathode strips.

### Voltage scan (in progress; NOT finished)

- Separate GitHub branch \`feature/p3-hv-scan-20261116\` scans 1900, 1950, 2000, 2050 V using fixed 1000 Geant4 crossings / fixed Heed seeds per point.
- Local pilot is being compiled/run. Do not assume complete results.
- Avalanche-size caps MUST be audited; larger voltage may invalidate a capped simulation and is sensitive to missing space-charge effects.
- Metrics: per-primary mean gain including zeros, positive-only gain, P(g=0), RMS(g)/mean(g), muon-effective gain fluctuations, signed A and station-conditioned checks; ROOT/C++ plots with per-muon bootstrap.

### Additional requested HV / wire-diameter pilot (prepared, not yet executed)

- New independent GitHub branch `feature/p3-hv-wire50-check-20261116`.
- Same Geant4/Heed inputs; 20 muons at 2000, 2050, 2100 V for **30 um** anodes; paired first five muons at 2000 V with **50 um** anodes.
- Separate local target, runner, ROOT/C++ analysis:
  `test/mwpc-stage-c/hv-diameter-study/`.
- Pilot=5 per configuration, full extends the 30 um set to 20; 50 um remains 5.
- Code prepared and source-reviewed only: no Garfield++ compilation, run or physics results have been confirmed yet. Abort on size-capped avalanches.
- Do not confuse this branch/run with the original 1900–2050 V scan.

### Future 30k per-wire sample

A new 30k wire-resolved Garfield++ campaign may be in progress on the HPC farm (status not independently confirmed here). Before use, require immutable manifests, exact configuration, per-wire endpoint counts, birth census, and discover whether raw ion birth positions/times or only bounded reservoirs are saved.

## D. Keep statistical units separate

- **Per primary Heed electron:** microscopic avalanche collected gain (can be zero), exact ion births, local cloud.
- **Per actual collecting wire per muon:** Q_k=e * collected_e[k], separate nearest-wire ion birth populations; not necessarily one-to-one associations.
- **Per muon:** effective gain G_mu=Sum(g_i)/N_Heed, signed A_mu and meaningful multi-wire fraction with explicit >5% definition.
- **Cathode strip:** induced current/charge waveform after electron/ion drift and weighting field; not equivalent to the anode charge Q_k.
- For signed A with mean close to 0, report RMS(A), mean|A|, full signed distribution, and correlations; NOT RMS(A)/mean(A).

## E. Open physics/design decisions

- Which anode HV to use for *real* cathode-response development after the scan?
- What is the correct mixture-dependent positive-ion mobility / ion composition? Existing Phase-C code explicitly labels mobility as provisional.
- Which strip layout on upper/lower cathodes per station, with stereo angle sign, pitch, width, gaps and boundary effects?
- Which electrostatic and weighting-field models provide validated strip currents (not whole cathode currents)?
- How to represent billions of ions: full ion-birth reference pilot vs unbiased weighted macroions vs transfer kernels; prove numerical convergence.
- What electronics response and thresholds should be implemented? Do not invent values; make them configurable and label provisional.
- What cathode-signal observables and hit-reconstruction performance are the acceptance metrics?
- Fast5 remains future development: consider conditional *whole-muon response-library* sampling rather than reconstructing local A_k from global A_mu. No validated 10x speed claim; prior Fast5-T prototype was only ~3.4x faster than Fast4 in response-only benchmark and had local asymmetry problems.

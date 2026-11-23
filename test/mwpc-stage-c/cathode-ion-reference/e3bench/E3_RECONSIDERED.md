# E3 reconsidered against the Fast4 / Fast5-T chain (2026-11-23)

Labels: **[identity]**, **[assumption-consequence]**, **[measured]** (with machine and sample), **[synthetic]**
(tool test, not NA60+ data), **[hypothesis]**, **[open]**. Nothing here is validated against detector data.

## 1. Where the cost is

### Existing measurements (not re-measured here)
* Fast4 matched Python benchmark (slide 40): 1.60 ± 0.02 s per 20,000 muons, about 80 µs per muon. Fast5-T:
  0.47 ± 0.06 s, about 23.5 µs per muon. Both use precomputed Heed inputs. Neither includes a cathode response.
  **[measured, Python, by the project]**
* Frozen RKF + spectral replay: 55.5 s per muon. **[measured, project]**

### New C++ stage benchmark: `digitizer_cost_bench.cc` **[measured, table contents synthetic]**
Machine: 2-core Xeon 2.3 GHz cloud VM, g++ -O2, single thread, 20,000 crossings, median of 3 repetitions.
Table sizes and loop counts are realistic: 52.78 primaries and 1.39 wires with ions per crossing (V0 means), LUT500
900 cells × 500 trials, 539,422 J/psi donor seeds, a 51,550-record joint 2100 V donor in 72 cells, and a
2 planes × 3 harmonics × 601 Δs × 400 t (2.5 ns) shaped-kernel table (7.7 MB). Physics content is random.

| Stage | ns per crossing (Tp=100 run / Tp=200 run) | Earlier runs |
|---|---|---|
| Heed **surrogate** (Poisson clusters + size table) | 410 / 433 | 410–600 |
| A: Fast4 draws (LUT500 + J/psi donor) + per-wire sums | 859 / 1,414 | 1,005–2,466 |
| A: joint 2100 V start-cell donor + per-wire sums | 557 / 987 | 554–1,257 |
| B: direct per-wire correlated draw | 34 / 35 | 36–106 |
| Kernel, per wire, one time (T4-like) | 3,870 / 3,812 | 3,815–4,659 |
| Kernel, per wire, mean + variance term (T3-like) | 4,839 / 4,540 | 4,411–5,541 |
| Kernel, per primary (T2-like) | 157,295 / 156,890 | 148,469–154,525 |
| Digitize (peak, parabolic refine, threshold time) | 412 / 453 | 389–546 |

Reproduce: `g++ -O2 -std=c++17 -o bench digitizer_cost_bench.cc && ./bench 20000 3 100 && ./bench 20000 3 200`
(arguments: crossings, repetitions, Tp in ns).

Findings:
* Evaluating the kernel per primary costs about 157 µs per crossing, about 30 times all the other stages combined
  (about 35 times the per-wire kernel).
  Aggregating per collecting wire before the kernel is the decisive saving, not the choice of generator. **[measured]**
* With per-wire aggregation the total is about 5.5–7.5 µs per crossing. Of that, the kernel is 3.8–4.8 µs. Strategy A
  generation costs 0.56–1.4 µs and Strategy B 0.03–0.1 µs, so B would save at most about 0.5–1.4 µs per crossing, or
  roughly 10–20% of the digitizer. That is not a large enough gain to give up the correlations A keeps by construction.
  **[measured]**
* Scale: 720 station hits per Pb–Pb event at about 6–7 µs each is about 4–5 ms of digitizer time per event per core,
  and about 1.2–1.4 core-hours per million events, **excluding Heed and Geant4**. **[arithmetic on the measurement]**

### Is Heed a relevant production cost? **[open: not measurable here]**
The Heed surrogate costs about 0.4 µs. Real TrackHeed is not measured in any document I found. In
`phaseB_native_muon_3p3_avalanche.cxx`, `HeedPrimaries` builds a new `TrackHeed` for every crossing, so initialization
is paid per crossing. If that cost is in the 0.1–1 ms range, it is 15–150 times the rest of the digitizer and would
need either a reused `TrackHeed` per geometry or a cluster-level surrogate. This is the largest unknown in the budget.
It is separate from microscopic avalanche multiplication, which the production chain does not run.

Measurement (WSL, existing binary, Heed only, no avalanches):
```bash
B=<build>/mwpc_phase_b_native_muon_3p3
head -n 101  NATIVE.csv > /tmp/h100.csv; head -n 1001 NATIVE.csv > /tmp/h1000.csv
/usr/bin/time -f "%e s" $B --validate-heed /tmp/h100.csv  > /dev/null
/usr/bin/time -f "%e s" $B --validate-heed /tmp/h1000.csv > /dev/null
# per-crossing Heed cost = (t1000 - t100) / 900 ; intercept = gas/startup
```

## 2. Keep LUT + per-primary generation? Yes, with a 2100 V input
**[measured cost] + [hypothesis on physics]** Keep Strategy A: one stochastic record per Heed primary carrying the
collecting wire, onset time, ion yield and cloud moments, then aggregation per wire. The generator is cheap, and the
correlations between collecting wire, time, yield and asymmetry come from the same draw. A direct per-wire generator
(B) must model those correlations explicitly. V0 says that is hard: the independent-primary variance ratio is 3.34
unconditioned and 2.21 after coarse conditioning, so the per-wire sums are not simple functions of the geometry.

Recommended input: a **joint 2100 V start-cell donor** built from the K256 archives. Each archive primary is a joint
sample of transport and avalanche given its start position: collected electrons per wire, ion births per wire, onset
time, delay and cloud moments. Drawing one such record per primary replaces LUT500 + J/psi donor at 2100 V. It also
avoids applying attachment twice, since the archive record already includes it.
Requirements still to check on the real V0 output: counts per start cell (72 cells, about 716 primaries per cell on
average), within-cell dependence on track angle, and the residual correlation between primaries (`rho_w` from
`v0_plots.py`). Correlation that start conditioning cannot explain must come from the muon-level context (angle,
cluster structure). The donor must then be stratified by that context, or B-like corrections added on top.

## 3. LUT500 parts that need 2100 V recalibration
All five LUT500 outputs depend on the field: collection probability (attachment against drift time), drift time,
longitudinal displacement, collecting wire near cell boundaries, and arrival azimuth. Gain, ion yield, A and cloud
moments come from the 1800 V 30k campaign. None can be transferred without a check **[assumption-consequence]**.
The joint 2100 V donor removes the need for a separate 2100 V LUT if cell coverage is adequate. A targeted 2100 V
`DriftElectron` LUT (same 900 cells, 500 trials, drift only) is cheap next to microscopic avalanches. It is worth
running only if (a) some cells are under-populated, or (b) the arrival azimuth is needed separately from the ion-cloud
moments **[hypothesis]**. The archives do not store electron drift lines or arrival azimuth.

Wire collection, three separate claims:
1. 0.121% (52 / 42,879): dominant collecting wire vs dominant ion-birth nearest wire. Both describe the actual
   avalanche outcome. This says nothing about transport. A bug that set the argmax to wire 0 for all-zero census rows
   was fixed today; re-check this count.
2. About 4%, from the old 1800 V nearest-start-position model. **Not recovered from the repository.** The closest
   controlled definition is `lut500-vs-geometry-100k` (branch `lutgeo`): geometric midline vs LUT500 with the same
   electron starts and gains, comparing >5% second-anode and third-anode rates per muon. That README warns that
   earlier cluster-nearest vs Fast4 comparisons also changed the gain model. The 4% figure needs its
   `comparison_table.md` or the slide source.
3. New plot `07_collecting_vs_start_anode.png`: dominant collecting anode minus the anode nearest the PRIMARY START,
   by start cell. This is the 2100 V counterpart of claim 2, per primary, not per muon.

## 4. What can be omitted for 100 / 200 ns shaping
* Shaping times: the VMM3a user guide gives 25/50/100/200 ns peaking times and a 3rd-order shaper (1 real + 2
  complex-conjugate poles) without pole values. I found no document that fixes the NA60+ operating value.
  **[open: needs the front-end configuration]**
* Within-avalanche birth-time spread: 0.0947 ns RMS (V0, 977 muons). The CR-RC^3 estimate of the relative peak loss
  is n σ² / (2 Tp²) ≈ 1.3e-6 at 100 ns. **[assumption-consequence]** Omit it; T1 keeps everything else.
* Between primaries on one wire: 10.15 ns residual RMS. The same formula gives about 1.5% (100 ns) and 0.4%
  (200 ns) peak loss. The larger effect is likely **threshold-time bias**: collapsing to the mean time moves the leading
  edge, which the earliest primaries set. **[hypothesis, synthetic support]** On 3 synthetic muons
  (`timing-ablation/`, real spectral LUT, straight-line paths), T1 changed every observable by < 3e-4 relative and
  < 0.2 ns. Per-wire time collapse (TW) moved the 2% leading-strip threshold time by +5.2 ns mean (Tp = 100 ns) and the
  lead peak by +1%. **[synthetic, n = 3; not physics]**
* Consequence for the kernel: keep per-primary times as a per-wire charge-time histogram (2.5–10 ns bins) or as a few
  time-ordered sub-pulses per wire, not one pulse per wire. Kernel cost is linear in the number of pulses per wire, so
  3–5 sub-pulses cost about 3–5 times the T4 kernel (≈ 12–20 µs per crossing). That is still about 8–13 times cheaper than
  per-primary evaluation. The real T-ladder run decides the number.
* Spatial moments: ⟨sin θ⟩ is the leading asymmetry carrier (Spearman with A ≈ +0.98, with ⟨dy⟩ ≈ 1.00). Keep the
  first harmonic (K0, Ks, Kc). Whether ⟨sin 2θ⟩, ⟨cos 2θ⟩ and the radial moment matter is a kernel-accuracy test against
  the exact replay **[open]**.

## 5. Minimal cathode-kernel design that attaches to the chain
Per plane p ∈ {U, V}, offline:
`G_p,h(Δs, t) = (shaper ∗ induced current of a unit ion + its electron prompt step, harmonic h of the birth angle)`,
with h ∈ {0, sin, cos}. Here Δs is the normalized stereo distance between strip centre and avalanche
(s = (w ∓ 0.1 u)/sqrt(1.01)), and the time axis is the shaper output at 2.5 ns. Sign convention: the frozen
chargeSign = +1. The kernel is built from the same spectral LUT and RKF paths as the replay. Size is about 7.7 MB per
shaper setting.

Online, per crossing: generate per primary with the joint donor, then accumulate per (wire, time bin) the charge N,
N⟨sin⟩ and N⟨cos⟩, plus the mean w. Then for each strip in ±6 of the avalanche:
`V_strip(t) = Σ_bins [N G_0 + N⟨sin⟩ G_s + N⟨cos⟩ G_c](Δs, t − t_bin)`; superpose wires; add noise; digitize.
Linear in charge, so the shaper is applied once offline. Replacing the shaper only needs a rebuilt table.

## 6. Implementation and validation, in order
1. **Done here:** the T-ladder tools (`timing-ablation/`). Next: run on 5, then 20 development muons with the real
   exporter and mobility (WSL; see its README). Gate: T1 vs T0 below the reconstruction tolerance at both Tp values.
   Then pick the number of time sub-pulses per wire from TW vs T1.
2. Build G_h from the replay machinery, without new physics. Compare the harmonic kernel on the frozen paths with the
   T1 replay, muon by muon (peak, threshold time, centroid, w/u). Development muons only; sealed muons once, at the end.
3. Build the joint 2100 V donor from the V0 exports. Validate the generated per-wire sums (yield, A, onset spread,
   multiplicity, second-anode rate) against held-out development archives, with the `rho_w` check.
4. Time the full A + kernel chain in C++ with real tables, and measure real Heed (section 1).

## 7. What is not done
* No real 2100 V archive or real path export is reachable from this environment. All plots and timing ablations here
  are tool tests on synthetic inputs. The real ones run on WSL with the commands in each README.
* ROOT is not available here. The V0 ROOT reader is still checked only against stub headers.

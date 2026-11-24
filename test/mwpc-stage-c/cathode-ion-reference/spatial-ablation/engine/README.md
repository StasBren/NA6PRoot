# Cathode U/V response engine (2100 V K-model, consolidated)

Inputs: archived ion birth records (`ion_births.root`) and an explicitly identified precomputed kernel directory.
Outputs: unshaped U/V strip currents on a 10 ns grid, the ion term and the APPROXIMATE avalanche-electron step returned
separately. No development-muon index, run directory or S0 file is used in normal operation; the validation runner
(`run_s_ladder.sh engine`) is the independent test.

## Interface
```python
import cathode_engine as ce
k   = ce.KernelTable(KERNEL_DIR, expected_id=SHA256)         # once per process; refuses a different kernel
eng = ce.CathodeResponseEngine(k, model="K2A", half_strips=5, window_ns=7000)
m   = ce.read_moments_csv(EVENT.moments.csv)                 # from p3_birth_moments (ROOT birth reader)
#   or ce.moments_from_births(primary, u_cm, v_cm, w_cm, t_ns, weight, patch_shift_cm, u_in_mm, w_in_mm)
r   = eng.evaluate(m, impl="compiled" | "vector" | "loop")
r.ion, r.electron                                            # (bins, 2 planes, strips) fC/ns, frozen sign convention
r.labels, r.t_center_ns, r.meta                              # meta: kernel_id, transport_id, weighting_id, model, impl
shaped = ce.Electronics.crrc_standin(100.0).shape(r.ion, r.dt_ns)   # or Electronics.from_table(impulse.txt)
```
* **Kernel identity**: `kernel_id` = sha256 over the sha256 of `K.npy`, `KPHI0.npy`, `kernel.json`
  (`engine_cli.py info KERNEL_DIR`). `transport_id` = sha256 of the Green grid-path vertices (gas, field, wires);
  `weighting_id` = sha256 of the weighting LUT (strip geometry). A strip-geometry study produces kernels with the same
  transport_id and a different weighting_id; a new shaper needs no new kernel.
* **Implementations** (same K-model basis, corners and time interpolation): `loop` = the validated s_kernel.py
  algorithm, kept as regression reference; `vector` = numpy, vectorized over clouds x strips x corners;
  `compiled` = `kernel_eval.cc` (C ABI via ctypes; CMake target `kernel_eval`, or
  `g++ -O3 -std=c++17 -shared -fPIC -o libkernel_eval.so kernel_eval.cc`; `NA60_KERNEL_EVAL_LIB` points to it).
  Synthetic agreement with `loop`: rel L2 <= 3e-15 (ion), <= 2e-16 (electron) for K1, K2, K1A, K2A.
* **Birth reading**: `p3_birth_moments --batch LIST OUTDIR` reads many archives in one process (frozen archive reader,
  all census checks) and writes one moments CSV per event + `batch_timing.json` (cold = first event).
* **Not in the engine**: Geant4, Heed, avalanche generation, digitization. Physics models unchanged.

## CLI
`engine_cli.py info | run | bench` (see `--help`). `run` keeps one kernel in memory for a whole batch and writes
`engine_timing.json` (import, kernel verify/load, per-event moments-read and evaluate, cold/warm, peak RSS, warm
throughput). `bench` measures latency and temporary memory vs clouds per event.

## Tests
`test_engine.sh BIRTH_TEST_DIR WORKDIR` (SYNTHETIC, after `../test_birth_moments.sh`): numpy moments_from_births ==
BirthMoments.hh (1.9e-16); wrong kernel id refused; loop == s_kernel.py apply (bit-identical); vector, compiled == loop
(< 1e-13 required; 2-3e-15 observed); Electronics stand-in == t_shaped_compare shaper.

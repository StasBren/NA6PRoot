# Inputs for the K2A signal-response presentation

`collect_presentation_inputs.py` gathers, read-only, the small files the presentation needs into one new upload
folder: an aggregate validation report (any report, e.g. a new aggregate of the running queue; `--also-report` adds the
report files of the frozen snapshot), the per-slot outcomes, comparison metrics, shaped observables and timings, the
overnight queue, the non-sealed preflight inventory and every p3_birth_moments product (relative paths preserved), the
probe outputs, the kernel axes plus a strictly verified < 1 MB kernel slice, and waveforms for a few example muons.

Safety: inputs are only read; the upload folder must be new and outside every input directory; sealed slots are never
copied and a report containing one is refused before anything is written; the preflight inventory must be readable
and contain the sealed rows (otherwise exclusion cannot be verified). Ion paths, numpy tables, ROOT files and group
exports are never copied.

Checks: `CHECKS.tsv` (inventory, sealed exclusion, preflight coverage, one kernel / commit across runs, failing steps
of FAIL slots, kernel-slice verification) and `MISSING_ESSENTIAL.tsv` (per completed slot: preflight moments / event /
frame files, run directory, provenance, outcome, steps, timing, and for accuracy-set slots the S0_T1 comparison and
shaped observables). Exit 0 = COMPLETE, 2 = INCOMPLETE (folder written, gaps listed), 1 = refused (nothing written).

Test (SYNTHETIC, no real data): `python3 tests/collector_test.py NEW_WORK_DIR`.

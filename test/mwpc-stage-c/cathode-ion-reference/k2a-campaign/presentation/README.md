# Inputs for the K2A signal-response presentation

`collect_presentation_inputs.py` gathers, read-only, the small files the presentation needs into one new upload
folder: the aggregate validation report, the per-slot outcomes, comparison metrics, shaped observables and timings,
the overnight queue, the non-sealed preflight features (for the station / channel / complexity join), the probe
outputs, the kernel axes plus a < 1 MB kernel slice, and waveforms for a few example muons.

It never copies ion paths (`vertices.f64`), numpy tables, ROOT files or group exports, and it refuses any sealed
slot (from the preflight inventory). `MANIFEST.tsv` lists every copied file with its sha256; `SKIPPED.tsv` lists what
was left out and why.

The synthetic check used during development: a fake two-run layout with a sealed preflight row, an ion-path file and a
kernel; the sealed moments and the ion path were skipped, examples picked as the median and largest model U shaped L2
plus 0637, and an existing output folder refused.

# Sources of the K2A signal-response deck (28 Nov 2026)

Scripts that produced `K2A_signal_response.pptx` from the collected upload `upload_20261128_143502`
(aggregate report `analysis_20261128_143452`, 195 slots; snapshot `analysis_20261128_115131` for comparison).

- `build_tables.py`: joins per_muon.csv, outcome.json, shaped observables, preflight features and the queue.
- `stats.py`: every number quoted on the slides (stats.json).
- `style.py`, `figs.py`: ATLAS-style figures; waveform shaping imports `timing-ablation/t_shaped_compare.py` unchanged
  (shaped L2 of the examples reproduces the reported values exactly, 0-6 us window).
- `deck.js` (pptxgenjs) and `embed_fonts.py` (embeds the Public Sans subset of the avalanche deck).

Paths are those of the build environment (upload unpacked under /home/claude/pres/data); adjust before rerunning.
Read-only with respect to all validation inputs.

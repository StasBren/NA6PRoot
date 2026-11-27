# Read-only K2A probes

These probes test hypotheses about the K2A model without changing the frozen kernel, Green table, S0 references or sealed set. They run no avalanche, no ion transport and no S0 replay.

| probe | question | inputs | cost |
|---|---|---|---|
| `k2a_wire_probe.py` | Where does an event's charge sit in the nine-wire patch? Which of its wires are outside the kernel? | preflight moments and their `.event.json` / `.frame.json`; `kernel.json`; optionally the domain scan | seconds |
| `k2a_kernel_probe.py` | How large are the tau- and spatial-interpolation errors? Is the wire-translation convention right? How long are the drift times (how much time a configurable horizon must cover)? | `kernel_TAG`; `green_paths_TAG` (`grid.csv`, `ions.csv`) | seconds to minutes |
| `k2a_closure_anatomy.py` | Which spatial simplification produces the early-time residual: Green interpolation (G), single radius (SR), angular truncation (S3/S3A) or longitudinal collapse (SW)? | existing S-ladder variants (`RUN/muon_*/<variant>_TAG`) or a layout from `k2a_group_synth.sh` | seconds |
| `k2a_group_synth.sh` | The same decomposition for one validation muon (for example 0637) | that muon's existing group exports and S0 groups; `green_TAG` | Python, per ion (minutes) |

The tests are in `tests/probes_test.py`, which part K of `validation/tests/test_validation_synthetic.sh` runs on SYNTHETIC data.

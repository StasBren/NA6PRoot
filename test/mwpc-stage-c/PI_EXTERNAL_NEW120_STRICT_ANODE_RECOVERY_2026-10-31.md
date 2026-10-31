# Independent 120-muon holdout: strict anode-metal recovery of completed physical chunk_07

**2026-10-31. This is a forensic software-validation exception, not a new physics model or a modified Garfield result.**

## User-observed WSL evidence

- New independent Geant4 parent identities: 120, disjoint from older 120, all 12 groups preflight-passed TrackHeed.
- New full Garfield native physics successfully SHA-sealed **seven** 10-parent groups, chunk_00 through chunk_06; 70 different muon parents, **15,016** original microscopic avalanche calls, zero capped/failed checkpoints. Original compiled C++ source/binary SHA `1a7d1e62f084524d12c4003c3876319dd20ebd83a38f0697da2895b943539eef`.
- Native chunk_07 (8th group) wrote its physics and C++ completion markers before an **old Python validator** aborted on `jpsi40-bore|2|8904|5|2|5|266`, `p3_noFW`, primary 7: `u=0.00237593484574722 mm`, `v=0.0107657602328344 mm`, `successful=0, capped=0, ne=1, ni=0, ion_births=0, collected_e=1, attached=0, other=0`, zero directional ion sectors. Distance from nearest anode axis **0.01102482 mm = 11.025 um**, strictly inside the anode wire whose diameter is 30 um, radius 15 um. This is consistent with Garfield's no-avalanche, collected-origin geometry outcome, NOT proof that all other 8th-chunk records are valid.
- Old `pi_run_muon_lut_validation.validate_reference` had a **narrowly defined** permitted-in-metal exception only for an original `75um` **field** wire, not a `30um` anode. **DO NOT modify this original validator** (frozen source chain). Never delete or rerun the unsealed chunk_07; need independent whole-chunk proof before sealing it.

## Separate strict recovery and continuation helper

New file `test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py` (GitHub commit `9b0145b`) runs the bytecode of the original trusted `validate_reference` in a private copied globals namespace with one additional check for original metal-origin no-ion, no-cap, no-attachment, `ne=1, ni=0, collected_e=1, other=0` anode starts within **radius 15um minus 1nm margin** at the nine true anode centers. The original `field_wire_blocked_zero_start` helper is untouched. All other original row and cluster conservation validators, plus `check_crossings`, remain unchanged. The exact former original status and charge census remain in physical CSV, not normalized to artificial successes.

The helper independently checks all 7 previous SHA-sealed chunks, original new parent selection/Geant4 SHA, C++ source/executable SHA, 10-parent original input bytes, all physical CSVs, the C++ full completion markers, absence of caps and all 4×Heed primary seeds. It explicitly counts physical **anode** and **field-wire** blocked starts separately; anything nonmetal/capped/nonzero-ion fails closed. `audit` reads and prints details, **does not write**. `recover` repeats audit and (only when successful) creates a **new, previously absent** `chunk_07/complete.json` with exact SHA of all original physics/log/input files. The original native physics timer was lost when prior validator aborted, and is **NOT fabricated**. Additional source code is separate; old `pi_external_muon120_guarded_lut_v1.py` and original 120 scripts remain unchanged (their hashes were frozen before seeing the new physics).

### First action: strict read-only forensic audit, no physics

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py selftest
python3 -u test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py audit
)
```

If audit fails: stop. Do **not** delete/modify chunk_07, do **not** rerun full Garfield; inspect report and original `chunk_07/physics.log`.

### Second action: only after full audit PASS

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py recover
python3 -u test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py status
)
```

If recovered, the SHA-verified 8th native physical checkpoint is present. This writes only missing `complete.json` and avoids recomputing native Garfield.

### Third action: resume from ninth native block using this helper

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py run --max-chunks 12 --hours 8
python3 test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py status
)
```

Run is resumable; original 8 completed blocks are reused and SHA checked. The physical C++ remains identical; only the native Python validation in the helper's **process-local isolated copy** is allowed to classify the same narrow, geometrically provable in-metal origin. The original runner's `run_config.json`, output CSVs, raw physics logs, original target library and original old120 muons are never changed.

If all 12 chunks seal and new reference is final, use **the same helper**, NOT direct original external wrapper, for blind evaluation; it installs the same strict geometry-aware validation in process memory before verifying the aggregated reference:

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py evaluate --reps 3000
python3 test/mwpc-stage-c/pi_external_muon120_wire_geometry_recovery_v1.py status
)
```

This keeps A20/A100 and guarded spatial-pooling 5x5 parameters **fixed** from before external truth. Final external accuracy is *pending* actual reference completion. The whole experiment is B=0 and ion-birth sector census only, not cathode induced U/V strip current.

**Audit status**: code compiled syntactically in the assistant sandbox; independent AST-extracted geometry/flags tests PASS for the exact user coordinates. Real new120 chunk_07 native audit **has not yet been executed** in user's WSL. No 8th-chunk SHA seal yet. Do not describe the 8th block as recovered until user confirms PASS.

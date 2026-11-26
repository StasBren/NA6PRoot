#!/usr/bin/env python3
"""Read-only, four-level data handoff from the verified 2100 V V0 population.

Output is a compact ZIP of CSV tables, NOT a drop-in wiretruth1800_full.root.
Original Garfield++ archives are never modified or simulated again.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from per_anode_census_2100 import read_per_anode


def need(flag, message):
    if not flag:
        raise RuntimeError("NA60+ Claude 2100 V handoff: " + message)


def collect_original(roots, filename, output):
    header = None
    count = 0
    with output.open("w", newline="") as target:
        writer = None
        for source in roots:
            slot_match = re.fullmatch(r"muon_(\d+)", source.parent.name)
            need(slot_match is not None, "bad muon source path: " + str(source))
            slot = int(slot_match.group(1))
            path = source.parent / filename
            need(path.is_file(), "original source file not found: " + str(path))
            with path.open(newline="") as f:
                rows = csv.DictReader(f)
                need(rows.fieldnames is not None, "empty original CSV: " + str(path))
                names = ["muon_slot"] + list(rows.fieldnames)
                if header is None:
                    header = names
                    writer = csv.DictWriter(target, fieldnames=names)
                    writer.writeheader()
                need(names == header, "source CSV schema drift: " + str(path))
                for row in rows:
                    writer.writerow({"muon_slot": slot, **row})
                    count += 1
    return count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, required=True,
                    help="Completed 2100_gain_asym/population_* analysis folder")
    ap.add_argument("--out", type=Path, required=True,
                    help="Fresh destination, e.g. ~/NA60_2100V_997_for_Claude.zip")
    args = ap.parse_args()
    run = args.run.expanduser().resolve()
    out = args.out.expanduser().resolve()
    need(not out.exists(), "output already exists: " + str(out))
    source = run / "rootfiles.txt"
    census = run / "census"
    need(source.is_file(), "missing V0 archive manifest " + str(source))
    files = [Path(line.strip()) for line in source.read_text().splitlines()
             if line.strip()]
    need(len(files) == 997 and len(set(files)) == 997,
         "requires exactly 997 distinct verified original muon paths")
    for name in ("v0_muons.csv", "v0_primaries.csv", "v0_fragments.csv",
                 "v0_wires.csv"):
        need((census / name).is_file(), "missing complete V0 export " + name)
    mu = pd.read_csv(census / "v0_muons.csv")
    primary = pd.read_csv(census / "v0_primaries.csv")
    need(len(mu) == 997 and mu["muon_slot"].is_unique,
         "V0 muon census does not contain 997 distinct events")
    need(set(mu["muon_slot"].astype(int)) ==
         {int(re.fullmatch(r"muon_(\d+)", p.parent.name).group(1))
          for p in files}, "V0 muon IDs differ from original ROOT list")
    need(not primary[["muon_slot", "primary"]].duplicated().any(),
         "duplicate primary in V0 census")

    # Existing audited loader recomputes original per-wire A and time
    # accounting and stops on ROOT-vs-CSV ion/electron mismatches.
    primary_wire, muon_wire = read_per_anode(run, primary, mu)
    need(int(primary_wire["ion_births"].sum()) ==
         int(mu["ion_births_exact"].sum()),
         "per-primary wire ion births do not close globally")
    need(int(primary_wire["collected_e"].sum()) ==
         int(mu["collected_electrons"].sum()),
         "per-primary wire collected electrons do not close globally")

    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="na60_claude_2100_") as tmp:
        d = Path(tmp)
        # Four levels correspond to the four 1800 V conceptual ROOT trees;
        # column names reflect 2100 V V0 schema, NOT identical 1800 branches.
        mu.to_csv(d / "muons.csv", index=False)
        primary.to_csv(d / "seeds.csv", index=False)
        muon_wire.to_csv(d / "wires.csv", index=False)
        primary_wire.to_csv(d / "seed_wires.csv", index=False)
        fragment = pd.read_csv(census / "v0_fragments.csv")
        fragment.to_csv(d / "k256_fragment_moments.csv", index=False)
        pd.read_csv(census / "v0_wires.csv").to_csv(
            d / "v0_wire_moments.csv", index=False)
        counters = {}
        for filename, output in (
            ("microscopic_seed.csv", "original_avalanche_status.csv"),
            ("microscopic_heed_primaries.csv", "original_heed_primaries.csv"),
            ("input.csv", "original_geant4_inputs.csv"),
        ):
            counters[output] = collect_original(files, filename, d / output)
        for label in ("v0_summary.txt", "v0_manifest.txt"):
            f = census / label
            if f.is_file():
                (d / label).write_bytes(f.read_bytes())
        (d / "README_for_Claude.txt").write_text(
            "NA60+ 2100 V / 997 original verified muon archives.\n\n"
            "This compact ZIP contains CSV tables, not a drop-in copy of\n"
            "the 1800 V wiretruth ROOT file. The four conceptual levels are\n"
            "muons, wires, seeds, seed_wires. Branch names in the 2100 V V0\n"
            "census differ; use the CSV headers and join using muon_slot,\n"
            "primary, and global_wire as appropriate.\n\n"
            "wires.csv and seed_wires.csv: exact Garfield ion-birth and true\n"
            "collecting-anode electron counters, exact A_UD, source timing.\n"
            "seeds.csv: all Heed primaries including zero-gain avalanches.\n"
            "original_avalanche_status.csv: Garfield success/capped flags,\n"
            "initial electron coordinates and loss counters.\n"
            "original_heed_primaries.csv: source primary information.\n"
            "original_geant4_inputs.csv: original muon crossing metadata.\n"
            "k256_fragment_moments.csv: sampled ion-birth moments, with\n"
            "N_hat, effective sample size and uncertainty; NOT full-ion\n"
            "positions. The ion-hemisphere census and total gains are exact.\n\n"
            "A_UD=(ion_up-ion_down)/(ion_up+ion_down); v>0 up, v<0 down;\n"
            "on-axis births excluded from A denominator only. Ion-to-wire\n"
            "assignment is NEAREST ANODE AT BIRTH, not ion-drift endpoint;\n"
            "collected electrons are assigned to physical collection wire.\n"
            "Mean production delay is mean ion birth time minus its own\n"
            "Heed primary start time, not subsequent ion drift time.\n"
            "Model: idealized 9-wire P3 noFW, 3mm pitch, 30um anode,\n"
            "cathodes +/-3mm, 2100V, Ar/CO2 70/30, Penning r=.57, B=0.\n"
            "Not the as-built prototype chamber. Descriptive 997 sample\n"
            "includes former holdouts; do not treat as independent validation.\n"
            "No 997-event ion drift waveforms / U/V strip signals here.\n"
        )
        manifest = {
            "voltage_V": 2100,
            "muons": len(mu),
            "primary_rows": len(primary),
            "primary_wire_rows": len(primary_wire),
            "muon_wire_rows": len(muon_wire),
            "k256_fragment_rows": len(fragment),
            "original_csv_counts": counters,
            "source_list_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "source_run": str(run),
            "format": "compact CSV level-equivalents (NOT 1800 ROOT branch-identical)",
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        (d / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        with zipfile.ZipFile(out, mode="x", compression=zipfile.ZIP_DEFLATED,
                             compresslevel=6, allowZip64=True) as z:
            for f in sorted(d.iterdir()):
                z.write(f, f.name)
    print("PASS 2100 V 997-muon source-validated compact handoff", flush=True)
    print("FILE=" + str(out), flush=True)
    print("SIZE_MB=" + str(round(out.stat().st_size / 1.e6, 2)), flush=True)
    print("MUONS=" + str(len(mu)) + " PRIMARIES=" + str(len(primary)),
          flush=True)


if __name__ == "__main__":
    main()

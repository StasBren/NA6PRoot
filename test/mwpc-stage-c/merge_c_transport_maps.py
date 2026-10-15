#!/usr/bin/env python3
"""
Merge a complete baseline Stage-C transport map with targeted high-stat boxes.

Rows from the override map replace the corresponding (iu,iv) boxes from the
baseline map.  This is intended for a hybrid response kernel such as:
  * 50 microscopic electrons/box over the full 20x30 grid;
  * 2000 electrons/box in selected boundary u-columns.

The output remains a single ordinary c_transport_box_electrons.csv file.
Downstream folding samples from whatever number of microscopic trials is
available in each box.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def key(row):
    return int(float(row["iu"])), int(float(row["iv"]))


def scan(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        if not r.fieldnames:
            raise SystemExit(f"No header in {path}")
        counts = Counter()
        centers = {}
        for row in r:
            k = key(row)
            counts[k] += 1
            centers[k] = (
                float(row["u_center_mm"]),
                float(row["v_center_mm"]),
            )
    if not counts:
        raise SystemExit(f"No rows in {path}")
    return list(r.fieldnames), counts, centers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True, type=Path)
    ap.add_argument("--override", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    base_fields, base_counts, base_centers = scan(args.baseline)
    over_fields, over_counts, over_centers = scan(args.override)

    if base_fields != over_fields:
        raise SystemExit("Baseline/override CSV schemas differ.")

    missing = sorted(set(over_counts) - set(base_counts))
    if missing:
        raise SystemExit(
            f"Override contains {len(missing)} boxes absent from baseline; "
            f"first={missing[0]}"
        )

    for k in over_counts:
        bu, bv = base_centers[k]
        ou, ov = over_centers[k]
        if abs(bu - ou) > 1e-10 or abs(bv - ov) > 1e-10:
            raise SystemExit(
                f"Box centre mismatch for {k}: "
                f"baseline=({bu},{bv}) override=({ou},{ov})"
            )

    override_keys = set(over_counts)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    base_rows_written = 0
    over_rows_written = 0

    with args.output.open("w", newline="", encoding="utf-8") as fout:
        w = csv.DictWriter(fout, fieldnames=base_fields)
        w.writeheader()

        with args.baseline.open(newline="", encoding="utf-8") as f:
            r = csv.DictReader(f)
            for row in r:
                if key(row) in override_keys:
                    continue
                w.writerow(row)
                base_rows_written += 1

        with args.override.open(newline="", encoding="utf-8") as f:
            r = csv.DictReader(f)
            for row in r:
                w.writerow(row)
                over_rows_written += 1

    hybrid_counts = Counter(base_counts)
    for k in override_keys:
        hybrid_counts[k] = over_counts[k]

    print("Hybrid transport map COMPLETE")
    print(f"  baseline boxes             : {len(base_counts)}")
    print(f"  override boxes             : {len(over_counts)}")
    print(f"  hybrid boxes               : {len(hybrid_counts)}")
    print(
        "  baseline trials/box       : "
        + ",".join(map(str, sorted(set(base_counts.values()))))
    )
    print(
        "  override trials/box       : "
        + ",".join(map(str, sorted(set(over_counts.values()))))
    )
    print(
        "  hybrid trial levels       : "
        + ",".join(map(str, sorted(set(hybrid_counts.values()))))
    )
    print(f"  baseline rows retained     : {base_rows_written}")
    print(f"  override rows inserted     : {over_rows_written}")
    print(f"  wrote                      : {args.output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Fold the PI geometric Poisson/Heed clusters through the Stage-C Garfield
(u,v) transport-box map.

This is the bridge between the two model layers:

  GEANT track
    -> Poisson clusters + Heed multiplicity
    -> geometric wire region (midpoint rule)
    -> local (u,v) box inside that wire cell
    -> empirical Garfield response kernel
         P(destination wire, drift time, collection | u,v)

The geometric wire assignment is retained as a reference.  Garfield is allowed
to move a primary electron to a neighbouring wire.  This lets us quantify
exactly where the midpoint/binomial approximation succeeds or fails.

The transport map is sampled with replacement.  Therefore a very large
synthetic event sample does not create more microscopic Garfield information
than the finite number of trials stored in each box.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


def as_int(row, key):
    return int(float(row[key]))


def as_float(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def load_transport(path: Path):
    boxes = defaultdict(list)
    u_centers = {}
    v_centers = {}

    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            iu = as_int(row, "iu")
            iv = as_int(row, "iv")
            boxes[(iu, iv)].append({
                "collected": as_int(row, "collected") == 1,
                "wire_index": as_int(row, "wire_index"),
                "drift_time_ns": as_float(row, "drift_time_ns"),
            })
            u_centers[iu] = as_float(row, "u_center_mm")
            v_centers[iv] = as_float(row, "v_center_mm")

    if not boxes:
        raise SystemExit("Transport map is empty.")

    iu_values = sorted(u_centers)
    iv_values = sorted(v_centers)

    def spacing(values, name):
        if len(values) < 2:
            raise SystemExit(f"Need at least two {name} bins.")
        ds = [values[i + 1] - values[i] for i in range(len(values) - 1)]
        d = sum(ds) / len(ds)
        if any(abs(x - d) > 1e-8 for x in ds):
            raise SystemExit(f"Non-uniform {name} grid.")
        return d

    us = [u_centers[i] for i in iu_values]
    vs = [v_centers[i] for i in iv_values]
    du = spacing(us, "u")
    dv = spacing(vs, "v")

    missing = [
        (iu, iv)
        for iu in iu_values
        for iv in iv_values
        if (iu, iv) not in boxes
    ]
    if missing:
        raise SystemExit(
            f"Incomplete transport map: {len(missing)} missing boxes."
        )

    trial_counts = [len(v) for v in boxes.values()]
    if any(n < 1 for n in trial_counts):
        raise SystemExit("Transport map contains an empty box.")

    return boxes, {
        "iu_min": iu_values[0],
        "iv_min": iv_values[0],
        "n_u": len(iu_values),
        "n_v": len(iv_values),
        "du": du,
        "dv": dv,
        "u_min_edge": us[0] - 0.5 * du,
        "v_min_edge": vs[0] - 0.5 * dv,
        "u_max_edge": us[-1] + 0.5 * du,
        "v_max_edge": vs[-1] + 0.5 * dv,
        "trials_per_box_min": min(trial_counts),
        "trials_per_box_max": max(trial_counts),
        "trials_per_box_unique": sorted(set(trial_counts)),
    }


def box_offset(x, min_edge, step, n, tol=1e-9):
    xmax = min_edge + n * step
    if x < min_edge - tol or x > xmax + tol:
        return None
    if x >= xmax:
        return n - 1
    i = int(math.floor((x - min_edge) / step))
    return i if 0 <= i < n else None


def mean(s, n):
    return s / n if n else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--transport-electrons", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    ap.add_argument("--wire-pitch-mm", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=20261015)
    args = ap.parse_args()

    if args.wire_pitch_mm <= 0:
        raise SystemExit("Wire pitch must be positive.")

    rng = random.Random(args.seed)
    boxes, meta = load_transport(args.transport_electrons)

    print("PI geometric clusters + Garfield box map")
    print(
        f"  box grid                    : {meta['n_u']} x {meta['n_v']} "
        f"({meta['du']:.3g} x {meta['dv']:.3g} mm)"
    )
    if meta["trials_per_box_min"] == meta["trials_per_box_max"]:
        print(
            f"  microscopic trials / box    : "
            f"{meta['trials_per_box_min']}"
        )
    else:
        print(
            f"  microscopic trials / box    : "
            f"{meta['trials_per_box_min']} .. "
            f"{meta['trials_per_box_max']} "
            f"(hybrid map; levels={meta['trials_per_box_unique']})"
        )
    print(
        f"  u coverage [mm]             : "
        f"{meta['u_min_edge']:.3g} .. {meta['u_max_edge']:.3g}"
    )
    print(
        f"  v coverage [mm]             : "
        f"{meta['v_min_edge']:.3g} .. {meta['v_max_edge']:.3g}"
    )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    event_path = Path(str(args.output_prefix) + "_crossings.csv")
    cluster_path = Path(str(args.output_prefix) + "_clusters.csv")
    wire_path = Path(str(args.output_prefix) + "_wires.csv")

    event_fields = [
        "crossing_id", "event", "track_id", "det_id", "station",
        "primary_electrons", "transported_electrons", "collected_electrons",
        "collection_fraction", "out_of_map_electrons",
        "geometric_active_wires", "garfield_active_wires",
        "geometric_dominant_wire", "garfield_dominant_wire",
        "dominant_wire_changed",
        "migrated_collected_electrons", "migration_fraction_collected",
        "geometric_wire_cog_mm", "garfield_wire_cog_mm",
        "mean_drift_time_ns",
        "fraction_collected_before_100ns",
        "fraction_collected_before_200ns",
    ]
    cluster_fields = [
        "crossing_id", "cluster_index",
        "event", "track_id", "det_id", "station",
        "geometric_wire_index", "u_cell_mm",
        "u_mm", "v_mm", "w_mm",
        "n_primary_electrons",
        "transported_electrons", "collected_electrons",
        "stayed_electrons", "migrated_electrons",
        "migration_fraction_collected",
        "mean_drift_time_ns",
        "fraction_collected_before_100ns",
        "fraction_collected_before_200ns",
        "out_of_map_electrons",
    ]
    wire_fields = [
        "crossing_id", "event", "track_id", "det_id", "station",
        "wire_index", "wire_u_mm",
        "geometric_primary_electrons", "geometric_primary_fraction",
        "garfield_collected_electrons", "garfield_collected_fraction",
        "mean_drift_time_ns",
        "fraction_before_100ns", "fraction_before_200ns",
    ]

    def new_state(row):
        return {
            "meta": {
                "crossing_id": as_int(row, "crossing_id"),
                "event": as_int(row, "event"),
                "track_id": as_int(row, "track_id"),
                "det_id": as_int(row, "det_id"),
                "station": as_int(row, "station"),
            },
            "primary": 0,
            "transported": 0,
            "collected": 0,
            "out": 0,
            "migrated": 0,
            "geom": defaultdict(int),
            "dest": defaultdict(
                lambda: {
                    "n": 0,
                    "sum_dt": 0.0,
                    "n_dt": 0,
                    "before100": 0,
                    "before200": 0,
                }
            ),
            "sum_dt": 0.0,
            "n_dt": 0,
            "before100": 0,
            "before200": 0,
        }

    total_events = 0
    total_primary = 0
    total_collected = 0
    total_migrated = 0
    total_out = 0

    with (
        args.clusters.open(newline="", encoding="utf-8") as fin,
        event_path.open("w", newline="", encoding="utf-8") as fe,
        cluster_path.open("w", newline="", encoding="utf-8") as fc,
        wire_path.open("w", newline="", encoding="utf-8") as fw,
    ):
        reader = csv.DictReader(fin)
        required = {
            "crossing_id", "cluster_index",
            "event", "track_id", "det_id", "station",
            "wire_index", "u_mm", "v_mm", "w_mm",
            "n_primary_electrons",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Cluster CSV missing: " + ", ".join(sorted(missing))
            )

        ew = csv.DictWriter(fe, fieldnames=event_fields)
        cw = csv.DictWriter(fc, fieldnames=cluster_fields)
        ww = csv.DictWriter(fw, fieldnames=wire_fields)
        ew.writeheader()
        cw.writeheader()
        ww.writeheader()

        current = None
        state = None

        def flush(s):
            nonlocal total_events, total_primary, total_collected
            nonlocal total_migrated, total_out
            if s is None:
                return

            m = s["meta"]
            primary = s["primary"]
            collected = s["collected"]

            geom_active = [k for k, n in s["geom"].items() if n > 0]
            dest_active = [k for k, d in s["dest"].items() if d["n"] > 0]

            if geom_active:
                geom_dom = max(
                    sorted(geom_active),
                    key=lambda k: s["geom"][k],
                )
                geom_cog = (
                    sum(k * args.wire_pitch_mm * n for k, n in s["geom"].items())
                    / primary
                    if primary
                    else math.nan
                )
            else:
                geom_dom = 999999
                geom_cog = math.nan

            if dest_active:
                dest_dom = max(
                    sorted(dest_active),
                    key=lambda k: s["dest"][k]["n"],
                )
                garf_cog = (
                    sum(
                        k * args.wire_pitch_mm * d["n"]
                        for k, d in s["dest"].items()
                    )
                    / collected
                )
            else:
                dest_dom = 999999
                garf_cog = math.nan

            ew.writerow({
                **m,
                "primary_electrons": primary,
                "transported_electrons": s["transported"],
                "collected_electrons": collected,
                "collection_fraction": (
                    collected / s["transported"]
                    if s["transported"] else math.nan
                ),
                "out_of_map_electrons": s["out"],
                "geometric_active_wires": len(geom_active),
                "garfield_active_wires": len(dest_active),
                "geometric_dominant_wire": geom_dom,
                "garfield_dominant_wire": dest_dom,
                "dominant_wire_changed": int(
                    geom_dom != 999999
                    and dest_dom != 999999
                    and geom_dom != dest_dom
                ),
                "migrated_collected_electrons": s["migrated"],
                "migration_fraction_collected": (
                    s["migrated"] / collected
                    if collected else math.nan
                ),
                "geometric_wire_cog_mm": geom_cog,
                "garfield_wire_cog_mm": garf_cog,
                "mean_drift_time_ns": mean(s["sum_dt"], s["n_dt"]),
                "fraction_collected_before_100ns": (
                    s["before100"] / collected if collected else math.nan
                ),
                "fraction_collected_before_200ns": (
                    s["before200"] / collected if collected else math.nan
                ),
            })

            all_wires = sorted(set(s["geom"]) | set(s["dest"]))
            for k in all_wires:
                d = s["dest"][k]
                ww.writerow({
                    **m,
                    "wire_index": k,
                    "wire_u_mm": k * args.wire_pitch_mm,
                    "geometric_primary_electrons": s["geom"][k],
                    "geometric_primary_fraction": (
                        s["geom"][k] / primary if primary else math.nan
                    ),
                    "garfield_collected_electrons": d["n"],
                    "garfield_collected_fraction": (
                        d["n"] / collected if collected else math.nan
                    ),
                    "mean_drift_time_ns": mean(d["sum_dt"], d["n_dt"]),
                    "fraction_before_100ns": (
                        d["before100"] / d["n"] if d["n"] else math.nan
                    ),
                    "fraction_before_200ns": (
                        d["before200"] / d["n"] if d["n"] else math.nan
                    ),
                })

            total_events += 1
            total_primary += primary
            total_collected += collected
            total_migrated += s["migrated"]
            total_out += s["out"]

        for row in reader:
            cid = as_int(row, "crossing_id")
            if current is None or cid != current:
                flush(state)
                current = cid
                state = new_state(row)

            geom_wire = as_int(row, "wire_index")
            u = as_float(row, "u_mm")
            v = as_float(row, "v_mm")
            wpos = as_float(row, "w_mm")
            ne = as_int(row, "n_primary_electrons")
            u_cell = u - geom_wire * args.wire_pitch_mm

            state["primary"] += ne
            state["geom"][geom_wire] += ne

            iu0 = box_offset(
                u_cell, meta["u_min_edge"], meta["du"], meta["n_u"]
            )
            iv0 = box_offset(
                v, meta["v_min_edge"], meta["dv"], meta["n_v"]
            )

            transported = 0
            collected = 0
            stayed = 0
            migrated = 0
            out = 0
            sum_dt = 0.0
            n_dt = 0
            before100 = 0
            before200 = 0

            if iu0 is None or iv0 is None:
                out = ne
                state["out"] += ne
            else:
                iu = meta["iu_min"] + iu0
                iv = meta["iv_min"] + iv0
                trials = boxes[(iu, iv)]

                for _ in range(ne):
                    transported += 1
                    state["transported"] += 1

                    response = rng.choice(trials)
                    if not response["collected"]:
                        continue

                    collected += 1
                    state["collected"] += 1

                    destination = geom_wire + response["wire_index"]
                    d = state["dest"][destination]
                    d["n"] += 1

                    if destination == geom_wire:
                        stayed += 1
                    else:
                        migrated += 1
                        state["migrated"] += 1

                    dt = response["drift_time_ns"]
                    if math.isfinite(dt) and dt >= 0:
                        sum_dt += dt
                        n_dt += 1
                        state["sum_dt"] += dt
                        state["n_dt"] += 1
                        d["sum_dt"] += dt
                        d["n_dt"] += 1

                        if dt <= 100:
                            before100 += 1
                            state["before100"] += 1
                            d["before100"] += 1
                        if dt <= 200:
                            before200 += 1
                            state["before200"] += 1
                            d["before200"] += 1

            cw.writerow({
                "crossing_id": cid,
                "cluster_index": as_int(row, "cluster_index"),
                "event": as_int(row, "event"),
                "track_id": as_int(row, "track_id"),
                "det_id": as_int(row, "det_id"),
                "station": as_int(row, "station"),
                "geometric_wire_index": geom_wire,
                "u_cell_mm": u_cell,
                "u_mm": u,
                "v_mm": v,
                "w_mm": wpos,
                "n_primary_electrons": ne,
                "transported_electrons": transported,
                "collected_electrons": collected,
                "stayed_electrons": stayed,
                "migrated_electrons": migrated,
                "migration_fraction_collected": (
                    migrated / collected if collected else math.nan
                ),
                "mean_drift_time_ns": mean(sum_dt, n_dt),
                "fraction_collected_before_100ns": (
                    before100 / collected if collected else math.nan
                ),
                "fraction_collected_before_200ns": (
                    before200 / collected if collected else math.nan
                ),
                "out_of_map_electrons": out,
            })

        flush(state)

    print()
    print("PI geometric + Garfield folding COMPLETE")
    print(f"  crossings                    : {total_events}")
    print(f"  primary electrons             : {total_primary}")
    print(f"  collected electrons           : {total_collected}")
    print(
        f"  collection fraction           : "
        f"{total_collected / (total_primary - total_out):.6f}"
        if total_primary > total_out else
        "  collection fraction           : nan"
    )
    print(f"  out-of-map electrons          : {total_out}")
    print(f"  migrated collected electrons  : {total_migrated}")
    print(
        f"  migration fraction            : "
        f"{total_migrated / total_collected:.6f}"
        if total_collected else
        "  migration fraction            : nan"
    )
    print(f"  crossing summary              : {event_path}")
    print(f"  cluster response              : {cluster_path}")
    print(f"  wire comparison               : {wire_path}")


if __name__ == "__main__":
    main()

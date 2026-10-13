#!/usr/bin/env python3
"""
Fold PI-style synthetic ionisation clusters through the empirical Stage-C
Garfield transport-box map.

Inputs
------
1) c_pi_clusters_heed2p4.csv
   Realistic track kinematics from the large MNP33 GEANT production, with
   Poisson-distributed cluster positions and empirical Heed cluster sizes.

2) c_transport_box_electrons.csv
   Single-electron Garfield responses for one 4-mm wire cell and the 2+4 mm
   gas gap.  Each box contains many microscopic trials, including collection,
   attachment/loss, destination wire, drift time and arrival angle.

For every primary electron in every synthetic cluster, sample one microscopic
response from the matching (u,v) box.  The local response wire index is shifted
by the cluster's wire_offset so the event can cross cell boundaries.

Outputs
-------
<prefix>_event_summary.csv
<prefix>_wire_summary.csv

This is a fast digitisation/response model.  It intentionally stops before
avalanche-gain fluctuations, cathode induction and electronics.
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
    except (KeyError, ValueError, TypeError):
        return math.nan


def load_transport(path: Path):
    boxes = defaultdict(list)
    u_centers = {}
    v_centers = {}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            iu = as_int(row, "iu")
            iv = as_int(row, "iv")
            boxes[(iu, iv)].append(
                {
                    "collected": as_int(row, "collected") == 1,
                    "status": as_int(row, "status"),
                    "wire_index": as_int(row, "wire_index"),
                    "drift_time_ns": as_float(row, "drift_time_ns"),
                    "cos_phi_hit": as_float(row, "cos_phi_hit"),
                    "hit_side_u": as_int(row, "hit_side_u"),
                    "incidence_normal_deg": as_float(
                        row, "incidence_normal_deg"
                    ),
                }
            )
            u_centers[iu] = as_float(row, "u_center_mm")
            v_centers[iv] = as_float(row, "v_center_mm")

    if not boxes:
        raise SystemExit("Transport map is empty.")

    iu_values = sorted(u_centers)
    iv_values = sorted(v_centers)
    if iu_values != list(range(min(iu_values), max(iu_values) + 1)):
        raise SystemExit("Transport map has missing u-bin indices.")
    if iv_values != list(range(min(iv_values), max(iv_values) + 1)):
        raise SystemExit("Transport map has missing v-bin indices.")

    # Require a rectangular complete map, not a partial checkpoint.
    missing = [
        (iu, iv)
        for iu in iu_values
        for iv in iv_values
        if (iu, iv) not in boxes
    ]
    if missing:
        raise SystemExit(
            f"Transport map is incomplete: {len(missing)} boxes missing; "
            f"first missing={missing[0]}"
        )

    u_sorted = [u_centers[i] for i in iu_values]
    v_sorted = [v_centers[i] for i in iv_values]

    def spacing(values, axis):
        if len(values) < 2:
            raise SystemExit(f"Need at least two {axis} bins.")
        ds = [values[i + 1] - values[i] for i in range(len(values) - 1)]
        d0 = sum(ds) / len(ds)
        if any(abs(d - d0) > 1e-8 for d in ds):
            raise SystemExit(f"Non-uniform {axis} grid is not supported.")
        return d0

    du = spacing(u_sorted, "u")
    dv = spacing(v_sorted, "v")
    u_min_edge = u_sorted[0] - 0.5 * du
    v_min_edge = v_sorted[0] - 0.5 * dv
    u_max_edge = u_sorted[-1] + 0.5 * du
    v_max_edge = v_sorted[-1] + 0.5 * dv

    trials = {len(v) for v in boxes.values()}
    if len(trials) != 1:
        raise SystemExit(
            "Transport boxes do not all have the same number of trials: "
            f"{sorted(trials)}"
        )

    meta = {
        "iu_min": iu_values[0],
        "iu_max": iu_values[-1],
        "iv_min": iv_values[0],
        "iv_max": iv_values[-1],
        "n_u": len(iu_values),
        "n_v": len(iv_values),
        "du": du,
        "dv": dv,
        "u_min_edge": u_min_edge,
        "u_max_edge": u_max_edge,
        "v_min_edge": v_min_edge,
        "v_max_edge": v_max_edge,
        "trials_per_box": next(iter(trials)),
    }
    return boxes, meta


def box_index(x, min_edge, step, n, tol=1e-9):
    # Values are expected inside [min_edge, min_edge+n*step].  Protect the
    # extremely unlikely exact upper endpoint by assigning it to the last bin.
    xmax = min_edge + n * step
    if x < min_edge - tol or x > xmax + tol:
        return None
    if x >= xmax:
        return n - 1
    i = int(math.floor((x - min_edge) / step))
    if i < 0 or i >= n:
        return None
    return i


def finite_mean(sum_value, count):
    return sum_value / count if count > 0 else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--transport-electrons", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    ap.add_argument("--wire-pitch-mm", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=181001)
    args = ap.parse_args()

    if args.wire_pitch_mm <= 0:
        raise SystemExit("--wire-pitch-mm must be positive.")

    rng = random.Random(args.seed)
    boxes, meta = load_transport(args.transport_electrons)

    print("PI fast response: transport map")
    print(
        f"  grid                 : {meta['n_u']} x {meta['n_v']} "
        f"(du={meta['du']:.6g} mm, dv={meta['dv']:.6g} mm)"
    )
    print(
        f"  u edges [mm]         : "
        f"{meta['u_min_edge']:.6g} .. {meta['u_max_edge']:.6g}"
    )
    print(
        f"  v edges [mm]         : "
        f"{meta['v_min_edge']:.6g} .. {meta['v_max_edge']:.6g}"
    )
    print(f"  trials / box         : {meta['trials_per_box']}")

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    event_path = Path(str(args.output_prefix) + "_event_summary.csv")
    wire_path = Path(str(args.output_prefix) + "_wire_summary.csv")

    event_fields = [
        "synthetic_track_id",
        "channel", "job", "event", "track", "pdg", "station",
        "p_GeV", "theta_u_deg", "theta_w_deg",
        "u0_cell_mm", "track_length_mm",
        "clusters", "primary_electrons",
        "transported_electrons", "collected_electrons",
        "collection_fraction", "active_wires",
        "dominant_wire_index", "dominant_wire_u_mm",
        "dominant_wire_fraction", "wire_cog_u_mm",
        "truth_electron_weighted_u_mm",
        "binary_residual_u0_mm", "binary_residual_charge_truth_mm",
        "wire_cog_residual_u0_mm", "wire_cog_residual_charge_truth_mm",
        "mean_drift_time_ns", "earliest_drift_time_ns",
        "dominant_wire_A_LR", "dominant_wire_mean_cos_phi",
        "dominant_wire_mean_incidence_deg",
        "out_of_map_electrons",
    ]
    wire_fields = [
        "synthetic_track_id", "channel", "station",
        "wire_index", "wire_u_mm",
        "collected_electrons", "collected_fraction",
        "A_LR", "mean_cos_phi", "mean_incidence_deg",
        "mean_drift_time_ns",
    ]

    total_tracks = 0
    total_primary = 0
    total_collected = 0
    total_out_of_map = 0

    with (
        args.clusters.open(newline="", encoding="utf-8") as fin,
        event_path.open("w", newline="", encoding="utf-8") as fevent,
        wire_path.open("w", newline="", encoding="utf-8") as fwire,
    ):
        reader = csv.DictReader(fin)
        ew = csv.DictWriter(fevent, fieldnames=event_fields)
        ww = csv.DictWriter(fwire, fieldnames=wire_fields)
        ew.writeheader()
        ww.writeheader()

        current_id = None
        state = None

        def new_state(row):
            return {
                "meta": {
                    "synthetic_track_id": as_int(row, "synthetic_track_id"),
                    "channel": row["channel"],
                    "job": row["job"],
                    "event": as_int(row, "event"),
                    "track": as_int(row, "track"),
                    "pdg": as_int(row, "pdg"),
                    "station": as_int(row, "station"),
                    "p_GeV": as_float(row, "p_GeV"),
                    "theta_u_deg": as_float(row, "theta_u_deg"),
                    "theta_w_deg": as_float(row, "theta_w_deg"),
                    "u0_cell_mm": as_float(row, "u0_cell_mm"),
                    "track_length_mm": as_float(row, "track_length_mm"),
                },
                "clusters": 0,
                "primary": 0,
                "transported": 0,
                "collected": 0,
                "out": 0,
                "sum_truth_u": 0.0,
                "wire": defaultdict(
                    lambda: {
                        "n": 0,
                        "side_plus": 0,
                        "side_minus": 0,
                        "sum_cos": 0.0,
                        "n_cos": 0,
                        "sum_inc": 0.0,
                        "n_inc": 0,
                        "sum_dt": 0.0,
                        "n_dt": 0,
                    }
                ),
                "sum_dt": 0.0,
                "n_dt": 0,
                "earliest_dt": math.inf,
            }

        def flush(s):
            nonlocal total_tracks, total_primary
            nonlocal total_collected, total_out_of_map
            if s is None:
                return

            m = s["meta"]
            primary = s["primary"]
            collected = s["collected"]
            wires = s["wire"]
            total_tracks += 1
            total_primary += primary
            total_collected += collected
            total_out_of_map += s["out"]

            truth_charge_u = (
                s["sum_truth_u"] / primary if primary > 0 else math.nan
            )

            if collected > 0 and wires:
                dominant = max(
                    sorted(wires),
                    key=lambda k: wires[k]["n"],
                )
                dominant_n = wires[dominant]["n"]
                dominant_fraction = dominant_n / collected
                dominant_u = dominant * args.wire_pitch_mm
                wire_cog = (
                    sum(
                        k * args.wire_pitch_mm * v["n"]
                        for k, v in wires.items()
                    )
                    / collected
                )
                wd = wires[dominant]
                lr_den = wd["side_plus"] + wd["side_minus"]
                dominant_a = (
                    (wd["side_plus"] - wd["side_minus"]) / lr_den
                    if lr_den > 0
                    else math.nan
                )
                dominant_cos = finite_mean(wd["sum_cos"], wd["n_cos"])
                dominant_inc = finite_mean(wd["sum_inc"], wd["n_inc"])
            else:
                dominant = 999999
                dominant_n = 0
                dominant_fraction = math.nan
                dominant_u = math.nan
                wire_cog = math.nan
                dominant_a = math.nan
                dominant_cos = math.nan
                dominant_inc = math.nan

            ew.writerow(
                {
                    **m,
                    "clusters": s["clusters"],
                    "primary_electrons": primary,
                    "transported_electrons": s["transported"],
                    "collected_electrons": collected,
                    "collection_fraction": (
                        collected / s["transported"]
                        if s["transported"] > 0
                        else math.nan
                    ),
                    "active_wires": len(wires),
                    "dominant_wire_index": dominant,
                    "dominant_wire_u_mm": dominant_u,
                    "dominant_wire_fraction": dominant_fraction,
                    "wire_cog_u_mm": wire_cog,
                    "truth_electron_weighted_u_mm": truth_charge_u,
                    "binary_residual_u0_mm": (
                        dominant_u - m["u0_cell_mm"]
                        if collected > 0
                        else math.nan
                    ),
                    "binary_residual_charge_truth_mm": (
                        dominant_u - truth_charge_u
                        if collected > 0 and math.isfinite(truth_charge_u)
                        else math.nan
                    ),
                    "wire_cog_residual_u0_mm": (
                        wire_cog - m["u0_cell_mm"]
                        if collected > 0
                        else math.nan
                    ),
                    "wire_cog_residual_charge_truth_mm": (
                        wire_cog - truth_charge_u
                        if collected > 0 and math.isfinite(truth_charge_u)
                        else math.nan
                    ),
                    "mean_drift_time_ns": finite_mean(
                        s["sum_dt"], s["n_dt"]
                    ),
                    "earliest_drift_time_ns": (
                        s["earliest_dt"]
                        if math.isfinite(s["earliest_dt"])
                        else math.nan
                    ),
                    "dominant_wire_A_LR": dominant_a,
                    "dominant_wire_mean_cos_phi": dominant_cos,
                    "dominant_wire_mean_incidence_deg": dominant_inc,
                    "out_of_map_electrons": s["out"],
                }
            )

            for k in sorted(wires):
                wd = wires[k]
                lr_den = wd["side_plus"] + wd["side_minus"]
                ww.writerow(
                    {
                        "synthetic_track_id": m["synthetic_track_id"],
                        "channel": m["channel"],
                        "station": m["station"],
                        "wire_index": k,
                        "wire_u_mm": k * args.wire_pitch_mm,
                        "collected_electrons": wd["n"],
                        "collected_fraction": wd["n"] / collected,
                        "A_LR": (
                            (wd["side_plus"] - wd["side_minus"]) / lr_den
                            if lr_den > 0
                            else math.nan
                        ),
                        "mean_cos_phi": finite_mean(
                            wd["sum_cos"], wd["n_cos"]
                        ),
                        "mean_incidence_deg": finite_mean(
                            wd["sum_inc"], wd["n_inc"]
                        ),
                        "mean_drift_time_ns": finite_mean(
                            wd["sum_dt"], wd["n_dt"]
                        ),
                    }
                )

        for row in reader:
            tid = as_int(row, "synthetic_track_id")
            if current_id is None or tid != current_id:
                flush(state)
                current_id = tid
                state = new_state(row)

            state["clusters"] += 1
            ne = as_int(row, "n_primary_electrons")
            state["primary"] += ne

            u_cell = as_float(row, "u_cell_mm")
            v_mm = as_float(row, "v_mm")
            wire_offset = as_int(row, "wire_offset")
            u_unfolded = as_float(row, "u_unfolded_mm")
            state["sum_truth_u"] += ne * u_unfolded

            iu0 = box_index(
                u_cell,
                meta["u_min_edge"],
                meta["du"],
                meta["n_u"],
            )
            iv0 = box_index(
                v_mm,
                meta["v_min_edge"],
                meta["dv"],
                meta["n_v"],
            )

            if iu0 is None or iv0 is None:
                state["out"] += ne
                continue

            iu = meta["iu_min"] + iu0
            iv = meta["iv_min"] + iv0
            trials = boxes[(iu, iv)]

            for _ in range(ne):
                state["transported"] += 1
                response = rng.choice(trials)
                if not response["collected"]:
                    continue

                wire_index = wire_offset + response["wire_index"]
                wd = state["wire"][wire_index]
                wd["n"] += 1
                state["collected"] += 1

                side = response["hit_side_u"]
                if side > 0:
                    wd["side_plus"] += 1
                elif side < 0:
                    wd["side_minus"] += 1

                cphi = response["cos_phi_hit"]
                if math.isfinite(cphi):
                    wd["sum_cos"] += cphi
                    wd["n_cos"] += 1

                inc = response["incidence_normal_deg"]
                if math.isfinite(inc):
                    wd["sum_inc"] += inc
                    wd["n_inc"] += 1

                dt = response["drift_time_ns"]
                if math.isfinite(dt) and dt >= 0:
                    wd["sum_dt"] += dt
                    wd["n_dt"] += 1
                    state["sum_dt"] += dt
                    state["n_dt"] += 1
                    state["earliest_dt"] = min(state["earliest_dt"], dt)

        flush(state)

    print()
    print("PI fast response COMPLETE")
    print(f"  tracks                 : {total_tracks}")
    print(f"  primary electrons      : {total_primary}")
    print(f"  collected electrons    : {total_collected}")
    print(
        f"  collection fraction    : "
        f"{total_collected / total_primary:.6f}"
        if total_primary
        else "  collection fraction    : nan"
    )
    print(f"  outside transport map  : {total_out_of_map}")
    print(f"  event summary           : {event_path}")
    print(f"  wire summary            : {wire_path}")


if __name__ == "__main__":
    main()

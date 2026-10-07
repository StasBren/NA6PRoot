#!/usr/bin/env python3
"""Stage C0: accepted local muon-incidence phase space.

Reads the persistent high-statistics MNP33-bore Geant4 production.  No Geant4
transport is rerun.

The nominal MWPC coordinate convention is the one already implemented and
unit-tested in ../mwpc-avalanche/MWPCCoordinateFrame.h:

  local u = detector-global Y  (across anode wires)
  local v = detector-global Z  (chamber normal)
  local w = detector-global X  (along anode wires)

Therefore, for the nominal station orientation,

  theta_u = atan2(p_u, p_v) = atan2(pY, pZ)
  theta_w = atan2(p_w, p_v) = atan2(pX, pZ)

The large downstream magnetic-bending tail is expected mainly in theta_w.
"""

from __future__ import annotations

import argparse
from array import array
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable

import ROOT as R

R.gROOT.SetBatch(True)
for _lib in ("libutilsLib", "libbaseLib", "libsimLib"):
    if R.gSystem.Load(_lib) < 0:
        raise RuntimeError(f"Cannot load {_lib}")

HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "stage_c_reference.json"
N_STUDY_STATIONS = 4
SELECTIONS = (
    "all",
    "geometry_bore",
    "useful_bore",
    "useful_bore_working_area",
)


class Bucket:
    def __init__(self) -> None:
        self.theta_u = array("d")
        self.theta_w = array("d")
        self.theta_3d = array("d")
        self.p = array("d")

    def fill(self, state: dict[str, float]) -> None:
        # Reuse the already established nominal chamber-frame mapping:
        #   p_u = pY, p_v = pZ, p_w = pX.
        p_u = state["py"]
        p_v = state["pz"]
        p_w = state["px"]
        self.theta_u.append(math.degrees(math.atan2(p_u, p_v)))
        self.theta_w.append(math.degrees(math.atan2(p_w, p_v)))
        self.theta_3d.append(
            math.degrees(math.atan2(math.hypot(p_u, p_w), p_v))
        )
        self.p.append(state["p"])

    def __len__(self) -> int:
        return len(self.theta_3d)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def request_for_run(run: Path) -> dict:
    for name in ("request.json", "input_request.json"):
        path = run / name
        if path.exists():
            return read_json(path)
    raise FileNotFoundError(f"No request.json or input_request.json in {run}")


def direct_daughter_ids(mckine, event: int, parent_pdg: int) -> set[int]:
    mckine.GetEntry(event)
    pdg = [int(t.GetPdgCode()) for t in mckine.tracks]
    mother = [int(t.GetFirstMother()) for t in mckine.tracks]
    return {
        i
        for i, p in enumerate(pdg)
        if abs(p) == 13
        and 0 <= mother[i] < len(pdg)
        and pdg[mother[i]] == int(parent_pdg)
    }


def quantile(values: Iterable[float], q: float, absolute: bool = False) -> float:
    data = [abs(float(v)) if absolute else float(v) for v in values]
    if not data:
        return float("nan")
    data.sort()
    if len(data) == 1:
        return data[0]
    pos = q * (len(data) - 1)
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    if lo == hi:
        return data[lo]
    f = pos - lo
    return data[lo] * (1.0 - f) + data[hi] * f


def mean_rms(values: Iterable[float]) -> tuple[float, float]:
    data = [float(v) for v in values]
    if not data:
        return float("nan"), float("nan")
    mean = sum(data) / len(data)
    rms = math.sqrt(
        max(0.0, sum(x * x for x in data) / len(data) - mean * mean)
    )
    return mean, rms


def fraction_abs_above(values: Iterable[float], threshold: float) -> float:
    data = [abs(float(v)) for v in values]
    if not data:
        return float("nan")
    return sum(v > threshold for v in data) / len(data)


def inside_bore(
    states: dict[float, dict[str, float]],
    bore_z_cm: tuple[float, ...],
    half_x_cm: float,
    half_y_cm: float,
) -> bool:
    for z in bore_z_cm:
        s = states.get(z)
        if s is None:
            return False
        if abs(s["x"]) > half_x_cm or abs(s["y"]) > half_y_cm:
            return False
    return True


def fill_muon(
    states: dict[float, dict[str, float]],
    nominal_z: list[float],
    buckets: dict[str, list[Bucket]],
    bore_z_cm: tuple[float, ...],
    geometry_half_x: float,
    geometry_half_y: float,
    useful_half_x: float,
    useful_half_y: float,
    working_area_cm: tuple[tuple[float, float], ...],
    counters: dict[str, int],
) -> None:
    if any(z in states for z in nominal_z[:N_STUDY_STATIONS]):
        counters["with_ms0_ms3_state"] += 1

    geometry_ok = inside_bore(
        states, bore_z_cm, geometry_half_x, geometry_half_y
    )
    useful_ok = inside_bore(
        states, bore_z_cm, useful_half_x, useful_half_y
    )

    if geometry_ok:
        counters["geometry_bore_muons"] += 1
    if useful_ok:
        counters["useful_bore_muons"] += 1

    for station, z in enumerate(nominal_z[:N_STUDY_STATIONS]):
        state = states.get(z)
        if state is None:
            continue

        working_x, working_y = working_area_cm[station]
        inside_working_area = (
            abs(state["x"]) <= 0.5 * working_x
            and abs(state["y"]) <= 0.5 * working_y
        )

        selected = {
            "all": True,
            "geometry_bore": geometry_ok,
            "useful_bore": useful_ok,
            "useful_bore_working_area": useful_ok and inside_working_area,
        }
        for selection in SELECTIONS:
            if selected[selection]:
                buckets[selection][station].fill(state)


def process_job(
    run: Path,
    job: dict,
    nominal_z: list[float],
    buckets: dict[str, list[Bucket]],
    bore_z_cm: tuple[float, ...],
    geometry_half_x: float,
    geometry_half_y: float,
    useful_half_x: float,
    useful_half_y: float,
    working_area_cm: tuple[tuple[float, float], ...],
    counters: dict[str, int],
) -> None:
    request = request_for_run(run)
    parent_pdg = int(request.get("parent_pdg", job["parent_pdg"]))
    plane_z = [float(z) for z in request["reference_z_cm"]]

    audit_path = run / "LocalTrajectoryAudit.root"
    kine_path = run / "MCKine.root"
    tf = R.TFile.Open(str(audit_path), "READ")
    kf = R.TFile.Open(str(kine_path), "READ")
    if not tf or tf.IsZombie() or not kf or kf.IsZombie():
        raise RuntimeError(f"Cannot open ROOT files in {run}")

    tree = tf.Get("Crossings")
    mckine = kf.Get("mckine")
    if not tree or not mckine:
        raise RuntimeError(f"Missing Crossings/mckine tree in {run}")

    counters["jobs"] += 1
    counters["generated_events"] += int(
        request.get("events", job.get("events", 0))
    )

    for event in range(int(mckine.GetEntries())):
        counters["generated_direct_muons"] += len(
            direct_daughter_ids(mckine, event, parent_pdg)
        )

    current_event = None
    current_ids: set[int] = set()
    event_states: Dict[int, dict[float, dict[str, float]]] = {}

    def flush() -> None:
        if current_event is None:
            return
        for track in current_ids:
            fill_muon(
                event_states.get(track, {}),
                nominal_z,
                buckets,
                bore_z_cm,
                geometry_half_x,
                geometry_half_y,
                useful_half_x,
                useful_half_y,
                working_area_cm,
                counters,
            )

    for row in tree:
        event = int(row.event)
        if current_event is None or event != current_event:
            flush()
            current_event = event
            current_ids = direct_daughter_ids(
                mckine, event, parent_pdg
            )
            event_states = {track: {} for track in current_ids}

        track = int(row.track)
        if track not in current_ids:
            continue

        plane = int(row.plane)
        if plane < 0 or plane >= len(plane_z):
            continue

        z = plane_z[plane]
        event_states[track][z] = {
            "x": float(row.x_cm),
            "y": float(row.y_cm),
            "px": float(row.px_GeV),
            "py": float(row.py_GeV),
            "pz": float(row.pz_GeV),
            "p": float(row.p_GeV),
        }

    flush()
    tf.Close()
    kf.Close()


def stats_row(
    source: str, station: int, selection: str, bucket: Bucket
) -> dict:
    tu_mean, tu_rms = mean_rms(bucket.theta_u)
    tw_mean, tw_rms = mean_rms(bucket.theta_w)
    t3_mean, t3_rms = mean_rms(bucket.theta_3d)
    p_mean, p_rms = mean_rms(bucket.p)

    return {
        "source": source,
        "station": station,
        "selection": selection,
        "N": len(bucket),
        "mean_theta_u_deg": tu_mean,
        "rms_theta_u_deg": tu_rms,
        "q50_abs_theta_u_deg": quantile(bucket.theta_u, 0.50, True),
        "q95_abs_theta_u_deg": quantile(bucket.theta_u, 0.95, True),
        "q99_abs_theta_u_deg": quantile(bucket.theta_u, 0.99, True),
        "mean_theta_w_deg": tw_mean,
        "rms_theta_w_deg": tw_rms,
        "q50_abs_theta_w_deg": quantile(bucket.theta_w, 0.50, True),
        "q95_abs_theta_w_deg": quantile(bucket.theta_w, 0.95, True),
        "q99_abs_theta_w_deg": quantile(bucket.theta_w, 0.99, True),
        "frac_abs_theta_w_gt20": fraction_abs_above(bucket.theta_w, 20.0),
        "frac_abs_theta_w_gt30": fraction_abs_above(bucket.theta_w, 30.0),
        "frac_abs_theta_w_gt40": fraction_abs_above(bucket.theta_w, 40.0),
        "frac_abs_theta_w_gt50": fraction_abs_above(bucket.theta_w, 50.0),
        "mean_theta_3d_deg": t3_mean,
        "rms_theta_3d_deg": t3_rms,
        "q50_theta_3d_deg": quantile(bucket.theta_3d, 0.50),
        "q95_theta_3d_deg": quantile(bucket.theta_3d, 0.95),
        "q99_theta_3d_deg": quantile(bucket.theta_3d, 0.99),
        "mean_p_GeV": p_mean,
        "rms_p_GeV": p_rms,
        "q05_p_GeV": quantile(bucket.p, 0.05),
        "q50_p_GeV": quantile(bucket.p, 0.50),
        "q95_p_GeV": quantile(bucket.p, 0.95),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "suite", type=Path, help="Saved high-statistics MNP33-bore suite"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Stage-C reference configuration",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("test_runs/mwpc_stage_c/c0_angles"),
        help="Output directory; the historical suite is never modified",
    )
    args = parser.parse_args()

    config = read_json(args.config)
    suite = args.suite.expanduser().resolve()
    manifest = read_json(suite / "manifest.json")

    bore_z_cm = tuple(
        float(z) for z in config["mnp33"]["reference_z_cm"]
    )
    useful_full_x, useful_full_y = (
        float(x)
        for x in config["mnp33"]["useful_aperture_full_cm_xy"]
    )
    useful_half_x = 0.5 * useful_full_x
    useful_half_y = 0.5 * useful_full_y

    working_area_cm = tuple(
        tuple(float(x) for x in config["working_area_cm_xy"][f"MS{i}"])
        for i in range(N_STUDY_STATIONS)
    )

    nominal_z = [float(z) for z in manifest["nominal_station_z_cm"]]
    if len(nominal_z) < N_STUDY_STATIONS:
        raise RuntimeError("Manifest has fewer than four nominal MS planes")

    ref_z = {float(z) for z in manifest["reference_z_cm"]}
    if not set(bore_z_cm).issubset(ref_z):
        raise RuntimeError(
            f"Suite does not contain all configured MNP33 planes {bore_z_cm}"
        )

    bore_size = manifest.get("analytic_mnp33_bore_cm", [320.0, 240.0])
    geometry_half_x = 0.5 * float(bore_size[0])
    geometry_half_y = 0.5 * float(bore_size[1])

    jobs_by_sample: dict[str, list[dict]] = {}
    for job in manifest["jobs"]:
        jobs_by_sample.setdefault(job["sample"], []).append(job)

    args.output.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    source_counts: dict[str, dict[str, int]] = {}

    print(f"Suite: {suite}")
    print(f"Manifest mode: {manifest.get('mode', 'unknown')}")
    print(f"Events/source: {manifest.get('events_per_source', 'unknown')}")
    print(
        "Local frame: u=global Y (across wires), "
        "v=global Z (normal), w=global X (along wires)"
    )
    print(
        "Fiducials: geometry "
        f"|x|<={geometry_half_x:g} cm, |y|<={geometry_half_y:g} cm; "
        f"useful |x|<={useful_half_x:g} cm, |y|<={useful_half_y:g} cm"
    )
    print(
        "Working areas MS0..MS3 [global x*y cm]: "
        + ", ".join(f"{x:g}x{y:g}" for x, y in working_area_cm)
    )

    for source, jobs in sorted(jobs_by_sample.items()):
        buckets = {
            selection: [Bucket() for _ in range(N_STUDY_STATIONS)]
            for selection in SELECTIONS
        }
        counters = {
            "jobs": 0,
            "generated_events": 0,
            "generated_direct_muons": 0,
            "with_ms0_ms3_state": 0,
            "geometry_bore_muons": 0,
            "useful_bore_muons": 0,
        }

        print(f"\nProcessing {source}: {len(jobs)} chunk(s)")
        for index, job in enumerate(
            sorted(jobs, key=lambda x: x["name"]), start=1
        ):
            run = suite / job["name"]
            print(f"  [{index}/{len(jobs)}] {run.name}", flush=True)
            process_job(
                run,
                job,
                nominal_z,
                buckets,
                bore_z_cm,
                geometry_half_x,
                geometry_half_y,
                useful_half_x,
                useful_half_y,
                working_area_cm,
                counters,
            )

        source_counts[source] = counters
        print(
            f"  generated direct muons={counters['generated_direct_muons']}, "
            f"geometry-bore={counters['geometry_bore_muons']}, "
            f"useful-bore={counters['useful_bore_muons']}"
        )
        print(
            "  station selection                    N      "
            "q95|theta_u| q99|theta_u|  "
            "q95|theta_w| q99|theta_w|  "
            "f(|theta_w|>30) f(>40) f(>50)"
        )

        for station in range(N_STUDY_STATIONS):
            for selection in SELECTIONS:
                row = stats_row(
                    source, station, selection,
                    buckets[selection][station]
                )
                all_rows.append(row)
                print(
                    f"  MS{station}    {selection:<27} "
                    f"{row['N']:8d}  "
                    f"{row['q95_abs_theta_u_deg']:12.3f} "
                    f"{row['q99_abs_theta_u_deg']:12.3f}  "
                    f"{row['q95_abs_theta_w_deg']:12.3f} "
                    f"{row['q99_abs_theta_w_deg']:12.3f}  "
                    f"{row['frac_abs_theta_w_gt30']:15.6f} "
                    f"{row['frac_abs_theta_w_gt40']:7.6f} "
                    f"{row['frac_abs_theta_w_gt50']:7.6f}"
                )

    summary_path = args.output / "incident_angle_summary.csv"
    if all_rows:
        with summary_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=list(all_rows[0].keys())
            )
            writer.writeheader()
            writer.writerows(all_rows)

    metadata = {
        "input_suite": str(suite),
        "config": str(args.config.resolve()),
        "manifest_mode": manifest.get("mode"),
        "events_per_source": manifest.get("events_per_source"),
        "nominal_station_z_cm": nominal_z,
        "mnp33_reference_z_cm": list(bore_z_cm),
        "geometry_bore_half_width_cm": [
            geometry_half_x, geometry_half_y
        ],
        "useful_bore_half_width_cm": [
            useful_half_x, useful_half_y
        ],
        "working_area_cm_xy": [
            list(x) for x in working_area_cm
        ],
        "coordinate_mapping": config["coordinate_frame"],
        "selection_scope": (
            "single direct daughter muon; no pair-level requirement"
        ),
        "source_counts": source_counts,
    }
    metadata_path = args.output / "metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    print(f"\nCSV:  {summary_path}")
    print(f"Meta: {metadata_path}")


if __name__ == "__main__":
    main()

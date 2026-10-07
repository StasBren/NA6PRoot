#!/usr/bin/env python3
"""Stage C0 incident-angle analysis from the saved high-statistics MNP33-bore suite.

This deliberately reuses persistent Geant4 trajectory data.  No transport is
rerun.  The analysis is single-muon based because Garfield response is evaluated
per charged-particle crossing: each direct daughter muon is classified by whether
its own transported trajectory crosses the three MNP33 reference planes inside

  geometry bore : half-widths from manifest (normally 160 x 120 cm)
  useful bore   : 122.5 x 120 cm (245 x 240 cm proposal useful aperture)

A fourth selection, useful_bore_working_area, additionally requires the muon
crossing position at the station plane to lie inside the same approximate
finite working rectangle already used by the offline stagger study.  This
directly tests whether the large downstream angle tails are still relevant to
finite station coverage.

Angles at MS0..MS3 are taken directly from the transported momentum stored at the
nominal reference plane:
  theta_x = atan2(px, pz)
  theta_y = atan2(py, pz)
  theta   = atan2(sqrt(px^2 + py^2), pz)
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

BORE_Z_CM = (365.0, 430.0, 495.0)
USEFUL_HALF_X_CM = 122.5
USEFUL_HALF_Y_CM = 120.0
N_STUDY_STATIONS = 4

# Approximate externally specified station working rectangles used by the
# existing offline stagger study, in detector-global x,y.  These are a
# physics/layout benchmark rather than a frozen construction envelope.
WORKING_AREA_CM = (
    (220.0, 220.0),  # MS0
    (230.0, 240.0),  # MS1
    (310.0, 310.0),  # MS2
    (320.0, 320.0),  # MS3
)

SELECTIONS = (
    "all",
    "geometry_bore",
    "useful_bore",
    "useful_bore_working_area",
)


class Bucket:
    def __init__(self) -> None:
        self.theta_x = array("d")
        self.theta_y = array("d")
        self.theta = array("d")
        self.p = array("d")

    def fill(self, state: dict[str, float]) -> None:
        px, py, pz = state["px"], state["py"], state["pz"]
        self.theta_x.append(math.degrees(math.atan2(px, pz)))
        self.theta_y.append(math.degrees(math.atan2(py, pz)))
        self.theta.append(math.degrees(math.atan2(math.hypot(px, py), pz)))
        self.p.append(state["p"])

    def __len__(self) -> int:
        return len(self.theta)


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
    n = 0
    s = 0.0
    s2 = 0.0
    for raw in values:
        x = float(raw)
        n += 1
        s += x
        s2 += x * x
    if n == 0:
        return float("nan"), float("nan")
    mean = s / n
    rms = math.sqrt(max(0.0, s2 / n - mean * mean))
    return mean, rms


def inside_bore(states: dict[float, dict[str, float]], half_x: float, half_y: float) -> bool:
    for z in BORE_Z_CM:
        s = states.get(z)
        if s is None:
            return False
        if abs(s["x"]) > half_x or abs(s["y"]) > half_y:
            return False
    return True


def fill_muon(
    states: dict[float, dict[str, float]],
    nominal_z: list[float],
    buckets: dict[str, list[Bucket]],
    geometry_half_x: float,
    geometry_half_y: float,
    counters: dict[str, int],
) -> None:
    has_any_station = any(z in states for z in nominal_z[:N_STUDY_STATIONS])
    if has_any_station:
        counters["with_ms0_ms3_state"] += 1

    geometry_ok = inside_bore(states, geometry_half_x, geometry_half_y)
    useful_ok = inside_bore(states, USEFUL_HALF_X_CM, USEFUL_HALF_Y_CM)
    if geometry_ok:
        counters["geometry_bore_muons"] += 1
    if useful_ok:
        counters["useful_bore_muons"] += 1

    for station, z in enumerate(nominal_z[:N_STUDY_STATIONS]):
        state = states.get(z)
        if state is None:
            continue

        working_x, working_y = WORKING_AREA_CM[station]
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
        for sel in SELECTIONS:
            if selected[sel]:
                buckets[sel][station].fill(state)


def process_job(
    run: Path,
    job: dict,
    nominal_z: list[float],
    buckets: dict[str, list[Bucket]],
    geometry_half_x: float,
    geometry_half_y: float,
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
    counters["generated_events"] += int(request.get("events", job.get("events", 0)))

    # Count generated direct daughters exactly from MCKine.
    for event in range(int(mckine.GetEntries())):
        counters["generated_direct_muons"] += len(direct_daughter_ids(mckine, event, parent_pdg))

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
                geometry_half_x,
                geometry_half_y,
                counters,
            )

    for row in tree:
        event = int(row.event)
        if current_event is None or event != current_event:
            flush()
            current_event = event
            current_ids = direct_daughter_ids(mckine, event, parent_pdg)
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


def fraction_abs_above(values: Iterable[float], threshold: float) -> float:
    data = [abs(float(v)) for v in values]
    if not data:
        return float("nan")
    return sum(v > threshold for v in data) / len(data)


def stats_row(source: str, station: int, selection: str, b: Bucket) -> dict:
    tx_mean, tx_rms = mean_rms(b.theta_x)
    ty_mean, ty_rms = mean_rms(b.theta_y)
    t_mean, t_rms = mean_rms(b.theta)
    p_mean, p_rms = mean_rms(b.p)
    return {
        "source": source,
        "station": station,
        "selection": selection,
        "N": len(b),
        "mean_theta_x_deg": tx_mean,
        "rms_theta_x_deg": tx_rms,
        "q50_abs_theta_x_deg": quantile(b.theta_x, 0.50, True),
        "q95_abs_theta_x_deg": quantile(b.theta_x, 0.95, True),
        "q99_abs_theta_x_deg": quantile(b.theta_x, 0.99, True),
        "frac_abs_theta_x_gt20": fraction_abs_above(b.theta_x, 20.0),
        "frac_abs_theta_x_gt30": fraction_abs_above(b.theta_x, 30.0),
        "frac_abs_theta_x_gt40": fraction_abs_above(b.theta_x, 40.0),
        "frac_abs_theta_x_gt50": fraction_abs_above(b.theta_x, 50.0),
        "mean_theta_y_deg": ty_mean,
        "rms_theta_y_deg": ty_rms,
        "q50_abs_theta_y_deg": quantile(b.theta_y, 0.50, True),
        "q95_abs_theta_y_deg": quantile(b.theta_y, 0.95, True),
        "q99_abs_theta_y_deg": quantile(b.theta_y, 0.99, True),
        "mean_theta_deg": t_mean,
        "rms_theta_deg": t_rms,
        "q50_theta_deg": quantile(b.theta, 0.50),
        "q95_theta_deg": quantile(b.theta, 0.95),
        "q99_theta_deg": quantile(b.theta, 0.99),
        "mean_p_GeV": p_mean,
        "rms_p_GeV": p_rms,
        "q50_p_GeV": quantile(b.p, 0.50),
        "q05_p_GeV": quantile(b.p, 0.05),
        "q95_p_GeV": quantile(b.p, 0.95),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("suite", type=Path, help="Saved mnp33-bore production suite")
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("test_runs/mwpc_incident_angles_highstat"),
        help="Output directory (historical suite is never modified)",
    )
    args = ap.parse_args()

    suite = args.suite.expanduser().resolve()
    manifest = read_json(suite / "manifest.json")
    nominal_z = [float(z) for z in manifest["nominal_station_z_cm"]]
    if len(nominal_z) < N_STUDY_STATIONS:
        raise RuntimeError("Manifest has fewer than four nominal MS planes")
    ref_z = {float(z) for z in manifest["reference_z_cm"]}
    if not set(BORE_Z_CM).issubset(ref_z):
        raise RuntimeError(f"Suite does not contain all MNP33 planes {BORE_Z_CM}")

    bore_size = manifest.get("analytic_mnp33_bore_cm", [320.0, 240.0])
    geometry_half_x = float(bore_size[0]) / 2.0
    geometry_half_y = float(bore_size[1]) / 2.0

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
        "Fiducials: geometry "
        f"|x|<={geometry_half_x:g} cm, |y|<={geometry_half_y:g} cm; "
        f"useful |x|<={USEFUL_HALF_X_CM:g} cm, |y|<={USEFUL_HALF_Y_CM:g} cm"
    )
    print(
        "Working areas MS0..MS3 [x*y cm]: "
        + ", ".join(f"{x:g}x{y:g}" for x, y in WORKING_AREA_CM)
    )

    for source, jobs in sorted(jobs_by_sample.items()):
        buckets = {sel: [Bucket() for _ in range(N_STUDY_STATIONS)] for sel in SELECTIONS}
        counters = {
            "jobs": 0,
            "generated_events": 0,
            "generated_direct_muons": 0,
            "with_ms0_ms3_state": 0,
            "geometry_bore_muons": 0,
            "useful_bore_muons": 0,
        }

        print(f"\nProcessing {source}: {len(jobs)} chunk(s)")
        for index, job in enumerate(sorted(jobs, key=lambda x: x["name"]), start=1):
            run = suite / job["name"]
            print(f"  [{index}/{len(jobs)}] {run.name}", flush=True)
            process_job(
                run,
                job,
                nominal_z,
                buckets,
                geometry_half_x,
                geometry_half_y,
                counters,
            )

        source_counts[source] = counters
        print(
            f"  generated direct muons={counters['generated_direct_muons']}, "
            f"geometry-bore={counters['geometry_bore_muons']}, "
            f"useful-bore={counters['useful_bore_muons']}"
        )
        print(
            "  station sel                 N      "
            "q50|tx| q95|tx| q99|tx|   q50|ty| q95|ty| q99|ty|   "
            "q95(theta)  p50"
        )

        for station in range(N_STUDY_STATIONS):
            for sel in SELECTIONS:
                row = stats_row(source, station, sel, buckets[sel][station])
                all_rows.append(row)
                print(
                    f"  MS{station:<1}    {sel:<15} "
                    f"{row['N']:8d}  "
                    f"{row['q50_abs_theta_x_deg']:8.3f} "
                    f"{row['q95_abs_theta_x_deg']:8.3f} "
                    f"{row['q99_abs_theta_x_deg']:8.3f}   "
                    f"{row['q50_abs_theta_y_deg']:8.3f} "
                    f"{row['q95_abs_theta_y_deg']:8.3f} "
                    f"{row['q99_abs_theta_y_deg']:8.3f}   "
                    f"{row['q95_theta_deg']:10.3f} "
                    f"{row['q50_p_GeV']:7.3f}"
                )

    summary_path = args.output / "incident_angle_summary.csv"
    if all_rows:
        with summary_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)

    metadata = {
        "input_suite": str(suite),
        "manifest_mode": manifest.get("mode"),
        "events_per_source": manifest.get("events_per_source"),
        "nominal_station_z_cm": nominal_z,
        "bore_z_cm": list(BORE_Z_CM),
        "geometry_bore_half_width_cm": [geometry_half_x, geometry_half_y],
        "useful_bore_half_width_cm": [USEFUL_HALF_X_CM, USEFUL_HALF_Y_CM],
        "selection_scope": "single direct daughter muon; no pair-level requirement",
        "working_area_cm_xy": [list(x) for x in WORKING_AREA_CM],
        "working_area_note": (
            "Approximate external working rectangles from the existing stagger study; "
            "not frozen construction envelopes."
        ),
        "source_counts": source_counts,
    }
    metadata_path = args.output / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print(f"\nCSV:  {summary_path}")
    print(f"Meta: {metadata_path}")


if __name__ == "__main__":
    main()

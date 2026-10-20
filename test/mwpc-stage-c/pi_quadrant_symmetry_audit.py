#!/usr/bin/env python3
"""Read-only GEANT4 quadrant preflight. No claim of exact fourfold symmetry.

Uses frozen full-campaign module hit counts and signed-angle histograms.
Segment counts are correlated within events; differences and TVDs are
DESCRIPTIVE, not significance tests or evidence that an event can be reflected.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from mwpc_checkerboard_projection import SOURCES, build_chambers, NX, NY
from pi_native_full_campaign_qa import ANGLE

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = REPO / "test_runs/mwpc_stage_c/pi_native_full_quality"
OUT = REPO / "test_runs/mwpc_stage_c/pi_quadrant_preflight"
FIELDS = ("theta_u_in_deg", "theta_w_in_deg")
N_BINS = len(ANGLE) - 1
QUADS = ("Q1", "Q2", "Q3", "Q4")


def quadrant(c):
    if abs(c.x) < 1.e-9 or abs(c.y) < 1.e-9:
        return "AXIS"
    if c.x > 0:
        return "Q1" if c.y > 0 else "Q4"
    return "Q2" if c.y > 0 else "Q3"


def geometry_summary():
    results = []
    geometry = build_chambers()
    for st, group in enumerate(geometry):
        coord = {(c.row, c.col): c for c in group}
        orbits = 0
        same_z = 0
        for c in group:
            if quadrant(c) != "Q1":
                continue
            mirrors = [
                coord[(c.row, c.col)],
                coord[(c.row, NX[st] - 1 - c.col)],
                coord[(NY[st] - 1 - c.row, NX[st] - 1 - c.col)],
                coord[(NY[st] - 1 - c.row, c.col)],
            ]
            assert [quadrant(v) for v in mirrors] == list(QUADS)
            orbits += 1
            same_z += all(v.z_offset == c.z_offset and v.layer == c.layer
                          for v in mirrors)
        results.append({
            "station": st, "modules": len(group),
            "modules_on_axes": sum(quadrant(c) == "AXIS" for c in group),
            "four_module_orbits": orbits,
            "orbits_identical_layer_and_z": same_z,
            "exact_quadrant_chamber_geometry":
                orbits > 0 and same_z == orbits,
        })
    return results


def total_variation(first, other):
    a, b = sum(first), sum(other)
    if not a or not b:
        return None
    return 0.5 * sum(abs(x / a - y / b) for x, y in zip(first, other))


def require_campaign(data):
    report = json.loads((data / "pi_native_all_QA_report.json")
                        .read_text(encoding="utf-8"))
    if not report.get("full_campaign_analyzed") or \
       report.get("processed_chunks_per_source") != 67:
        raise ValueError("Quadrant study requires validated full 201-chunk QA")


def analyse(data):
    require_campaign(data)
    geometry = {c.det_id: c for group in build_chambers() for c in group}
    exposures = {}
    with (data / "pi_native_all_chamber_quality.csv").open(
            newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = row["source"], int(row["det_id"])
            if key[0] not in SOURCES or key[1] not in geometry or key in exposures:
                raise ValueError(f"Invalid/duplicated QA chamber {key}")
            if int(row["station"]) != geometry[key[1]].station:
                raise ValueError("Inconsistent station in QA chamber")
            n = int(row["gas_hit_rows"])
            if n < 0:
                raise ValueError("Negative gas-hit exposure")
            exposures[key] = n
    if len(exposures) != len(SOURCES) * len(geometry):
        raise ValueError("Expected all 894 source/module QA counts")
    totals = defaultdict(int)
    nmodules = defaultdict(int)
    for (source, det), n in exposures.items():
        c = geometry[det]
        key = source, c.station, quadrant(c)
        totals[key] += n
        nmodules[key] += 1

    # Full signed-angle PDF bin edges are symmetric around zero. Mirroring
    # global Y flips the local-u angle; mirroring X flips local-w angle.
    hist = defaultdict(lambda: [0] * N_BINS)
    seen = set()
    with (data / "pi_native_all_histograms.csv").open(
            newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            field = row["variable"]
            if field not in FIELDS:
                continue
            source, det = row["source"], int(row["det_id"])
            if (source, det) not in exposures:
                raise ValueError("Angular PDF refers to unknown module")
            c = geometry[det]
            if int(row["station"]) != c.station or \
               int(row["gas_hit_exposure"]) != exposures[source, det]:
                raise ValueError("Angular PDF exposure/station mismatch")
            lo, hi = float(row["bin_lower"]), float(row["bin_upper"])
            idx = next((i for i in range(N_BINS)
                        if abs(ANGLE[i] - lo) < 1.e-8 and
                           abs(ANGLE[i+1] - hi) < 1.e-8), None)
            if idx is None:
                raise ValueError("Unknown exact signed-angle bin")
            unique = (source, det, field, idx)
            if unique in seen:
                raise ValueError("Duplicate signed-angle QA bin")
            seen.add(unique)
            q = quadrant(c)
            # Bring the marginal distributions into Q1 sign conventions.
            flip = (field == "theta_u_in_deg" and c.y < 0) or \
                   (field == "theta_w_in_deg" and c.x < 0)
            dest = N_BINS - 1 - idx if flip else idx
            count = int(row["hit_count"])
            if count < 0:
                raise ValueError("Negative angle-bin count")
            hist[source, c.station, q, field][dest] += count
    if len(seen) != len(exposures) * len(FIELDS) * N_BINS:
        raise ValueError("Incomplete 894 x 2 signed angle QA grid")

    rows = []
    for source in SOURCES:
        for st in range(6):
            for q in QUADS + ("AXIS",):
                n = totals[source, st, q]
                modules = nmodules[source, st, q]
                if q in QUADS and not modules:
                    raise ValueError("Missing quadrants")
                for field in FIELDS:
                    if sum(hist[source, st, q, field]) != n:
                        raise ValueError("Angle/segment histogram conservation failed")
                ref = totals[source, st, "Q1"]
                ref_modules = nmodules[source, st, "Q1"]
                rate = n / modules if modules else None
                ref_rate = ref / ref_modules if ref_modules else None
                rows.append({
                    "source": source, "station": st, "quadrant": q,
                    "chambers": modules, "gas_hit_segments": n,
                    "segments_per_chamber": rate,
                    "ratio_to_Q1_per_chamber":
                        rate / ref_rate if rate is not None and ref_rate else None,
                    "theta_u_marginal_TVD_to_Q1":
                        total_variation(hist[source, st, "Q1", FIELDS[0]],
                                        hist[source, st, q, FIELDS[0]]),
                    "theta_w_marginal_TVD_to_Q1":
                        total_variation(hist[source, st, "Q1", FIELDS[1]],
                                        hist[source, st, q, FIELDS[1]]),
                })
            if sum(totals[source, st, q] for q in QUADS + ("AXIS",)) != \
               sum(exposures[source, c.det_id] for c in geometry.values()
                   if c.station == st):
                raise ValueError("Quadrants do not conserve station segment counts")
    return rows


def run(data, out, geometry_only=False):
    out.mkdir(parents=True, exist_ok=True)
    geom = geometry_summary()
    summary = {
        "geometry": geom,
        "exact_four_quadrant_geometry": all(
            row["exact_quadrant_chamber_geometry"] for row in geom),
        "interpretation": (
            "Different ABCD layer/z placement prevents exact four-quadrant "
            "mapping of the complete installed chambers; full spectrometer "
            "materials, fields, beam and dimuon pair correlations also need "
            "validation before any quadrant-only generator is allowed."),
        "statistics": (
            "Geant4 gas-hit SEGMENTS, not generated muons or independent "
            "observations. Counts and signed marginal TVD are descriptive "
            "only. No sqrt(N) inference: segment and event correlations."),
        "orientation": (
            "Global X reflection reverses signed theta_w; Y reflection "
            "reverses signed theta_u, under the current TGeo-local conventions."),
        "full_native_QA_checked": not geometry_only,
    }
    rows = [] if geometry_only else analyse(data)
    if rows:
        with (out / "pi_quadrant_native_counts.csv").open(
                "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (out / "pi_quadrant_geometry_preflight.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    for r in geom:
        print(f'MS{r["station"]}: exact fourfold layer/z orbits '
              f'{r["orbits_identical_layer_and_z"]}/{r["four_module_orbits"]}; '
              f'axis modules={r["modules_on_axes"]}')
    if rows:
        print(f"PASS: {len(rows)} source/station/quadrant QA summaries, "
              "including signed entrance-angle marginal comparisons.")
        print("Per-chamber hit-density ratios and signed-angle TVD to Q1 "
              "(DESCRIPTIVE only; no significance test):")
        for source in SOURCES:
            for st in range(6):
                chosen = {row["quadrant"]:row for row in rows
                          if row["source"] == source and row["station"] == st}
                def val(q, name):
                    number = chosen[q][name]
                    return "NA" if number is None else f"{number:.3f}"
                comparisons = []
                for q in ("Q2","Q3","Q4"):
                    comparisons.append(
                        f"{q}/Q1={val(q,'ratio_to_Q1_per_chamber')} "
                        f"TVDu={val(q,'theta_u_marginal_TVD_to_Q1')} "
                        f"TVDw={val(q,'theta_w_marginal_TVD_to_Q1')}")
                print(f"{source} MS{st}: " + "; ".join(comparisons))
    print("RESTRICTION: no quadrant-only Geant4 generator approved by this test.")
    return summary, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa", type=Path, default=DATA)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--geometry-only", action="store_true")
    args = parser.parse_args()
    run(args.qa, args.out, args.geometry_only)


if __name__ == "__main__":
    main()

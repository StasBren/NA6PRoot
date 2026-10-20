#!/usr/bin/env python3
"""Honest chunk timing summary. Never estimate runtime from file mtimes.

New native Geant4 chunks store real perf_counter wall seconds. Historical
completed.json without clocks is reported as NOT MEASURED, not assigned times.
Mean and RMS use per-chunk seconds per generated event; errors are resampled
over independent chunk units, not over their correlated muon/electron rows.
"""
from __future__ import annotations
import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_SUITE = REPO / "test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1"
STAGES = ("geant4", "root_validation", "chunk_total")


def stats(values, seed=271828, replicates=2000):
    values = [float(x) for x in values]
    if any(not math.isfinite(x) or x <= 0 for x in values):
        raise ValueError("Invalid positive real timing")
    n = len(values)
    if not n:
        return {"n": 0, "mean": None, "rms": None,
                "mean_error": None, "rms_error": None}
    def moments(a):
        mean = sum(a) / len(a)
        rms = math.sqrt(sum((x - mean)**2 for x in a) / len(a))
        return mean, rms
    mean, rms = moments(values)
    if n < 2:
        return {"n": n, "mean": mean, "rms": rms,
                "mean_error": None, "rms_error": None}
    rng = random.Random(seed)
    boot = [moments([values[rng.randrange(n)] for _ in range(n)])
            for _ in range(replicates)]
    err = lambda v: math.sqrt(sum((x - sum(v)/len(v))**2 for x in v)
                              / (len(v)-1))
    return {"n": n, "mean": mean, "rms": rms,
            "mean_error": err([v[0] for v in boot]),
            "rms_error": err([v[1] for v in boot])}


def analyse(suite):
    frozen = json.loads((suite / "campaign.json").read_text(
        encoding="utf-8"))["frozen_inputs"]
    channels = frozen["channels"]
    target = int(frozen["events_per_source"])
    chunk_size = int(frozen["chunk_size"])
    if target <= 0 or chunk_size <= 0:
        raise ValueError("Invalid frozen campaign event/chunk counts")
    results = {}
    for channel in channels:
        rows = []
        completed = 0
        missing_timing = 0
        for i in range(math.ceil(target/chunk_size)):
            done = suite / channel.lower() / f"chunk{i+1:04d}" / "completed.json"
            if not done.is_file():
                continue
            info = json.loads(done.read_text(encoding="utf-8"))
            n = min(chunk_size, target-i*chunk_size)
            if info["events"] != n or info["channel"] != channel or \
               info["chunk"] != i+1:
                raise ValueError(f"Inconsistent completed checkpoint: {done}")
            completed += n
            timers = info.get("wall_seconds")
            if timers is None:
                missing_timing += 1
                continue
            if info.get("timing_clock") != \
                    "time.perf_counter, elapsed wall seconds":
                raise ValueError(f"Unknown wall-clock provenance: {done}")
            values = {stage: float(timers[stage]) for stage in STAGES}
            if not (0 < values["geant4"] <= values["chunk_total"] and
                    0 < values["root_validation"] <= values["chunk_total"]):
                raise ValueError(f"Invalid process timing: {done}")
            rows.append({"chunk": i+1, "events": n,
                         **{stage: values[stage]/n for stage in STAGES}})
        report = {
            "expected_events": target, "completed_events": completed,
            "timed_chunks": len(rows),
            "legacy_chunks_without_timings": missing_timing,
            "timing_basis": "measured wall seconds / actual generated events",
            "units": "seconds per generated parent event",
            "stages": {stage: stats([r[stage] for r in rows],
                                     seed=271828 + j)
                       for j, stage in enumerate(STAGES)}
        }
        results[channel] = report
    return {
        "source_suite": str(suite),
        "status": "MEASURED" if any(
            r["timed_chunks"] for r in results.values()) else "NOT_MEASURED",
        "bootstrap_unit": "completed Geant4 chunk",
        "note": (
            "Mean/RMS are across per-chunk per-event runtimes; bootstrap "
            "uncertainties only quantify chunk-to-chunk variation. "
            "Never equate throughput of Geant4 with Heed/Garfield folding. "
            "No CPU seconds/RAM claim unless separately instrumented. "
            "Legacy chunks without measured timings stay unavailable."),
        "channels": results,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suite", type=Path, required=True,
                        help="Frozen campaign directory (old or new)")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    report = analyse(args.suite)
    out = args.out or args.suite / "pi_stagec_wall_performance.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"Status: {report['status']}; result: {out}")
    for channel, item in report["channels"].items():
        x = item["stages"]["geant4"]
        if x["n"]:
            print(f"{channel}: n={x['n']} chunks, "
                  f"mean={x['mean']:.6g} s/event, "
                  f"RMS={x['rms']:.6g} s/event, "
                  f"mean_error={x['mean_error']}, "
                  f"RMS_error={x['rms_error']}")
        else:
            print(f"{channel}: no measured timing in completed checkpoints; "
                  f"legacy_without_timing={item['legacy_chunks_without_timings']}")


if __name__ == "__main__":
    main()

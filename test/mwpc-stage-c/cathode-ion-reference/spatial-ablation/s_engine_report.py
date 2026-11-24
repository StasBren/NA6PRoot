#!/usr/bin/env python3
"""Summarize an engine validation/benchmark directory written by run_s_ladder.sh engine-validate / engine-bench.

  s_engine_report.py ENGINE_DIR --out REPORT.md
"""
import argparse
import json
import re
from pathlib import Path


def l2(md):
    t = md.read_text()
    return ("PASS" if "PASS waveforms" in t else "FAIL"), max((float(x) for x in re.findall(r"\| ([0-9.eE+-]+) \|\n", t)), default=float("nan"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("d", type=Path); ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    L = [f"# Response engine validation and benchmark: `{a.d}`", ""]
    cm = sorted(a.d.glob("compare/*_moments.md"))
    if cm:
        L += ["## Moments: batch birth reader vs export-path cloud_moments.csv", "", "| muon | result |", "|---|---|"]
        L += [f"| {p.name.split('_moments')[0]} | {'PASS' if 'PASS cloud moments' in p.read_text() else 'FAIL'} |" for p in cm]
        L.append("")
    cw = sorted(a.d.glob("compare/*_wave_*.md"))
    if cw:
        L += ["## Waveforms (ion term and electron step, full precision)", "", "| muon | implementation | reference | result | max rel L2 |",
              "|---|---|---|---|---|"]
        for p in cw:
            m = re.match(r"(muon_\d+)_wave_(\w+?)_vs_(\w+)\.md", p.name)
            r, v = l2(p)
            L.append(f"| {m.group(1)} | {m.group(2)} | {m.group(3)} | {r} | {v:.2e} |")
        L.append("")
    bt = a.d / "moments" / "batch_timing.json"
    if bt.exists():
        b = json.loads(bt.read_text())
        L += ["## Birth reading + moment construction (C++, one process)", "",
              f"process wall {b['process_wall_s']:.3f} s for {b['n_events']} events, peak RSS {b['peak_rss_kib'] / 1024:.0f} MiB", "",
              "| event | rows | clouds | archive [MiB] | read+validate [ms] | moments [ms] | cold/warm |", "|---|---|---|---|---|---|---|"]
        for e in b["events"]:
            L.append(f"| {e['id']} | {e['rows']} | {e['clouds']} | {e['archive_bytes'] / 2**20:.2f} | {e['read_validate_s'] * 1e3:.2f} | "
                     f"{e['moments_s'] * 1e3:.3f} | {'warm' if e['warm'] else 'cold'} |")
        L.append("")
    for impl in ("vector", "compiled", "loop"):
        tj = a.d / f"bench_{impl}" / "engine_timing.json"
        if not tj.exists():
            continue
        t = json.loads(tj.read_text())
        k = t["kernel"]
        L += [f"## Engine evaluate, implementation `{impl}` (one process, kernel opened once)", "",
              f"import {t['import_s']:.3f} s; kernel verify (sha256) {k['verify_s']:.2f} s; kernel load {k['load_s']:.2f} s; "
              f"resident kernel {k['resident_bytes'] / 2**20:.0f} MiB; process peak RSS {t['peak_rss_mib']:.0f} MiB; "
              f"warm throughput {t['throughput_events_per_s_warm'] or float('nan'):.1f} events/s", "",
              "| event | clouds | read moments [ms] | evaluate [ms] (min of repeats) | cold/warm |", "|---|---|---|---|---|"]
        for e in t["events"]:
            L.append(f"| {e['event']} | {e['clouds']} | {e['read_moments_s'] * 1e3:.2f} | {e['evaluate_s'] * 1e3:.2f} | {'warm' if e['warm'] else 'cold'} |")
        L.append("")
    if bt.exists() and (a.d / "bench_compiled" / "engine_timing.json").exists():
        b = {e["id"]: e for e in json.loads(bt.read_text())["events"]}
        t = {e["event"]: e for e in json.loads((a.d / "bench_compiled" / "engine_timing.json").read_text())["events"]}
        L += ["## Per-event fast path (birth file -> unshaped U/V ion + electron currents), compiled evaluator", "",
              "| event | read+validate | moments | moments CSV read | evaluate | total [ms] |", "|---|---|---|---|---|---|"]
        for ev in b:
            if ev in t:
                parts = [b[ev]["read_validate_s"], b[ev]["moments_s"], t[ev]["read_moments_s"], t[ev]["evaluate_s"]]
                L.append(f"| {ev} | " + " | ".join(f"{x * 1e3:.2f}" for x in parts) + f" | {sum(parts) * 1e3:.2f} |")
        L += ["", "One-time costs (kernel generation, Green table, grid paths) are not included; see the fast-path report and "
              "timing.txt. The moments CSV hand-off between the C++ reader and the Python engine is a file in this test "
              "setup and is timed separately."]
    sj = a.d / "scaling.json"
    if sj.exists():
        s = json.loads(sj.read_text())
        L += ["", "## Scaling with clouds per event (largest event replicated)", "",
              "| clouds | vector [ms] | vector temp [MiB] | compiled [ms] | loop [ms] |", "|---|---|---|---|---|"]
        g = lambda r, i: f"{r[i]['median_ms']:.2f}" if i in r else "-"
        for r in s["rows"]:
            L.append(f"| {r['clouds']} | {g(r, 'vector')} | {r['vector']['temp_peak_mib']:.1f} | {g(r, 'compiled')} | {g(r, 'loop')} |")
    a.out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

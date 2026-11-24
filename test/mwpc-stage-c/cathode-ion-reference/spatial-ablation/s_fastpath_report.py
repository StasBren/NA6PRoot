#!/usr/bin/env python3
"""Fast-path validation and benchmark summary for birth-record moments (run_s_ladder.sh births).

  s_fastpath_report.py RUN TAG --out REPORT.md

Per muon (each RUN/muon_*/births_*/ directory, latest first): moment agreement (compare_moments.md), K2A/K2 waveform
agreement vs the export-path kernel outputs (compare_waveforms_*.md), timings and memory:
  per event (measured): ROOT read + validation, moment construction (p3_birth_moments, min over repeats),
                        kernel application (s_kernel.py apply, median), and their sum = birth-to-waveform fast path
  per event (reference): S0 frozen replay of the T1 export (summary.json wall_s); the path export itself is not timed
                        by prepare and is listed as "not measured"
  one time (measured where recorded): Green grid paths, Green table, harmonic kernels
File sizes: ion_births.root, frozen/T1 exports, S0 outputs, Green table, kernels, moments.
The speed ratio quoted is fast path vs S0 replay only, labelled as such (export and shaping excluded on both sides).
"""
import argparse
import json
import re
from pathlib import Path


def size(p):
    p = Path(p)
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0


def mib(b):
    return f"{b / 2**20:.2f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path); ap.add_argument("tag"); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--archive", type=Path, default=None, help="archive root (for ion_births.root sizes)")
    a = ap.parse_args()
    T = a.tag
    timing = {}
    tf = a.run / "timing.txt"
    if tf.exists():
        for l in tf.read_text().split():
            if "=" in l:
                k, v = l.split("=", 1); timing[k] = float(v)
    kj = a.run / f"kernel_{T}" / "kernel.json"
    kb = json.loads(kj.read_text())["build_s"] if kj.exists() else None
    gj = a.run / f"green_{T}" / "axes.json"
    gb = json.loads(gj.read_text()).get("build_s") if gj.exists() else None
    L = [f"# Birth-record fast path: validation and benchmark (tag {T})", "",
         "## One-time precomputation (shared by all events of this transport/geometry domain)", "",
         f"| step | wall [s] | size [MiB] |", "|---|---|---|",
         f"| Green grid paths (RKF) | {timing.get(f'green_paths_{T}_s', 'not recorded')} | {mib(size(a.run / f'green_paths_{T}'))} |",
         f"| Green table | {gb if gb is not None else 'n/a'} | {mib(size(a.run / f'green_{T}'))} |",
         f"| harmonic kernels | {kb if kb is not None else 'n/a'} | {mib(size(a.run / f'kernel_{T}'))} |", "",
         "## Per event", "",
         "| muon | rows | clouds | moments vs export | K2A waveforms | K2 waveforms | read+validate [ms] | moments [ms] | "
         "kernel K2A [ms] | fast path [ms] | S0 replay [s] | S0 / fast path | peak RSS converter / apply [MiB] |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    files = ["", "## File sizes per event [MiB]", "",
             "| muon | ion_births.root | frozen_export (no vertices) | T1_export | S0 | births moments |", "|---|---|---|---|---|---|"]
    for R in sorted(a.run.glob("muon_*")):
        bds = sorted(R.glob("births_*"), reverse=True)
        if not bds:
            continue
        B = bds[0]
        cj = next(B.glob("*.csv.json"), None)
        conv = json.loads(cj.read_text()) if cj else {}
        aj = B / f"K2Akern_{T}_births" / "kernel_apply.json"
        app = json.loads(aj.read_text()) if aj.exists() else {}
        s0 = json.loads((R / "S0" / "replay" / "summary.json").read_text()) if (R / "S0" / "replay" / "summary.json").exists() else {}
        cm = B / "compare_moments.md"
        mm = ("PASS" if "PASS cloud moments" in cm.read_text() else "FAIL") if cm.exists() else "n/a"
        def wf(m):
            p = B / f"compare_waveforms_{m}.md"
            if not p.exists():
                return "n/a"
            t = p.read_text()
            l2 = max((float(x) for x in re.findall(r"\| ([0-9.eE+-]+) \|\n", t)), default=float("nan"))
            return ("PASS" if "PASS waveforms" in t else "FAIL") + f" (max rel L2 {l2:.2g})"
        rd, mo, ka = conv.get("read_validate_s_min"), conv.get("moments_s_min"), app.get("wall_s_median")
        fp = (rd + mo + ka) if None not in (rd, mo, ka) else None
        g = lambda x, sc=1.0: "n/a" if x is None else f"{x * sc:.4g}"
        ratio = s0["wall_s"] / fp if fp and s0.get("wall_s") else None
        L.append(f"| {R.name} | {conv.get('rows')} | {conv.get('clouds')} | {mm} | {wf('K2A')} | {wf('K2')} | "
                 f"{g(rd, 1e3)} | {g(mo, 1e3)} | {g(ka, 1e3)} | {g(fp, 1e3)} | {g(s0.get('wall_s'))} | {g(ratio)} | "
                 f"{conv.get('peak_rss_kib', 0) / 1024:.0f} / {app.get('peak_rss_kib', 0) / 1024:.0f} |")
        arch = (a.archive / R.name / "ion_births.root") if a.archive else None
        files.append(f"| {R.name} | {mib(conv.get('archive_bytes', size(arch) if arch else 0))} | {mib(size(R / 'frozen_export'))} | "
                     f"{mib(size(R / 'T1_export'))} | {mib(size(R / 'S0'))} | {mib(conv.get('moments_bytes', 0))} |")
    L += files
    L += ["", "## Notes", "",
          "* fast path = ROOT read + validation + moment construction (C++, p3_birth_moments) + K2A kernel application "
          "(numpy, unoptimized per-cloud loop). Shaping and digitization are excluded, as they are for S0.",
          "* S0 replay excludes the RKF path export that it needs; the fast path needs no export. The ratio above is "
          "therefore conservative for the fast path, but it is NOT an end-to-end detector-simulation speedup: Geant4, "
          "Heed, avalanche generation, shaping and digitization are outside both numbers.",
          "* Intermediates: see the README section 'Which intermediates could later be removed'. Nothing is deleted by "
          "this tool or by run_s_ladder.sh births."]
    a.out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

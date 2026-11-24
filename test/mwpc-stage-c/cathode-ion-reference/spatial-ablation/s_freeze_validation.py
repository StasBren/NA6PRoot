#!/usr/bin/env python3
"""Freeze a validated S-ladder / birth-path state: write a NEW manifest directory, modify nothing.

  s_freeze_validation.py RUN TAG --archive ARCHIVE_ROOT --repo REPO --out NEW_DIR

Records: git commit and dirty state; per muon the ion_births.root sha256 (and the archive's own SHA256SUMS check),
prepare provenance, the latest births_* results (moments and waveform comparisons, timings) with sha256 of every
moments/waveform file; the kernel identity used by the engine (sha256 of K.npy, KPHI0.npy, kernel.json), the Green
table axes (LUT, grid-path and grid-ion sha256), the Green grid-path seal; copies (not moves) of the small reports.
Ends with NEW_DIR/SHA256SUMS over the manifest files.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


def sha(p, chunk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def tree_sha(d):
    return {str(p.relative_to(d)): sha(p) for p in sorted(Path(d).rglob("*")) if p.is_file()} if Path(d).exists() else {}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path); ap.add_argument("tag"); ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--repo", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    a.out.mkdir(parents=True)
    git = lambda *x: subprocess.run(["git", "-C", str(a.repo), *x], capture_output=True, text=True).stdout.strip()
    T = a.tag
    kd = a.run / f"kernel_{T}"
    kid = hashlib.sha256("".join(sha(kd / f) for f in ("K.npy", "KPHI0.npy", "kernel.json")).encode()).hexdigest()
    gax = json.loads((a.run / f"green_{T}" / "axes.json").read_text())
    man = {"git_commit": git("rev-parse", "HEAD"), "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
           "git_dirty_files": len([l for l in git("status", "--porcelain").splitlines() if l.strip()]),
           "run": str(a.run.resolve()), "tag": T,
           "kernel": {"dir": str(kd), "engine_kernel_id": kid, "K.npy_sha256": sha(kd / "K.npy"),
                      "KPHI0.npy_sha256": sha(kd / "KPHI0.npy"), "kernel.json_sha256": sha(kd / "kernel.json"),
                      "basis": json.loads((kd / "kernel.json").read_text())["basis"]},
           "green": {"axes.json_sha256": sha(a.run / f"green_{T}" / "axes.json"), "G.npy_sha256": sha(a.run / f"green_{T}" / "G.npy"),
                     "lut_sha256": gax["lut_sha256"], "grid_ions_sha256": gax["grid_ions_sha256"],
                     "grid_vertices_sha256": gax["grid_vertices_sha256"], "shape": gax["shape"],
                     "grid_seal": (a.run / f"green_paths_{T}" / "complete.txt").read_text() if (a.run / f"green_paths_{T}" / "complete.txt").exists() else None},
           "muons": {}}
    for R in sorted(a.run.glob("muon_*")):
        bds = sorted((d for d in R.glob("births_*") if (d / f"K2Akern_{T}_births").exists()), reverse=True)
        if not bds:
            continue
        B = bds[0]
        arch = a.archive / R.name / "ion_births.root"
        sums = subprocess.run(["sha256sum", "-c", "SHA256SUMS"], cwd=a.archive / R.name, capture_output=True, text=True)
        cmp = {}
        for f in ("compare_moments.md", "compare_waveforms_K2A.md", "compare_waveforms_K2.md"):
            t = (B / f).read_text() if (B / f).exists() else ""
            cmp[f] = {"pass": "PASS" in t.splitlines()[-1] if t.strip() else False,
                      "max_rel_L2": max((float(x) for x in re.findall(r"\| ([0-9.eE+-]+) \|\n", t)), default=None) if "waveforms" in f else None}
        man["muons"][R.name] = {
            "archive": str(arch), "archive_sha256": sha(arch), "archive_SHA256SUMS_ok": sums.returncode == 0,
            "births_dir": str(B), "comparisons": cmp,
            "files_sha256": {**{f"{B.name}/{k}": v for k, v in tree_sha(B).items()},
                             **{f"K2Akern_{T}/{k}": v for k, v in tree_sha(R / f"K2Akern_{T}").items()},
                             **{f"S0/{k}": v for k, v in tree_sha(R / "S0").items()},
                             "cloud_moments.csv": sha(R / "cloud_moments.csv") if (R / "cloud_moments.csv").exists() else None}}
        dst = a.out / R.name; dst.mkdir()
        for f in list(B.glob("*.md")) + list(B.glob("*.json")) + list(B.glob("*.txt")) + list(B.glob("*.log")):
            shutil.copy2(f, dst / f.name)
    for f in sorted(a.run.glob(f"fastpath_report_{T}_*.md")) + sorted(a.run.glob(f"plan_{T}*.txt")) + [a.run / "timing.txt"]:
        if f.exists():
            shutil.copy2(f, a.out / f.name)
    (a.out / "manifest.json").write_text(json.dumps(man, indent=1) + "\n")
    with open(a.out / "SHA256SUMS", "w") as f:
        for p in sorted(a.out.rglob("*")):
            if p.is_file() and p.name != "SHA256SUMS":
                f.write(f"{sha(p)}  {p.relative_to(a.out)}\n")
    ok = all(v["archive_SHA256SUMS_ok"] and all(c["pass"] for c in v["comparisons"].values()) for v in man["muons"].values())
    print(f"{'PASS' if ok else 'CHECK'} frozen manifest {a.out}: commit {man['git_commit'][:12]}, kernel id {kid[:16]}, "
          f"{len(man['muons'])} muons, all comparisons PASS and archives intact: {ok}")


if __name__ == "__main__":
    main()

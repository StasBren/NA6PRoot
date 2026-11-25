#!/usr/bin/env python3
"""Allowed-diff check: the K2A-inline generator copy may differ from the FROZEN generator only inside marked blocks.

  check_inline_generator_diff.py [--frozen FROZEN.cxx] [--inline INLINE.cxx] [--report OUT.json]

Rules (exit 0 = PASS, 1 = FAIL):
  * the frozen source has its recorded sha256 (1436a715...f874f) - the frozen file itself is never modified;
  * blocks are '// K2A-INLINE BEGIN <name>' ... '// K2A-INLINE END <name>', not nested, each name exactly once, and the
    set of names is exactly ALLOWED (header, include, sink, end_primary, options, frame, subset, exit);
  * removing every block and putting back the block's '//@frozen <line>' lines (verbatim, in order) gives a file that is
    BYTE-IDENTICAL to the frozen source; so every frozen line is either untouched or kept verbatim inside the block that
    replaces it;
  * no '//@frozen' line and no K2A identifier outside a block; the replaced frozen lines are listed in the report.
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AV = HERE.parents[2] / "mwpc-avalanche"
FROZEN_SHA = "1436a71575cc34c002b49130a017409364e551ed11a35d1ee31df17224f5874f"
ALLOWED = {
    "header": "comment describing the copy",
    "include": "MomentAccumulator.hh include and the hook state (accumulator pointer, options)",
    "sink": "RunGeometry, per Heed primary: BeginPrimary + compose the existing fullBirthSink with the accumulator",
    "end_primary": "RunGeometry, after the frozen outside-patch abort: EndPrimary(births.total)",
    "options": "main: strip --k2a-moments-out / --max-heed-primaries before the frozen argument parser",
    "frame": "main, per crossing: muon frame (frozen Archive::Muon expressions) and accumulator creation",
    "subset": "main: RunGeometry call with onlyIndices = primaries 0..N-1 (frozen call kept as //@frozen), Finish + write",
    "exit": "main: exit code 4 when the moments carry a blocking status",
}
BEGIN = re.compile(r"^// K2A-INLINE BEGIN (\S+)$")
END = re.compile(r"^// K2A-INLINE END (\S+)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen", type=Path, default=AV / "phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx")
    ap.add_argument("--inline", dest="inl", type=Path, default=AV / "phaseB_native_muon_p3_noFW_penning057_k2a_inline_v1.cxx")
    ap.add_argument("--report", type=Path)
    a = ap.parse_args()
    problems = []
    fb = a.frozen.read_bytes()
    sha = hashlib.sha256(fb).hexdigest()
    if sha != FROZEN_SHA:
        problems.append(f"frozen source sha256 {sha} != recorded {FROZEN_SHA}")
    lines = a.inl.read_text().split("\n")
    out, blocks, cur, replaced = [], {}, None, []
    for no, l in enumerate(lines, 1):
        mb, me = BEGIN.match(l), END.match(l)
        if mb:
            if cur:
                problems.append(f"line {no}: nested block {mb.group(1)} inside {cur}")
            cur = mb.group(1)
            if cur in blocks:
                problems.append(f"line {no}: block {cur} appears twice")
            blocks[cur] = {"begin_line": no, "frozen_lines_kept": 0, "new_lines": 0}
            continue
        if me:
            if me.group(1) != cur:
                problems.append(f"line {no}: END {me.group(1)} does not close {cur}")
            else:
                blocks[cur]["end_line"] = no
            cur = None
            continue
        if cur:
            if l.startswith("//@frozen "):
                out.append(l[len("//@frozen "):]); blocks[cur]["frozen_lines_kept"] += 1
                replaced.append({"block": cur, "frozen_line": l[len("//@frozen "):]})
            else:
                blocks[cur]["new_lines"] += 1
            continue
        if "//@frozen" in l or "k2a" in l.lower() or "K2AInline" in l:
            problems.append(f"line {no}: K2A text outside a marked block: {l.strip()[:100]}")
        out.append(l)
    if cur:
        problems.append(f"unterminated block {cur}")
    if set(blocks) != set(ALLOWED):
        problems.append(f"block names {sorted(blocks)} != allowed {sorted(ALLOWED)}")
    restored = "\n".join(out).encode()
    identical = restored == fb
    if not identical:
        fl = fb.decode().split("\n")
        for i, (x, y) in enumerate(zip(out, fl), 1):
            if x != y:
                problems.append(f"restored line {i} differs from frozen: {x!r} vs {y!r}"); break
        else:
            problems.append(f"restored length {len(out)} lines vs frozen {len(fl)} lines")
    # a replaced frozen line must not survive unchanged next to its replacement (it is restored only via //@frozen)
    rep = {
        "tool": "check_inline_generator_diff.py", "status": "PASS" if not problems else "FAIL", "frozen": str(a.frozen),
        "frozen_sha256": sha, "inline": str(a.inl), "inline_sha256": hashlib.sha256(a.inl.read_bytes()).hexdigest(),
        "restored_identical_to_frozen": identical, "blocks": {k: {**v, "purpose": ALLOWED.get(k, "NOT ALLOWED")} for k, v in blocks.items()},
        "replaced_frozen_lines": replaced, "problems": problems}
    if a.report:
        a.report.write_text(json.dumps(rep, indent=1) + "\n")
    for k, v in blocks.items():
        print(f"  block {k:<12} lines {v['begin_line']}-{v.get('end_line', '?')}: {v['new_lines']} new, "
              f"{v['frozen_lines_kept']} frozen kept verbatim")
    for p in problems:
        print("  PROBLEM:", p)
    print(("PASS" if not problems else "FAIL") + f" inline generator diff: {len(blocks)} marked blocks; "
          f"restored source {'IS' if identical else 'is NOT'} byte-identical to the frozen generator (sha256 {sha[:12]}...)")
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())

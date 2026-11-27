#!/usr/bin/env python3
"""Put several replay-format variants on ONE absolute-strip label set for t_shaped_compare (which requires it).

  k2a_label_union.py OUT_MUON_DIR NAME=DIR [NAME=DIR ...]

Each DIR has replay/waveforms.csv and prompt/prompt_waveforms.csv (first column time, then strip labels U_i / V_j).
If every variant has the same labels in the same order, OUT_MUON_DIR/NAME is a symlink to DIR (nothing copied).
Otherwise the union label set (U ascending, then V ascending) is used: a variant lacking a strip gets a column of
exact zeros there (the strip was not computed by that variant: its window did not include it), other values are
copied as TEXT (bit-exact). OUT_MUON_DIR/label_union.json records the union and, per variant, the padded labels.
The charge the other variants carry on padded strips is reported by k2a_recenter_compare (missing strip tails).
Exit 0 written; 1 a variant's replay and prompt label sets differ, or an input is missing.
"""
import argparse
import json
import os
import sys
from pathlib import Path

FILES = ("replay/waveforms.csv", "prompt/prompt_waveforms.csv")


def header(p):
    with open(p) as f:
        return f.readline().rstrip("\n").split(",")


def key(l):
    s, i = l.split("_", 1)
    return (0 if s == "U" else 1, int(i))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("variants", nargs="+")
    a = ap.parse_args()
    vs = []
    for x in a.variants:
        n, d = x.split("=", 1); d = Path(d)
        hs = []
        for f in FILES:
            if not (d / f).is_file():
                print(f"FAIL label union: missing {d / f}"); return 1
            hs.append(header(d / f))
        if hs[0][1:] != hs[1][1:]:
            print(f"FAIL label union: {n} replay and prompt label sets differ"); return 1
        vs.append((n, d, hs[0][1:]))
    a.out.mkdir(parents=True, exist_ok=True)
    union = sorted({l for _, _, ls in vs for l in ls}, key=key)
    same = all(ls == vs[0][2] for _, _, ls in vs)
    rec = {"tool": "k2a_label_union.py", "identical": same, "union": union, "variants": {}}
    for n, d, ls in vs:
        tgt = a.out / n
        if same or ls == union:
            os.symlink(os.path.abspath(d), tgt)
            rec["variants"][n] = {"source": str(d), "padded_labels": [], "copied": False}
            continue
        miss = [l for l in union if l not in ls]
        for f in FILES:
            (tgt / f).parent.mkdir(parents=True, exist_ok=True)
            with open(d / f) as fi, open(tgt / f, "w") as fo:
                h = fi.readline().rstrip("\n").split(",")
                col = {l: i for i, l in enumerate(h)}
                fo.write(",".join([h[0]] + union) + "\n")
                for line in fi:
                    p = line.rstrip("\n").split(",")
                    fo.write(",".join([p[0]] + [p[col[l]] if l in col else "0" for l in union]) + "\n")
        rec["variants"][n] = {"source": str(d), "padded_labels": miss, "copied": True}
    (a.out / "label_union.json").write_text(json.dumps(rec, indent=1) + "\n")
    pad = {n: v["padded_labels"] for n, v in rec["variants"].items() if v["padded_labels"]}
    print("OK label union: " + ("identical label sets (symlinks)" if same else
          f"{len(union)} union strips; zero-padded: " + "; ".join(f"{n} {len(p)} ({','.join(p[:6])}{'...' if len(p) > 6 else ''})" for n, p in pad.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())

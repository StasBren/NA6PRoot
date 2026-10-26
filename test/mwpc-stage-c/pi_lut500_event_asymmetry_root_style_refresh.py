#!/usr/bin/env python3
"""Regenerate only the ROOT/ATLAS-inspired HTML from a verified LUT500 summary.

Conservatively refuse to touch unverifiable data. Never re-run Heed,
Geant4 or Garfield. Preserve summary.json and raw event_counts.csv.gz.
Update HTML and manifest atomically after verifying the prior SHA
fingerprints; do not touch the old full 2M Jpsi LUT20 dashboard.
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

from pi_lut500_event_asymmetry_root_dashboard import SCHEMA,TEMPLATE,sha


def refresh(folder:Path,template:Path=TEMPLATE)->None:
    folder=folder.resolve()
    summary=folder/"summary.json"
    raw=folder/"event_counts.csv.gz"
    html=folder/"dashboard.html"
    manifest=folder/"manifest.json"
    for path in (summary,raw,html,manifest,template):
        if not path.is_file():
            raise FileNotFoundError(f"Required file absent: {path}")
    info=json.loads(manifest.read_text(encoding="utf-8"))
    data=json.loads(summary.read_text(encoding="utf-8"))
    if data.get("schema")!=SCHEMA or info.get("schema")!=SCHEMA:
        raise ValueError("Not a p3 LUT500 parent-asymmetry dashboard")
    if set(info.get("files_sha256",{}))!={
            "summary.json","event_counts.csv.gz","dashboard.html"}:
        raise ValueError("Manifest SHA index incomplete")
    for f in (summary,raw,html):
        expected=info["files_sha256"][f.name]
        if sha(f)!=expected:
            raise ValueError(f"SHA mismatch, refusing style refresh: {f}")
    source=template.read_text(encoding="utf-8")
    if source.count("__PAYLOAD__")!=1:
        raise ValueError("HTML template payload marker count !=1")
    page=source.replace("__PAYLOAD__",
      json.dumps(data,ensure_ascii=False,separators=(",",":"),
                 allow_nan=False).replace("</","<\\/"))
    # Safe intermediate files in SAME directory ensure atomic os.replace.
    with tempfile.NamedTemporaryFile(mode="w",encoding="utf-8",
        prefix=".lut500-style-",suffix=".html.tmp",
        dir=folder,delete=False) as f:
        f.write(page)
        htmltmp=Path(f.name)
    try:
        newsha=sha(htmltmp)
        info["files_sha256"]["dashboard.html"]=newsha
        with tempfile.NamedTemporaryFile(mode="w",encoding="utf-8",
              prefix=".lut500-manifest-",suffix=".json.tmp",
              dir=folder,delete=False) as f:
            f.write(json.dumps(info,indent=2)+"\n")
            manitmp=Path(f.name)
        try:
            # Both replace operations are atomic individually. Re-running
            # the command after a crash requires verifying any mismatch.
            os.replace(htmltmp,html)
            os.replace(manitmp,manifest)
        finally:
            manitmp.unlink(missing_ok=True)
    finally:
        htmltmp.unlink(missing_ok=True)
    if sha(summary)!=info["files_sha256"]["summary.json"] or \
            sha(raw)!=info["files_sha256"]["event_counts.csv.gz"]:
        raise AssertionError("Frozen calculation files changed during refresh")
    print("PASS ROOT/ATLAS-inspired style-only refresh:",html)
    print("PASS preserved physics summary and raw event gzip SHA256")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder",type=Path,required=True)
    args=p.parse_args()
    refresh(args.folder)


if __name__=="__main__":
    main()

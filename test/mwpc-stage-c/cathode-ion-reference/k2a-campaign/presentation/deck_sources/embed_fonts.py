#!/usr/bin/env python3
"""Embed the Public Sans fonts of the avalanche deck (same embedded subset) into the K2A deck, as that deck did."""
import re
import shutil
import sys
import zipfile
from pathlib import Path

SRC = Path("/tmp/claude-0/-home-claude/6e19f127-ea1b-521a-80a5-970a9287006d/scratchpad/aval/x/ppt/fonts")
deck = Path(sys.argv[1])
tmp = deck.with_suffix(".tmp.pptx")
zin = zipfile.ZipFile(deck)
zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
for item in zin.infolist():
    data = zin.read(item.filename)
    if item.filename == "[Content_Types].xml":
        t = data.decode()
        if 'Extension="fntdata"' not in t:
            t = t.replace("<Default ", '<Default Extension="fntdata" ContentType="application/x-fontdata"/><Default ', 1)
        data = t.encode()
    elif item.filename == "ppt/_rels/presentation.xml.rels":
        t = data.decode()
        t = t.replace("</Relationships>",
                      '<Relationship Id="rIdFont1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/font" Target="fonts/font1.fntdata"/>'
                      '<Relationship Id="rIdFont2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/font" Target="fonts/font2.fntdata"/></Relationships>')
        data = t.encode()
    elif item.filename == "ppt/presentation.xml":
        t = data.decode()
        def attrs(m):
            a = m.group(1)
            if "embedTrueTypeFonts" not in a:
                a += ' embedTrueTypeFonts="1"'
            if "saveSubsetFonts" not in a:
                a += ' saveSubsetFonts="1"'
            return "<p:presentation " + a + ">"
        t = re.sub(r"<p:presentation ([^>]*)>", attrs, t, count=1)
        lst = ('<p:embeddedFontLst><p:embeddedFont><p:font typeface="Public Sans" pitchFamily="34" charset="0"/>'
               '<p:regular r:id="rIdFont1"/><p:bold r:id="rIdFont2"/></p:embeddedFont></p:embeddedFontLst>')
        t = re.sub(r"(<p:notesSz[^>]*/>)", r"\1" + lst, t, count=1)
        data = t.encode()
    zout.writestr(item, data)
for f in ("font1.fntdata", "font2.fntdata"):
    zout.write(SRC / f, f"ppt/fonts/{f}")
zout.close()
shutil.move(tmp, deck)
print("embedded Public Sans into", deck)

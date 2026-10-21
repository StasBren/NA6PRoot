#!/usr/bin/env python3
"""Split 496k-group native-first-chunk six-geometry JSON to usable local HTML.

READ-ONLY on original validated pi_six_geometry_native_timing.json.
Streams individual JSON objects from the huge top-level groups array via
Python's stdlib incremental JSONDecoder (no ijson/pandas dependency), so
the entire ~half-million-group document is not loaded into RAM at once.
Makes one small all-source/station overview, plus 18 Jpsi/Phi/Omega x
MS0..MS5 single-station HTML pages with full per-chamber/per-anode and
signed-muon-angle histograms. Temporary shard JSONs are removed afterward.
No Garfield/Heed/Geant4 computation and NO physics-data modification.
"""
from __future__ import annotations
import argparse
import html
import json
import tempfile
from collections import Counter
from pathlib import Path

from pi_six_geometry_native_timing import OUT, SOURCES
from pi_six_geometry_native_timing_html import build

PAGES=OUT/"pages"
INPUT=OUT/"pi_six_geometry_native_timing.json"


def read_groups_stream(path):
    """Yield header then individual group objects with bounded memory.

    Expected writer format: {"meta":...,"group_count":N,"groups":[...]}.
    decode only header (small) and individual group dictionaries.
    """
    decoder=json.JSONDecoder()
    with path.open("r",encoding="utf-8") as stream:
        buffer=stream.read(131072)
        mark=buffer.find(',"groups":[')
        while mark<0:
            more=stream.read(131072)
            if not more or len(buffer)>8_000_000:
                raise ValueError("Missing or excessively large JSON header")
            buffer+=more
            mark=buffer.find(',"groups":[')
        header=json.loads(buffer[:mark]+"}")
        if "meta" not in header or "group_count" not in header:
            raise ValueError("Incomplete native full-JSON metadata")
        yield header
        index=mark+len(',"groups":[')
        decoded=0
        while True:
            while index<len(buffer) and buffer[index] in " \n\r\t,":
                index+=1
            if index>=len(buffer):
                extra=stream.read(131072)
                if not extra:
                    raise ValueError("Premature end of JSON group array")
                buffer+=extra
                continue
            if buffer[index]=="]":
                index+=1
                while index<len(buffer) and buffer[index].isspace():
                    index+=1
                if index>=len(buffer):
                    buffer+=stream.read(32)
                if buffer[index]!="}":
                    raise ValueError("Group array missing enclosing JSON brace")
                if stream.read().strip():
                    raise ValueError("Unexpected trailing JSON data")
                break
            try:
                group,last=decoder.raw_decode(buffer,index)
            except json.JSONDecodeError as error:
                extra=stream.read(131072)
                if not extra:
                    raise ValueError(
                        f"Unable to decode native group {decoded}: {error}") from error
                buffer+=extra
                continue
            if not isinstance(group,dict) or "source" not in group \
                    or "station" not in group:
                raise ValueError("Invalid native group")
            decoded+=1
            yield group
            index=last
            if index>262144:
                buffer=buffer[index:]
                index=0
        if decoded!=header["group_count"]:
            raise ValueError(f"Missing native groups: {decoded}/"
                             f"{header['group_count']}")


def compact_index(records,root,source_names):
    """Give direct per-source/station navigation and compact real data table."""
    rows=[]
    for src in source_names:
        for station in range(6):
            key=(src,station,"p4p0_anode_only")
            g=records.get(key)
            path=f"{src}/MS{station}.html"
            name=f"{src} — MS{station}"
            if g is None:
                raise ValueError(f"Absent native station overview: {name}")
            mean="—" if g["mean_ns"] is None else f"{g['mean_ns']:.2f}"
            rms="—" if g["rms_ns"] is None else f"{g['rms_ns']:.2f}"
            rows.append(
                "<tr><td><a href='"+html.escape(path,quote=True)+"'>"+
                html.escape(name)+"</a></td><td>"+
                str(g["segments"])+"</td><td>"+str(g["collected"])+
                "</td><td>"+mean+"</td><td>"+rms+"</td></tr>")
    out=root/"index.html"
    page="""<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ — навигация по 18 временным и угловым отчётам</title>
<style>body{background:#112940;color:#eaf3ff;max-width:1020px;margin:30px auto;
font:15px/1.5 system-ui,sans-serif;padding:0 20px}
a{color:#87d3f4}p{color:#c6d7e4}
table{width:100%;border-collapse:collapse}
td,th{padding:10px;border-bottom:1px solid #385d75;text-align:left}
tr:hover{background:#1e3e56}
small{color:#b7ccda}</style></head><body>
<h1>NA60+ — реальные временные и угловые спектры</h1>
<p>Первый нативный чанк Geant4 для каждого J/ψ/ϕ/ω; одни и те же
Heed-электроны, шесть вариантов Garfield: шаг 3 / 3.5 / 4 мм,
с промежуточной 75-мкм заземлённой проволокой и без неё.
Выбери источник и станцию, затем в отчёте открой камеру или анод.
Время дрейфа [0,300) нс показано шагом 10 нс.</p>
<p>Сводные Mean/RMS ниже приведены для контрольного варианта p=4 мм,
без полевой проволоки. Они характеризуют собранные электроны с известным
временем, а не всю популяцию мюонов. Статистические ошибки родительских
событий будут рассчитаны отдельно.</p>
<table><thead><tr><th>Открыть</th><th>Газовых сегментов</th>
<th>Электронов на анодах</th><th>Mean, нс</th><th>RMS, нс</th></tr></thead>
<tbody>__ROWS__</tbody></table>
<p><small>Все полные машинно-читаемые распределения и проверки сохранения
числа электронов остаются в исходном JSON. HTML разделён исключительно
для снижения нагрузки на браузер.</small></p></body></html>"""
    out.write_text(page.replace("__ROWS__","\n".join(rows)),
                   encoding="utf-8")
    return out


def generate(source_path,output):
    source_names=list(SOURCES)
    output.mkdir(parents=True,exist_ok=True)
    counts=Counter()
    station_summary={}
    with tempfile.TemporaryDirectory(prefix="six_native_shards_") as temp:
        tmp=Path(temp)
        streams={}
        first={}
        try:
            iterator=read_groups_stream(source_path)
            header=next(iterator)
            meta=header["meta"]
            if len(meta.get("garfield_luts",{}))!=6:
                raise ValueError("Expected exactly six pinned Garfield LUTs")
            for src in source_names:
                folder=output/src
                folder.mkdir(parents=True,exist_ok=True)
                for station in range(6):
                    key=(src,station)
                    path=tmp/f"{src}_MS{station}.json"
                    f=path.open("w",encoding="utf-8")
                    f.write('{"meta":')
                    json.dump(meta,f,ensure_ascii=False,separators=(",",":"))
                    f.write(',"groups":[')
                    streams[key]=f
                    first[key]=True
            for group in iterator:
                src=group["source"]
                station=int(group["station"])
                if src in source_names and 0<=station<6:
                    k=(src,station)
                    f=streams[k]
                    if not first[k]:
                        f.write(",")
                    first[k]=False
                    json.dump(group,f,ensure_ascii=False,
                              separators=(",",":"),allow_nan=False)
                    counts[k]+=1
                    if (int(group["det_id"])==-1 and
                        int(group["anode"])==-1 and
                        group["axis"]=="all" and
                        int(group["anglebin"])==-1):
                        station_summary[(src,station,group["config"])]=group
        finally:
            for f in streams.values():
                if not f.closed:
                    f.write("]}")
                    f.close()
        for (src,station),n in sorted(counts.items()):
            if n==0:
                raise ValueError(f"No native groups in {src} MS{station}")
            inp=tmp/f"{src}_MS{station}.json"
            dest=output/src/f"MS{station}.html"
            build(inp,dest)
            print(f"{src} MS{station}: {n:,} groups; "
                  f"{dest.stat().st_size/(1024*1024):.1f} MiB HTML",
                  flush=True)
    index=compact_index(station_summary,output,source_names)
    print("PASS: real-native six-geometry station HTML pages:",sum(counts.values()),
          "source-specific groups across",len(counts),"pages")
    print("Open overview:",index)
    print("Full original JSON stays unchanged:",source_path)
    return index


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,default=INPUT)
    p.add_argument("--out",type=Path,default=PAGES)
    args=p.parse_args()
    if not args.input.is_file():
        p.error(f"Missing complete native first-chunk JSON: {args.input}")
    generate(args.input,args.out)


if __name__=="__main__":
    main()

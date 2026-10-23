#!/usr/bin/env python3
"""Scientific read-only analysis of the 100k-generated-parents/source MWPC data.

Only uses the 18 already MERGED, SHA-validated native six-geometry
Heed/Garfield station JSON files. Never runs Geant4, Garfield or Heed.
Uses a bounded-memory streaming reader of group rows.

Correct estimands:
- Station collection fraction = collected anode electrons / all born Heed
  primary electrons (including no-hit and no-cluster gas segments).
- Station fast collection fraction = collected anode electrons with known
  t<100ns / all born primary electrons; *not* P(>=2 wires/event).
- Drift mean and centered RMS are conditional on collected electrons with
  known drift time (including overflow t>=300 ns).
- Per-anode data HAVE NO primary-electron denominator. Only electron-
  weighted, conditional-on-collected timing spectra can be compared there.
- A wire's integer index is pitch-dependent. Compare NORMALISED physical
  position within each same-area chamber, NOT anode # i across pitches.
- Parent-event statistical errors and LUT finite-trial systematics have
  NOT been calculated here.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from pathlib import Path

from pi_six_geometry_native_dashboard_pages import read_groups_stream
from pi_six_geometry_native_timing import n_anodes
from pi_2m_native_six_geometry_fold import OUTPUT
from pi_pitch_fieldwire_pilot import tag

SOURCES={"Jpsi":"jpsi40-bore","Phi":"phi40-bore","Omega":"omega40-bore"}
CONFIGS=[tag(p,field) for p in (3.,3.5,4.) for field in (0.,75.)]
NUM_POSITION_BINS=20
DATA_ROOT=OUTPUT/"full"/"merged"
OUT=OUTPUT/"full"/"anode_timing_analysis_en"
COLS=["source","station","configuration","segments","primary_electrons",
      "collected_electrons","known_times","unknown_times",
      "collection_fraction","fast_collected_fraction_all_primaries",
      "fast_fraction_known_times","mean_ns","centered_rms_ns",
      "fraction_ge300_of_known","center_mean_ns","outer_mean_ns",
      "outer_minus_center_mean_ns","center_known_electrons",
      "outer_known_electrons"]


def blank():
    return {"segments":0,"primaries":0,"collected":0,"known":0,
            "unknown":0,"edge":0,"hist":[0]*31,
            "sum_drift_ns":0.,"sum_drift_ns2":0.}


def add(b,row):
    for k in ("segments","primaries","collected","known","unknown","edge"):
        b[k]+=int(row[k])
    h=row["hist"]
    if len(h)!=31:
        raise ValueError("Expected thirty 10-ns bins plus >=300 ns overflow")
    for i,x in enumerate(h):b["hist"][i]+=int(x)
    b["sum_drift_ns"]+=float(row["sum_drift_ns"])
    b["sum_drift_ns2"]+=float(row["sum_drift_ns2"])


def validate_bucket(b):
    if b["collected"]!=b["known"]+b["unknown"]:
        raise AssertionError("Collection = known + unknown violated")
    if sum(b["hist"])!=b["known"]:
        raise AssertionError("31-bin full drift-time conservation violated")
    if b["known"] and (b["sum_drift_ns"]<0 or b["sum_drift_ns2"]<0):
        raise AssertionError("Nonphysical negative drift moments")
    if b["primaries"] and b["collected"]>b["primaries"]:
        raise AssertionError("More collected electrons than born primaries")


def mean_rms(b):
    n=b["known"]
    if not n:return None,None
    mean=b["sum_drift_ns"]/n
    variance=b["sum_drift_ns2"]/n-mean*mean
    if variance< -1e-6:raise AssertionError("Negative physical drift variance")
    return mean,math.sqrt(max(0.,variance))


def stats(b):
    validate_bucket(b)
    mu,rms=mean_rms(b)
    n=b["primaries"];k=b["known"]
    return {
        "segments":b["segments"],"primary_electrons":n,
        "collected_electrons":b["collected"],"known_times":k,
        "unknown_times":b["unknown"],
        "collection_fraction":b["collected"]/n if n else None,
        "fast_collected_fraction_all_primaries":
            sum(b["hist"][:10])/n if n else None,
        "fast_fraction_known_times":
            sum(b["hist"][:10])/k if k else None,
        "mean_ns":mu,"centered_rms_ns":rms,
        "fraction_ge300_of_known":b["hist"][30]/k if k else None,
        "hist":list(b["hist"]),
    }


def analyze_file(source,station,root=DATA_ROOT):
    file=root/source.lower()/f"MS{station}.json"
    if not file.is_file():raise FileNotFoundError(file)
    iterator=iter(read_groups_stream(file))
    header=next(iterator)
    meta=header["meta"]
    if meta.get("n_generated_parent_events_per_source")!=100_000:
        raise ValueError(f"Not the 100k-generated-parents dataset: {file}")
    if meta.get("geometries")!=6 or set(meta["garfield_luts"])!=set(CONFIGS):
        raise ValueError(f"Missing/mixed six frozen Garfield configurations: {file}")
    if meta.get("source_station")!=f"{SOURCES[source]} MS{station}":
        raise ValueError(f"Wrong source/station identity: {file}")
    sums={}
    anodes={}
    chambers={}
    for cfg in CONFIGS:
        anodes[cfg]=[blank() for _ in range(NUM_POSITION_BINS)]
        chambers[cfg]=blank()
    seen=0
    for row in iterator:
        seen+=1
        if row["source"]!=SOURCES[source] or row["station"]!=station:
            raise ValueError(f"Foreign source/station data inside {file}")
        cfg=row["config"]
        if cfg not in anodes:raise ValueError(f"Foreign Garfield geometry: {cfg}")
        axis=row["axis"];angle=int(row["anglebin"])
        if axis!="all" or angle!=-1:
            continue
        det=int(row["det_id"]);wire=int(row["anode"])
        if det==-1 and wire==-1:
            if cfg in sums:raise ValueError(f"Duplicate station totals: {file} {cfg}")
            sums[cfg]=row
        elif det>=0 and wire==-1:
            add(chambers[cfg],row)
        elif det>=0 and wire>=0:
            pitch=3. if cfg.startswith("p3p0_") else 3.5 if cfg.startswith("p3p5_") else 4.
            n=n_anodes(station,pitch)
            if wire>=n:
                raise ValueError(f"Wire out of physical centred finite lattice: {wire}/{n}")
            binindex=min(NUM_POSITION_BINS-1, (NUM_POSITION_BINS*wire)//n)
            add(anodes[cfg][binindex],row)
    if seen!=header["group_count"]:
        raise ValueError("Truncated merged native timing groups")
    if set(sums)!=set(CONFIGS):
        raise ValueError(f"Missing station complete summaries: {file}")
    analyzed={}
    for cfg in CONFIGS:
        station_b=blank()
        add(station_b,sums[cfg])
        ch=chambers[cfg]
        wires=blank()
        for region in anodes[cfg]:add(wires,region)
        validate_bucket(station_b);validate_bucket(ch);validate_bucket(wires)
        fields=("collected","known","unknown","sum_drift_ns","sum_drift_ns2")
        for name in fields:
            for alias,aggregated in (("chambers",ch),("anodes",wires)):
                delta=aggregated[name]-station_b[name]
                if name.startswith("sum_"):
                    if abs(delta)>1.e-7*max(1.,abs(station_b[name])):
                        raise AssertionError(f"Moment conservation {cfg} {alias} {name}: {delta}")
                elif delta!=0:
                    raise AssertionError(f"Electron conservation {cfg} {alias} {name}: {delta}")
        if wires["hist"]!=station_b["hist"] or ch["hist"]!=station_b["hist"]:
            raise AssertionError(f"Per-anode/chamber 10-ns histogram mismatch: {cfg}")
        if station_b["segments"]!=meta["n_muon_gas_crossings_station"]:
            raise AssertionError(f"Real gas-hit exposure mismatch: {file}")
        parts=[stats(region) for region in anodes[cfg]]
        # centre = central 20% of anode positions; outer = outer 20%.
        centre=blank();outer=blank()
        for j,region in enumerate(anodes[cfg]):
            r=2*((j+.5)/NUM_POSITION_BINS)-1
            if abs(r)<.2:add(centre,region)
            if abs(r)>.8:add(outer,region)
        center_time=mean_rms(centre)[0];outer_time=mean_rms(outer)[0]
        overall=stats(station_b)
        overall.update({
            "source":source,"station":station,"configuration":cfg,
            "per_position_bin":parts,
            "center_mean_ns":center_time,"outer_mean_ns":outer_time,
            "outer_minus_center_mean_ns":(
                outer_time-center_time if outer_time is not None
                and center_time is not None else None),
            "center_known_electrons":centre["known"],
            "outer_known_electrons":outer["known"],
            "n_position_bins":NUM_POSITION_BINS,
        })
        analyzed[cfg]=overall
    return analyzed


def sum_across_stations(records):
    # Summation is valid for count-weighted electron observables; NOT
    # unweighted averaging of station means/standard deviations.
    combined={}
    for source in SOURCES:
        for cfg in CONFIGS:
            b=blank()
            for st in range(6):
                r=records[(source,st)][cfg]
                b["segments"]+=r["segments"]
                b["primaries"]+=r["primary_electrons"]
                b["collected"]+=r["collected_electrons"]
                b["known"]+=r["known_times"]
                b["unknown"]+=r["unknown_times"]
                for j,n in enumerate(r["hist"]):b["hist"][j]+=n
                # RMS cannot be recovered from the station-level rounded
                # moments alone. Instead preserve exact mean and RMS by
                # reconstructing E[t^2] from the unrounded station values.
                if r["known"]:
                    b["sum_drift_ns"]+=r["known_times"]*r["mean_ns"]
                    b["sum_drift_ns2"]+=r["known_times"]*(
                        r["centered_rms_ns"]**2+r["mean_ns"]**2)
            combined[(source,cfg)]=dict(
                source=source,station="ALL",configuration=cfg,**stats(b))
    return combined


def safe_percent(x):
    return "—" if x is None else f"{100*x:6.2f}%"


def render_text(records,totals):
    lines=[
        "NA60+ MWPC: VERIFIED 100,000 generated parent events/source",
        "REAL NATIVE Geant4 -> shared Poisson/Heed -> six Garfield LUTs",
        "OBSERVABLES: electron drift since birth (NOT signal/ASIC trigger time).",
        "Mean/RMS: only anode-collected electrons with known drift time.",
        "P(t<100 & anode-collected) / ALL Heed primaries separately shown.",
        "Outer/centre: electron-weighted conditional means in outermost /",
        "central 20% of physical anode indices pooled across chambers.",
        "NO parent-event bootstrap error bars or Garfield-LUT systematic",
        "uncertainty reported; source and station populations are different.",
        "",
        "Source | Config | Gas hits | Anode collection | Mean ns | RMS ns"
        " | Fast/all primaries | P(t>=300|known)",
    ]
    for source in SOURCES:
        for cfg in CONFIGS:
            r=totals[(source,cfg)]
            lines.append(
                f"{source:5s} | {cfg:20s} | {r['segments']:8,d} | "
                f"{safe_percent(r['collection_fraction'])} | "
                f"{r['mean_ns']:7.2f} | {r['centered_rms_ns']:7.2f} | "
                f"{safe_percent(r['fast_collected_fraction_all_primaries'])} | "
                f"{safe_percent(r['fraction_ge300_of_known'])}")
    lines+=["","Paired FIELD WIRE effect at FIXED anode pitch (fw minus no-fw):"]
    for source in SOURCES:
        for st in range(6):
            cfgs=records[(source,st)]
            for p in (3.,3.5,4.):
                x=cfgs[tag(p,0.)];y=cfgs[tag(p,75.)]
                lines.append(f"  {source:5s} MS{st} p={p:g}mm: "
                    f"Δcollection = {100*(y['collection_fraction']-x['collection_fraction']):+.2f} pp; "
                    f"Δmean = {y['mean_ns']-x['mean_ns']:+.2f} ns; "
                    f"ΔRMS = {y['centered_rms_ns']-x['centered_rms_ns']:+.2f} ns; "
                    f"Δfast/all = {100*(y['fast_collected_fraction_all_primaries']-x['fast_collected_fraction_all_primaries']):+.2f} pp")
    lines+=["","Station MS0 and MS5 selected anode-position time differences "
            "(outer - centre, known-collected electrons; no error bars):"]
    for source in SOURCES:
        for st in (0,5):
            for cfg in CONFIGS:
                r=records[(source,st)][cfg]
                v=r["outer_minus_center_mean_ns"]
                lines.append(f"  {source:5s} MS{st} {cfg:20s}: "
                     f"Δt = {('—' if v is None else f'{v:+.2f} ns')}; "
                     f"N_center={r['center_known_electrons']:,}, "
                     f"N_outer={r['outer_known_electrons']:,}")
    return "\n".join(lines)+"\n"


HTML=r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ MWPC · Anode timing & six-geometry analysis</title>
<style>
:root{color-scheme:dark}*{box-sizing:border-box}
body{margin:0;background:#0b1e2d;color:#eaf5ff;font:14px/1.55 system-ui,Arial,sans-serif}
header,main,footer{max-width:1260px;margin:auto;padding:20px 24px}
header{padding-top:30px}.top{font-size:12px;color:#60dcc5;letter-spacing:.13em;font-weight:700}
h1{font-size:28px;line-height:1.2;margin:10px 0}h2{font-size:19px;margin:4px 0 12px}
p{color:#c0d6e5;margin:9px 0}.controls{display:flex;flex-wrap:wrap;gap:12px;margin:20px 0}
label{color:#cae1f2;font-size:12px}select{display:block;background:#193e56;border:1px solid #7492a8;color:white;
padding:9px;border-radius:7px;min-width:130px;margin-top:5px}
.panel{border:1px solid #365774;background:#132e43;border-radius:11px;margin:15px 0;padding:20px}
.chart{width:100%;height:auto;border-radius:6px;background:#10283b}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media(max-width:950px){.grid{grid-template-columns:1fr}}
table{width:100%;border-collapse:collapse;font-size:12px}
th,td{padding:9px 5px;border-bottom:1px solid #33536a;text-align:right;font-variant-numeric:tabular-nums}
th:first-child,td:first-child{text-align:left}
th{color:#9edce3}.hint{font-size:12px;color:#b0c9d8}.warning{padding-left:10px;border-left:3px solid #edba6f}
#legend{display:flex;flex-wrap:wrap;gap:8px;margin-top:12px}
.legend{font-size:11px;border:1px solid #395d72;padding:4px 9px;border-radius:15px}
</style></head><body>
<header><div class="top">NA60+ / MUON SPECTROMETER / MWPC</div>
<h1>Primary-electron drift vs physical anode position</h1>
<p>Native Geant4 (100,000 generated parents per source) → shared Heed births →
six frozen Garfield response maps. All 18 source/station datasets, no new simulation.</p>
<div class="controls"><label>Muon source<select id="src"></select></label>
<label>Muon station<select id="st"></select></label></div>
<p class="warning">These are electron drift times since birth, not ASIC timing resolution
or induced current. Mean and RMS are conditional on anode collection and a known
time; collection efficiency uses all primary electrons. No parent-event error bars yet.
Anode number is pitch dependent, so spatial comparisons use normalised positions
in the same physical chamber aperture. Event-level multiwire probability is not
available in this merged timing dataset.</p></header>
<main>
<section class="panel"><h2>Six-configuration station summary</h2>
<table><thead><tr><th>Pitch / field wire</th><th>Primary e⁻</th><th>Collected</th>
<th>Collection</th><th>Mean [ns]</th><th>RMS [ns]</th><th>Collected &lt;100 ns / all primary</th>
<th>≥300 ns / known</th></tr></thead><tbody id="rows"></tbody></table>
</section>
<section class="panel"><h2>Drift-time spectrum — six geometries</h2>
<p>Fraction of anode-collected electrons with a known drift time per 10 ns bin;
the ≥300 ns overflow is included in the denominator, but not drawn.</p>
<svg class="chart" id="time" viewBox="0 0 1000 350" aria-label="Electron time PDF"></svg>
<div id="legend"></div></section>
<div class="grid">
<section class="panel"><h2>Mean drift time vs anode position</h2>
<p>Electron-weighted mean, conditional on collected and known drift times.
Each point pools physically similar anode positions across station chambers.</p>
<svg class="chart" id="anode" viewBox="0 0 750 360" aria-label="Mean time vs normalized wire position"></svg></section>
<section class="panel"><h2>Fast electrons vs anode position</h2>
<p>P(t &lt; 100 ns | anode-collected, known drift time).
This is NOT per-anode collection efficiency.</p>
<svg class="chart" id="fast" viewBox="0 0 750 360" aria-label="Fast time fraction per anode position"></svg></section>
</div>
<section class="panel"><h2>What the plots can and cannot tell us</h2>
<p><strong>Field-wire effect:</strong> compare neighbouring rows at fixed anode pitch.
A smaller mean conditioned on collected electrons is not by itself better
if collection efficiency falls. Use the collected-within-100-ns / all-primaries
column as a joint speed-and-collection metric.</p>
<p><strong>Spatial effects:</strong> a centre-edge difference can reflect the
native muon illumination pattern as well as finite anode ends, gas boundaries
and electron transport. Sparse outer-anode bins must not be overinterpreted.</p>
<p><strong>Limitations:</strong> B=0, analytic Garfield field, 20 microscopic LUT
trials per position box, no avalanche current, threshold or strip electronics.
Parent-event bootstrap uncertainties and finite-LUT uncertainty are future work.</p>
</section></main>
<footer>Source data: verified two-chunk native 100k/source merge. Original result
files and Russian/English V3 pages remain unchanged.</footer>
<script type="application/json" id="dataset">__DATA__</script>
<script>
"use strict";
const dataset=JSON.parse(document.getElementById("dataset").textContent);
const $=id=>document.getElementById(id);
const configs=dataset.configs,sources=dataset.sources,stations=[0,1,2,3,4,5];
const colors=["#50d7c4","#84b7ff","#eab268","#c19afb","#ff8d98","#add374"];
function select(id,opts){for(const [k,v] of opts){
let x=document.createElement("option");x.value=k;x.textContent=v;$(id).append(x)}}
select("src",sources.map(s=>[s,s==="Jpsi"?"J/ψ":s==="Phi"?"ϕ":"ω"]));
select("st",stations.map(s=>[s,"MS"+s]));
function svg(name,attrs,host,txt){
let node=document.createElementNS("http://www.w3.org/2000/svg",name);
for(const [key,v] of Object.entries(attrs))node.setAttribute(key,v);
if(txt!==undefined)node.textContent=txt;host.append(node);return node;
}
function plot(id,values,xmin,xmax,ymax,titleX,titleY,percent=false){
const host=$(id);host.replaceChildren();
const W=Number(host.getAttribute("viewBox").split(" ")[2]);
const H=Number(host.getAttribute("viewBox").split(" ")[3]);
const margin={l:64,r:17,t:20,b:58};
const xx=x=>margin.l+(x-xmin)/(xmax-xmin)*(W-margin.l-margin.r);
const yy=y=>H-margin.b-y/ymax*(H-margin.t-margin.b);
for(let i=0;i<=5;i++){
 const yyval=ymax*i/5,py=yy(yyval);
 svg("line",{x1:margin.l,x2:W-margin.r,y1:py,y2:py,stroke:"#39566a"},host);
 svg("text",{x:margin.l-6,y:py+5,fill:"#cde6f8","text-anchor":"end","font-size":12},host,
   (percent?100*yyval:yyval).toFixed(percent?0:1));
}
for(let i=0;i<=5;i++){
 const x=xmin+(xmax-xmin)*i/5;
 svg("text",{x:xx(x),y:H-29,fill:"#cde6f8","text-anchor":"middle","font-size":12},host,
   x.toFixed(id==="time"?0:1));
}
svg("text",{x:W/2,y:H-7,fill:"#cde6f8","text-anchor":"middle","font-size":13},host,titleX);
svg("text",{x:margin.l,y:15,fill:"#cde6f8","font-size":12},host,titleY);
values.forEach((item,ci)=>{
  let segment=[];
  function flush(){
    if(segment.length>1)svg("polyline",{points:segment.join(" "),stroke:colors[ci],
      fill:"none","stroke-width":2.3},host);
    else if(segment.length===1){let [x,y]=segment[0].split(",");
      svg("circle",{cx:x,cy:y,r:2.5,fill:colors[ci]},host);}
    segment=[];
  }
  item.forEach(d=>{
   if(d===null||!Number.isFinite(d[1])){flush();return;}
   segment.push(xx(d[0])+","+yy(Math.max(0,d[1])));
  });flush();
});
}
const pct=n=>n===null?"—":(100*n).toFixed(2)+"%";
const num=n=>n===null?"—":n.toFixed(2);
const integer=n=>n.toLocaleString("en-US");
function draw(){
 const source=$("src").value,station=Number($("st").value);
 const key=source+"|"+station,entries=dataset.records[key];
 if(!entries){throw Error("No verified station: "+key)}
 $("rows").replaceChildren();$("legend").replaceChildren();
 const time=[],mean=[],fast=[];
 configs.forEach((cfg,i)=>{
  const r=entries[cfg];if(!r)throw Error("Missing six-config sample: "+cfg);
  let row=document.createElement("tr");
  const cols=[cfg,integer(r.primary_electrons),integer(r.collected_electrons),
    pct(r.collection_fraction),num(r.mean_ns),num(r.centered_rms_ns),
    pct(r.fast_collected_fraction_all_primaries),pct(r.fraction_ge300_of_known)];
  cols.forEach(c=>{let td=document.createElement("td");td.textContent=c;row.append(td)});
  $("rows").append(row);
  time.push(r.hist.slice(0,30).map((v,j)=>[j*10+5,
    r.known_times?100*v/r.known_times:0]));
  mean.push(r.per_position_bin.map((b,j)=>[2*(j+.5)/20-1,
    b.known_times>=100?b.mean_ns:null]));
  fast.push(r.per_position_bin.map((b,j)=>[2*(j+.5)/20-1,
    b.known_times>=100?b.fast_fraction_known_times:null]));
  const leg=document.createElement("span");leg.className="legend";
  leg.style.borderLeft="5px solid "+colors[i];leg.textContent=cfg;
  $("legend").append(leg);
 });
 plot("time",time,0,300,Math.max(1,...time.flat().map(x=>x[1]))*1.15,
   "Drift time [ns]","% per 10 ns bin");
 plot("anode",mean,-1,1,Math.max(1,...mean.flat().filter(x=>x&&x[1]!==null).map(x=>x[1]))*1.16,
   "Normalised anode position: −1 left, +1 right","Mean [ns]");
 plot("fast",fast,-1,1,1,
   "Normalised anode position: −1 left, +1 right","Fraction [%]",true);
}
$("src").addEventListener("change",draw);
$("st").addEventListener("change",draw);
draw();
</script></body></html>'''


def write_reports(records,totals,out):
    out.mkdir(parents=True,exist_ok=True)
    compact={f"{src}|{st}":records[(src,st)]
             for src in SOURCES for st in range(6)}
    payload={"sources":list(SOURCES),"configs":CONFIGS,"records":compact,
             "exposure":"100000 generated parent events per source",
             "interpretation":"electron times since birth; conditional on collection"}
    datafile=out/"station_anode_timing.json"
    text=json.dumps(payload,ensure_ascii=False,separators=(",",":"),allow_nan=False)
    datafile.write_text(text+"\n",encoding="utf-8")
    page=HTML.replace("__DATA__",text.replace("</","<\\/"))
    (out/"anode_timing_report_en.html").write_text(page,encoding="utf-8")
    with (out/"station_summary.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=COLS);w.writeheader()
        for source in SOURCES:
            for st in range(6):
                for cfg in CONFIGS:
                    r=records[(source,st)][cfg]
                    w.writerow({k:r.get(k) for k in COLS})
    summary=render_text(records,totals)
    (out/"summary.txt").write_text(summary,encoding="utf-8")
    print(summary,flush=True)
    print("PASS: 18 station datasets, six frozen geometries, "
          "native anode/chamber/station conservation")
    print("English interactive anode timing:",out/"anode_timing_report_en.html")
    print("Compact stats for physics discussion:",out/"summary.txt")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-root",type=Path,default=DATA_ROOT)
    p.add_argument("--output",type=Path,default=OUT)
    args=p.parse_args()
    records={}
    for source in SOURCES:
        for st in range(6):
            records[(source,st)]=analyze_file(source,st,args.data_root)
            print(f"PASS streamed/validated {source} MS{st}: "
                  "six configurations and full anode timing conservation",
                  flush=True)
    totals=sum_across_stations(records)
    write_reports(records,totals,args.output)


if __name__=="__main__":
    main()

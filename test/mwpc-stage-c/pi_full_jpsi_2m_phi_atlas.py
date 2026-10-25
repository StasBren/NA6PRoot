#!/usr/bin/env python3
"""Offline J/psi 2M near-anode impact-position azimuth atlas for EVERY chamber.

READ-ONLY input: six already merged full J/psi station JSON files. Extract
only real, conserved station/chamber Garfield phi_hit histograms, and the
signed Geant4 muon entrance-angle slices. No drift or fold is re-run.

phi_hit = atan2(v_end, u_end-u_anode) around the wire circumference;
NOT the electron's velocity azimuth or last-step incidence to the normal.
The 20 Garfield LUT trials/spatial cell were reused statistically, so
electron histogram counts are not independent Garfield trajectories.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

from mwpc_checkerboard_projection import build_chambers

CONFIGS=(
    "p3p0_anode_only","p3p0_field75p0um",
    "p3p5_anode_only","p3p5_field75p0um",
    "p4p0_anode_only","p4p0_field75p0um",
)
PITCHES=(3.0,3.5,4.0)
SCHEMA="NA60_JPSI_2M_PHI_HIT_ALL_CHAMBERS_V1"

HTML=r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ / J/ψ / wire-surface electron impact azimuth</title>
<style>
:root{color-scheme:dark;--bg:#091626;--card:#142d44;--muted:#b9cfdd;--fg:#edf8ff}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.48 system-ui,Arial,sans-serif}main{max-width:1300px;margin:auto;padding:22px}
h1{font-size:25px;margin:7px 0}h2{font-size:17px}h3{font-size:14px;margin:8px 0}
a{color:#8de9e6}.nav{display:flex;flex-wrap:wrap;gap:9px;margin:15px 0}
.nav a{border:1px solid #446681;border-radius:7px;padding:7px 12px;
text-decoration:none}.nav a.active{background:#2b6770;color:#fff}
.panel{border:1px solid #36556a;background:var(--card);padding:15px;
border-radius:11px;margin:14px 0}.controls{display:flex;gap:14px;flex-wrap:wrap}
label{display:grid;gap:4px;color:var(--muted);font-size:13px}
select{font:inherit;padding:7px;background:#203f55;color:var(--fg);
border:1px solid #547389;border-radius:6px;min-width:170px}
.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:13px}
.plot{background:#19364b;border:1px solid #385e74;border-radius:9px;
padding:8px;text-align:center}.plot svg{max-width:295px;width:100%;
display:block;margin:auto}.plot small{display:block;color:var(--muted)}
.meta{color:var(--muted);font-size:13px}.compare{display:flex;
gap:12px;flex-wrap:wrap}.compare span{padding:7px 11px;
background:#24485b;border-radius:6px;font-variant-numeric:tabular-nums}
.note{color:#bfcedb;font-size:13px;border-left:3px solid #d5aa77;
padding:8px 13px}strong{color:#fff}
@media(max-width:820px){.grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media(max-width:520px){.grid{grid-template-columns:1fr}}
</style></head><body><main>
<p class="meta">NA60+ / STAGE C / FROZEN GEANT4 + HEED + GARFIELD LUT</p>
<h1>Electron impact position around the anode · J/ψ</h1>
<p class="meta">2,000,000 <strong>generated J/ψ parents</strong>; real station
gas-crossing exposure and collected-electron counts displayed below.
All six geometries share the same primary Heed births. This is
<em>φ<sub>hit</sub></em>, the endpoint <strong>position</strong> around the
wire circumference, not the electron's arrival direction.</p>
<nav class="nav" id="stations" aria-label="Muon spectrometer station"></nav>
<p><a href="../index.html">Full J/ψ dashboard</a> ·
<a id="xy">Physical XY map</a> · <a id="timing">Drift-time/angle dashboard</a> ·
<a href="../../phi_atlas_stats.csv">All chamber numerical comparisons (CSV)</a></p>
<section class="panel"><div class="controls">
<label>Scope / physical detector<select id="det"></select></label>
<label>Geant4 muon entrance angle<select id="axis">
<option value="all">All entrance angles</option>
<option value="u">Signed θu</option>
<option value="w">Signed θw</option></select></label>
<label>Muon angle interval<select id="angle"></select></label>
</div><p id="exposure" class="meta"></p></section>
<section class="panel">
<h2>Six matched pitch / field-wire configurations</h2>
<p class="meta">Equal 10° azimuth bins and a common radial percentage scale
across the six plots. The gray disk represents the anode cross-section;
each sector extends outwards in proportion to the fraction of
collected electrons with a finite φ<sub>hit</sub>.
0° = +u, +90° = +v, 180° = −u, −90° = −v.</p>
<div id="plots" class="grid"></div>
<h2>Field-wire effect at fixed anode pitch</h2>
<div id="comparison" class="compare"></div>
</section>
<section class="panel note">
<strong>Interpretation and limitations.</strong>
Conditional distributions are P(φ<sub>hit</sub> | anode collected, entrance
angle bin), not absolute collection probabilities. Each plot has its own
finite-φ collected-electron denominator, and that sample can change
when field wires are present. Total variation (TV) is 0.5 × sum over
10° bins of the absolute difference between two normalized histograms.
No confidence bands: electrons sampled from the reused 20-trial/cell
Garfield lookup are not independent events. B=0; the grounded 75 µm
field wire is a design-study assumption. No avalanche, electronics
response, or inference of detector efficiency. Unpopulated chambers
remain explicitly empty; no simulated values are interpolated.
</section>
</main>
<script type="application/json" id="payload">__DATA__</script>
<script>
"use strict";
const D=JSON.parse(document.getElementById("payload").textContent);
const C=D.configs, $=id=>document.getElementById(id);
function opt(select,values,preferred) {
  select.replaceChildren();
  for(const [value,label] of values) {
    const o=document.createElement("option");o.value=String(value);
    o.textContent=label;select.appendChild(o);
  }
  if(values.some(x=>String(x[0])===String(preferred)))select.value=String(preferred);
}
const nav=$("stations");
for(let s=0;s<6;s++){
 const a=document.createElement("a");a.href="MS"+s+"_phi_atlas.html";
 a.textContent="MS"+s;if(s===D.station)a.className="active";nav.appendChild(a);
}
$("xy").href="MS"+D.station+"_layout_v3_en.html";
$("timing").href="MS"+D.station+"_en.html";
opt($("det"),[[-1,"Entire MS"+D.station+" station"],
  ...D.detectors.map(r=>[r.id,"DET"+r.id+" · row "+r.row+
  " / col "+r.col+" / layer "+r.layer])],-1);
function angleOptions(){
 const axis=$("axis").value,b=$("angle"),old=b.value;
 if(axis==="all"){opt(b,[[-1,"All angles"]],-1);b.disabled=true;return;}
 b.disabled=false;
 const intervals=D.muonAngleEdges.slice(0,-1).map((v,i)=>
  [i,"["+v+", "+D.muonAngleEdges[i+1]+")°"]);
 opt(b,intervals,old==="-1"?0:old);
}
$("axis").addEventListener("change",()=>{angleOptions();draw()});
$("det").addEventListener("change",draw);
$("angle").addEventListener("change",draw);
function cell(config){
 const det=$("det").value,axis=$("axis").value;
 const a=D.hist[det][config];
 return axis==="all"?a.all:a[axis][Number($("angle").value)];
}
const NS="http://www.w3.org/2000/svg";
function node(parent,tag,attrs,content){
 const e=document.createElementNS(NS,tag);
 for(const [k,v] of Object.entries(attrs||{}))e.setAttribute(k,String(v));
 if(content!==undefined)e.textContent=content;
 parent.appendChild(e);return e;
}
function polar(radius,angle){return [radius*Math.cos(angle),-radius*Math.sin(angle)]}
function pathSlice(ri,ro,a,b){
 const [xi,yi]=polar(ri,a),[xo,yo]=polar(ro,a);
 const [xe,ye]=polar(ro,b),[xj,yj]=polar(ri,b);
 return "M"+xi+" "+yi+" L"+xo+" "+yo+" A"+ro+" "+ro+" 0 0 0 "+
  xe+" "+ye+" L"+xj+" "+yj+" A"+ri+" "+ri+" 0 0 1 "+xi+" "+yi+" Z";
}
function plot(svg,r,maxFraction,color){
 svg.setAttribute("viewBox","-169 -169 338 338");
 for(let a=0;a<360;a+=90){
   const [x,y]=polar(145,a*Math.PI/180);
   node(svg,"line",{x1:0,y1:0,x2:x,y2:y,stroke:"#46657a",
     "stroke-width":0.9,"stroke-dasharray":"3 5"});
 }
 const names=[[0,"0° +u"],[90,"+90° +v"],[180,"180° −u"],[270,"−90° −v"]];
 for(const [angle,label] of names){
   const [x,y]=polar(156,angle*Math.PI/180);
   node(svg,"text",{x,y,fill:"#c5d7e7","font-size":10,
     "text-anchor":"middle","dominant-baseline":"central"},label);
 }
 const n=r.n;
 if(n>0){
  for(let i=0;i<36;i++){
   const z=r.hist[i]/n;
   const ro=34+(maxFraction>0?Math.max(0,z/maxFraction)*105:0);
   if(z>0){
    const a=(-180+i*10+0.5)*Math.PI/180;
    const b=(-180+(i+1)*10-0.5)*Math.PI/180;
    node(svg,"path",{d:pathSlice(34,ro,a,b),fill:color,opacity:0.89});
   }
  }
 }
 node(svg,"circle",{cx:0,cy:0,r:32,fill:"#949ba5",stroke:"#e1e8ed",
   "stroke-width":1.5});
 node(svg,"text",{x:0,y:-2,fill:"#132637","font-size":12,
   "font-weight":"bold","text-anchor":"middle"},"ANODE");
 node(svg,"text",{x:0,y:10,fill:"#1e384b","font-size":8,
   "text-anchor":"middle"},"cross-section");
}
function tv(a,b) {
 if(!a.n||!b.n)return null;
 let s=0;for(let i=0;i<36;i++)
  s+=Math.abs(a.hist[i]/a.n-b.hist[i]/b.n);
 return s/2;
}
function fmt(x,d=1){return Number.isFinite(x)?x.toFixed(d):"—"}
function draw(){
 const series=C.map(k=>[k,cell(k)]);
 const peak=Math.max(0,...series.flatMap(x=>
   x[1].n?x[1].hist.map(v=>v/x[1].n):[0]));
 const plots=$("plots");plots.replaceChildren();
 const axis=$("axis").value;
 const active=series[0][1];
 const interval=axis==="all"?"all muon angles":
  "θ"+axis+" ["+D.muonAngleEdges[Number($("angle").value)]+", "+
   D.muonAngleEdges[Number($("angle").value)+1]+")°";
 const det=$("det").value;
 $("exposure").textContent=(det==="-1"?"MS"+D.station:"DET"+det)+
  " · "+interval+" · muon gas-crossing segments: "+
  active.segments.toLocaleString("en-US")+
  " · common max azimuth-bin probability: "+fmt(100*peak,2)+"%";
 series.forEach(([k,r],i)=>{
   const div=document.createElement("div");div.className="plot";
   const heading=document.createElement("h3");
   const p=i<2?"3.0":i<4?"3.5":"4.0";
   heading.textContent="p = "+p+" mm · "+
    (i%2?"field wire 75 µm":"no field wire");
   div.appendChild(heading);
   const svg=node(div,"svg",{"role":"img",
      "aria-label":"Polar histogram of phi hit for "+k});
   plot(svg,r,peak,i%2?"#f1ae69":"#53d3c3");
   const foot=document.createElement("small");
   foot.textContent="Known φ endpoints: "+r.n.toLocaleString("en-US")+
     " · anode-collected: "+r.collected.toLocaleString("en-US")+
     " / Heed births: "+r.primaries.toLocaleString("en-US")+
     (r.n===0?" · no finite-φ electrons":"");
   div.appendChild(foot);plots.appendChild(div);
 });
 const comparison=$("comparison");comparison.replaceChildren();
 for(let j=0;j<3;j++){
  const a=series[j*2][1],b=series[j*2+1][1],dist=tv(a,b);
  const span=document.createElement("span");
  span.textContent="p="+[3,3.5,4][j].toFixed(1)+" mm · FW vs no-FW TV = "+
     (dist===null?"n/a":fmt(dist,3))+" · collection "+
     (a.primaries?fmt(100*a.collected/a.primaries,1)+"%":"n/a")+
     " → "+(b.primaries?fmt(100*b.collected/b.primaries,1)+"%":"n/a");
  comparison.appendChild(span);
 }
}
angleOptions();draw();
</script></body></html>'''


def source_index(blob,station,source,parents,configs=CONFIGS):
    """Validate ALL physical chamber phi projections, including zero-hit modules."""
    m=blob["meta"]
    if (m.get("n_generated_parent_events_per_source")!=parents or
        m.get("source_channel")!="Jpsi" or
        int(m.get("source_station","").split("MS")[-1])!=station or
        tuple(m.get("garfield_luts",{}))!=tuple(configs) or
        m.get("electron_phi_hit_edges_deg")!=list(range(-180,181,10))):
        raise ValueError("Wrong J/psi exposure, station, or six frozen LUTs")
    edges=m["angle_bins_deg"]
    if len(edges)<3 or any(edges[i]>=edges[i+1] for i in range(len(edges)-1)):
        raise ValueError("Bad signed-muon angle edges")
    shape=[
        {"id":c.det_id,"row":c.row,"col":c.col,"layer":c.layer}
        for group in build_chambers() for c in group if c.station==station
    ]
    shape.sort(key=lambda c:c["id"])
    ids={x["id"] for x in shape}
    if not ids or len(ids)!=len(shape):
        raise ValueError("Physical chamber IDs missing/duplicated")
    lookup={}
    for r in blob["groups"]:
        if (r["station"]!=station or r["source"]!=source or
            r["config"] not in configs):
            raise ValueError("Mixed J/psi station or config")
        if r["anode"]!=-1:
            continue   # Anode sums are NOT chamber sums. Avoid double counting.
        det=int(r["det_id"]);axis=r["axis"];bi=int(r["anglebin"])
        if det!=-1 and det not in ids:
            raise ValueError("Nonphysical detector ID")
        if not ((axis=="all" and bi==-1) or
                (axis in ("u","w") and 0<=bi<len(edges)-1)):
            raise ValueError("Unknown angular projection")
        h=r["phi_hist"]
        if (len(h)!=36 or any(not isinstance(v,int) or v<0 for v in h) or
            sum(h)!=r["phi_n"] or r["phi_n"]>r["collected"] or
            r["segments"]<0 or r["primaries"]<0 or r["collected"]<0):
            raise ValueError("Invalid collected-electron phi histogram")
        key=(det,r["config"],axis,bi)
        if key in lookup:raise ValueError("Duplicate station/chamber phi group")
        lookup[key]={"hist":h,"n":r["phi_n"],
                     "segments":r["segments"],"primaries":r["primaries"],
                     "collected":r["collected"]}
    zero=lambda:{"hist":[0]*36,"n":0,"segments":0,"primaries":0,"collected":0}
    hist={}
    for det in [-1]+sorted(ids):
        groups={}
        for cfg in configs:
            all_row=lookup.get((det,cfg,"all",-1),zero())
            if det==-1 and (det,cfg,"all",-1) not in lookup:
                raise ValueError("Missing station-wide physics row")
            value={"all":all_row}
            for axis in ("u","w"):
                rows=[lookup.get((det,cfg,axis,i),zero())
                      for i in range(len(edges)-1)]
                for k in ("n","segments","primaries","collected"):
                    if sum(r[k] for r in rows)!=all_row[k]:
                        raise ValueError(f"Muon {axis} conditional {k} != all: {det} {cfg}")
                if any(sum(r["hist"][b] for r in rows)!=all_row["hist"][b]
                       for b in range(36)):
                    raise ValueError(f"Muon {axis} phi bin != all: {det} {cfg}")
                value[axis]=rows
            groups[cfg]=value
        hist[str(det)]=groups
    # Both full and angular station PDFs must equal the sum of independent
    # physical DET chamber rows, with zero-hit DETs retained in the dropdown.
    for cfg in configs:
        for axis in ("all","u","w"):
            bins=[-1] if axis=="all" else range(len(edges)-1)
            for bi in bins:
                station_row=(hist["-1"][cfg]["all"] if axis=="all" else
                             hist["-1"][cfg][bi if isinstance(bi,int) else 0])
                if axis=="w":station_row=hist["-1"][cfg]["w"][bi]
                if axis=="u":station_row=hist["-1"][cfg]["u"][bi]
                module_rows=[hist[str(d)][cfg]["all"] if axis=="all" else
                             hist[str(d)][cfg][axis][bi] for d in sorted(ids)]
                for k in ("n","segments","primaries","collected"):
                    if sum(r[k] for r in module_rows)!=station_row[k]:
                        raise ValueError(f"Chambers don't sum to station {axis} {bi} {k}")
                if any(sum(r["hist"][b] for r in module_rows)!=
                       station_row["hist"][b] for b in range(36)):
                    raise ValueError(f"Chamber phi bins don't sum to station: {axis} {bi}")
    return {
        "schema":SCHEMA,"station":station,"source":source,"parents":parents,
        "muonAngleEdges":edges,"phiEdges":list(range(-180,181,10)),
        "configs":list(configs),"detectors":shape,"hist":hist,
    }


def render_station(source_file,dest,source,station,parents=2_000_000):
    blob=json.loads(Path(source_file).read_text(encoding="utf-8"))
    payload=source_index(blob,station,source,parents)
    encoded=json.dumps(payload,ensure_ascii=False,separators=(",",":"),
                       allow_nan=False).replace("</","<\\/")
    if HTML.count("__DATA__")!=1:raise AssertionError("Missing HTML payload marker")
    page=HTML.replace("__DATA__",encoded)
    path=Path(dest);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".incomplete")
    tmp.write_text(page,encoding="utf-8")
    os.replace(tmp,path)
    stats=[]
    for det in [-1]+[r["id"] for r in payload["detectors"]]:
        for p,first in zip(PITCHES,(0,2,4)):
            k0,kf=CONFIGS[first:first+2]
            a=payload["hist"][str(det)][k0]["all"]
            b=payload["hist"][str(det)][kf]["all"]
            tv=(sum(abs(x/a["n"]-y/b["n"]) for x,y in
                    zip(a["hist"],b["hist"]))/2
                if a["n"] and b["n"] else None)
            stats.append({
                "station":station,"det_id":det,"pitch_mm":p,
                "generated_parent_events":parents,
                "muon_gas_crossings":a["segments"],
                "heed_primary_births":a["primaries"],
                "anode_collected_without_FW":a["collected"],
                "anode_collected_with_FW":b["collected"],
                "finite_phi_without_FW":a["n"],"finite_phi_with_FW":b["n"],
                "collection_fraction_without_FW":
                    a["collected"]/a["primaries"] if a["primaries"] else None,
                "collection_fraction_with_FW":
                    b["collected"]/b["primaries"] if b["primaries"] else None,
                "phi_hit_total_variation_FW_vs_noFW":tv,
            })
    return stats


def render_all(merged_folder,pages_folder,source,parents=2_000_000):
    all_stats=[]
    result={}
    for st in range(6):
        inp=Path(merged_folder)/f"MS{st}.json"
        out=Path(pages_folder)/source/f"MS{st}_phi_atlas.html"
        all_stats.extend(render_station(inp,out,source,st,parents))
        result[st]=out
        print(f"PASS φ_hit J/ψ MS{st}: all physical DETs, six geometries "
              f"and signed θu/θw bins -> {out}",flush=True)
    output=Path(pages_folder).parent/"phi_atlas_stats.csv"
    cols=list(all_stats[0])
    tmp=output.with_name(output.name+".incomplete")
    with tmp.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader();w.writerows(all_stats)
    os.replace(tmp,output)
    return result,output


if __name__=="__main__":
    from pi_2m_native_six_geometry_fold import OUTPUT, SOURCE
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--merged",type=Path,default=
                   OUTPUT.parent/"pi_full_jpsi_2m_dashboard_v1/merged/jpsi")
    p.add_argument("--pages",type=Path,default=
                   OUTPUT.parent/"pi_full_jpsi_2m_dashboard_v1/pages")
    args=p.parse_args()
    render_all(args.merged,args.pages,SOURCE["Jpsi"])

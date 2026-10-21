#!/usr/bin/env python3
"""Build an entirely offline PI-facing local Garfield pitch comparison page."""
from __future__ import annotations

import html
import json
from pathlib import Path


def render(report, profile, out):
    groups = report["groups"]
    if len(groups) != 3:
        raise ValueError("Expected THREE Garfield pitch pairs")
    payload = {
        "summary":[{
            "pitch":g["pitch_mm"],
            "metrics":g["point_estimate"],
            "errors":g["bootstrap_standard_errors"],
            "transitions":g["transitions"],
            "endpoints":g["field_wire_endpoints"],
        } for g in groups],
        "profiles":profile,
        "field_wire_um": report["field_wire_diameter_um"],
    }
    encoded=json.dumps(payload,allow_nan=False,separators=(",",":")).replace(
        "</","<\\/")
    content=r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NA60+ MWPC: Garfield pitch / field-wire paired study</title>
<style>
:root {color-scheme:light;--fg:#1c2932;--muted:#526676;--line:#d7dfe6;--surface:#f7f9fa;--a:#2363a1;--b:#c26821}
body{font:15px/1.5 system-ui,sans-serif;max-width:1120px;margin:26px auto;padding:0 22px;color:var(--fg)}
h1{font-size:1.7rem;line-height:1.2} h2{font-size:1.1rem;margin-top:24px}
label{font-weight:600} select{font:inherit;background:#fff;border:1px solid #a0adb9;border-radius:6px;padding:7px}
.controls{display:flex;gap:20px;flex-wrap:wrap;align-items:center;margin:22px 0}
.controls label{display:flex;gap:9px;align-items:center}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:10px}
.stat{border:1px solid var(--line);border-radius:8px;padding:12px}
.num{font-size:1.5rem;font-variant-numeric:tabular-nums}
small,.foot{color:var(--muted)} table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}
td,th{text-align:left;border-bottom:1px solid var(--line);padding:9px}
svg{display:block;width:100%;height:auto;overflow:visible}
.legend{display:flex;gap:20px;margin:12px 0}.legend span:before{content:"";display:inline-block;width:18px;height:3px;background:var(--a);vertical-align:middle;margin-right:6px}.legend .field:before{background:var(--b)}
.note{background:var(--surface);padding:14px;border-radius:8px;border-left:3px solid var(--line)}
@media(prefers-color-scheme:dark){:root{color-scheme:dark;--fg:#e6eef3;--muted:#b0c0cb;--line:#405260;--surface:#1d2a34;--a:#69aef2;--b:#ffb678}body{background:#141d25}select{background:#243442;color:var(--fg)}}
</style></head><body>
<h1>NA60+ MWPC — Garfield drift and field-wire comparison</h1>
<p>True Garfield microscopic single-electron endpoints. Same initial electron
birth sites are paired <em>within each anode pitch</em>. Uniform local u-v box
prior; 3+3 mm gas, B=0, anodes +1800 V, grounded field wires __WIRE_UM__ µm.
Results are <b>not</b> weighted by GEANT/Heed muon incidence angles.</p>
<div class="controls">
<label>Wire pitch <select id="pitch">
<option value="3">3.0 mm</option><option value="3.5">3.5 mm</option>
<option value="4">4.0 mm</option></select></label>
<label>Spatial response <select id="metric">
<option value="mean_time">Mean drift time [ns]</option>
<option value="rms_time">RMS drift time [ns]</option>
<option value="anode_collection">Anode collection probability</option>
<option value="mean_angle">Mean incidence-to-normal angle [deg]</option>
<option value="rms_angle">RMS incidence-to-normal angle [deg]</option>
<option value="p50">P(anode collected in ≤50 ns)</option>
<option value="p100">P(anode collected in ≤100 ns)</option>
</select></label>
</div>
<div class="grid" id="stats"></div>
<h2>Response versus birth position within one wire pitch</h2>
<p class="foot">Horizontal axis: u / pitch, −0.5 to +0.5.
Each point pools the full v=−3..+3 mm range for that u-column.
Pitches cover different *physical* u coordinates, so cross-pitch
comparisons here are not paired electron-by-electron.</p>
<div class="legend"><span>Anode only</span><span class="field">+ grounded field wires</span></div>
<svg id="plot" viewBox="0 0 980 410" role="img" aria-label="Microphysical response across a wire cell">
<rect x="88" y="20" width="862" height="325" fill="none" stroke="currentColor" opacity=".18"/></svg>
<h2>Collection outcomes for matched electrons</h2>
<table><thead><tr><th>Without field wire</th><th>With field wire</th><th>Number</th></tr></thead>
<tbody id="transitions"></tbody></table>
<p id="field-ends" class="foot"></p>
<div class="note">
<b>Reading the statistics:</b> The mean/RMS describe ANODE-collected
electrons only, so the conditional mean can decrease if slower electrons
become uncollected. The <b>same-collected-electron</b> time difference
compares identical starting electrons that reached an anode in BOTH setups.
Uncertainties use a paired, u-v-cell-stratified bootstrap of Garfield
electron trials; they are not parent-Geant4-event uncertainties.
No field-wire endpoint does not establish that the wire is irrelevant.
Final-velocity incidence angle uses the last microscopic transport step.
</div>
<p class="foot">Data and provenance: frozen Garfield pilot CSVs, SHA-pinned
in the paired_summary.json; generated locally without a new physics run.</p>
<script type="application/json" id="payload">__PAYLOAD__</script>
<script>
"use strict";
const raw=JSON.parse(document.getElementById("payload").textContent);
const psel=document.getElementById("pitch"), msel=document.getElementById("metric");
const mappings={
mean_time:["mean_time_no_ns","mean_time_field_ns","ns"],
rms_time:["rms_time_no_ns","rms_time_field_ns","ns"],
anode_collection:["anode_fraction_no","anode_fraction_field",""],
mean_angle:["angle_mean_no_deg","angle_mean_field_deg","deg"],
rms_angle:["angle_rms_no_deg","angle_rms_field_deg","deg"],
p50:["anode_le50ns_fraction_no","anode_le50ns_fraction_field",""],
p100:["anode_le100ns_fraction_no","anode_le100ns_fraction_field",""]};
const ns="http://www.w3.org/2000/svg";
function E(tag,attrs,text,parent){const e=document.createElementNS(ns,tag);
 for(const [k,v] of Object.entries(attrs)) e.setAttribute(k,v);
 if(text!==null&&text!==undefined) e.textContent=text;
 (parent||document.getElementById("plot")).appendChild(e);return e;}
function num(x,d=2){return x===null||!Number.isFinite(x)?"N/A":x.toFixed(d);}
function textStat(label,value,sub){return '<div class="stat"><small>'+label+
 '</small><div class="num">'+value+'</div><small>'+sub+'</small></div>';}
function redraw(){
 const pitch=Number(psel.value), key=msel.value, conf=mappings[key],
 group=raw.summary.find(x=>x.pitch===pitch), data=raw.profiles.filter(
 x=>Number(x.pitch_mm)===pitch);
 const m=group.metrics,err=group.errors;
 const a=conf[0],b=conf[1],unit=conf[2];
 const diff=m[b]-m[a], dkey=key==="mean_time"?"delta_conditional_mean_time_ns":
 key==="rms_time"?"delta_conditional_rms_time_ns":
 key==="mean_angle"?"delta_conditional_mean_angle_deg":
 key==="rms_angle"?"delta_conditional_rms_angle_deg":
 key==="anode_collection"?"delta_anode_fraction_pp":null;
 const dec=x=>num(x,3);
 const delta=key==="anode_collection"?100*diff:diff;
 const errmsg=dkey&&err[dkey]!=null?" ± "+dec(err[dkey]):"";
 document.getElementById("stats").innerHTML=
  textStat("Anode only",dec(m[a])+" "+unit,"Garfield conditional / launched prior")+
  textStat("With "+raw.field_wire_um+" µm field wire",dec(m[b])+" "+unit,
           "Garfield conditional / launched prior")+
  textStat("Field minus no field",(delta>=0?"+":"")+dec(delta)+errmsg+
      (key==="anode_collection"?" pp":" "+unit),
      "Paired stratified error when available")+
  textStat("Same-collected Δmean time",
      dec(m.paired_mean_time_delta_both_ns)+" ± "+
      dec(err.paired_mean_time_delta_both_ns)+" ns",
      "Only electrons that collected in BOTH configs; N="+
      m.n_paired_valid_times);
 const tbody=document.getElementById("transitions");tbody.replaceChildren();
 for(const pair of [["0->0","Not collected","Not collected"],
 ["0->1","Not collected","Collected"],["1->0","Collected","Not collected"],
 ["1->1","Collected","Collected"]]){
 const tr=document.createElement("tr");
 for(const cell of [pair[1],pair[2],group.transitions[pair[0]]]){
 const td=document.createElement("td");td.textContent=cell;tr.appendChild(td);}
 tbody.appendChild(tr);}
 document.getElementById("field-ends").textContent=
 "Endpoints classified on grounded field wires: "+group.endpoints+
 " (among "+(group.metrics.n_both_collected+
 group.transitions["0->0"]+group.transitions["0->1"]+
 group.transitions["1->0"])+" launched electron pairs).";
 const svg=document.getElementById("plot");
 while(svg.children.length>1)svg.lastChild.remove();
 const ys=data.flatMap(d=>[Number(d[a]),Number(d[b])]).filter(Number.isFinite);
 if(!ys.length)return;
 let mn=Math.min(...ys),mx=Math.max(...ys);
 let pad=Math.max((mx-mn)*.10,.005);
 mn-=pad;mx+=pad;
 const xx=u=>88+(u+.5)*862, yy=y=>345-(y-mn)/(mx-mn)*325;
 for(let i=0;i<=4;i++){
 const y=mn+(mx-mn)*i/4,Y=yy(y);
 E("line",{x1:88,y1:Y,x2:950,y2:Y,stroke:"currentColor",opacity:.12});
 E("text",{x:78,y:Y+4,"text-anchor":"end","font-size":13,fill:"currentColor"},num(y,2));}
 for(let i=0;i<=4;i++){
 const x=-.5+i/4,X=xx(x);
 E("line",{x1:X,y1:20,x2:X,y2:345,stroke:"currentColor",opacity:.1});
 E("text",{x:X,y:368,"text-anchor":"middle","font-size":13,fill:"currentColor"},num(x,2));}
 E("text",{x:520,y:398,"text-anchor":"middle","font-size":14,fill:"currentColor"},
   "u / anode pitch");
 for(const [name,col] of [[a,"var(--a)"],[b,"var(--b)"]]){
 let d="";
 for(const item of data){
 const value=Number(item[name]);
 if(!Number.isFinite(value))continue;
 d+=(d?" L":"M")+xx(Number(item.u_over_pitch)).toFixed(2)+
 " "+yy(value).toFixed(2);
 }
 E("path",{d,stroke:col,"stroke-width":3,fill:"none"});
 for(const item of data){
 const value=Number(item[name]);
 if(!Number.isFinite(value))continue;
 const cir=E("circle",{cx:xx(Number(item.u_over_pitch)),
 cy:yy(value),r:3,fill:col});
 E("title",{},name+"; u/p="+num(Number(item.u_over_pitch),3)+
 "; response="+num(value,4),cir);}
 }
}
psel.addEventListener("change",redraw);
msel.addEventListener("change",redraw);
redraw();
</script></body></html>'''
    content=content.replace("__PAYLOAD__",encoded).replace(
        "__WIRE_UM__",html.escape(str(report["field_wire_diameter_um"])))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(content,encoding="utf-8")
    return out

#!/usr/bin/env python3
"""Self-contained station-level parent-bootstrap error bars, first native chunk."""
import argparse
import json
from pathlib import Path
from pi_six_geometry_native_timing import OUT


PAGE=r'''<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>NA60+ — ошибки первого чанка по родительским событиям</title>
<style>body{background:#112940;color:#e4f3ff;max-width:1050px;margin:26px auto;
padding:0 18px;font:15px/1.5 system-ui,sans-serif}
h1{font-size:1.5rem} h2{font-size:1.1rem}
p,small{color:#bbd0df}.panel{padding:14px;background:#18364d;
border:1px solid #294861;border-radius:8px;margin-top:12px}
.controls{display:flex;flex-wrap:wrap;gap:16px}
label{display:flex;flex-direction:column;gap:5px;color:#b7ccdb}
select{padding:6px;background:#294d67;color:#e4f3ff;
border:1px solid #567b98;border-radius:5px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.card{padding:10px;background:#23445b;border-radius:6px}
.card b{display:block;font-variant-numeric:tabular-nums;font-size:1.3rem}
.card span{color:#b7ccdb;font-size:.85rem}
svg{width:100%;height:auto;display:block}
</style></head><body>
<h1>NA60+ — статистические ошибки времени дрейфа</h1>
<p>Первый нативный чанк: 3000 сгенерированных родительских событий
J/ψ, ϕ, ω на источник. В каждом bootstrap-повторении все электроны,
газовые пересечения и камеры одного родительского события повторно
выбираются вместе. Это не независимая выборка электронов.
Полная проверка распределений выполнена против исходного Geant4/Heed/
Garfield first-chunk fold, без повторного запуска Garfield.</p>
<div class="panel controls">
<label>Источник<select id="source"></select></label>
<label>Станция<select id="station"></select></label>
<label>Шаг, полевые проволоки<select id="config"></select></label>
<label>Спектр<select id="which">
<option value="time">Время дрейфа [10 нс]</option>
<option value="approach">Угол прихода электрона [10°]</option>
<option value="phi">Угол попадания φ [10°]</option>
</select></label>
</div>
<div id="stats" class="panel cards"></div>
<div class="panel">
<h2 id="title">Доля электронов на бин, % с ±1 bootstrap SD</h2>
<svg id="plot" viewBox="0 0 970 425"></svg>
<p id="subnote"></p>
</div>
<div class="panel">
<p><b>Статистический смысл:</b> столбцы показывают оценённую долю
собранных электронов; вертикальные отрезки — стандартное отклонение
bootstrap распределения оценки для 3000 исходных сгенерированных событий.
Mean ± error и centered RMS ± error относятся только к электронам с
известным временем. Ошибка от конечных 20 Garfield-траекторий в каждой
ячейке LUT и неопределённость физической модели сюда <b>не входят</b>.
Отдельных parent-errorbars по физическому аноду пока нет.</p>
<p id="scope"></p>
</div>
<script type="application/json" id="payload">__PAYLOAD__</script>
<script>
"use strict";
const data=JSON.parse(document.getElementById("payload").textContent);
const groups=data.groups;
const byKey=new Map(groups.map(x=>[[x.source,x.station,x.config].join("|"),x]));
const els=Object.fromEntries(["source","station","config","which"].map(k=>
 [k,document.getElementById(k)]));
const sourceNames=[...new Set(groups.map(g=>g.source))];
const configs=[...new Set(groups.map(g=>g.config))];
function options(el,values){
 for(const [v,label] of values){
  let opt=document.createElement("option");opt.value=v;opt.textContent=label;
  el.appendChild(opt);
 }}
options(els.source,sourceNames.map(s=>[s,s]));
options(els.station,Array.from({length:6},(_,i)=>[String(i),"MS"+i]));
options(els.config,configs.map(s=>[s,s]));
const NS="http://www.w3.org/2000/svg";
function svg(type,attrs,text){
 const node=document.createElementNS(NS,type);
 for(const [k,v] of Object.entries(attrs))node.setAttribute(k,v);
 if(text!=null)node.textContent=text;
 document.getElementById("plot").append(node);
 return node;
}
function fmt(x,n=3){return x==null||!Number.isFinite(x)?"—":x.toFixed(n);}
function draw(){
 const key=[els.source.value,els.station.value,els.config.value].join("|");
 const g=byKey.get(key);if(!g)return;
 const m=g.statistic,e=g.bootstrap_sd,w=els.which.value;
 const p=document.getElementById("stats");p.replaceChildren();
 function stat(label,x,err,suffix=""){
  let d=document.createElement("div");d.className="card";
  let n=document.createElement("b");
  n.textContent=fmt(x)+" ± "+fmt(err)+" "+suffix;
  let t=document.createElement("span");t.textContent=label;
  d.append(n,t);p.append(d);
 }
 stat("Mean drift (conditional)",m.mean_drift_ns,e.mean_drift_ns,"нс");
 stat("Centered RMS drift",m.centered_rms_drift_ns,e.centered_rms_drift_ns,"нс");
 stat("Анодный сбор",m.anode_collection_fraction==null?null:100*m.anode_collection_fraction,
      e.anode_collection_fraction==null?null:100*e.anode_collection_fraction,"%");
 stat("Mean arrival angle",m.mean_approach_angle_deg,
      e.mean_approach_angle_deg,"°");
 const entries={
  time:["time_pdf_fraction","time_pdf_fraction",0,300,30],
  approach:["approach_angle_pdf_fraction","approach_angle_pdf_fraction",0,180,18],
  phi:["impact_phi_pdf_fraction","impact_phi_pdf_fraction",-180,180,36],
 };
 const [field,efield,x0,x1,n]=entries[w];
 const y=m[field]||Array(n).fill(0),sd=e[efield]||Array(n).fill(0);
 const ymax=Math.max(1,...y.map((v,i)=>100*(v||0)+100*(sd[i]||0)))*1.14;
 const sx=v=>78+v*(840)/(x1-x0)-x0*(840)/(x1-x0);
 const sy=v=>350-v/ymax*303;
 const chart=document.getElementById("plot");chart.replaceChildren();
 for(let k=0;k<=4;k++){
  const Y=sy(ymax*k/4);
  svg("line",{x1:78,y1:Y,x2:918,y2:Y,stroke:"#3c6178"});
  svg("text",{x:65,y:Y+4,"text-anchor":"end","font-size":12,fill:"#e5f3ff"},
     fmt(ymax*k/4,1));
 }
 for(let k=0;k<=6;k++){
  const X=sx(x0+(x1-x0)*k/6);
  svg("text",{x:X,y:377,"text-anchor":"middle","font-size":12,fill:"#e5f3ff"},
      fmt(x0+(x1-x0)*k/6,0));
 }
 const dx=(x1-x0)/n;
 y.forEach((v,i)=>{
  const value=100*(v||0),se=100*(sd[i]||0);
  const x=sx(x0+(i+.5)*dx),left=sx(x0+i*dx)+2;
  svg("rect",{x:left,y:sy(value),width:Math.max(1,840/n-4),
              height:350-sy(value),fill:"#52d5c3"});
  if(Number.isFinite(se)&&se>0){
   svg("line",{x1:x,x2:x,y1:sy(value+se),y2:sy(Math.max(0,value-se)),
       stroke:"#ffda8a","stroke-width":2});
   for(const level of [value+se,Math.max(0,value-se)]){
    svg("line",{x1:x-4,x2:x+4,y1:sy(level),y2:sy(level),
        stroke:"#ffda8a","stroke-width":2});
   }
  }
 });
 const cap=w==="time"?"Время дрейфа, нс":w==="approach"?
 "Угол к нормали анода, °":"Знаковый угол на поверхности анода φ, °";
 svg("text",{x:500,y:406,"text-anchor":"middle","font-size":15,fill:"#e4f3ff"},cap);
 document.getElementById("title").textContent=
  "Доля электронов, %, с ±1 bootstrap SD — "+els.config.value;
 document.getElementById("subnote").textContent=
  w==="time"?"Для исходного времени имеется отдельный overflow ≥300 нс: "+
  "его доля учтена в нормировке и ошибке, но на оси 0–300 не показана.":
  "Угол вычислен для собранных электронов с известным углом; "
  +"нормировка отдельна от гистограммы времени.";
}
for(const el of Object.values(els))el.addEventListener("change",draw);
document.getElementById("scope").textContent=data.scope+
 " Bootstrap repetitions: "+data.replicas+".";
draw();
</script></body></html>'''


def render(source,out):
    data=json.loads(source.read_text(encoding="utf-8"))
    if len(data["groups"])!=3*6*6:
        raise ValueError("Expected 108 source/station/config bootstrap groups")
    output=PAGE.replace("__PAYLOAD__",json.dumps(
        data,ensure_ascii=False,separators=(",",":"),allow_nan=False
    ).replace("</","<\\/"))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(output,encoding="utf-8")
    print("Interactive parent-event bootstrap errorbars:",out)
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,
         default=OUT/"parent_bootstrap/pi_six_geometry_parent_bootstrap.json")
    ap.add_argument("--out",type=Path,
         default=OUT/"parent_bootstrap/pi_six_geometry_parent_bootstrap.html")
    args=ap.parse_args()
    render(args.input,args.out)


if __name__=="__main__":
    main()

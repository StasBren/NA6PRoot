#!/usr/bin/env python3
"""Create self-contained offline interactive six-geometry native timing spectra."""
from __future__ import annotations
import argparse
import html
import json
from pathlib import Path
from pi_six_geometry_native_timing import OUT


HTML=r'''<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NA60+ — время дрейфа для разных анодов</title>
<style>
:root{--bg:#112940;--bg2:#1d3c55;--fg:#e9f4ff;--quiet:#bbcfda;--accent:#4bd5c3;--grid:#355775}
*{box-sizing:border-box} body{max-width:1150px;margin:auto;padding:22px;
font:15px/1.5 system-ui,Arial,sans-serif;background:var(--bg);color:var(--fg)}
h1{font-size:1.6rem;margin-bottom:5px}h2{font-size:1.15rem}
small,p{color:var(--quiet)}
[hidden]{display:none!important}
.panel{background:#142d44;border:1px solid #29465e;border-radius:9px;padding:14px;margin-top:14px}
.controls{display:flex;flex-wrap:wrap;gap:13px;margin-top:12px}
label{font-size:13px;color:var(--quiet);display:flex;flex-direction:column;gap:5px}
select{background:#25465e;color:var(--fg);border:1px solid #446b89;border-radius:5px;padding:7px;font:inherit;min-width:125px}
.statgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:9px;margin-top:13px}
.stat{padding:10px;background:var(--bg2);border-radius:6px}
.stat b{display:block;font-size:1.2rem;font-variant-numeric:tabular-nums}
.stat span{font-size:12px;color:var(--quiet)}
svg{display:block;width:100%;height:auto}
.legend{display:flex;flex-wrap:wrap;gap:10px;font-size:12px}
.circ{display:inline-block;width:9px;height:9px;border-radius:100%;margin-right:4px}
.info{font-size:13px;color:var(--quiet)}
a{color:#89d0ff} .note{padding:12px;border-left:3px solid #7890a0}
</style></head><body>
<h1>NA60+ · Дрейф электронов к анодным проволокам</h1>
<p class="info">Версия интерфейса: <b>MWPC-NATIVE-STATIC-FALLBACK-20261021-v2</b></p>
<p>Реальные локальные газовые пересечения Geant4 → общие для шести вариантов
первичные электроны Heed → Garfield LUT. Время отсчитывается от рождения
электрона, а не от пролёта мюона или порога ASIC.</p>
<div class="panel">
<div class="controls">
<label>Источник<select id="source"></select></label>
<label>Станция<select id="station"></select></label>
<label>Камера<select id="det"></select></label>
<label>Физический анод<select id="wire"></select></label>
<label>Шаг и полевая проволока<select id="config"></select></label>
<label>Угол входа мюона<select id="axis"></select></label>
<label>Интервал угла<select id="angle"></select></label>
<label>Наблюдаемая величина<select id="observable">
<option value="time">Время дрейфа [нс]</option>
<option value="incidence">Угол прихода относительно нормали [°]</option>
<option value="phi">Угол положения попадания на аноде [°]</option>
</select></label>
<label>Отображение<select id="view"><option value="bar">Гистограмма выбранной конфигурации</option><option value="overlay">Наложение всех шести конфигураций</option></select></label>
</div>
<div id="coverage" class="info">Ожидание JavaScript. Если график не появился, статический спектр ниже доступен без скриптов.</div>
<div id="stats" class="statgrid"></div>
</div>
<div class="panel" id="staticFallback">__STATIC_FALLBACK__</div>
<div class="panel">
<h2 id="chartTitle">Доля электронов в интервале 10 нс, %</h2>
<svg id="chart" viewBox="0 0 960 405" role="img" aria-label="Спектр времени дрейфа с шагом 10 наносекунд"></svg>
<div id="legend" class="legend"></div>
<div class="info">Временная нормировка включает ≥300 нс; они вынесены
в отдельный счётчик. Угловые гистограммы используют электроны с известным
углом, что не обязательно совпадает с выборкой известных времен.
Штриховая линия 100 нс показывается только на временном графике.
Угол попадания φ — геометрический угол на поверхности проволоки, а
угол к нормали — направление последнего микроскопического шага Garfield.
Пустые группы не заменяются данными других камер.</div>
</div>
<div class="panel note"><b>Границы физической интерпретации</b>
<p id="meta"></p>
<p>В режиме «анод» число газовых сегментов означает сегменты с ≥1
собранным электроном на выбранном аноде; это не количество входящих
мюонов в его геометрическую область. Статистические ошибки по родительским
Geant4-событиям ещё не вычислены. Mean и centered RMS ниже — описательные
оценки данного первого чанка, не приборная временная разрешающая способность.</p>
</div>
<script type="application/json" id="data">__PAYLOAD__</script>
<script>
"use strict";
const data=JSON.parse(document.getElementById("data").textContent);
const groups=data.groups,meta=data.meta;
const map=new Map(groups.map(g=>[[g.source,g.config,g.station,g.det_id,g.anode,g.axis,g.anglebin].join("|"),g]));
const confs=Object.keys(meta.garfield_luts);
const ids=["source","station","det","wire","config","axis","angle","observable","view"];
const el={
 ...Object.fromEntries(ids.map(x=>[x,document.getElementById(x)])),
 stats:document.getElementById("stats"),
 coverage:document.getElementById("coverage")
};
const fills=["#50d5c4","#80b5ff","#efb86a","#b4a2ff","#ff8f8f","#9ad27c"];
function sel(element,options,wanted){
 const old=wanted===undefined?element.value:wanted;
 element.replaceChildren();
 for(const [value,label] of options){
   const o=document.createElement("option");o.value=String(value);
   o.textContent=label;element.appendChild(o);
 }
 if(options.some(([value])=>String(value)===String(old)))element.value=String(old);
}
const presentSources=[...new Set(groups.map(g=>g.source))];
const sourceOrder=["all",...Object.keys(meta.gas_hit_counts)];
sel(el.source,sourceOrder.filter(k=>presentSources.includes(k))
    .map(k=>[k,k==="all"?"Все источники":k]),
    presentSources.includes("all")?"all":presentSources[0]);
const presentStations=[...new Set(groups.map(g=>g.station))].sort((a,b)=>a-b);
sel(el.station,presentStations.map(i=>[i,"MS"+i]),presentStations[0]);
sel(el.config,confs.map((k,i)=>[k,k.replace("anode_only","без field wire").replace(/field75p0um/,"field wire 75 µm")]),confs[0]);
sel(el.axis,[["all","Без отбора"],["u","θu (со знаком)"],["w","θw (со знаком)"]],"all");
function updateChoices(){
 let st=Number(el.station.value),source=el.source.value;
 let ds=groups.filter(g=>g.source===source&&g.station===st&&g.det_id>=0&&g.anode===-1&&g.axis==="all");
 const detectors=[...new Set(ds.map(g=>g.det_id))].sort((a,b)=>a-b);
 sel(el.det,[[-1,"Все камеры"],...detectors.map(k=>[k,"DET "+k])]);
 const det=Number(el.det.value),config=el.config.value;
 let wires=groups.filter(g=>g.source===source&&g.station===st&&
     g.det_id===det&&g.config===config&&g.anode>=0&&g.axis==="all")
 const wireIds=[...new Set(wires.map(g=>g.anode))].sort((a,b)=>a-b);
 sel(el.wire,[[-1,"Все аноды"],...wireIds.map(i=>[i,"анод "+i])]);
 el.wire.disabled=det===-1;
 el.axis.disabled=Number(el.wire.value)>=0;
 if(el.axis.disabled)el.axis.value="all";
 const axis=el.axis.value;
 const bounds=meta.angle_bins_deg;
 sel(el.angle,axis==="all"?[[-1,"Все углы"]]:
    Array.from({length:bounds.length-1},(_,i)=>[i,
      "["+bounds[i]+", "+bounds[i+1]+")°"]));
 el.angle.disabled=axis==="all";
}
function svg(tag,attr,text){
 const ns="http://www.w3.org/2000/svg",e=document.createElementNS(ns,tag);
 for(const [k,v] of Object.entries(attr||{}))e.setAttribute(k,v);
 if(text!==undefined)e.textContent=text;
 document.getElementById("chart").appendChild(e);return e;
}
function fmt(x,n=2){return x===null||!Number.isFinite(x)?"—":x.toFixed(n)}
function find(conf){
 const k=[el.source.value,conf,Number(el.station.value),
          Number(el.det.value),Number(el.wire.value),el.axis.value,
          el.axis.value==="all"?-1:Number(el.angle.value)].join("|");
 return map.get(k);
}
function draw(){
 updateChoices();
 let rec=find(el.config.value);
 let active=el.view.value==="overlay"?confs.map(c=>[c,find(c)]):
 [[el.config.value,rec]];
 el.stats.replaceChildren();
 function stat(label,value){
 const div=document.createElement("div");div.className="stat";
 const v=document.createElement("b");v.textContent=value;const l=document.createElement("span");l.textContent=label;
 div.append(v,l);el.stats.appendChild(div);
 }
 stat("Газовых сегментов",rec?rec.segments.toLocaleString("ru-RU"):"—");
 stat("Первичных электронов",rec&&Number(el.wire.value)<0?rec.primaries.toLocaleString("ru-RU"):"—");
 stat("Собрано на анод",rec?rec.collected.toLocaleString("ru-RU"):"—");
 stat("Время известно",rec?rec.known.toLocaleString("ru-RU"):"—");
 stat("Неизвестное время",rec?rec.unknown.toLocaleString("ru-RU"):"—");
 stat("Потери на конечных краях",rec&&Number(el.wire.value)<0?rec.edge.toLocaleString("ru-RU"):"—");
 stat("Mean [нс]",rec?fmt(rec.mean_ns):"—");
 stat("Centered RMS [нс]",rec?fmt(rec.rms_ns):"—");
 stat("≥ 300 нс",rec?rec.hist[30].toLocaleString("ru-RU"):"—");
 stat("Mean angle [deg]",rec?fmt(rec.mean_angle_deg):"—");
 let title="MS"+el.station.value+", "+el.source.value+
     " — "+el.config.value;
 el.coverage.textContent=title+(rec?"":" — выбранная группа отсутствует");
 const chart=document.getElementById("chart");chart.replaceChildren();
 const quantity=el.observable.value;
 const field=quantity==="time"?"hist":quantity==="incidence"?"angle_hist":"phi_hist";
 const denominator=quantity==="time"?"known":quantity==="incidence"?"angle_n":"phi_n";
 const nBins=quantity==="time"?30:quantity==="incidence"?18:36;
 const xMin=quantity==="phi"?-180:0;
 const xMax=quantity==="time"?300:180;
 const width=(xMax-xMin)/nBins;
 const xlabel=quantity==="time"?"Время дрейфа от рождения электрона, нс":
   quantity==="incidence"?"Угол подхода электрона к нормали, градусы":
   "Геометрический угол попадания φ, градусы";
 document.getElementById("chartTitle").textContent=
   quantity==="time"?"Доля электронов в интервале 10 нс, %":
   "Доля электронов в интервале 10°, %";
 const valid=active.filter(([c,g])=>g&&g[denominator]>0);
 const yMax=Math.max(1,...valid.flatMap(([c,g])=>
       g[field].slice(0,nBins).map(x=>100*x/g[denominator])))*1.12;
 const xx=v=>75+(v-xMin)/(xMax-xMin)*844, yy=v=>342-v/yMax*295;
 for(let i=0;i<=4;i++){let value=i*yMax/4;let Y=yy(value);
 svg("line",{x1:75,y1:Y,x2:919,y2:Y,stroke:"#385b73","stroke-width":1});
 svg("text",{x:62,y:Y+4,"text-anchor":"end","font-size":12,fill:"#ddeaf0"},fmt(value,1));}
 for(let i=0;i<=6;i++){let x=xMin+(xMax-xMin)*i/6,X=xx(x);
 svg("text",{x:X,y:366,"text-anchor":"middle","font-size":12,fill:"#ddeaf0"},x);}
 svg("text",{x:497,y:394,"text-anchor":"middle","font-size":14,fill:"#daeafa"},xlabel);
 if(quantity==="time"){
  svg("line",{x1:xx(100),y1:42,x2:xx(100),y2:342,stroke:"#fac778",
       "stroke-dasharray":"4,5","stroke-width":2});
  svg("text",{x:xx(100)+4,y:34,"font-size":12,fill:"#fac778"},"100 нс");
 }
 const legend=document.getElementById("legend");legend.replaceChildren();
 active.forEach(([conf,g])=>{
 if(!g||g[denominator]===0)return;
 const color=fills[confs.indexOf(conf)];
 if(el.view.value==="bar"){
  g[field].slice(0,nBins).forEach((count,i)=>{
   const y=100*count/g[denominator];
   svg("rect",{x:xx(xMin+i*width)+1,y:yy(y),
       width:Math.max(1,844/nBins-2),height:342-yy(y),fill:color});
  });
 }else{
 let path=g[field].slice(0,nBins).map((count,i)=>
   (i?"L":"M")+xx(xMin+(i+.5)*width)+" "+
   yy(100*count/g[denominator])).join(" ");
 svg("path",{d:path,stroke:color,"stroke-width":3,fill:"none"});
 }
 let x=document.createElement("span");
 x.innerHTML="<i class='circ' style='background:"+color+"'></i>";
 x.append(document.createTextNode(conf));
 legend.append(x);
 });
}
function safeDraw(){
 try{
  draw();
  if(el.stats.children.length<10 ||
     document.getElementById("chart").children.length<5){
    throw new Error("Недостаточно элементов отрисовано: stats "+
      el.stats.children.length+", svg "+
      document.getElementById("chart").children.length);
  }
  document.getElementById("staticFallback").hidden=true;
 }
 catch(error){
  console.error("Native six-geometry dashboard rendering failed",error);
  const banner=document.getElementById("coverage");
  banner.textContent="Ошибка построения графика: "+String(error)+
    ". Исходные физические данные сохранены; сообщи текст ошибки.";
  banner.style.color="#ffca91";
 }
}
for(const k of ids.filter(k=>k!=="view"))el[k].addEventListener("change",safeDraw);
el.view.addEventListener("change",safeDraw);
document.getElementById("meta").textContent=
 meta.data_scope+" · "+meta.input_model+
 " · "+meta.statistical_errors+
 " · физическая ширина газового объёма сохраняется; анодов на станциях: "+
 JSON.stringify(meta.finite_anodes_by_pitch)+".";
safeDraw();
</script></body></html>'''




UI_REVISION = "MWPC-NATIVE-STATIC-FALLBACK-20261021-v2"


def static_native_spectrum(payload):
    """Server-render a REAL selected station spectrum, with zero JavaScript.

    This must remain visible if an older browser aborts its interactive
    payload processing, and is hidden only after all chart cards render.
    """
    items=[g for g in payload["groups"]
           if g.get("det_id")==-1 and g.get("anode")==-1
           and g.get("axis")=="all" and g.get("anglebin")==-1]
    if not items:
        raise ValueError("No station reference to render static histogram")
    conf="p3p0_anode_only"
    first=next((g for g in items if g["config"]==conf),items[0])
    count=int(first["known"])
    hist=first["hist"]
    if len(hist)!=31 or count<0 or sum(hist)!=count:
        raise ValueError("Static page baseline does not conserve time bins")
    maximum=max([100.*c/count for c in hist[:30]]+[1.]) if count else 1.
    def xx(v):
        return 70+v/300*850
    def yy(v):
        return 310-v/(maximum*1.12)*270
    parts=['<h2>Резервная гистограмма — построена Python без JavaScript</h2>',
           '<p class="info">Версия интерфейса: '+UI_REVISION+
           '. Если интерактивный график остаётся пустым, этот результат '
           'отображает действительные данные для выбранной станции и '
           'конфигурации p=3 мм без полевой проволоки.</p>',
           '<div class="statgrid">']
    for label,value in [
        ("Источник",first["source"]),("Станция","MS"+str(first["station"])),
        ("Конфигурация",first["config"]),
        ("Газовых сегментов",str(first["segments"])),
        ("Первичных электронов",str(first["primaries"])),
        ("Собрано на анод",str(first["collected"])),
        ("Известное время",str(count)),
        ("Mean, нс",format(first["mean_ns"],".3f") if first["mean_ns"] is not None else "—"),
        ("Centered RMS, нс",format(first["rms_ns"],".3f") if first["rms_ns"] is not None else "—"),
        ("≥300 нс",str(hist[30])),
    ]:
        parts+=['<div class="stat"><b>',html.escape(str(value)),
                '</b><span>',html.escape(label),'</span></div>']
    parts+=['</div>',
            '<svg viewBox="0 0 970 365" aria-label="Гистограмма времени дрейфа, резервный график">']
    for i in range(5):
        y=yy(maximum*1.12*i/4)
        parts.append(f'<line x1="70" x2="920" y1="{y:.2f}" y2="{y:.2f}" '
                     f'stroke="#395775"/>')
        parts.append(f'<text x="59" y="{y+4:.2f}" fill="#e4f4ff" '
                     f'font-size="12" text-anchor="end">{maximum*1.12*i/4:.1f}</text>')
    for i in range(7):
        x=xx(i*50.)
        parts.append(f'<text x="{x:.2f}" y="338" fill="#e4f4ff" '
                     f'font-size="12" text-anchor="middle">{i*50}</text>')
    if count:
        for i,n in enumerate(hist[:30]):
            value=100.*n/count
            y=yy(value)
            parts.append(f'<rect x="{xx(i*10)+1:.2f}" y="{y:.2f}" '
                         f'width="{850/30-2:.2f}" height="{310-y:.2f}" '
                         f'fill="#4bd5c3"/>')
    x=xx(100.)
    parts+=[
        f'<line x1="{x:.2f}" x2="{x:.2f}" y1="28" y2="310" '
        'stroke="#fac778" stroke-dasharray="4 5"/>',
        f'<text x="{x+5:.2f}" y="20" fill="#fac778" font-size="12">100 нс</text>',
        '<text x="500" y="361" fill="#e4f4ff" font-size="13" '
        'text-anchor="middle">Время дрейфа электрона, нс</text>',
        '</svg>']
    return "".join(parts)


def build(source,dest):
    payload=json.loads(source.read_text(encoding="utf-8"))
    if len(payload["meta"]["garfield_luts"])!=6:
        raise ValueError("Six Garfield map identities required")
    static=static_native_spectrum(payload)
    page=HTML.replace("__PAYLOAD__",json.dumps(payload,ensure_ascii=False,
                separators=(",",":"),allow_nan=False).replace("</","<\\/"))
    page=page.replace("__STATIC_FALLBACK__",static)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(page,encoding="utf-8")
    print("Open offline interactive native timing dashboard:",dest)
    return dest


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,default=OUT/"pi_six_geometry_native_timing.json")
    ap.add_argument("--out",type=Path,default=OUT/"pi_six_geometry_native_timing.html")
    args=ap.parse_args()
    build(args.input,args.out)


if __name__=="__main__":
    main()

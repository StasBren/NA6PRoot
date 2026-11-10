#!/usr/bin/env python3
"""Read-only clickable (u,v) real-muon cross-section of Penning-ON MWPC.

Only *existing* Geant4 source and TrackHeed cluster centroids are shown.
The diagram does NOT invent drift trajectories, electron arrival points,
ion clouds, or cathode signals. Anodes are the actual p3_noFW analytic
wire lattice translated back into the original Geant4 u coordinates.

Returns the clicked original Heed cluster index to the Streamlit muon
chronology view. Selection uses built-in st.plotly_chart(on_select="rerun").
"""
from __future__ import annotations
import csv
import hashlib
import math
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
PITCH_MM=3.0
CATHODE_V_MM=3.0
WIRE_DIAMETER_MM=.03
CROSSING_COLUMNS=("source","chunk","event","hit_index","track_id",
                  "station","det_id","u_in_mm","v_in_mm","w_in_mm",
                  "u_out_mm","v_out_mm","w_out_mm")
CLUSTER_COLUMNS=("crossing_id","cluster_index","u_rel_cm","v_abs_cm",
                 "w_rel_cm","t_ns","n_original","n_accepted")

def need(c,reason):
    if not c:raise RuntimeError(reason)
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            digest.update(block)
    return digest.hexdigest()
def csvrows(path,required):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        need(r.fieldnames is not None and set(required).issubset(r.fieldnames)
             and len(set(r.fieldnames))==len(r.fieldnames),
             "Physical muon / TrackHeed cluster CSV schema changed: "+str(path))
        yield from r
def cpp_round(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)
def entrance_anode_base(u_in,station):
    """Original real station parity and p3_noFW 3mm wire lattice."""
    full_width=307.2 if station==0 else 512.
    n=math.floor(full_width/PITCH_MM+.5)
    phase=PITCH_MM/2 if n%2==0 else 0.
    nearest=(cpp_round((u_in-phase)/PITCH_MM)*PITCH_MM+phase)
    offset=u_in-nearest
    need(abs(offset)<=1.5+1e-8,"Unphysical nearest-anode phase")
    return nearest,offset

def get_real_muon(index,meta):
    """All return coordinates in true global Geant4 u,v,w [mm]."""
    need(0<=index<len(meta["muons"])==1000,
         "Selected muon is not in the sealed real 1000-muon campaign")
    expected=meta["muons"][index]
    input_path=ROOT/"input/native_muon_crossings_input.csv"
    input_manifest=ROOT/"input/input_manifest.json"
    need(input_path.is_file() and input_manifest.is_file(),
         "Missing original selected native real Geant4 muon input")
    im=json.loads(input_manifest.read_text())
    need(im.get("n_parents")==1000 and
         sha(input_path)==im["selected_input_sha256"],
         "Original real Geant4 muon input SHA differs from frozen input")
    with input_path.open(newline="") as f:
        rs=csv.DictReader(f)
        need(rs.fieldnames is not None and set(CROSSING_COLUMNS).issubset(rs.fieldnames),
             "Original muon source CSV missing spatial fields")
        native=None
        for j,r in enumerate(rs):
            if j==index:native=r;break
    need(native is not None,"Selected Geant4 muon not in 1000 frozen input")
    identifier="|".join(native[f] for f in
        ("source","chunk","event","hit_index","track_id","station","det_id"))
    need(identifier==expected["id"],
         "Native Geant4 spatial input does not match timing-cache muon identity")
    station=int(native["station"])
    a=np.array([float(native[f"{c}_in_mm"]) for c in ("u","v","w")])
    b=np.array([float(native[f"{c}_out_mm"]) for c in ("u","v","w")])
    delta=b-a
    need(np.isfinite(a).all() and np.isfinite(b).all()
         and np.linalg.norm(delta)>.05
         and max(abs(a[1]),abs(b[1]))<=3.05,
         "Bad real Geant4 muon geometry source")
    nearest,phase_offset=entrance_anode_base(float(a[0]),station)
    chunk=index//20
    d=ROOT/"native_chunks"/f"chunk_{chunk:04d}"
    rec=json.loads((d/"complete.json").read_text())
    need(rec["chunk"]==chunk and rec["n_parents"]==20
         and rec["stats"]["n_capped_seeds"]==0,
         "Real muon chunk is not a complete immutable Garfield checkpoint")
    name="muon_reference_heed_clusters.csv"
    need(sha(d/name)==rec["file_sha256"][name],
         "Original TrackHeed cluster geometry failed SHA")
    clusters=[]
    for r in csvrows(d/name,CLUSTER_COLUMNS):
        if r["crossing_id"]!=identifier:continue
        point=np.array([a[0]+10*float(r["u_rel_cm"]),
                        10*float(r["v_abs_cm"]),
                        a[2]+10*float(r["w_rel_cm"])])
        need(np.isfinite(point).all(),
             "Non-finite native TrackHeed cluster 3D coordinates")
        parallel=float(np.dot(point-a,delta)/np.dot(delta,delta))
        transverse=float(np.linalg.norm(point-(a+parallel*delta)))
        need(-.02<=parallel<=1.02 and transverse<1.,
             "TrackHeed cluster is not located on the original Geant4 muon chord")
        cl={
            "index":int(r["cluster_index"]),
            "u":float(point[0]),"v":float(point[1]),"w":float(point[2]),
            "t_ns":float(r["t_ns"]),
            "n_original":int(r["n_original"]),
            "n_accepted":int(r["n_accepted"]),
            "anode_u":nearest+cpp_round((point[0]-nearest)/PITCH_MM)*PITCH_MM,
            "fraction_path":parallel,
        }
        need(cl["n_original"]>=cl["n_accepted"]>=0,
             "Invalid actual TrackHeed cluster electron census")
        clusters.append(cl)
    clusters.sort(key=lambda cl:cl["index"])
    need(len(clusters)==expected["n_clusters"] and
         [c["index"] for c in clusters]==list(range(len(clusters))) and
         sum(c["n_accepted"] for c in clusters)==expected["n_Heed"],
         "Physical TrackHeed clusters or accepted-primary counts inconsistent")
    return {
        "muon_index":index,
        "a":a.tolist(),"b":b.tolist(),
        "wire_origin":nearest,"station":station,
        "phase_offset_mm":phase_offset,
        "clusters":clusters,
        "id":identifier,
        "n_Heed":expected["n_Heed"]
    }

def build_figure(model):
    """Return Plotly figure and the EXACT index of the selectable cluster trace."""
    import plotly.graph_objects as go
    a,b=np.asarray(model["a"]),np.asarray(model["b"])
    cls=model["clusters"]
    coordinates=[a[0],b[0]]+[c["u"] for c in cls]
    low=min(coordinates)-PITCH_MM*1.2
    high=max(coordinates)+PITCH_MM*1.2
    wire_base=model["wire_origin"]
    kmin=math.floor((low-wire_base)/PITCH_MM)
    kmax=math.ceil((high-wire_base)/PITCH_MM)
    # For a near-horizontal grazing muon, do not fabricate or
    # silently crop wires. Use actual full relevant lattice.
    wires=wire_base+np.arange(kmin,kmax+1)*PITCH_MM
    fig=go.Figure()
    for v,name,col in [
        (+CATHODE_V_MM,"Upper cathode (+3 mm)","#777777"),
        (-CATHODE_V_MM,"Lower cathode (−3 mm)","#777777"),
    ]:
        fig.add_hline(y=v,line_color=col,line_width=2,
                      annotation_text=name,annotation_position="top left",
                      annotation_font_size=11)
    # Nearest-wire midpoint boundaries are displayed as geometric
    # guides, not simulated electron drift streamlines or a measured LUT.
    for w in wires[:-1]+PITCH_MM*.5:
        fig.add_vline(x=float(w),line_color="#B7BFC7",
                      line_dash="dot",line_width=.8,opacity=.35)
    fig.add_trace(go.Scatter(
        x=wires,y=np.zeros(len(wires)),mode="markers",
        marker=dict(symbol="circle",size=11,color="#111111"),
        name="Anode wires (3 mm pitch; marker enlarged)",
        hovertemplate="Anode u=%{x:.3f} mm<br>v=0 mm"
                      "<br>Wire diameter=30 μm<extra></extra>"))
    fig.add_trace(go.Scatter(
        x=[a[0],b[0]],y=[a[1],b[1]],mode="lines+markers",
        line=dict(color="#3D82A6",width=3),
        marker=dict(size=9,symbol=["triangle-up","triangle-down"]),
        name="Geant4 muon entry → exit",
        text=["Gas entry","Gas exit"],
        hovertemplate="%{text}: u=%{x:.3f} mm, v=%{y:.3f} mm<extra></extra>"))
    if cls:
        u=[c["u"] for c in cls]
        v=[c["v"] for c in cls]
        sizes=[max(10,min(24,10+2.6*math.sqrt(c["n_accepted"])))
               for c in cls]
        cd=[[c["index"],c["n_original"],c["n_accepted"],
             c["t_ns"],c["w"],c["anode_u"]] for c in cls]
        cluster_trace_index=len(fig.data)
        fig.add_trace(go.Scatter(
            x=u,y=v,mode="markers",
            customdata=cd,
            marker=dict(size=sizes,color="#E46B08",
                        line=dict(color="#FFFFFF",width=1.5),opacity=.95),
            name="TrackHeed cluster — CLICK to highlight timeline",
            hovertemplate="Cluster #%{customdata[0]:.0f}"
              "<br>u=%{x:.3f} mm; v=%{y:.3f} mm"
              "<br>w=%{customdata[4]:.3f} mm"
              "<br>Original electrons=%{customdata[1]:.0f}"
              "<br>Accepted electrons=%{customdata[2]:.0f}"
              "<br>Heed t=%{customdata[3]:.5g} ns"
              "<br>Nearest anode u≈%{customdata[5]:.3f} mm"
              "<extra></extra>"))
    else:
        cluster_trace_index=None
    fig.update_layout(
        title="Physical (u,v) projection of the real muon crossing",
        xaxis_title="Global Geant4 u [mm]",
        yaxis_title="Global Geant4 v [mm]",
        template="plotly_white",height=510,
        font=dict(family="Arial,Helvetica,sans-serif",size=13),
        legend=dict(orientation="h",x=0,y=-.16),
        margin=dict(l=65,r=40,t=70,b=110))
    fig.update_xaxes(range=[low,high],showgrid=False,
                     zeroline=False)
    fig.update_yaxes(range=[-3.48,3.48],showgrid=False,
                     zeroline=False,scaleanchor="x",scaleratio=1)
    return fig,cluster_trace_index

def select_cluster(event,cluster_trace,clusters):
    """Use Streamlit point selection, rejecting non-cluster traces."""
    if cluster_trace is None or event is None:return None
    try:
        points=event["selection"]["points"]
    except (KeyError,TypeError):
        return None
    for p in reversed(points):
        if int(p.get("curve_number",-1))!=cluster_trace:continue
        data=p.get("customdata")
        if isinstance(data,(list,tuple)) and data:
            j=int(data[0])
            if any(c["index"]==j for c in clusters):return j
        index=p.get("point_index",p.get("point_number"))
        if index is not None:
            i=int(index)
            if 0<=i<len(clusters):return int(clusters[i]["index"])
    return None

def render(st,index,meta):
    """Draw true cross-section and return selected physical cluster index."""
    geom=get_real_muon(index,meta)
    fig,trace=build_figure(geom)
    st.markdown("#### Muon path and primary ionisation — interactive geometry")
    st.caption(
        "True Geant4 entrance/exit projected onto u–v, real TrackHeed "
        "cluster positions, and p3_noFW 3 mm anode lattice with cathodes "
        "at v = ±3 mm. Anode circles are enlarged for visibility. "
        "Dotted lines show closest-anode midpoint boundaries, **not** "
        "electron trajectories. Select an **orange cluster** to highlight "
        "its individual avalanches in the chronology below.")
    event=st.plotly_chart(fig,use_container_width=True,
                          key=f"penning-real-muon-geom-{index}",
                          on_select="rerun",selection_mode="points",
                          theme=None,
                          config={"displaylogo":False})
    selected=select_cluster(event,trace,geom["clusters"])
    choices=["Use plot click"]+[cl["index"] for cl in geom["clusters"]]
    manual=st.selectbox("Alternatively, select a cluster by index",
                        choices,key=f"penning-geometry-cluster-choice-{index}")
    if manual!="Use plot click":selected=int(manual)
    if selected is not None:
        c=next(v for v in geom["clusters"] if v["index"]==selected)
        st.info(
            f"**Selected TrackHeed cluster #{selected}:** "
            f"(u,v,w)=({c['u']:.3f}, {c['v']:.3f}, {c['w']:.3f}) mm; "
            f"{c['n_accepted']} accepted electrons "
            f"({c['n_original']} originally produced); "
            f"closest anode **u≈{c['anode_u']:.3f} mm** "
            "by the geometric midpoint rule. "
            "Its timing/avalanche row is highlighted below. "
            "A cluster with zero accepted electrons has no Garfield avalanche.")
    else:
        st.caption("No cluster selected yet; the full chronology is shown.")
    return selected

def selftest():
    # Grounded reproduction of C++ EntrancePhase and actual MWPC pitch.
    assert abs(entrance_anode_base(1.0,0)[0]-1.5)<1e-10
    assert abs(entrance_anode_base(1.0,1)[0]-0.)<1e-10
    assert abs(entrance_anode_base(4.4,0)[1]-(-.1))<1e-10
    assert cpp_round(-.5)==-1 and cpp_round(.5)==1
    sample=[{"index":0},{"index":2}]
    event={"selection":{"points":[
        {"curve_number":2,"point_index":1,"customdata":[2,3,2,0,0,1]}
    ]}}
    assert select_cluster(event,2,sample)==2
    assert select_cluster(event,1,sample) is None
    print("PASS physical 3-mm wire parity, coordinate phase and "
          "Plotly linked cluster-selection mapping; no output.")

if __name__=="__main__":
    selftest()

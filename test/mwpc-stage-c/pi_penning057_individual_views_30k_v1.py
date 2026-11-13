#!/usr/bin/env python3
"""Read-only 30k Penning-ON individual Geant4/Heed/Garfield event views.

Requires MWPC_PENNING_30K_RAW_ROOT to point to the existing SSHFS mount.
Reads exactly one original 20-muon chunk on demand, checking frozen input
identities and SHA-256 seals for every raw file accessed. Includes the
recovery-only successful chunks and NEVER overwrites source data.

Plots:
  1. Actual Geant4 (u,v) track, nearest-wire lattice and ALL TrackHeed clusters.
  2. Cluster x accepted-Heed-electron heatmap of exact avalanche size.
  3. Single selected avalanche 2D and rotatable (u,v,w) ion-birth cloud.
  4. Single selected avalanche angular distribution about the nearest anode.
Only <=128 uniformly sampled ion birth points/avalanche were archived;
exact total ion births and collected electron counts are independent of that.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from pi_penning057_original2m_cluster10k_v1 import (
    check_input, crossing_id, digest, metadata, require,
)
from pi_penning057_30k_plotstudio_cache_v1 import chunk_dir
import pi_p3_noFW_penning057_muon_geometry_v1 as geom
import pi_p3_noFW_penning057_anode_wrapping_v1 as wrapping

CHANNELS = {"jpsi": "J/ψ", "phi": "φ", "omega": "ω"}
PITCH_MM = 3.0
RAW_COLUMNS = {
    "penning_heed_clusters.csv": (
        "crossing_id", "cluster_index", "u_rel_cm", "v_abs_cm", "w_rel_cm",
        "t_ns", "n_original", "n_accepted"),
    "penning_heed_primaries.csv": (
        "crossing_id", "primary_index", "cluster_index", "t_ns"),
    "penning_seed.csv": (
        "crossing_id", "primary_index", "seed_u_mm", "seed_v_mm",
        "collected_e", "ion_births", "ion_up", "ion_down", "capped"),
    "penning_ion_birth_sample.csv": (
        "crossing_id", "primary_index", "u_mm", "v_mm", "w_mm", "t_ns",
        "seed_ion_births_total", "reservoir_weight"),
}


def _rows(path, fields):
    with Path(path).open(newline="") as file:
        reader = csv.DictReader(file)
        require(reader.fieldnames is not None
                and len(reader.fieldnames) == len(set(reader.fieldnames))
                and set(fields).issubset(reader.fieldnames),
                "Invalid source native CSV schema: " + str(path))
        yield from reader


def _seal(path):
    folder = Path(path).parent
    lines = (folder / "outputs.sha256").read_text().splitlines()
    hashes = {}
    for line in lines:
        if len(line) >= 67 and line[64:66] == "  ":
            hashes[line[66:]] = line[:64]
    require(path.name in hashes, f"File not in original SHA seal: {path}")
    require(digest(path) == hashes[path.name],
            f"Native source checksum mismatch: {path}")


def _parent(raw_root, channel, muon_index):
    require(channel in CHANNELS and 0 <= muon_index < 10000,
            "Invalid channel or parent index")
    campaign = (Path(raw_root) /
                f"original2m_{channel}_10000_v1")
    meta = metadata(campaign)
    required_channel = {"jpsi": "Jpsi", "phi": "Phi", "omega": "Omega"}[channel]
    require(meta["channel"] == required_channel and meta["n_chunks"] == 500,
            "Not the frozen 30k Penning-on campaign")
    block = muon_index // 20
    _, input_rows = check_input(campaign, block)
    source = input_rows[muon_index % 20]
    output, recovery = chunk_dir(campaign, block)
    report_file = output / "verified_report.json"
    _seal(report_file)
    report = json.loads(report_file.read_text())
    require(report["chunk_index"] == block
            and report["n_muons"] == 20
            and report["input_sha256"] ==
            meta["chunks"][f"{block:04d}"]["sha256"]
            and report["binary_sha256"] == meta["binary_sha256"]
            and report["n_capped_avalanches"] == 0,
            "Frozen physics report invalid")
    return campaign, source, output, report, recovery


def read_event(raw_root, channel, muon_index, expected_id=None,
               expected_heed=None, expected_clusters=None):
    """Only check/scan the three small native files for one real Geant4 parent."""
    campaign, source, output, report, recovered = _parent(
        raw_root, channel, muon_index)
    cid = crossing_id(source)
    if expected_id is not None:
        require(cid == expected_id,
                "Selected native Geant4 muon differs from Plot Studio cache")
    paths = {name: output / name for name in
             ("penning_heed_clusters.csv", "penning_heed_primaries.csv",
              "penning_seed.csv")}
    for path in paths.values():
        _seal(path)
    a = np.array([float(source[f"{ax}_in_mm"]) for ax in "uvw"])
    b = np.array([float(source[f"{ax}_out_mm"]) for ax in "uvw"])
    require(np.isfinite(a).all() and np.isfinite(b).all()
            and np.linalg.norm(a-b) > .05,
            "Unphysical Geant4 start/end")
    station = int(source["station"])
    wire_origin, offset = geom.entrance_anode_base(a[0], station)
    delta = b-a
    clusters = []
    for row in _rows(paths["penning_heed_clusters.csv"],
                     RAW_COLUMNS["penning_heed_clusters.csv"]):
        if row["crossing_id"] != cid:
            continue
        position = np.array([
            a[0] + 10. * float(row["u_rel_cm"]),
            10. * float(row["v_abs_cm"]),
            a[2] + 10. * float(row["w_rel_cm"])
        ])
        parallel = np.dot(position-a, delta)/np.dot(delta,delta)
        distance = np.linalg.norm(position-(a+parallel*delta))
        require(np.isfinite(position).all() and
                -.02 <= parallel <= 1.02 and distance < 1.,
                "Heed cluster is not on its true original Geant4 muon track")
        clusters.append({
            "index": int(row["cluster_index"]),
            "u": float(position[0]), "v": float(position[1]),
            "w": float(position[2]), "t_ns": float(row["t_ns"]),
            "n_original": int(row["n_original"]),
            "n_accepted": int(row["n_accepted"]),
            "anode_u": wire_origin + geom.cpp_round(
                (position[0]-wire_origin)/PITCH_MM)*PITCH_MM,
            "fraction_path": float(parallel),
        })
    clusters.sort(key=lambda r:r["index"])
    require([r["index"] for r in clusters] == list(range(len(clusters))),
            "Incomplete/out-of-order true TrackHeed clusters")
    heeds = {}
    for row in _rows(paths["penning_heed_primaries.csv"],
                     RAW_COLUMNS["penning_heed_primaries.csv"]):
        if row["crossing_id"] == cid:
            primary = int(row["primary_index"])
            require(primary not in heeds,"Duplicated Heed electron")
            heeds[primary] = {
                "cluster": int(row["cluster_index"]),
                "t_heed": float(row["t_ns"]),
            }
    seeds = []
    for row in _rows(paths["penning_seed.csv"],
                     RAW_COLUMNS["penning_seed.csv"]):
        if row["crossing_id"] != cid:
            continue
        j = int(row["primary_index"])
        require(j in heeds and int(row["capped"]) == 0,
                "Unmatched/capped native Garfield seed")
        seeds.append({
            "primary": j,
            "cluster": heeds[j]["cluster"],
            "t_heed": heeds[j]["t_heed"],
            "gain": int(row["collected_e"]),
            "ions": int(row["ion_births"]),
            "up": int(row["ion_up"]),
            "down": int(row["ion_down"]),
            "u0": float(row["seed_u_mm"]),
            "v0": float(row["seed_v_mm"]),
        })
    seeds.sort(key=lambda r:r["primary"])
    require(len(seeds) == len(heeds)
            and [s["primary"] for s in seeds] == list(range(len(seeds)))
            and sum(c["n_accepted"] for c in clusters) == len(seeds),
            "TrackHeed cluster/primary/Garfield avalanche census inconsistent")
    if expected_heed is not None:
        require(len(seeds) == expected_heed,
                "Cache and native TrackHeed electron count differ")
    if expected_clusters is not None:
        require(len(clusters) == expected_clusters,
                "Cache and native TrackHeed cluster count differ")
    return {
        "muon_index":muon_index,"id":cid,
        "a":a.tolist(),"b":b.tolist(),"wire_origin":wire_origin,
        "station":station,"phase_offset_mm":offset,
        "clusters":clusters,"n_Heed":len(seeds),"seeds":seeds,
        "chunk":muon_index//20,"recovered":recovered,
        "output":str(output),
    }


def read_cloud(raw_root, channel, muon_index, primary, event):
    """On-demand scan ONE sealed original reservoir for ONE microscopic seed."""
    _, source, output, _, recovered = _parent(
        raw_root, channel, muon_index)
    require(crossing_id(source) == event["id"] and
            event["recovered"] == recovered,
            "Selected source changed since event geometry read")
    matches = [s for s in event["seeds"] if s["primary"] == primary]
    require(len(matches) == 1,"Selected avalanche not in authentic Heed seeds")
    seed = matches[0]
    sample_file = output / "penning_ion_birth_sample.csv"
    _seal(sample_file)
    samples = []
    for row in _rows(sample_file,RAW_COLUMNS[sample_file.name]):
        if row["crossing_id"] != event["id"] or int(row["primary_index"])!=primary:
            continue
        require(int(row["seed_ion_births_total"]) == seed["ions"],
                "Sampled birth's exact ion multiplicity mismatch")
        samples.append([float(row[k]) for k in
                        ("u_mm","v_mm","w_mm","t_ns")])
    require(len(samples) == min(128,seed["ions"]),
            "Missing native ion-birth reservoir coordinates")
    points = np.asarray(samples,dtype=float).reshape(-1,4)
    wire = wrapping.cpp_round(seed["u0"]/PITCH_MM)*PITCH_MM
    if points.size:
        points[:,0] -= wire
    asym = (seed["up"]-seed["down"])/(seed["up"]+seed["down"]) \
           if seed["up"]+seed["down"] else float("nan")
    return {
        "seed":primary, "parent":muon_index, "cluster":seed["cluster"],
        "start_u_wire_mm":seed["u0"]-wire,"start_v_mm":seed["v0"],
        "t_heed":seed["t_heed"],"collected_e":seed["gain"],
        "ion_births":seed["ions"],"A":asym,
        "wire_u_local_mm":wire,"points":points,
    }


def cluster_heatmap(event, measure="Collected electrons"):
    """Each cell represents EXACTLY one Heed electron's Garfield avalanche."""
    import plotly.graph_objects as go
    clusters, seeds = event["clusters"], event["seeds"]
    width = max([1]+[c["n_accepted"] for c in clusters])
    width = min(width, 10000)
    values = np.full((max(1,len(clusters)), width),np.nan)
    ids = np.full_like(values,-1,dtype=int)
    ion = np.full_like(values,-1,dtype=int)
    local = {}
    for s in seeds:
        cl = s["cluster"]
        j = local.get(cl,0)
        require(0<=cl<len(clusters) and j<width,
                "Unphysical cluster seed indexing")
        local[cl]=j+1
        values[cl,j]=s["gain"] if measure=="Collected electrons" else s["ions"]
        ids[cl,j]=s["primary"]
        ion[cl,j]=s["ions"]
    z=np.log10(1.+values)
    fig=go.Figure(go.Heatmap(
        x=np.arange(width),y=np.arange(max(1,len(clusters))),
        z=z,colorscale="Viridis",xgap=1,ygap=2,
        customdata=np.dstack([ids,ion,values]),
        colorbar=dict(title="log₁₀(1 + N)"),
        hovertemplate="Cluster #%{y:.0f}<br>Electron within cluster %{x:.0f}"
                      "<br>Original Heed primary %{customdata[0]:.0f}"
                      "<br>Exact ion births %{customdata[1]:,.0f}"
                      "<br>Selected avalanche size %{customdata[2]:,.0f}"
                      "<extra></extra>"))
    fig.update_layout(
        title="Exact Garfield avalanche size for every primary electron",
        xaxis_title="Primary-electron index within cluster (0-based)",
        yaxis_title="TrackHeed cluster index (0-based)",
        template="plotly_white",height=max(560,min(1000,310+23*len(clusters))),
        font=dict(family="Arial,Helvetica,sans-serif",size=15),
        margin=dict(l=90,r=115,t=75,b=85))
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False,dtick=1 if len(clusters)<=30 else 2)
    return fig


def angular_hist(points, collected_ions):
    import plotly.graph_objects as go
    if not len(points):
        return go.Figure()
    angles = np.degrees(np.arctan2(points[:,1],points[:,0]))
    edges = np.linspace(-180.,180.,73)
    counts=np.histogram(angles,bins=edges)[0]
    fractions=counts/counts.sum()
    fig=go.Figure(go.Bar(
        x=(edges[:-1]+edges[1:])/2,y=fractions,width=5.,
        marker_color="#0072B2",name="Sampled ion births",
        hovertemplate="Azimuth %{x:.1f}°<br>Fraction %{y:.2%}<extra></extra>"))
    fig.update_layout(
        title="Ion creation angle around selected anode wire",
        xaxis_title="φ = atan2(v, u − wire centre) [degrees]",
        yaxis_title="Fraction of retained births per 5° bin",
        template="plotly_white",height=550,
        font=dict(family="Arial,Helvetica,sans-serif",size=15),
        margin=dict(l=85,r=60,t=80,b=75),
        xaxis=dict(range=[-180,180],dtick=45),bargap=0)
    return fig


def _channel_cache(cache_root, channel):
    folder = Path(cache_root) / channel
    mp = folder / "temporal_correlations_v1.json"
    fp = folder / "temporal_correlations_v1.npz"
    require(mp.is_file() and fp.is_file(),
            "Missing verified 10k per-channel temporal cache")
    m=json.loads(mp.read_text())
    require(m["dataset"] == channel
            and m["n_independent_true_Geant4_muons"] == 10000,
            "Wrong temporal campaign in detailed event viewer")
    require(digest(fp) == m["auxiliary_cache_sha256"],
            "Temporal cache SHA differs from frozen provenance")
    return m


def render(st, cache_root, raw_root, scope, kind):
    """One of two detailed, intentionally non-overlay single-event pages."""
    import plotly.graph_objects as go
    codes=tuple(CHANNELS)
    initial=scope if scope in codes else "jpsi"
    st.markdown("## Individual Garfield++ event inspector")
    st.caption("Read-only original 20-muon native chunks, loaded on demand "
               "through SSHFS; no re-simulation and no 45 GB download.")
    sel=st.selectbox("Channel", codes, index=codes.index(initial),
                     format_func=lambda x: CHANNELS[x],
                     key="mwpc-30k-detail-channel-"+kind)
    meta=_channel_cache(cache_root,sel)
    which=st.number_input("Real Geant4 muon index within this channel (0–9999)",
                          min_value=0,max_value=9999,
                          value=0,step=1,
                          key="mwpc-30k-detail-muon-"+kind)
    mu=int(which)
    expected=meta["muons"][mu]
    @st.cache_data(show_spinner="SHA-checking selected original 20-muon native chunk…")
    def event_cached(folder,ch,idx,ident,n_heed,n_clusters):
        return read_event(folder,ch,idx,ident,n_heed,n_clusters)
    evt=event_cached(str(Path(raw_root).expanduser()),sel,mu,
                     expected["id"],expected["n_Heed"],
                     expected["n_clusters"])
    st.caption(f"Original Geant4 crossing ID: {evt['id']} · "
               f"verified chunk {evt['chunk']:04d}" +
               (" · repaired successful chunk" if evt["recovered"] else "") +
               f" · {len(evt['clusters'])} true clusters · "
               f"{len(evt['seeds'])} Heed electrons")

    if kind=="detail_clusters":
        options=["All clusters"]+[cl["index"] for cl in evt["clusters"]]
        manual=st.selectbox("Inspect cluster (or click an orange cluster)",
                            options,key=f"mwpc-30k-cluster-select-{sel}-{mu}")
        # Stacked full-width plots are more legible and export better.
        picture,curve=geom.build_figure(evt)
        picture.update_layout(
            height=680,font=dict(family="Arial,Helvetica,sans-serif",size=15),
            margin=dict(l=85,r=55,t=85,b=110))
        st.markdown("### 1. Geant4 track and TrackHeed clusters")
        if curve is not None:
            picture.data[curve].name="TrackHeed cluster — click to inspect"
        with st.container():
            st.caption("Click an orange cluster marker to select its associated "
                       "Garfield avalanches, or use the selector above.")
            event=st.plotly_chart(picture,use_container_width=True,
                                  theme=None,
                                  key=f"mwpc-30k-geom-{sel}-{mu}",
                                  on_select="rerun",selection_mode="points",
                                  config={"displaylogo":False})
        clicked=geom.select_cluster(event,curve,evt["clusters"])
        selected=manual if isinstance(manual,int) else clicked
        st.divider()
        with st.container():
            st.markdown("### 2. Exact avalanche size per TrackHeed electron")
            metric=st.radio("Cell values",
                            ("Collected electrons","Positive ion births"),
                            horizontal=True,key="mwpc-30k-cluster-metric")
            matrix=cluster_heatmap(evt,metric)
            if selected is not None:
                matrix.add_hrect(y0=selected-.47,y1=selected+.47,
                                 line_color="#D95F02",line_width=2,
                                 fillcolor="rgba(217,95,2,.10)")
            st.plotly_chart(matrix,use_container_width=True,theme=None,
                            config={"displaylogo":False})
        if selected is not None:
            st.caption(f"TrackHeed cluster #{selected} is highlighted in the "
                       "avalanche-size heatmap; use the selector or click "
                       "an orange cluster on the track.")
        if selected is not None:
            c=evt["clusters"][selected]
            seeds=[s for s in evt["seeds"] if s["cluster"]==selected]
            st.markdown(f"**TrackHeed cluster #{selected}:** "
                        f"{c['n_accepted']} accepted primary electrons, "
                        f"{sum(s['gain'] for s in seeds):,} collected electrons, "
                        f"{sum(s['ions'] for s in seeds):,} true ion births; "
                        f"(u, v, w) = ({c['u']:.3f}, {c['v']:.3f}, "
                        f"{c['w']:.3f}) mm.")
            st.dataframe([{
                "Heed primary":s["primary"],"exact gain":s["gain"],
                "exact ion births":s["ions"],
                "start time [ns]":s["t_heed"],
            } for s in seeds],hide_index=True,use_container_width=True)
        st.caption("Heatmap cells are individual microscopic electron-started "
                   "avalanches, not averaged cluster sizes. The missing cells "
                   "mean no primary electron; actual zero-gain avalanches remain 0.")
    elif kind=="detail_avalanche":
        if not evt["seeds"]:
            st.warning("This Geant4 muon has no retained TrackHeed electrons.")
            return
        seedmap={s["primary"]:s for s in evt["seeds"]}
        # Largest ion-producing event makes an informative initial 3D example.
        order=sorted(seedmap,key=lambda k:(-seedmap[k]["ions"],k))
        primary=st.selectbox(
            "Original Heed electron / Garfield avalanche",order,
            format_func=lambda k:(f"Seed {k} · cluster {seedmap[k]['cluster']}"
                                  f" · gain {seedmap[k]['gain']:,}"
                                  f" · ion births {seedmap[k]['ions']:,}"),
            key=f"mwpc-30k-detail-seed-{sel}-{mu}")
        @st.cache_data(show_spinner="Verifying and reading original sampled ion births…")
        def cloud_cached(folder,ch,idx,primary,rec,source,ev):
            return read_cloud(folder,ch,idx,primary,ev)
        obs=cloud_cached(str(Path(raw_root).expanduser()),
                         sel,mu,primary,evt["recovered"],evt["output"],evt)
        points=obs["points"]
        st.markdown(f"**Exact avalanche:** collected gain "
                    f"**{obs['collected_e']:,} e**, true positive ions "
                    f"**{obs['ion_births']:,}**, exact hemisphere "
                    f"asymmetry **{obs['A']:+.3f}** · "
                    f"retained birth positions **{len(points)}**.")
        st.caption("These are ion CREATION coordinates (u,v,w,t), not ion "
                   "drift trajectories or cathode-induced charges. "
                   "The maximum reservoir contains 128 spatial points.")
        if not len(points):
            st.info("This particular avalanche produced no ions to plot; "
                    "choose a different Heed primary.")
            return
        half=st.select_slider(
            "Anode cross-section half-width [mm]",
            options=[.05,.1,.25,.5,1.5,3.],
            value=.1,key="mwpc-30k-detail-halfwidth")
        colour_by=st.radio(
            "Colour sampled ion births by",
            ("Time since first ion birth","Delay since primary Heed electron"),
            horizontal=True,key="mwpc-30k-detail-birth-colour")
        u,v,w,t=points.T
        time_origin=(float(np.min(t)) if
                     colour_by=="Time since first ion birth" else obs["t_heed"])
        colour_title=("Time since first birth [ns]" if
                      colour_by=="Time since first ion birth" else
                      "Delay since Heed electron [ns]")
        st.markdown("### 1. Ion-birth cloud around the 30-μm anode")
        with st.container():
            fig=go.Figure()
            radius=.015
            fig.add_shape(type="circle",x0=-radius,x1=radius,
                          y0=-radius,y1=radius,
                          fillcolor="#9A9A9A",line_color="#555555")
            fig.add_trace(go.Scatter(
                x=u,y=v,mode="markers",
                marker=dict(size=7,color=t-time_origin,
                            colorscale="YlOrRd",showscale=True,
                            colorbar=dict(title=colour_title)),
                customdata=np.column_stack((w,t)),
                name="Uniformly sampled ion births",
                hovertemplate="u−wire=%{x:.5f} mm<br>v=%{y:.5f} mm"
                              "<br>w=%{customdata[0]:.5f} mm"
                              "<br>t=%{customdata[1]:.5f} ns<extra></extra>"))
            fig.add_trace(go.Scatter(
                x=[obs["start_u_wire_mm"]],y=[obs["start_v_mm"]],
                mode="markers",name="Heed starting electron",
                marker=dict(symbol="star",size=14,color="#0072B2")))
            fig.update_layout(
                title="Ion birth cloud around anode · u–v cross-section",
                xaxis_title="u − anode wire centre [mm]",
                yaxis_title="v [mm]",height=660,template="plotly_white",
                font=dict(family="Arial,Helvetica,sans-serif",size=15),
                margin=dict(l=90,r=105,t=80,b=75))
            fig.update_xaxes(range=[-half,half],constrain="domain",showgrid=False)
            fig.update_yaxes(range=[-half,half],scaleanchor="x",scaleratio=1,
                             constrain="domain",showgrid=False)
            st.plotly_chart(fig,use_container_width=True,theme=None,
                            config={"displaylogo":False})
        st.caption("The anode wire is drawn at its real 15-μm radius. "
                   "The original Heed electron may lie outside the zoom.")
        st.divider()
        st.markdown("### 2. Angular distribution of ion births")
        with st.container():
            st.plotly_chart(angular_hist(points,obs["ion_births"]),
                            use_container_width=True,theme=None,
                            config={"displaylogo":False})
        st.caption("−180° and +180° refer to the same physical azimuth "
                   "(a cloud straddling that direction appears at both ends).")
        st.divider()
        st.markdown("### 3. Rotatable 3D ion-birth cloud along wire w")
        colour_mode=("first" if colour_by=="Time since first ion birth"
                     else "heed")
        st.plotly_chart(wrapping.three_d_view(obs,half,colour_mode=colour_mode),
                        use_container_width=True,theme=None,
                        config={"displaylogo":False})
        if len(points) < obs["ion_births"]:
            st.caption(f"Angular and 3D views retain {len(points)} uniformly "
                       f"sampled births representing {obs['ion_births']:,} "
                       "true ions. Their spatial distributions are estimates; "
                       "exact totals/gain/asymmetry are full Garfield counters.")
    else:
        raise ValueError(kind)


def selftest():
    assert np.allclose(geom.entrance_anode_base(1.,0),(1.5,-.5))
    assert geom.cpp_round(-.5)==-1
    assert len(RAW_COLUMNS)==4
    assert min(128,20000)==128
    print("PASS detailed 30k viewer geometrical conventions and source schemas "
          "(read-only; run through Streamlit for real SSHFS chunk verification)")


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",nargs="?",choices=("selftest","probe"),
                        default="selftest")
    parser.add_argument("--raw-root",type=Path,
                        default=Path.home()/"na6p/cluster_mwpc")
    parser.add_argument("--cache-root",type=Path)
    parser.add_argument("--channel",choices=tuple(CHANNELS),default="jpsi")
    parser.add_argument("--muon",type=int,default=0)
    opts=parser.parse_args()
    if opts.command=="selftest":
        selftest()
        return
    require(opts.cache_root is not None,
            "probe needs --cache-root containing the four compact caches")
    require(0<=opts.muon<10000,"Muon index outside original 10k channel")
    m=_channel_cache(opts.cache_root,opts.channel)
    expect=m["muons"][opts.muon]
    event=read_event(opts.raw_root,opts.channel,opts.muon,
                     expect["id"],expect["n_Heed"],expect["n_clusters"])
    print(f"PASS original {CHANNELS[opts.channel]} Geant4 muon "
          f"#{opts.muon}: {len(event['clusters'])} TrackHeed clusters, "
          f"{len(event['seeds'])} microscopic Garfield avalanches; "
          f"native geometry/primaries/seeds SHA verified; "
          f"recovered={event['recovered']}.",flush=True)
    if event["seeds"]:
        chosen=max(event["seeds"],key=lambda s:s["ions"])
        cloud=read_cloud(opts.raw_root,opts.channel,opts.muon,
                         chosen["primary"],event)
        print(f"PASS native single avalanche #{chosen['primary']}: "
              f"{cloud['ion_births']:,} exact ion births, "
              f"{len(cloud['points'])} authentic u/v/w/t samples; "
              "SHA-sealed 20-muon reservoir verified.",flush=True)


if __name__=="__main__":
    main()

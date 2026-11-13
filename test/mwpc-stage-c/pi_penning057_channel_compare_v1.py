#!/usr/bin/env python3
"""Read-only, independent-channel comparison for the 30k MWPC Plot Studio.

Show 1--3 separate Geant4 parent samples together. Never treat the pooled ALL
sample as a fourth independent physics population; raw-count normalization is
explicit, while probability normalizes each selected channel independently.
2D quantities are faceted with a shared colour scale (not overlaid).
"""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import numpy as np

CODES = ("jpsi", "phi", "omega")
LABELS = {"jpsi": "J/ψ", "phi": "φ", "omega": "ω"}
COLORS = {"jpsi": "#D95F02", "phi": "#0072B2", "omega": "#009E73"}
HISTOGRAMS = {
    "gain_all", "gain_pos", "angle_eq", "angle_ion",
    "asym_seed", "asym_muon", "time", "span",
}
HEATMAPS = {
    "ion_xy", "gain_xy", "time_phi", "t10_start",
    "t10_gain_position", "wrap_v0_angle", "wrap_u0_angle",
}
CORRELATIONS = {
    "gain_asym", "t10_gain", "t10_asym", "t10_angles",
    "wrap_v0_A", "wrap_v0_opposite", "wrap_dispersion",
}
SUPPORTED = HISTOGRAMS | HEATMAPS | CORRELATIONS


def _normalise(counts, edges, mode):
    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(edges, dtype=float)
    if mode == "Raw counts / weighted counts":
        return counts.copy()
    if counts.sum() <= 0:
        return counts.copy()
    p = counts / counts.sum()
    if mode == "Probability density":
        return p / np.diff(edges)
    return p


def _shared_edges(kind, datasets, opts):
    if kind in {"gain_all", "gain_pos"}:
        maximum = max(int(np.max(d["gain"])) for d in datasets.values())
        w = int(opts["gain_width"])
        edges = np.arange(1, max(1, maximum) + w + 1, w, dtype=float)
        return np.r_[0., edges] if kind == "gain_all" else edges
    if kind in {"angle_eq", "angle_ion"}:
        factor = int(opts["angular_deg"]) // 5
        return next(iter(datasets.values()))["phi_edges"][::factor]
    if kind in {"asym_seed", "asym_muon"}:
        return np.linspace(-1., 1., int(opts["asym_bins"]) + 1)
    if kind == "time":
        return next(iter(datasets.values()))["t_edges"]
    if kind == "span":
        maximum = max(float(np.nanmax(d["seed_span_ns"])) for d in datasets.values())
        return np.linspace(0., max(1.e-6, maximum), 201)
    raise ValueError("Not a histogram " + kind)


def _hist_values(kind, d, edges, opts):
    if kind in {"gain_all", "gain_pos"}:
        values = d["gain"]
        if kind == "gain_pos":
            values = values[values > 0]
        counts = np.histogram(values, bins=edges)[0].astype(float)
        if int(counts.sum()) != len(values):
            raise ValueError("Gain histogram silently discarded observations")
        return counts, values
    if kind in {"angle_eq", "angle_ion"}:
        raw = d["angle_equal"] if kind == "angle_eq" else d["angle_ion"]
        factor = int(opts["angular_deg"]) // 5
        return raw.reshape(-1, factor).sum(axis=1), None
    if kind in {"asym_seed", "asym_muon"}:
        raw = d["seed_A"] if kind == "asym_seed" else d["muon_A"]
        v = raw[np.isfinite(raw)]
        return np.histogram(v, bins=edges)[0].astype(float), v
    if kind == "time":
        return d["t_weights"].astype(float), None
    if kind == "span":
        v = d["seed_span_ns"]
        v = v[np.isfinite(v)]
        return np.histogram(v, bins=edges)[0].astype(float), v
    raise ValueError(kind)


def _group_time(counts, edges, group, logx):
    # Exactly the same cuts and zero-time treatment for EVERY channel.
    if logx:
        counts, edges = counts[1:], edges[1:]
    cuts = list(range(0, len(counts), int(group)))
    if cuts[-1] != len(counts):
        cuts.append(len(counts))
    cuts = np.asarray(cuts, dtype=int)
    return np.add.reduceat(counts, cuts[:-1]), edges[cuts]


def _rgb_alpha(hex_color, alpha=.16):
    return ("rgba(" + ",".join(str(int(hex_color[i:i+2], 16))
            for i in (1, 3, 5)) + f",{alpha:g})")


def _write_csv(st, rows, filename):
    if not rows:
        st.info("For this plot, download the interactive HTML or original cache.")
        return
    columns = list(dict.fromkeys(k for row in rows for k in row))
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    st.download_button("Download all selected channels (.csv)",
                       stream.getvalue(), file_name=filename,
                       mime="text/csv", use_container_width=True)


def show(st, cache_root, selected, kind, opts, load_cache, temporal, wrapping,
         plot_explanation):
    """Render a complete independent-channel comparison with shared controls."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    codes = [c for c in CODES if c in selected]
    if not codes or kind not in SUPPORTED:
        st.warning("Select 1–3 channels and an aggregate plot.")
        return
    root = Path(cache_root).expanduser()
    datasets = {}
    manifests = {}
    auxiliary = {}
    for code in codes:
        folder = root / code
        npz, manifest = folder / "dashboard_data.npz", folder / "provenance.json"
        if not npz.is_file() or not manifest.is_file():
            st.error(f"Missing {code} cache; all three channels must be prepared.")
            return
        d, m = load_cache(str(npz), npz.stat().st_mtime_ns,
                          str(manifest), manifest.stat().st_mtime_ns)
        if m.get("dataset") != code or m.get("n_muons") != 10000:
            st.error(f"Unexpected {code} cache provenance; abort comparison.")
            return
        datasets[code], manifests[code] = d, m
        if kind.startswith("t10_") or kind.startswith("wrap_"):
            tp, mp = (folder / "temporal_correlations_v1.npz",
                      folder / "temporal_correlations_v1.json")
            if not tp.is_file() or not mp.is_file():
                st.error(f"Missing {code} temporal cache.")
                return
            # Read from this channel's explicit paths, not the legacy module
            # globals (which are overwritten when the selected scope changes).
            @st.cache_data(show_spinner="Loading channel-specific time correlations…")
            def read_temporal(path, mtime, meta_path, meta_mtime):
                with np.load(path, allow_pickle=False) as z:
                    arrays = {k: z[k] for k in z.files}
                return arrays, json.loads(Path(meta_path).read_text())
            auxiliary[code] = read_temporal(
                str(tp), tp.stat().st_mtime_ns, str(mp), mp.stat().st_mtime_ns)

    labels = [LABELS[c] for c in codes]
    st.markdown("#### Compare " + " vs ".join(labels))
    st.caption("Independent 10,000-muon samples · curves are not averaged "
               "across channels · ALL is not plotted as a fourth sample. "
               "Probability and density normalise each channel separately.")
    for box, code in zip(st.columns(len(codes)), codes):
        g = datasets[code]["gain"]
        box.metric(f"{LABELS[code]} · full avalanches", f"{len(g):,}",
                   f"mean gain {np.mean(g):,.1f}", delta_color="off")

    rows = []
    notes = []
    if kind in HEATMAPS:
        panels = []
        for code in codes:
            d, m = datasets[code], manifests[code]
            fig = go.Figure()
            if kind == "ion_xy":
                z = d["xy_ion"].astype(float)
                x = (d["u_edges"][1:] + d["u_edges"][:-1]) / 2
                y = (d["v_edges"][1:] + d["v_edges"][:-1]) / 2
                caption = "Ion-weighted creation counts, not an electrode signal"
            elif kind == "gain_xy":
                count = d["start_n"]
                z = np.divide(d["start_sum"], count,
                              out=np.full(count.shape, np.nan, dtype=float),
                              where=count >= 5)
                x = (d["u_edges"][1:] + d["u_edges"][:-1]) / 2
                y = (d["v_edges"][1:] + d["v_edges"][:-1]) / 2
                caption = "Mean gain per initial electron; cells with <5 entries hidden"
            elif kind == "time_phi":
                z = d["t_phi"].astype(float)
                x = (d["phi_edges"][1:] + d["phi_edges"][:-1]) / 2
                y = (d["t_phi_edges"][1:] + d["t_phi_edges"][:-1]) / 2
                caption = "Ion-weighted creation counts per time/angle cell"
            else:
                td, tm = auxiliary[code]
                controls = dict(opts["time_controls"] if kind.startswith("t10_")
                                else opts["wrap_controls"])
                controls.update(color=COLORS[code], legend=LABELS[code],
                                line_width=opts["line_width"],
                                zscale=opts["zscale"])
                if kind.startswith("t10_"):
                    _, note = temporal.render(kind, fig, td, tm, controls, st)
                else:
                    wrapping.N_BOOTSTRAP_PARENTS = int(m["n_muons"])
                    _, note, _ = wrapping.render(kind, fig, td, tm, controls, st)
                if not fig.data:
                    st.error(f"No {code} heatmap was produced.")
                    return
                trace = fig.data[0]
                z = np.asarray(trace.z, dtype=float)
                x, y = np.asarray(trace.x), np.asarray(trace.y)
                caption = note
                # The specialized renderer already applies log10 scaling.
            if opts["zscale"] == "Log" and kind in {
                "ion_xy", "gain_xy", "time_phi"
            }:
                with np.errstate(divide="ignore", invalid="ignore"):
                    z = np.where(z > 0, np.log10(z), np.nan)
            panels.append((code, x, y, z))
            notes.append(f"{LABELS[code]}: {caption}")

        # Shared colour limits are essential for meaningful side-by-side maps.
        finite = [z[np.isfinite(z)] for _, _, _, z in panels]
        nonempty=[v for v in finite if v.size]
        if not nonempty:
            st.warning("All selected heatmaps are empty under these cuts.")
            return
        values = np.concatenate(nonempty)
        zmin, zmax = float(values.min()), float(values.max())
        fig = make_subplots(rows=1, cols=len(codes),
                            subplot_titles=[LABELS[c] for c in codes],
                            shared_xaxes=True, shared_yaxes=True,
                            horizontal_spacing=.055)
        for index, (code, x, y, z) in enumerate(panels, 1):
            fig.add_trace(go.Heatmap(
                x=x, y=y, z=z, zmin=zmin, zmax=zmax,
                colorscale="YlOrBr", showscale=index == len(panels),
                colorbar=dict(title=opts["label"]) if index == len(panels)
                else None,
                hovertemplate=f"{LABELS[code]}<br>x=%{{x:.4g}}"
                              "<br>y=%{y:.4g}<br>z=%{z:.5g}<extra></extra>"
            ), row=1, col=index)
        fig.update_layout(height=650, showlegend=False)
        st.caption("2D quantities are displayed side by side with the SAME "
                   "colour scale; overlaying heatmaps would conceal structure.")
    else:
        fig = go.Figure()
        if kind in HISTOGRAMS:
            edges = _shared_edges(kind, datasets, opts)
            for code in codes:
                counts, observations = _hist_values(
                    kind, datasets[code], edges, opts)
                use_edges = edges
                if kind == "time":
                    counts, use_edges = _group_time(
                        counts, edges, opts["time_group"],
                        opts["xscale"] == "Log")
                    if opts["xscale"] == "Log":
                        mass = float(datasets[code]["t_weights"][0])
                        total = float(np.sum(datasets[code]["t_weights"]))
                        notes.append(f"{LABELS[code]}: excluded exactly-zero "
                                     f"time bin ({mass/total:.3%}); renormalised")
                values = _normalise(counts, use_edges, opts["norm"])
                if opts["norm"] == "Raw counts / weighted counts":
                    notes.append(f"{LABELS[code]}: raw count/ion-weight scale; "
                                 "amplitude differences also reflect sample size")
                color = COLORS[code]
                if opts["plot_style"] == "Bars":
                    fig.add_trace(go.Bar(
                        x=(use_edges[:-1]+use_edges[1:])/2,
                        y=values, width=np.diff(use_edges), name=LABELS[code],
                        marker_color=color, opacity=.43))
                else:
                    yy=np.repeat(values,2).astype(float)
                    if opts["yscale"] == "Log":
                        # Missing/empty bins must not be joined by misleading
                        # lines on a logarithmic probability axis.
                        yy[yy<=0]=np.nan
                    fig.add_trace(go.Scatter(
                        x=np.repeat(use_edges,2)[1:-1],
                        y=yy,
                        mode="lines", name=LABELS[code],
                        line=dict(color=color, width=opts["line_width"]),
                        hovertemplate=f"{LABELS[code]}<br>x=%{{x:.4g}}"
                                      "<br>height=%{y:.5g}<extra></extra>"))
                if opts["add_stats"] and observations is not None:
                    o = observations[np.isfinite(observations)]
                    if o.size:
                        fig.data[-1].name += (
                            f" · mean={np.mean(o):,.2f}, RMS={np.std(o):,.2f}")
                for lo, hi, v in zip(use_edges[:-1], use_edges[1:], values):
                    rows.append(dict(channel=LABELS[code], bin_lo=lo,
                                     bin_hi=hi, value=v))
            if opts["plot_style"] == "Bars":
                fig.update_layout(barmode="overlay")
            if kind == "span":
                notes.append("Time spans are lower bounds from at most 128 "
                             "sampled ion births per avalanche.")
        elif kind == "gain_asym":
            max_gain = max(int(np.max(d["gain"])) for d in datasets.values())
            shared = np.linspace(0., max(1., np.log10(max_gain)), 31)
            edges = 10. ** shared
            for code in codes:
                d = datasets[code]
                g = d["gain"].astype(float)
                A = d["seed_absA"].astype(float)
                mask = (g > 0) & np.isfinite(A)
                lx = np.log10(g[mask])
                n = np.histogram(lx, bins=shared)[0]
                total = np.histogram(lx, bins=shared, weights=A[mask])[0]
                avg = np.divide(total,n,out=np.full_like(total,np.nan),
                                where=n>=20)
                fig.add_trace(go.Scatter(
                    x=np.sqrt(edges[:-1]*edges[1:]), y=avg,
                    mode="lines+markers", name=LABELS[code],
                    line=dict(color=COLORS[code], width=opts["line_width"]),
                    customdata=n,
                    hovertemplate="Gain=%{x:,.1f}<br>mean |A|=%{y:.4f}"
                                  "<br>N=%{customdata}<extra></extra>"))
                for lo, hi, v, nn in zip(edges[:-1],edges[1:],avg,n):
                    rows.append(dict(channel=LABELS[code], bin_lo=lo,
                                     bin_hi=hi, mean_absA=v, n_avalanches=nn))
            notes.append("Shared gain bin edges; bins with fewer than 20 "
                         "avalanches in each channel are hidden.")
        else:
            for code in codes:
                td, tm = auxiliary[code]
                controls = dict(opts["time_controls"] if kind.startswith("t10_")
                                else opts["wrap_controls"])
                controls.update(color=COLORS[code], line_width=opts["line_width"],
                                legend=LABELS[code], add_stats=opts["add_stats"])
                before = len(fig.data)
                if kind.startswith("t10_"):
                    output, note = temporal.render(kind,fig,td,tm,controls,st)
                else:
                    wrapping.N_BOOTSTRAP_PARENTS = int(manifests[code]["n_muons"])
                    output, note, _ = wrapping.render(kind,fig,td,tm,controls,st)
                if len(fig.data) == before:
                    notes.append(f"{LABELS[code]}: {note}")
                    continue
                for j in range(before, len(fig.data)):
                    trace = fig.data[j]
                    trace.legendgroup = code
                    if trace.showlegend is not False:
                        trace.name = (LABELS[code] +
                            (" · early" if kind == "t10_angles" and j == before
                             else " · late" if kind == "t10_angles"
                             else ""))
                    if trace.type == "scatter":
                        if trace.line is not None and trace.line.width != 0:
                            trace.line.color = COLORS[code]
                            if kind == "t10_angles" and j > before:
                                trace.line.dash = "dash"
                        if trace.fill:
                            trace.fillcolor = _rgb_alpha(COLORS[code])
                if output:
                    keys = list(output)
                    for values in zip(*(output[k] for k in keys)):
                        rows.append(dict(channel=LABELS[code],
                                         **dict(zip(keys,values))))
                notes.append(f"{LABELS[code]}: {note}")
            if kind in ("wrap_v0_A", "wrap_v0_opposite","wrap_dispersion"):
                notes.append("95% bands bootstrap independent parent muons "
                             "WITHIN each channel, not individual electrons.")
            if kind in ("t10_gain","t10_asym"):
                notes.append("These time-correlation curves use separate "
                             "equal-occupancy time bins for each channel.")
            if kind == "t10_angles":
                notes.append("Solid lines: early births; dashed lines: late "
                             "births. Quantile groups are selected per channel.")

    x_log = opts["xscale"] == "Log"
    y_log = opts["yscale"] == "Log"
    if x_log and kind in {"gain_all","angle_eq","angle_ion",
                           "asym_seed","asym_muon"}:
        st.warning("X-log scale not valid for signed angles or the zero-gain bin.")
        x_log = False
    if x_log and kind not in HEATMAPS:
        fig.update_xaxes(type="log")
    if y_log:
        fig.update_yaxes(type="log")
    if (kind == "time" and not x_log
        and opts["time_range"] == "Central 99.9% of ion births"):
        limits = []
        for d in datasets.values():
            count = d["t_weights"]
            e = d["t_edges"]
            idx = min(len(count)-1, int(np.searchsorted(
                np.cumsum(count), .999*float(np.sum(count)))))
            limits.append(e[idx+1])
        fig.update_xaxes(range=[0., max(limits)])
        st.caption("Time X-axis includes at least the central 99.9% of "
                   "weighted births in EVERY selected channel.")
    for axis, mini, maxi, log in (
        ("x", opts["xmin"], opts["xmax"], x_log),
        ("y", opts["ymin"], opts["ymax"], y_log),
    ):
        try:
            lo = float(mini) if mini.strip() else None
            hi = float(maxi) if maxi.strip() else None
        except ValueError:
            st.warning(f"Invalid {axis.upper()} axis limit; ignored.")
            continue
        if lo is None and hi is None:
            continue
        if lo is None or hi is None:
            st.warning(f"Enter BOTH {axis.upper()} axis bounds for comparison.")
            continue
        if hi <= lo or (log and lo <= 0):
            st.warning(f"Invalid {axis.upper()} range, ignoring custom limits.")
            continue
        target = [float(np.log10(lo)),float(np.log10(hi))] if log else [lo,hi]
        if axis == "x":
            fig.update_xaxes(range=target)
        else:
            fig.update_yaxes(range=target)

    positions = {"Top right":(.99,.99,"right","top"),
                 "Top left":(.01,.99,"left","top"),
                 "Bottom right":(.99,.01,"right","bottom"),
                 "Bottom left":(.01,.01,"left","bottom")}
    lx,ly,lxa,lya = positions[opts["legpos"]]
    fig.update_layout(
        title=dict(text=opts["title"], font=dict(size=opts["title_font"])),
        xaxis_title=opts["xlab"], yaxis_title=opts["ylab"],
        template="plotly_white", height=650,
        font=dict(family="Arial, Helvetica, sans-serif",
                  size=opts["label_font"]),
        margin=dict(l=80,r=45,t=85,b=76),
        showlegend=opts["show_legend"] if kind not in HEATMAPS else False,
        legend=dict(x=lx,y=ly,xanchor=lxa,yanchor=lya,
                    bgcolor="rgba(255,255,255,.85)"),
        hovermode="closest",
    )
    if opts["show_interpretation"]:
        st.markdown("**How to read the comparison:** Each colour is one "
                    "independent original-Geant4 channel. The 10,000-muon "
                    "samples must not be interpreted as equal physical "
                    "production rates; no weighting by meson yield is implied.")
        if kind in HISTOGRAMS:
            st.caption("Probability or density uses the same physical bins "
                       "and is normalized independently for each channel.")
        elif kind in HEATMAPS:
            st.caption("Heatmaps have a common numerical colour scale; "
                       "their coordinates and units are those of the selected plot.")
        else:
            st.caption("Means, bootstrap intervals and time windows are "
                       "calculated using each channel's original events.")
        if kind not in HEATMAPS and kind in HISTOGRAMS | {"gain_asym"}:
            with st.expander("Channel-by-channel physical interpretation"):
                for code in codes:
                    st.markdown(f"**{LABELS[code]}**")
                    what, how, observed = plot_explanation(
                        kind,datasets[code],manifests[code],opts["angular_deg"])
                    st.markdown(what)
                    st.markdown(how)
                    st.markdown(observed)
    for note in dict.fromkeys(notes):
        st.caption(note)
    st.plotly_chart(
        fig, use_container_width=True,
        config={"displaylogo":False,
                "toImageButtonOptions":{
                    "format":opts["export_format"],
                    "filename":"compare_"+"_".join(codes)+"_"+kind,
                    "scale":3,
                }})
    col_html, col_csv = st.columns(2)
    name="compare_"+"_".join(codes)+"_"+kind
    with col_html:
        st.download_button(
            "Download interactive comparison (.html)",
            fig.to_html(full_html=True,include_plotlyjs=True),
            file_name=name+".html",
            mime="text/html", use_container_width=True)
    with col_csv:
        _write_csv(st,rows,name+".csv")


def selftest():
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    d = {
        "jpsi":{"gain":np.array([0,0,1,5,12]),"angle_equal":np.ones(72),
                "phi_edges":np.linspace(-180,180,73)},
        "phi":{"gain":np.array([0,2,9,11]),"angle_equal":np.ones(72)*2,
               "phi_edges":np.linspace(-180,180,73)},
    }
    opts={"gain_width":10, "angular_deg":10, "asym_bins":50}
    edges = _shared_edges("gain_all",d,opts)
    assert np.array_equal(edges,np.array([0.,1.,11.,21.]))
    for code in d:
        n, _ = _hist_values("gain_all",d[code],edges,opts)
        assert sum(n) == len(d[code]["gain"])
        assert abs(_normalise(n,edges,"Probability per bin").sum()-1.)<1e-12
        assert abs(np.sum(_normalise(n,edges,"Probability density") *
                          np.diff(edges))-1.)<1e-12
    t = np.arange(0,10,dtype=float)
    n,edge=_group_time(np.ones(9),t,4,False)
    assert len(n)==3 and len(edge)==4 and n.sum()==9
    assert len(_shared_edges("angle_eq",d,opts))==37
    assert _rgb_alpha("#0072B2")=="rgba(0,114,178,0.16)"
    fig=make_subplots(rows=1,cols=2)
    fig.add_trace(go.Heatmap(z=[[1,2]]),row=1,col=1)
    fig.add_trace(go.Heatmap(z=[[3,4]]),row=1,col=2)
    assert len(fig.data)==2
    print("PASS channel comparison: shared bins, zero-gain bin, "
          "independent normalization, grouping and facets")


if __name__ == "__main__":
    selftest()

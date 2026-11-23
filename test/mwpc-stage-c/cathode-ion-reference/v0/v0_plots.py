#!/usr/bin/env python3
"""Physics plots from the V0 CSV exports (2100 V K256 archives; development population, sealed slots excluded upstream).

  python3 v0_plots.py V0_OUT_DIR PLOT_DIR [--lut500-map]   (V0_OUT_DIR: v0_primaries.csv, v0_fragments.csv, ...)

Every panel states N (primaries / fragments / muons) and draws 68% cluster-bootstrap bands, resampling MUONS (primaries
of one muon are not independent). Exact-census columns (microscopic_wire_seed.csv) are used where they exist; reservoir
(K256) moments only where the census has no equivalent (angular / radial cloud moments).
Coordinates: start_du_mm = primary START u minus nearest anode (signed), start_v_mm = start normal coordinate.
If v0_manifest.txt says SYNTHETIC, every figure is stamped SYNTHETIC.

Also reports the between-primary residual correlation (v0_pair_correlation.py) in v0_plots_summary.txt.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from v0_pair_correlation import format_result, pair_correlation  # noqa: E402

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#2b2b29", "#8a8984"
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#ec9a7a", "#b5332b"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "lines.linewidth": 2})
RNG = np.random.default_rng(20261123)
NB = 400
DU_EDGES = np.linspace(-1.5, 1.5, 13)
ADU_EDGES = np.linspace(0, 1.5, 7)
V_EDGES = np.linspace(-3, 3, 13)


def stamp(fig, tag, extra=""):
    fig.text(0.005, 0.005, f"{tag}  |  2100 V K256 V0 archives{extra}", fontsize=7, color=MUTED)


def wmean(v, w):
    s = np.sum(w)
    return np.sum(w * v) / s if s > 0 else np.nan


def binned(df, xcol, ycol, edges, stat="mean", wcol=None, thr=None):
    """Per-bin statistic with muon-cluster bootstrap. stat: mean | rms | frac (fraction ycol>thr)."""
    d = df[[xcol, ycol, "muon_slot"] + ([wcol] if wcol else [])].dropna()
    b = np.digitize(d[xcol].to_numpy(), edges) - 1
    ok = (b >= 0) & (b < len(edges) - 1)
    d = d[ok].assign(_b=b[ok])
    y = d[ycol].to_numpy()
    if stat == "frac":
        y = (y > thr).astype(float)
    w = d[wcol].to_numpy() if wcol else np.ones(len(d))
    mu = d["muon_slot"].to_numpy()
    um, inv = np.unique(mu, return_inverse=True)
    nb = len(edges) - 1
    # per-muon per-bin sums so the bootstrap is fast
    S0 = np.zeros((len(um), nb)); S1 = np.zeros_like(S0); S2 = np.zeros_like(S0); NN = np.zeros_like(S0)
    np.add.at(S0, (inv, d["_b"].to_numpy()), w)
    np.add.at(S1, (inv, d["_b"].to_numpy()), w * y)
    np.add.at(S2, (inv, d["_b"].to_numpy()), w * y * y)
    np.add.at(NN, (inv, d["_b"].to_numpy()), 1)

    def f(s0, s1, s2):
        with np.errstate(invalid="ignore", divide="ignore"):
            m = s1 / s0
            return np.sqrt(np.maximum(s2 / s0 - m * m, 0)) if stat == "rms" else m
    est = f(S0.sum(0), S1.sum(0), S2.sum(0))
    bs = []
    for _ in range(NB):
        k = RNG.integers(0, len(um), len(um))
        bs.append(f(S0[k].sum(0), S1[k].sum(0), S2[k].sum(0)))
    bs = np.array(bs)
    lo, hi = np.nanpercentile(bs, 16, 0), np.nanpercentile(bs, 84, 0)
    n = NN.sum(0)
    est[n < 20] = np.nan; lo[n < 20] = np.nan; hi[n < 20] = np.nan
    return 0.5 * (edges[1:] + edges[:-1]), est, lo, hi, n


def map2d(df, xcol, ycol, zcol, xe, ye, stat="mean", thr=None, wcol=None):
    d = df[[xcol, ycol, zcol] + ([wcol] if wcol else [])].dropna()
    ix = np.digitize(d[xcol], xe) - 1; iy = np.digitize(d[ycol], ye) - 1
    ok = (ix >= 0) & (ix < len(xe) - 1) & (iy >= 0) & (iy < len(ye) - 1)
    z = d[zcol].to_numpy()[ok]
    if stat == "frac":
        z = (z > thr).astype(float)
    w = d[wcol].to_numpy()[ok] if wcol else np.ones(ok.sum())
    S0 = np.zeros((len(ye) - 1, len(xe) - 1)); S1 = np.zeros_like(S0); S2 = np.zeros_like(S0); N = np.zeros_like(S0)
    np.add.at(S0, (iy[ok], ix[ok]), w); np.add.at(S1, (iy[ok], ix[ok]), w * z); np.add.at(S2, (iy[ok], ix[ok]), w * z * z)
    np.add.at(N, (iy[ok], ix[ok]), 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = S1 / S0
        r = np.sqrt(np.maximum(S2 / S0 - m * m, 0))
    out = r if stat == "rms" else m
    out[N < 20] = np.nan
    return out, N


def show_map(ax, Z, N, xe, ye, cmap, label, vmin=None, vmax=None, fmt="{:.2g}"):
    im = ax.pcolormesh(xe, ye, Z, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat")
    plt.colorbar(im, ax=ax, label=label, shrink=0.85)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if N[i, j] >= 20:
                ax.text(0.5 * (xe[j] + xe[j + 1]), 0.5 * (ye[i] + ye[i + 1]), f"{int(N[i, j])}", fontsize=5,
                        ha="center", va="center", color=INK, alpha=0.6)
    ax.grid(False)


def band(ax, x, est, lo, hi, color, label):
    ax.fill_between(x, lo, hi, color=color, alpha=0.22, lw=0)
    ax.plot(x, est, color=color, label=label, marker="o", ms=3.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("v0dir", type=Path)
    ap.add_argument("plotdir", type=Path)
    a = ap.parse_args()
    a.plotdir.mkdir(parents=True, exist_ok=True)
    man = (a.v0dir / "v0_manifest.txt").read_text() if (a.v0dir / "v0_manifest.txt").exists() else ""
    tag = "SYNTHETIC (tool test, not NA60+ data)" if "SYNTHETIC" in man else "REAL development population"
    P = pd.read_csv(a.v0dir / "v0_primaries.csv")
    F = pd.read_csv(a.v0dir / "v0_fragments.csv")
    P["abs_du"] = P["start_du_mm"].abs(); F["abs_du"] = F["start_du_mm"].abs()
    P["absA"] = P["exact_A"].abs()
    P["Nw"] = P["ion_births_exact"].astype(float)
    nmu = P["muon_slot"].nunique()
    hdr = f"{nmu} muons, {len(P)} primaries"
    txt = [f"# V0 plots: {tag}; {hdr}; {len(F)} (primary, nearest-birth wire) fragments"]

    # 1. delay mean / RMS vs start position (exact census, per primary, ion-weighted mean over its wires)
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    for k, (stat, lab) in enumerate((("mean", "mean ion-birth delay [ns]"), ("rms", "RMS across primaries [ns]"))):
        Z, N = map2d(P, "start_du_mm", "start_v_mm", "exact_delay_ns", DU_EDGES, V_EDGES, stat)
        show_map(ax[0, k], Z, N, DU_EDGES, V_EDGES, SEQ, lab)
        ax[0, k].set(xlabel="primary start u - nearest anode [mm]", ylabel="primary start v [mm]",
                     title=f"{lab} (unweighted over primaries); cell N shown")
    for cut, col, lab in ((P["start_v_mm"] < 0, C1, "v0 < 0"), (P["start_v_mm"] > 0, C2, "v0 > 0")):
        x, e, lo, hi, n = binned(P[cut], "abs_du", "exact_delay_ns", ADU_EDGES)
        band(ax[1, 0], x, e, lo, hi, col, f"{lab} (N={int(n.sum())})")
        x, e, lo, hi, n = binned(P[cut], "abs_du", "exact_delay_ns", ADU_EDGES, "rms")
        band(ax[1, 1], x, e, lo, hi, col, f"{lab} (N={int(n.sum())})")
    ax[1, 0].set(xlabel="|start u - nearest anode| [mm]", ylabel="mean delay [ns]")
    ax[1, 1].set(xlabel="|start u - nearest anode| [mm]", ylabel="delay RMS [ns]")
    for x_ in ax[1]:
        x_.legend(frameon=False)
    fig.suptitle(f"1. Ion-birth delay after primary creation vs start position — exact census; {hdr}; 68% muon-bootstrap bands")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.02, 1, 0.96)); fig.savefig(a.plotdir / "01_delay_vs_start.png", dpi=150); plt.close(fig)

    # 2. hemisphere asymmetry A vs signed v0 and |du|
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    Z, N = map2d(P, "start_du_mm", "start_v_mm", "exact_A", DU_EDGES, V_EDGES)
    show_map(ax[0, 0], Z, N, DU_EDGES, V_EDGES, DIV, "mean exact A", -1, 1)
    ax[0, 0].set(xlabel="start u - nearest anode [mm]", ylabel="start v [mm]", title="mean A per primary (exact census)")
    Z, N = map2d(P, "start_du_mm", "start_v_mm", "exact_A", DU_EDGES, V_EDGES, "rms")
    show_map(ax[0, 1], Z, N, DU_EDGES, V_EDGES, SEQ, "RMS of A", 0, 1)
    ax[0, 1].set(xlabel="start u - nearest anode [mm]", ylabel="start v [mm]", title="RMS of A within cell")
    terc = np.quantile(P["abs_du"].dropna(), [0, 1 / 3, 2 / 3, 1])
    for i, col in enumerate((C1, C2, C3)):
        cut = (P["abs_du"] >= terc[i]) & (P["abs_du"] <= terc[i + 1])
        x, e, lo, hi, n = binned(P[cut], "start_v_mm", "exact_A", V_EDGES)
        band(ax[1, 0], x, e, lo, hi, col, f"|du| {terc[i]:.2f}-{terc[i+1]:.2f} mm (N={int(n.sum())})")
    ax[1, 0].set(xlabel="signed start v [mm]", ylabel="mean A", ylim=(-1.05, 1.05)); ax[1, 0].legend(frameon=False, fontsize=7)
    for cut, col, lab in ((P["start_v_mm"] < 0, C1, "v0 < 0"), (P["start_v_mm"] > 0, C2, "v0 > 0")):
        x, e, lo, hi, n = binned(P[cut], "abs_du", "exact_A", ADU_EDGES)
        band(ax[1, 1], x, e, lo, hi, col, f"{lab} (N={int(n.sum())})")
    ax[1, 1].set(xlabel="|start u - nearest anode| [mm]", ylabel="mean A", ylim=(-1.05, 1.05)); ax[1, 1].legend(frameon=False)
    fig.suptitle(f"2. Avalanche hemisphere asymmetry A vs start position; {hdr}")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.02, 1, 0.96)); fig.savefig(a.plotdir / "02_A_vs_start.png", dpi=150); plt.close(fig)

    # 3. frac(|A| > 0.9)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    Z, N = map2d(P, "start_du_mm", "start_v_mm", "absA", DU_EDGES, V_EDGES, "frac", 0.9)
    show_map(ax[0], Z, N, DU_EDGES, V_EDGES, SEQ, "P(|A| > 0.9)", 0, 1)
    ax[0].set(xlabel="start u - nearest anode [mm]", ylabel="start v [mm]")
    for i, col in enumerate((C1, C2, C3)):
        cut = (P["abs_du"] >= terc[i]) & (P["abs_du"] <= terc[i + 1])
        x, e, lo, hi, n = binned(P[cut], "start_v_mm", "absA", V_EDGES, "frac", thr=0.9)
        band(ax[1], x, e, lo, hi, col, f"|du| tercile {i+1} (N={int(n.sum())})")
    ax[1].set(xlabel="signed start v [mm]", ylabel="P(|A| > 0.9)", ylim=(0, 1)); ax[1].legend(frameon=False, fontsize=7)
    fig.suptitle(f"3. One-sided avalanches, fraction |A|>0.9 (exact census); {hdr}")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.03, 1, 0.94)); fig.savefig(a.plotdir / "03_frac_absA_gt_0p9.png", dpi=150); plt.close(fig)
    txt.append(f"pooled P(|A|>0.9) = {np.mean(P['absA'].dropna() > 0.9):.4f} over {P['absA'].notna().sum()} primaries")

    # 4. correlations as conditional profiles (fragments: reservoir cloud moments; exact A and yield)
    Fd = F.copy(); Fd["logN"] = np.log10(Fd["exact_births"].where(Fd["exact_births"] > 0))
    pairs = [("sin", "exact_A", "<sin theta> (reservoir)", "exact A"), ("sin", "dy_um", "<sin theta>", "<dy> [um]"),
             ("rbar_um", "exact_A", "<r> [um]", "exact A"), ("exact_delay_ns", "exact_A", "delay [ns]", "exact A"),
             ("logN", "exact_A", "log10 ion yield", "exact A"), ("logN", "exact_delay_ns", "log10 ion yield", "delay [ns]")]
    fig, ax = plt.subplots(2, 3, figsize=(12, 7))
    for k, (xc, yc, xl, yl) in enumerate(pairs):
        aa = ax.flat[k]
        x = Fd[xc].dropna()
        e = np.unique(np.quantile(x, np.linspace(0, 1, 13)))
        xm, est, lo, hi, n = binned(Fd, xc, yc, e)
        _, rm, _, _, _ = binned(Fd, xc, yc, e, "rms")
        aa.fill_between(xm, est - rm, est + rm, color=MUTED, alpha=0.15, lw=0, label="+-1 conditional RMS")
        band(aa, xm, est, lo, hi, C1, "conditional mean (68% band)")
        dd = Fd[[xc, yc]].dropna()
        rho = dd.rank().corr().iloc[0, 1]
        aa.set(xlabel=xl, ylabel=yl, title=f"Spearman {rho:+.3f}, N={len(dd)} fragments")
        txt.append(f"Spearman({xc},{yc}) = {rho:+.4f}  N={len(dd)} fragments")
    ax.flat[0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"4. Conditional profiles (quantile bins) — not only pooled correlation; {nmu} muons")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.02, 1, 0.95)); fig.savefig(a.plotdir / "04_correlation_profiles.png", dpi=150); plt.close(fig)

    # 5. within-avalanche ion-birth times
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    ion_path = a.v0dir / "v0_ion_examples.csv"
    if ion_path.exists():
        I = pd.read_csv(ion_path)
        y = I.groupby(["muon_slot", "primary"])["weight"].agg(["sum", "size"]).reset_index()
        y = y[y["size"] >= 64]
        lo_r = y.iloc[(y["sum"] - y["sum"].quantile(0.1)).abs().argmin()]
        hi_r = y.iloc[(y["sum"] - y["sum"].quantile(0.9)).abs().argmin()]
        for aa, r, col, nm in ((ax[0], lo_r, C1, "low-gain (10th pct yield)"), (ax[1], hi_r, C2, "high-gain (90th pct yield)")):
            s = I[(I.muon_slot == r.muon_slot) & (I.primary == r.primary)]
            t = s["t_minus_tprim_ns"].to_numpy(); w = s["weight"].to_numpy()
            tm = wmean(t, w); sd = np.sqrt(wmean((t - tm) ** 2, w))
            aa.hist(t - tm, bins=40, weights=w, color=col, alpha=0.85)
            aa.set(xlabel="ion birth time - avalanche mean [ns]", ylabel="ions (weighted)",
                   title=f"{nm}\nmuon {int(r.muon_slot)} prim {int(r.primary)}: N={r['sum']:.0f} ions, {len(s)} saved, RMS {sd:.3f} ns", )
            aa.title.set_fontsize(8)
    else:
        for aa in ax[:2]:
            aa.text(0.5, 0.5, "run V0 with --ion-dump-muons N", ha="center", transform=aa.transAxes)
    st = F["exact_sigt_ns"].dropna(); st = st[st > 0]
    ax[2].hist(st, bins=np.logspace(np.log10(max(st.min(), 1e-3)), np.log10(st.max()), 50), color=C1)
    ax[2].set_xscale("log")
    ax[2].set(xlabel="exact within-fragment birth-time sigma [ns]", ylabel="fragments",
              title=f"all fragments: median {st.median():.3f} ns, N={len(st)}")
    fig.suptitle("5. Ion-birth time spread inside single avalanches (reservoir ions for the examples; exact census for the distribution)")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.03, 1, 0.93)); fig.savefig(a.plotdir / "05_within_avalanche_times.png", dpi=150); plt.close(fig)

    # 6. onset distributions within a muon and within (muon, collecting wire)
    Q = P.dropna(subset=["exact_onset_ns"]).copy()
    Q = Q[Q["collected_electrons"] > 0]
    Q["w"] = Q["collected_electrons"].astype(float)
    # Use grouped sums/transform rather than groupby.apply(include_groups=False).
    # The latter requires pandas >=2.2 and causes plotting to stop after
    # figure 5 on older WSL pandas installations.  This is also faster.
    Q["_weighted_onset"] = Q["exact_onset_ns"] * Q["w"]
    for key, nm in ((["muon_slot"], "mu"), (["muon_slot", "argmax_collected_wire"], "wire")):
        sw = Q.groupby(key, sort=False)["w"].transform("sum")
        st = Q.groupby(key, sort=False)["_weighted_onset"].transform("sum")
        Q[f"m_{nm}"] = st / sw
        Q[f"res_{nm}"] = Q["exact_onset_ns"] - Q[f"m_{nm}"]
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    rng_ = np.linspace(-50, 50, 101)
    for col, nm, lab in ((C1, "mu", "minus muon mean"), (C2, "wire", "minus (muon, collecting wire) mean")):
        r = Q[f"res_{nm}"]; w = Q["w"]
        rms = np.sqrt(wmean(r ** 2, w))
        ax[0].hist(r, bins=rng_, weights=w / w.sum(), histtype="step", lw=2, color=col, label=f"{lab}: RMS {rms:.2f} ns")
        txt.append(f"onset residual {lab}: electron-weighted RMS {rms:.3f} ns over {len(r)} primaries")
    ax[0].set(xlabel="primary onset (mean collected-electron time) residual [ns]", ylabel="fraction (electron-weighted)")
    ax[0].legend(frameon=False, fontsize=7)
    Q["_weighted_res_mu2"] = Q["w"] * Q["res_mu"] ** 2
    by_mu = Q.groupby("muon_slot")
    per_mu = np.sqrt(by_mu["_weighted_res_mu2"].sum() / by_mu["w"].sum())
    ax[1].hist(per_mu, bins=30, color=C1)
    ax[1].set(xlabel="per-muon onset RMS [ns]", ylabel="muons", title=f"N={len(per_mu)} muons, median {per_mu.median():.2f} ns")
    x, e, lo, hi, n = binned(Q.assign(ab=Q["start_v_mm"].abs()), "ab", "exact_onset_ns", np.linspace(0, 3, 7), wcol="w")
    band(ax[2], x, e, lo, hi, C1, "mean onset")
    ax[2].set(xlabel="|start v| [mm]", ylabel="onset after primary creation [ns]", title="onset vs drift distance")
    fig.suptitle(f"6. Effective onset times across primaries of one muon / one collecting wire; {len(Q)} primaries with collected electrons")
    stamp(fig, tag); fig.tight_layout(rect=(0, 0.03, 1, 0.93)); fig.savefig(a.plotdir / "06_onset_spread.png", dpi=150); plt.close(fig)

    # 7. collecting anode vs anode nearest the primary START
    C = P[(P["argmax_collected_wire"] != -999) & (P["start_nearest_wire"] != -999)].copy()
    C["off"] = C["argmax_collected_wire"] - C["start_nearest_wire"]
    C["mis"] = (C["off"] != 0).astype(float)
    C["multi"] = (C["collected_frac_max_wire"] < 0.95).astype(float)
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    Z, N = map2d(C, "start_du_mm", "start_v_mm", "mis", DU_EDGES, V_EDGES)
    show_map(ax[0, 0], Z, N, DU_EDGES, V_EDGES, SEQ, "P(collecting != nearest-start anode)", 0, None)
    ax[0, 0].set(xlabel="start u - nearest anode [mm]", ylabel="start v [mm]")
    for cut, col, lab in ((C["start_v_mm"].abs() < 1.5, C1, "|v0| < 1.5 mm"), (C["start_v_mm"].abs() >= 1.5, C2, "|v0| >= 1.5 mm")):
        x, e, lo, hi, n = binned(C[cut], "abs_du", "mis", ADU_EDGES)
        band(ax[0, 1], x, e, lo, hi, col, f"{lab} (N={int(n.sum())})")
    ax[0, 1].set(xlabel="|start u - nearest anode| [mm]", ylabel="mismatch probability"); ax[0, 1].legend(frameon=False)
    vals, cnt = np.unique(C["off"], return_counts=True)
    ax[1, 0].bar(vals, cnt / cnt.sum(), color=C1, width=0.6)
    ax[1, 0].set_yscale("log")
    ax[1, 0].set(xlabel="collecting wire - nearest-start wire", ylabel="fraction of primaries",
                 title=f"N={len(C)}; mismatch {C['mis'].mean():.4f}")
    x, e, lo, hi, n = binned(C, "abs_du", "multi", ADU_EDGES)
    band(ax[1, 1], x, e, lo, hi, C3, f"collected on >1 wire (max share < 0.95), N={int(n.sum())}")
    ax[1, 1].set(xlabel="|start u - nearest anode| [mm]", ylabel="probability"); ax[1, 1].legend(frameon=False, fontsize=7)
    fig.suptitle(f"7. Dominant collecting anode vs anode nearest the PRIMARY START (2100 V exact census); {len(C)} primaries")
    stamp(fig, tag, "  |  NOT the 0.121% collecting-vs-ion-birth-wire comparison")
    fig.tight_layout(rect=(0, 0.03, 1, 0.95)); fig.savefig(a.plotdir / "07_collecting_vs_start_anode.png", dpi=150); plt.close(fig)
    txt.append(f"collecting != nearest-START anode: {int(C['mis'].sum())}/{len(C)} = {C['mis'].mean():.5f} (dominant collecting wire)")
    for v_, c_ in zip(vals, cnt):
        txt.append(f"   offset {int(v_):+d}: {c_}")

    # 8. LUT500 1800 V side-by-side (transport map) if present
    T = pd.read_csv(a.v0dir / "v0_transport_map.csv")
    if T["lut1800_n"].notna().any() and (T["lut1800_n"] > 0).any():
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for vb, col in ((-0.25, C1), (-1.75, C2), (0.25, C3)):
            t = T[(T.v_lo_mm <= vb) & (T.v_hi_mm > vb)]
            xm = 0.5 * (t.du_lo_mm + t.du_hi_mm)
            ax[0].plot(xm, t.v2100_same_wire_frac, color=col, marker="o", label=f"2100 V archive, v in [{t.v_lo_mm.iloc[0]},{t.v_hi_mm.iloc[0]})")
            ax[0].plot(xm, t.lut1800_same_wire_frac, color=col, ls="--", marker="s", label="1800 V LUT500 (single e-)")
            ax[1].plot(xm, t.v2100_onset_mean_ns, color=col, marker="o")
            ax[1].plot(xm, t.lut1800_drift_mean_ns, color=col, ls="--", marker="s")
        ax[0].set(xlabel="|start u - nearest anode| [mm]", ylabel="P(collected on nearest-start anode)")
        ax[1].set(xlabel="|start u - nearest anode| [mm]", ylabel="ns", title="solid: 2100 V mean collected-e time;  dashed: 1800 V single-e drift time")
        ax[1].title.set_fontsize(8)
        ax[0].legend(frameon=False, fontsize=6)
        fig.suptitle("8. 1800 V LUT500 vs 2100 V archives — DIFFERENT voltages and observables; comparison only, not calibration")
        stamp(fig, tag); fig.tight_layout(rect=(0, 0.03, 1, 0.93)); fig.savefig(a.plotdir / "08_lut500_1800V_vs_2100V.png", dpi=150); plt.close(fig)
    else:
        txt.append("plot 8 skipped: v0_transport_map.csv has no LUT500 columns (run V0 with MWPC_LUT500_DIR)")

    # repaired correlation estimator
    txt.append("\n# Between-primary residual correlation within muons (assumptions and bounds: v0_pair_correlation.py)")
    txt.append("# context = 12 signed-du x 12 v start cells x yield tercile; weight = exact ion yield; leave-muon-out means")
    P["_cu"] = np.digitize(P["start_du_mm"], DU_EDGES); P["_cv"] = np.digitize(P["start_v_mm"], V_EDGES)
    P["_wN"] = P["ion_births_exact"].astype(float)
    # exact-census quantities have no reservoir noise; reservoir A carries its own se_A. Reservoir <sin> etc. are
    # omitted because v0_primaries.csv has no per-primary standard error for them (denominator would be biased).
    for yc, sc in (("exact_A", None), ("A", "se_A"), ("exact_delay_ns", None), ("exact_onset_ns", None)):
        if yc in P and P[yc].notna().sum() > 100:
            txt.append(format_result(pair_correlation(P, yc, sc, ["_cu", "_cv"], "_wN", rng=RNG)))
    (a.plotdir / "v0_plots_summary.txt").write_text("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()

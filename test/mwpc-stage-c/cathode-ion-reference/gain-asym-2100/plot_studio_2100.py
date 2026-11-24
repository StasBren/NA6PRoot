#!/usr/bin/env python3
"""NA60+ 2100-V Plot Studio: read-only visualization of v0_archive_stats outputs.

No Garfield/Heed/ion-drift jobs. Computations for trends and figures use
the native C++/ROOT exported exact counters or labelled K256 estimators.
"""
from __future__ import annotations
import io
import math
import os
from pathlib import Path
from dataclasses import dataclass
from typing import Callable
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from per_anode_census_2100 import read_per_anode

st.set_page_config(page_title="NA60+ | 2100 V Plot Studio", layout="wide",
                   page_icon="📈", initial_sidebar_state="expanded")

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE.parents[3] / "test_runs" / "mwpc_stage_c" / "2100_gain_asym"
ROOT = Path(os.environ.get("NA60_PLOTSTUDIO_DATA_ROOT", str(DEFAULT_ROOT))).expanduser()

@dataclass(frozen=True)
class Spec:
    key: str
    section: str
    title: str
    dataset: str
    x: str
    y: str | None
    chart: str
    xlabel: str
    ylabel: str
    meaning: str
    calculation: str
    limitations: str
    logx: bool = False
    logy: bool = False

SPECS = [
 Spec("mu_ion_gain","Gas gain","Muon ionization gain","mu","G_ion",None,"hist",
      "Ion births / Heed primaries","Muons",
      "Distribution of ionization gain per muon.",
      "Exact Garfield ion-birth count divided by accepted Heed primary count for each muon. The population mean is a ratio of sums.",
      "A muon with zero Heed primaries has undefined gain. Zero-gain muons with primaries are retained."),
 Spec("mu_col_gain","Gas gain","Muon collected-electron gain","mu","G_collected",None,"hist",
      "Collected electrons / Heed primaries","Muons",
      "Actual electron-collection gain per muon.",
      "Exact collected-anode electron endpoints / Heed primaries; distinct from ion-birth yield.",
      "A zero-Heed muon has undefined gain. Zeros with nonzero primaries are retained."),
 Spec("primary_ions","Gas gain","Ion births per primary avalanche","primary","G_ion",None,"hist",
      "Ion births per Heed primary","Heed primary avalanches",
      "Single-primary microscopic avalanche ion-yield fluctuations.",
      "One primary Heed electron starts each avalanche; exact integer birth count is its ionization yield.",
      "Zero-ion avalanches are kept; this is not the muon-level gain."),
 Spec("gain_corr","Gas gain","Collected versus ionization gain","mu","G_ion","G_collected","density",
      "Muon ionization gain","Muon collected-electron gain",
      "Relation between positive-ion production and collected electrons.",
      "Both axes are independent exact muon counters divided by the same accepted primary-electron count.",
      "A correlation is not proof that every ion creates a collected anode electron.",True,True),
 Spec("total_ions","Gas gain","Total ion births per muon","mu","ion_births_exact",None,"hist",
      "Ion births / muon","Muons","Absolute ionization burden of each muon.",
      "Sum of all exact ion-birth callbacks across accepted primary avalanches.",
      "Do not divide by Heed count for this figure."),
 Spec("total_collected","Gas gain","Total collected electrons per muon","mu","collected_electrons",None,"hist",
      "Collected electrons / muon","Muons","Absolute collected electron charge as an electron count.",
      "Sum of physical-anode collecting electron endpoints, not induced cathode charge.",
      "Electron counts are not electronics output."),

 Spec("wire_mult","Anode charge sharing","Collecting anodes per muon","mu","n_collecting_wires",None,"hist",
      "Collecting anode multiplicity","Muons",
      "How many physical anodes receive electrons from one muon.",
      "Count wires with exact collected_e>0. Include muons collecting zero electrons.",
      "Derived from the existing V0 per-wire census and guarded by sum of collected electrons matching the muon total."),
 Spec("second_frac","Anode charge sharing","Second-wire collected fraction","mu","f2",None,"hist",
      "Second-largest wire charge / total","Muons with Q_total > 0",
      "Amplitude of multiwire charge sharing.",
      "f2 = Q2 / sum Q, where Q2 is the second-largest physical collecting-anode electron count.",
      "Undefined when total collection is zero; those events are excluded from this ratio."),
 Spec("sharing_angle","Anode charge sharing","P(f2 > 0.05) versus |theta_u|","mu","abs_theta_u","share5","profile",
      "|Muon incidence angle theta_u| [deg]","P(f2 > 0.05)",
      "Frequency of appreciable second-wire electron sharing versus incidence.",
      "For each muon, indicator 1 if f2>0.05, else 0. Empty-collection muons remain 0 in the denominator.",
      "Muon-level grouped estimates, not a cathode-strip efficiency curve."),
 Spec("sharing_gain","Anode charge sharing","Second-wire fraction versus ionization gain","mu","G_ion","f2","density",
      "Muon ionization gain","f2 = Q2/Qtotal",
      "How physical-anode charge sharing correlates with gas ionization gain.",
      "f2 from the actual collecting-wire census; ionization gain from exact birth/Heed counts.",
      "Muon-level only. Muons without collected electrons have undefined f2.",True,False),

 Spec("angular_spread_gain","Avalanche angular structure","Circular angular spread versus ion gain","primary","G_ion","sigma_phi_deg","profile",
      "Ion births / primary","Circular spread sigma_phi [deg]",
      "How broadly an individual avalanche produces ions around its collecting-anode neighborhood.",
      "R = sqrt(<sin phi>²+<cos phi>²), sigma_phi=sqrt(-2 ln R). <sin>,<cos> are weighted K256 ion-birth reservoir estimates aggregated per Heed primary.",
      "Approximate reservoir estimator, NOT an exact full-ion circular spread. Undefined if resultant length is zero or fewer than 3 sampled births.",True),
 Spec("angular_density","Avalanche angular structure","Angular spread population density","primary","G_ion","sigma_phi_deg","density",
      "Ion births / primary","Circular spread sigma_phi [deg]",
      "Population density of avalanche angular dispersion across gain.",
      "2D histogram of per-primary circular standard deviation estimates versus exact primary ion births.",
      "Each primary contributes one observation; this is NOT an ion-weighted angular histogram.",True),
 Spec("r_vs_gain","Avalanche angular structure","Angular resultant length versus ion gain","primary","G_ion","R_phi","profile",
      "Ion births / primary","Mean circular resultant R",
      "Degree of angular localization, with R near 1 indicating concentrated azimuths.",
      "R from sample-estimated mean sine and cosine for each Heed primary.",
      "A small R can indicate a broad or multimodal distribution; R alone cannot distinguish those.",True),
 Spec("angular_A","Avalanche angular structure","Angular spread versus exact hemispheric A","primary","A_exact","sigma_phi_deg","density",
      "Exact A_UD","Circular spread sigma_phi [deg]",
      "Association of signed hemisphere preference with circular angular spread.",
      "Exact up/down A from original counters; circular spread from K256 sampled positions.",
      "The two observables differ in sampling precision; both describe ion births, not drift."),

 Spec("long_spread_gain","Along-wire ion structure","Along-wire sigma_w versus ion gain","fragment","G_ion_primary","exact_sigw_mm","profile",
      "Ion births / primary","sigma_w [mm]",
      "Longitudinal extent of ions within one (Heed primary, nearest ion-birth anode) fragment.",
      "Exact per-wire ion-birth longitudinal RMS from frozen Garfield sums of w and w², in millimeters.",
      "This is per-fragment, not full-muon width. Multiple fragments from one primary are correlated; uncertainty resamples muons.",True),
 Spec("long_density","Along-wire ion structure","Along-wire sigma_w population density","fragment","G_ion_primary","exact_sigw_mm","density",
      "Ion births / primary","sigma_w [mm]",
      "Population density of within-fragment along-wire spread.",
      "Each nonzero-birth fragment contributes one exact RMS-versus-gain point.",
      "This is not the distribution of delta-w individual ion births; raw primary-relative longitudinal coordinates are not in the current compact export.",True),
 Spec("long_angle","Along-wire ion structure","Along-wire sigma_w versus |theta_u|","fragment","abs_theta_u","exact_sigw_mm","profile",
      "|Muon incidence angle theta_u| [deg]","sigma_w [mm]",
      "Dependence of within-fragment avalanche ion-cloud length on the muon incidence angle.",
      "Exact per-fragment RMS matched to the original parent-muon theta_u angle.",
      "No causal inference from the plot without controlling geometry and gain."),

 Spec("A_primary","Up/down hemisphere","Avalanche A_UD distribution","primary","A_exact",None,"hist",
      "A_UD","Primary avalanches with defined A",
      "Up/down ion-birth hemispheric asymmetry for individual Heed primary avalanches.",
      "A=(N_up-N_down)/(N_up+N_down) using exact Garfield per-wire sector census summed within primary.",
      "Births exactly on the dividing axis are excluded from A denominator, not ion gain; zero-denominator avalanches have undefined A."),
 Spec("A_muon","Up/down hemisphere","Muon A_UD distribution","mu","A_exact",None,"hist",
      "A_UD","Muons with defined A",
      "Net up/down birth asymmetry after summing a muon's avalanches.",
      "V0 exports the ion-birth-weighted mean of primary exact hemisphere asymmetries, A_exact = sum(Nbirth_primary*A_primary)/sum(Nbirth_primary).",
      "This is not cathode-induced charge. With births on the dividing boundary, the V0 weighted mean can differ from a direct ratio of globally summed N_up and N_down; inspect boundary counts before treating them as strictly identical."),
 Spec("A_gain","Up/down hemisphere","Mean signed A versus ion gain","primary","G_ion","A_exact","profile",
      "Ion births / primary","Mean A_UD",
      "Does avalanche ion yield correlate with a preferred up/down birth hemisphere?",
      "Each positive-ion Heed primary contributes one exact A value. Points are equal-primary binned means.",
      "Sign can cancel across opposite hemispheres; muon-block bootstrap handles correlated primaries.",True),
 Spec("absA_gain","Up/down hemisphere","Mean |A| versus ion gain","primary","G_ion","absA_exact","profile",
      "Ion births / primary","Mean |A_UD|",
      "Does a larger avalanche have stronger localization within one up/down hemisphere?",
      "Mean of absolute per-primary exact A in logarithmic gain bins, equal weight per primary; not |mean A|.",
      "Zero-ion avalanches do not have defined A; variation by starting position can confound gain trends.",True),
 Spec("A_muon_gain","Up/down hemisphere","Muon |A| versus effective ion gain","mu","G_ion","absA_exact","profile",
      "Ionization gain / Heed primary","Mean |A_UD| per muon",
      "How muon-level aggregation changes the observed hemispheric preference.",
      "Absolute value of each muon's exact ion-yield-weighted hemispheric A; bins are equal-muon means.",
      "Separate from the avalanche-level distribution and from induced cathode charge.",True),
 Spec("A_gain_density","Up/down hemisphere","A versus ion gain: population density","primary","G_ion","A_exact","density",
      "Ion births / primary","A_UD",
      "Joint population distribution of exact avalanche ion yield and signed hemispheric asymmetry.",
      "2D histogram with one entry per primary with defined exact A.",
      "Bin population density is not an ion-weighted count.",True),
 Spec("A_startV","Up/down hemisphere","Mean A versus primary start v","primary","start_v_mm","A_exact","profile",
      "Primary start coordinate v [mm]","Mean A_UD",
      "Effect of which side of the anode the initiating electron begins on.",
      "Binned exact per-primary A conditioned on original Heed electron start v coordinate.",
      "Source geometry and impact distributions can drive the asymmetry trend."),
 Spec("A_startDU","Up/down hemisphere","Mean A versus signed transverse wire offset","primary","start_du_mm","A_exact","profile",
      "Primary start u relative to nearest anode [mm]","Mean A_UD",
      "How the transverse starting position relates to ion-birth hemisphere.",
      "Binned exact A versus signed nearest-anode offset of the Heed electron start.",
      "Nearest collecting wire and nearest initial position are different possible assignments."),
 Spec("A_theta","Up/down hemisphere","Mean A versus signed incidence theta_u","primary","theta_u_deg","A_exact","profile",
      "Signed muon incidence theta_u [deg]","Mean A_UD",
      "Relation of the original muon crossing direction to avalanche hemisphere preference.",
      "Primary avalanches inherit their parent muon's signed theta_u; uncertainties resample parent muons.",
      "Many primaries share each angle: do not treat them as independent muon observations."),
 Spec("long_A","Up/down hemisphere","Along-wire sigma_w versus A","fragment","exact_A","exact_sigw_mm","density",
      "Exact fragment A_UD","Exact fragment sigma_w [mm]",
      "Connection between up/down ion-birth localization and the along-wire cloud width.",
      "Both exact from original per-(primary, nearest ion-birth wire) census.",
      "One entry per physical ion-birth wire fragment; the sample is not an ion-weighted angular distribution."),
 # One exact original Garfield CSV row per (Heed primary, local anode).
 Spec("anode_fragment_A","Individual anodes","Ion-birth A_UD — each primary × anode",
      "anode_fragment","A_UD_exact",None,"hist",
      "A_UD (exact up/down births)","Primary–anode fragments",
      "Hemispheric asymmetry of ion births attributed to one physical anode in a single primary avalanche.",
      "Original microscopic_wire_seed.csv exact ion_up/ion_down/ion_boundary counters; A=(up-down)/(up+down), with dividing-axis births excluded from A only.",
      "Anode assignment is the wire nearest ion births, not necessarily the wire collecting avalanche electrons."),
 Spec("anode_muon_A","Individual anodes","Ion-birth A_UD — each muon × anode",
      "anode","A_UD_exact",None,"hist",
      "A_UD (exact sum of up/down births)","Muon–anode records",
      "Net ion-birth hemisphere preference on a given anode over all primary avalanches from one muon.",
      "For each (muon, global anode), sum exact N_up and N_down across its primaries, then form (up-down)/(up+down).",
      "Exact ratio-of-sums, not the mean of fragment A values; ion-birth assignment, not collected cathode current."),
 Spec("anode_fragment_t","Individual anodes","Primary → mean ion-birth time (per avalanche–anode)",
      "anode_fragment","ion_birth_delay_ns",None,"hist",
      "Mean ion-birth delay from primary [ns]","Primary–anode fragments",
      "Time elapsed between creation of a Heed primary electron and the mean formation time of its avalanche ions near an anode.",
      "The original exact per-wire mean Garfield ion birth timestamp minus this primary's t_primary_ns.",
      "This is a mean ion-production delay, not the first ion birth, not ion drift to the cathode and not electronics latency."),
 Spec("anode_fragment_sigma_t","Individual anodes","Ion-production time spread within an avalanche–anode",
      "anode_fragment","sigma_t_birth_ns",None,"hist",
      "Within-fragment sigma_t [ns]","Primary–anode fragments",
      "Microscopic temporal extent of positive-ion production during one avalanche.",
      "Original exact sigma_t_birth_ns from the microscopic per-primary/per-wire counters.",
      "Do not confuse this sub-avalanche spread with across-primary arrival-time spread on an anode."),
 Spec("anode_fragment_coll_t","Individual anodes","Primary → mean collected-electron time",
      "anode_fragment","electron_collection_delay_ns",None,"hist",
      "Mean electron collection delay [ns]","Primary–anode fragments with collected electrons",
      "Mean time for electrons produced by a primary avalanche to reach this physical collecting anode.",
      "Original exact mean_t_collected_ns minus the original primary's t_primary_ns.",
      "Only wires actually collecting electrons contribute; electron collection is not positive-ion collection."),
 Spec("anode_fragment_time_compare","Individual anodes","Ion formation versus electron collection delay",
      "anode_fragment","ion_birth_delay_ns","electron_collection_delay_ns","density",
      "Mean ion-birth delay from primary [ns]","Mean collected-electron delay [ns]",
      "Compare the two physically distinct microscopic avalanche timing observables.",
      "Both delays are computed from original exact source means using the same Heed primary timestamp.",
      "Requires a wire with both positive ion births and actual collected electrons; does not include ion drift."),
 Spec("anode_fragment_t_gain","Individual anodes","Ion-birth delay versus primary ionization gain",
      "anode_fragment","G_ion_primary","ion_birth_delay_ns","profile",
      "Ion births per Heed primary","Mean ion-birth delay [ns]",
      "Tests whether larger avalanche ion yields are associated with later ion formation.",
      "Exact per-primary ionization gain, exact per-(primary,wire) mean ion-birth delay. One entry per primary–anode fragment; muon bootstrap.",
      "If a primary's births span multiple anodes, the full primary yield is used on the x-axis for each fragment.",True),
 Spec("anode_fragment_A_gain","Individual anodes","Per-anode |A_UD| versus primary ionization gain",
      "anode_fragment","G_ion_primary","absA_UD_exact","profile",
      "Ion births per Heed primary","Mean |A_UD|",
      "How ion production hemisphere localization varies with individual avalanche yield.",
      "Absolute exact up/down ratio for each primary–anode record, binned by exact ion births across that primary.",
      "Equal fragment weighting, not ion weighting; muon bootstrap accounts for shared parents.",True),
 Spec("anode_fragment_A_time","Individual anodes","Anode ion-birth asymmetry versus production delay",
      "anode_fragment","ion_birth_delay_ns","A_UD_exact","density",
      "Mean ion-birth delay [ns]","A_UD",
      "Joint distribution of signed ion-birth hemisphere and elapsed avalanche formation time.",
      "Exact per-primary–anode A and exact mean ion production delay from original Garfield counters.",
      "Does not resolve a within-avalanche time-dependent A(t); this is a correlation of summaries."),
 Spec("anode_fragment_w_gain","Individual anodes","Within-avalanche longitudinal spread versus gain",
      "anode_fragment","G_ion_primary","sigma_w_birth_mm","profile",
      "Ion births per Heed primary","Within-fragment sigma_w [mm]",
      "Longitudinal ion-birth extent along the anode, rather than across the chamber gap.",
      "Exact sigma_w_birth_mm, original source per (primary, nearest-birth wire), against exact primary ion yield.",
      "Not the spatial width of the entire muon or the collection charge sharing.",True),
 # One exact aggregation per (muon, physical/nearest-birth anode).
 Spec("anode_muon_births","Individual anodes","Total ion births on one anode per muon",
      "anode","ion_births",None,"hist",
      "Total ion births on anode","Muon–anode records",
      "How many positive ions in a muon were produced nearest to a given anode.",
      "Exact ion_births summed over all contributing primaries on (muon, global_wire).",
      "Counts are birth locations, not ions arriving at that wire or at the cathode.",True),
 Spec("anode_muon_collected","Individual anodes","Total electrons collected on one anode per muon",
      "anode","collected_e",None,"hist",
      "Electrons collected at anode","Muon–anode records",
      "Physical electron charge collecting on one anode for a given muon.",
      "Exact collected_e counts summed over all Heed primary avalanches.",
      "Not induced signal charge. Zero-electron wires with ions remain in the denominator.",True),
 Spec("anode_muon_t","Individual anodes","Mean ion-birth delay on each anode",
      "anode","ion_birth_delay_ns",None,"hist",
      "Ion-weighted mean delay from primary [ns]","Muon–anode records",
      "Typical delay between each ion's originating Heed primary and the ion creation, pooled on one anode.",
      "Birth-count-weighted mean of per-primary–anode (t_ion_mean-t_primary). Different primary time origins are individually subtracted.",
      "Not a common-event-clock arrival time; it excludes all positive-ion drift to cathodes."),
 Spec("anode_muon_sigma_t","Individual anodes","Ion-birth timing RMS on each anode",
      "anode","sigma_t_birth_relative_ns",None,"hist",
      "RMS of (t_ion−t_primary) [ns]","Muon–anode records",
      "Overall distribution width of ion production delays relative to each ion's own Heed primary.",
      "Combine exact within-fragment sigma_t and between-fragment mean-delay variance with ion-birth count weights.",
      "Unlike within-avalanche sigma_t, this includes differences between primaries contributing to the same anode."),
 Spec("anode_muon_share","Individual anodes","Muon ion-birth fraction on one anode",
      "anode","ion_share_muon",None,"hist",
      "Births on anode / births in muon","Muon–anode records",
      "Ionization charge-sharing fraction among wires defined by nearest ion-birth positions.",
      "Exact total ion births near a given wire divided by the muon's total exact ion births.",
      "Not equal to collected-electron anode sharing, nor induced U/V cathode-strip sharing."),
 Spec("anode_muon_e_share","Individual anodes","Muon collected-electron fraction on one anode",
      "anode","electron_share_muon",None,"hist",
      "Collected electrons on anode / muon total","Muon–anode records",
      "Fraction of actual anode-collected avalanche electrons on a given physical wire.",
      "Exact per-wire collecting electron counts normalized by total collected electron count for this muon.",
      "Zero-collected-electron muons have undefined fractions."),
 Spec("anode_muon_A_t","Individual anodes","Muon–anode A_UD versus ion-birth delay",
      "anode","ion_birth_delay_ns","A_UD_exact","density",
      "Ion-weighted mean birth delay [ns]","A_UD for muon–anode",
      "How whole-muon ion-birth timing and hemispheric preference vary across individual anodes.",
      "Both derived from exact original per-wire ion-birth counters, after aggregating over Heed primaries.",
      "Temporal mean is relative to each primary; this is not a strip response."),
 Spec("anode_muon_A_gain","Individual anodes","Muon–anode |A_UD| versus ion yield",
      "anode","ion_births","absA_UD_exact","profile",
      "Ion births on anode in one muon","Mean |A_UD| across muon–anode records",
      "Does a heavily ionizing anode have more or less hemispheric birth localization?",
      "Exact per-(muon,wire) ratio of summed ion-up/down counts; binned by exact anode ion births.",
      "Equal muon–anode record weighting; uncertainty resamples whole muons.",True),
]
BY_KEY={s.key:s for s in SPECS}
SECTIONS=list(dict.fromkeys(s.section for s in SPECS))
COLOR="#1768AC"

@st.cache_data(show_spinner="Loading read-only 2100 V census exports…")
def csv_frame(file: str, modified: int):
    return pd.read_csv(file, low_memory=False)

def discover():
    if not ROOT.is_dir(): return []
    results=[]
    for child in ROOT.iterdir():
        if not child.is_dir() or not (child.name.startswith("population_") or child.name.startswith("pilot_")):
            continue
        required=("v0_muons.csv","v0_primaries.csv","v0_fragments.csv","v0_wires.csv")
        if all((child/"census"/n).is_file() for n in required):
            results.append(child)
    return sorted(results,key=lambda p:(not p.name.startswith("population_"),-p.stat().st_mtime))

def read_one(root:Path, name:str):
    p=root/"census"/name
    return csv_frame(str(p),p.stat().st_mtime_ns).copy()

def ratio(n,d):
    return np.divide(n,d,out=np.full(len(n),np.nan,dtype=float),where=np.asarray(d,dtype=float)>0)

def add_circular(frame):
    s=pd.to_numeric(frame.get("sin"),errors="coerce").to_numpy(dtype=float)
    c=pd.to_numeric(frame.get("cos"),errors="coerce").to_numpy(dtype=float)
    R=np.clip(np.hypot(s,c),0.,1.)
    frame["R_phi"]=np.where(np.isfinite(R)&(R>0),R,np.nan)
    frame["sigma_phi_deg"]=np.degrees(np.sqrt(-2.*np.log(np.maximum(R,1e-15))))
    frame.loc[~np.isfinite(frame["R_phi"]),"sigma_phi_deg"]=np.nan
    if "n_saved" in frame: frame.loc[frame["n_saved"]<3,["sigma_phi_deg","R_phi"]]=np.nan
    if "n_reservoir" in frame: frame.loc[frame["n_reservoir"]<3,["sigma_phi_deg","R_phi"]]=np.nan
    return frame

def prepared(root:Path):
    mu=read_one(root,"v0_muons.csv")
    pr=read_one(root,"v0_primaries.csv")
    fr=read_one(root,"v0_fragments.csv")
    wi=read_one(root,"v0_wires.csv")
    for name,d,required in [
        ("muons",mu,["muon_slot","n_primaries","ion_births_exact","collected_electrons","theta_u_deg","A_exact","channel","station"]),
        ("primaries",pr,["muon_slot","ion_births_exact","collected_electrons","exact_A","sin","cos","start_v_mm","start_du_mm"]),
        ("fragments",fr,["muon_slot","exact_births","exact_A","exact_sigw_mm","primary_ion_births_exact","sin","cos"]),
        ("wires",wi,["muon_slot","global_wire","collected_e_exact"])
    ]:
        missing=set(required)-set(d.columns)
        if missing: raise ValueError(f"Invalid {name} export; missing: {sorted(missing)}")
    if mu["muon_slot"].duplicated().any(): raise ValueError("Duplicate muon IDs in V0 export")
    mu["G_ion"]=ratio(mu["ion_births_exact"].to_numpy(float),mu["n_primaries"].to_numpy(float))
    mu["G_collected"]=ratio(mu["collected_electrons"].to_numpy(float),mu["n_primaries"].to_numpy(float))
    mu["absA_exact"]=mu["A_exact"].abs()
    mu["abs_theta_u"]=mu["theta_u_deg"].abs()
    pr["G_ion"]=pr["ion_births_exact"]
    pr["G_collected"]=pr["collected_electrons"]
    pr["A_exact"]=pd.to_numeric(pr["exact_A"],errors="coerce")
    pr["absA_exact"]=pr["A_exact"].abs()
    fr["G_ion_primary"]=fr["primary_ion_births_exact"]
    fr["absA_exact"]=fr["exact_A"].abs()
    pr=add_circular(pr)
    fr=add_circular(fr)
    wire_stats=[]
    for slot,g in wi.groupby("muon_slot",sort=False):
        q=np.sort(g.loc[g["collected_e_exact"]>0,"collected_e_exact"].to_numpy(float))[::-1]
        wire_stats.append((slot,float(q.sum()),len(q),float(q[1]/q.sum()) if len(q)>1 else (0. if len(q)==1 else np.nan)))
    ws=pd.DataFrame(wire_stats,columns=["muon_slot","wire_collected","n_collecting_wires","f2"])
    mu=mu.merge(ws,on="muon_slot",how="left")
    mu["n_collecting_wires"]=mu["n_collecting_wires"].fillna(0)
    mu["wire_collected"]=mu["wire_collected"].fillna(0)
    mu["wire_census_ok"]=np.isclose(mu["wire_collected"],mu["collected_electrons"],atol=.1,rtol=1e-12)
    mu.loc[(mu["collected_electrons"]==0)&mu["wire_census_ok"],"f2"]=np.nan
    mu["share5"]=np.where(mu["wire_census_ok"],(mu["f2"].fillna(0)>.05).astype(float),np.nan)
    mu.loc[~mu["wire_census_ok"],["f2","n_collecting_wires"]]=np.nan
    # Restore exact parent-muon labels on both hierarchical levels.
    # Keep the two dataframes separate when reassigning after merge.
    for level,d in (("primary",pr),("fragment",fr)):
        d=d.drop(columns=[c for c in ("channel","station","theta_u_deg","abs_theta_u") if c in d])
        d=d.merge(mu[["muon_slot","channel","station","theta_u_deg","abs_theta_u"]],
                  on="muon_slot",how="left",validate="many_to_one")
        if level=="primary": pr=d
        else: fr=d
    # Anode-specific observables are read from the *original* exact Garfield
    # per-wire census, with strict ROOT/CSV closure (not from K256 samples).
    anode_frag,anode=read_per_anode(root,pr,mu)
    return {"mu":mu,"primary":pr,"fragment":fr,"wire":wi,
            "anode_fragment":anode_frag,"anode":anode}

@st.cache_data(show_spinner="Joining the original 2100 V census levels…")
def prepared_cached(directory: str, signatures: tuple):
    # Signatures make data changes explicit, so cached results are not stale.
    return prepared(Path(directory))

def st_filter(datasets,channels,stations):
    mu=datasets["mu"]
    ids=set(mu.loc[mu["channel"].isin(channels)&mu["station"].isin(stations),"muon_slot"])
    return {k:v.loc[v["muon_slot"].isin(ids)].copy() if "muon_slot" in v else v.copy()
            for k,v in datasets.items()}

def numeric_subset(df,spec):
    cols=[spec.x]+([spec.y] if spec.y else [])
    ok=np.ones(len(df),dtype=bool)
    for c in cols:
        vals=pd.to_numeric(df[c],errors="coerce").to_numpy(float)
        ok &= np.isfinite(vals)
    return df.loc[ok,["muon_slot"]+cols].copy()

def bootstrap_profile(frame,spec,bins,log_bins,min_entries,seed=2100):
    x=frame[spec.x].to_numpy(float)
    y=frame[spec.y].to_numpy(float)
    n=len(x)
    if n==0:return pd.DataFrame(),0
    if log_bins: x_ok=x>0
    else:x_ok=np.isfinite(x)
    x=x[x_ok];y=y[x_ok];slots=frame["muon_slot"].to_numpy()[x_ok]
    if not len(x):return pd.DataFrame(),0
    lo,hi=float(x.min()),float(x.max())
    if hi<=lo:hi=lo+max(1e-6,abs(lo)*.01)
    edges=np.geomspace(lo,hi*(1+1e-9),bins+1) if log_bins else np.linspace(lo,hi+max(1e-9,abs(hi)*1e-9),bins+1)
    k=np.searchsorted(edges,x,side="right")-1
    k=np.clip(k,0,bins-1)
    codes,uniq=pd.factorize(slots,sort=False)
    count=np.zeros((len(uniq),bins),dtype=np.float64)
    total=np.zeros_like(count)
    np.add.at(count,(codes,k),1.)
    np.add.at(total,(codes,k),y)
    n_by=count.sum(axis=0);s_by=total.sum(axis=0)
    allowed=(n_by>=min_entries)&((count>0).sum(axis=0)>=min(4,len(uniq)))
    # Cluster bootstrap: resample independent MUONS, retaining their whole set of primaries/fragments.
    rng=np.random.default_rng(seed)
    replicates=np.empty((160,bins),dtype=float)
    for b in range(160):
        ids=rng.integers(0,len(uniq),len(uniq))
        n_boot=count[ids].sum(axis=0)
        s_boot=total[ids].sum(axis=0)
        replicates[b]=np.divide(s_boot,n_boot,out=np.full(bins,np.nan),where=n_boot>0)
    err=np.nanstd(replicates,axis=0,ddof=1)
    centers=np.sqrt(edges[:-1]*edges[1:]) if log_bins else (edges[:-1]+edges[1:])/2
    result=pd.DataFrame({"x":centers,"y":np.divide(s_by,n_by,out=np.full(bins,np.nan),where=n_by>0),
        "error":err,"n":n_by.astype(int),"n_muons":(count>0).sum(axis=0),
        "low":edges[:-1],"high":edges[1:]})
    return result.loc[allowed].copy(),int(n)

def fig_style(fig,title,xlabel,ylabel,custom):
    # Explicit legend typography is required: plotly_white's legend font
    # otherwise overrides the black global font for a single-trace figure.
    legend_positions = {
        "Inside upper right": dict(x=.98, y=.98, xanchor="right", yanchor="top"),
        "Inside upper left": dict(x=.02, y=.98, xanchor="left", yanchor="top"),
        "Inside lower right": dict(x=.98, y=.02, xanchor="right", yanchor="bottom"),
        "Inside lower left": dict(x=.02, y=.02, xanchor="left", yanchor="bottom"),
        "Outside right": dict(x=1.02, y=1., xanchor="left", yanchor="top"),
        "Above plot": dict(x=.5, y=1.02, xanchor="center", yanchor="bottom"),
    }
    legend_position = custom.get("legend_position", "Inside upper right")
    legend_layout = legend_positions[legend_position]
    fig.update_layout(template="plotly_white",title={"text":title,"x":.02,"font":{"size":custom["title_font"],"color":"#000000"}},
        font={"family":"Arial, Helvetica, sans-serif","size":custom["axis_font"],"color":"#000000"},
        paper_bgcolor="white",plot_bgcolor="white",height=630,
        margin=dict(l=75,r=165 if legend_position == "Outside right" and custom["legend"] else 40,t=95,b=78),
        showlegend=custom["legend"],
        legend=dict(bgcolor="rgba(255,255,255,0.92)",borderwidth=0,
                    font=dict(family="Arial, Helvetica, sans-serif",
                              size=custom["axis_font"],color="#000000"),
                    **legend_layout),
        # Force black axis TITLE, numeric TICK LABELS, tick marks, and frame.
        # plotly_white's axis tick/title font defaults can override layout.font.
        xaxis=dict(title=dict(text=xlabel,font=dict(color="#000000")),
                   tickfont=dict(color="#000000"),tickcolor="#000000",
                   showline=True,linewidth=1.5,linecolor="#000000",
                   mirror=True,ticks="inside",showgrid=False,zeroline=False),
        yaxis=dict(title=dict(text=ylabel,font=dict(color="#000000")),
                   tickfont=dict(color="#000000"),tickcolor="#000000",
                   showline=True,linewidth=1.5,linecolor="#000000",
                   mirror=True,ticks="inside",showgrid=False,zeroline=False))
    fig.update_xaxes(type="log" if custom["xscale"]=="Log" else "linear")
    fig.update_yaxes(type="log" if custom["yscale"]=="Log" else "linear")
    if custom["xmin"] is not None and custom["xmax"] is not None:
        a,b=custom["xmin"],custom["xmax"]
        if b>a and (custom["xscale"]!="Log" or a>0):
            fig.update_xaxes(range=[np.log10(a),np.log10(b)] if custom["xscale"]=="Log" else [a,b])
    if custom["ymin"] is not None and custom["ymax"] is not None:
        a,b=custom["ymin"],custom["ymax"]
        if b>a and (custom["yscale"]!="Log" or a>0):
            fig.update_yaxes(range=[np.log10(a),np.log10(b)] if custom["yscale"]=="Log" else [a,b])
    fig.add_annotation(xref="paper",yref="paper",x=1,y=1.095,showarrow=False,
       text="NA60+ Simulation · 2100 V",font=dict(size=11,color="#454545"),xanchor="right")
    return fig

def draw(frame,spec,custom):
    df=numeric_subset(frame,spec)
    omitted=0
    if custom["xscale"]=="Log":
        omitted=int((df[spec.x]<=0).sum())
        df=df.loc[df[spec.x]>0]
    if custom["yscale"]=="Log" and spec.y:
        omitted+=int((df[spec.y]<=0).sum())
        df=df.loc[df[spec.y]>0]
    if df.empty: return None,pd.DataFrame(),0,omitted
    x=df[spec.x].to_numpy(float)
    col=custom["color"];fig=go.Figure(); export=df
    if spec.chart=="hist":
        # Bin *in physical coordinates*; geometric edges give correct log-scale bins.
        lo,hi=float(x.min()),float(x.max())
        if hi<=lo:hi=lo+max(1e-6,abs(lo)*.01)
        edges=(np.geomspace(lo,hi*(1+1e-9),custom["bins"]+1)
               if custom["xscale"]=="Log" else
               np.linspace(lo,hi+max(1e-9,abs(hi)*1e-9),custom["bins"]+1))
        counts,_=np.histogram(x,bins=edges)
        heights=counts.astype(float)
        if custom["normalization"]=="Probability":heights/=max(1.,heights.sum())
        if custom["normalization"]=="Density":heights/=np.maximum(1.,heights.sum())*np.diff(edges)
        centers=np.sqrt(edges[:-1]*edges[1:]) if custom["xscale"]=="Log" else (edges[:-1]+edges[1:])/2
        fig.add_trace(go.Bar(x=centers,y=heights,width=np.diff(edges),
              marker_color=col,marker_line_color="black",marker_line_width=.4,name="Entries"))
        export=pd.DataFrame({"bin_lo":edges[:-1],"bin_hi":edges[1:],
                             "n":counts,"height":heights})
    elif spec.chart=="density":
        yy=df[spec.y].to_numpy(float)
        fig.add_trace(go.Histogram2d(x=x,y=yy,nbinsx=custom["bins"],nbinsy=custom["bins"],
            colorscale="Greys" if custom["mono"] else "Blues",
            colorbar=dict(title="Entries"),name="Population density"))
    elif spec.chart=="profile":
        prof,full=bootstrap_profile(df,spec,custom["bins"],
           custom["xscale"]=="Log",custom["min_bin"])
        export=prof
        if not prof.empty:
            fig.add_trace(go.Scatter(x=prof.x,y=prof.y,
                 error_y=dict(type="data",array=prof.error,visible=custom["errors"],
                              color=col,thickness=1.25,width=3),
                 mode="lines+markers" if custom["lines"] else "markers",
                 line=dict(color=col,width=custom["line_width"]),
                 marker=dict(color=col,size=custom["marker_size"],symbol="circle"),
                 customdata=np.stack([prof.n,prof.n_muons],axis=-1),
                 hovertemplate="%{x:.4g}, %{y:.4g}<br>entries=%{customdata[0]}"
                               "<br>muons=%{customdata[1]}<extra></extra>",
                 name="Muon-bootstrap mean"))
    fig=fig_style(fig,custom["title"],custom["xlabel"],custom["ylabel"],custom)
    return fig,export,int(len(df)),omitted

def plotted_values(frame,spec,custom):
    """The same finite, scale-eligible ORIGINAL observations used by draw().
    These are unbinned physical values, never histogram centers or profile means.
    """
    df=numeric_subset(frame,spec)
    if custom["xscale"]=="Log":
        df=df.loc[df[spec.x]>0]
    if custom["yscale"]=="Log" and spec.y:
        df=df.loc[df[spec.y]>0]
    return df


def mean_and_rms(values):
    """ROOT TH1 GetMean/GetRMS convention: RMS is about the mean (ddof=0)."""
    a=np.asarray(values,dtype=np.float64)
    a=a[np.isfinite(a)]
    if not len(a): return float("nan"),float("nan")
    mean=float(np.mean(a))
    # Centre before squaring to avoid catastrophic cancellation.
    rms=float(np.sqrt(np.mean(np.square(a-mean))))
    return mean,rms


def format_plot_stats(frame,spec,custom):
    """Text appended to a legend entry, computed before plotting/binning.

    For histograms, stats describe x; for profile/scatter plots, stats describe
    y (original, unbinned observations); for 2D densities both x and y.
    Numeric moments stay in the original physical units even on logarithmic axes.
    """
    df=plotted_values(frame,spec,custom)
    n=len(df)
    if not n:return ""
    if spec.chart=="hist":
        axes=[("x",spec.x)]
    elif spec.chart=="density":
        axes=[("x",spec.x),("y",spec.y)]
    else:
        axes=[("y",spec.y)]
    lines=[f"N = {n:,}"]
    for axis,col in axes:
        mean,rms=mean_and_rms(df[col].to_numpy())
        prefix=f"{axis}: " if len(axes)>1 else ""
        lines.append(f"{prefix}Mean = {mean:.4g}")
        lines.append(f"{prefix}RMS = {rms:.4g}")
    return "<br>".join(lines)


def as_html(fig):
    return fig.to_html(full_html=True,include_plotlyjs=True).encode("utf-8")

def main():
    st.title("NA60+ MWPC · 2100 V Plot Studio")
    st.caption("Interactive plots from immutable Garfield++ microscopic archives · ATLAS/ROOT-style defaults · no re-simulation")
    found=discover()
    if not found:
        st.error(f"No completed ROOT/C++ census outputs found in {ROOT}. Run the five-muon pilot or 997-muon population analysis first.")
        st.code("bash test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_gain_asym_2100.sh pilot")
        st.stop()
    st.sidebar.header("Dataset and scope")
    paths={p.name:p for p in found}
    selected=st.sidebar.selectbox("Saved analysis",list(paths),index=0)
    folder=paths[selected]
    required_files=("v0_muons.csv","v0_primaries.csv","v0_fragments.csv","v0_wires.csv")
    signatures=tuple((folder/"census"/n).stat().st_mtime_ns for n in required_files)
    try: datasets=prepared_cached(str(folder),signatures)
    except Exception as exc:
        st.error(f"Could not load this analysis: {exc}")
        st.stop()
    allmu=datasets["mu"]
    channels=sorted(allmu["channel"].dropna().unique().tolist())
    stations=sorted(allmu["station"].dropna().unique().tolist())
    selected_channels=st.sidebar.multiselect("Geant4 source channel",channels,default=channels)
    selected_stations=st.sidebar.multiselect("Detector station",stations,default=stations)
    if not selected_channels or not selected_stations:st.info("Select at least one channel and station.");st.stop()
    data=st_filter(datasets,selected_channels,selected_stations)
    sec=st.selectbox("Physics section",SECTIONS,index=SECTIONS.index("Individual anodes"))
    if sec=="Individual anodes":
        st.sidebar.subheader("Individual anode selection")
        detector_ids=sorted(data["mu"]["det_id"].dropna().unique().tolist())
        selected_det=st.sidebar.multiselect("Detector identifier",detector_ids,
                                            default=detector_ids)
        allowed_ids=set(data["mu"].loc[data["mu"]["det_id"].isin(selected_det),
                                         "muon_slot"])
        data={k:(d.loc[d["muon_slot"].isin(allowed_ids)].copy()
                 if "muon_slot" in d else d.copy()) for k,d in data.items()}
        candidates=sorted(data["anode"]["global_wire"].dropna().astype(int).unique().tolist())
        anode_idx=st.sidebar.selectbox("Global anode ID",
                        ["All anodes"]+candidates,format_func=str)
        if anode_idx!="All anodes":
            for lev in ("anode_fragment","anode"):
                data[lev]=data[lev].loc[data[lev]["global_wire"]==anode_idx].copy()
        st.sidebar.caption("The global wire index must be interpreted within the "
                           "chosen station/detector. Source ions are assigned to "
                           "their nearest anode at birth; collected electrons refer "
                           "to their actual collecting wire.")
    mu=data["mu"];pr=data["primary"];fr=data["fragment"]
    count_mu=len(mu);count_pr=len(pr)
    if not count_mu:st.info("No muons match these filters.");st.stop()
    total_heed=float(mu["n_primaries"].sum())
    avg_ion=float(mu["ion_births_exact"].sum()/total_heed) if total_heed else np.nan
    avg_col=float(mu["collected_electrons"].sum()/total_heed) if total_heed else np.nan
    st.sidebar.caption(f"Source: {folder}")
    st.sidebar.info("Descriptive data, including previous holdout slots. NOT FOR MODEL CALIBRATION. Exact ion-birth and anode census; no 997-muon ion-drift replay.")
    if "pilot_" in selected:st.warning("PILOT ONLY: five-muon figures are software checks, not population conclusions. Select a completed population_* run for final results.")
    k1,k2,k3,k4=st.columns(4)
    k1.metric("Muons",f"{count_mu:,}")
    k2.metric("Heed primaries",f"{count_pr:,}")
    k3.metric("Global ion gain",f"{avg_ion:,.1f}")
    k4.metric("Global collected gain",f"{avg_col:,.1f}")
    bad=int((~mu["wire_census_ok"]).sum())
    if bad:st.warning(f"Physical collecting-wire totals do not close for {bad} muons. Those entries are excluded from wire-sharing observables.")
    in_section=[s for s in SPECS if s.section==sec]
    spec=BY_KEY[st.selectbox("Observable",[s.key for s in in_section],
                            format_func=lambda key:BY_KEY[key].title,
                            index=next((i for i,s in enumerate(in_section) if s.key=="absA_gain"),0))]
    st.sidebar.divider()
    st.sidebar.header("Plot editor")
    name=spec.key
    title=st.sidebar.text_input("Plot title",spec.title,key=f"title_{name}")
    xlabel=st.sidebar.text_input("Horizontal axis title",spec.xlabel,key=f"xlabel_{name}")
    ylabel=st.sidebar.text_input("Vertical axis title",spec.ylabel,key=f"ylabel_{name}")
    cx,cy=st.sidebar.columns(2)
    xscale=cx.selectbox("X scale",["Linear","Log"],index=int(spec.logx),key=f"xscale_{name}")
    yscale=cy.selectbox("Y scale",["Linear","Log"],index=int(spec.logy),key=f"yscale_{name}")
    bins=st.sidebar.slider("Histogram / profile bins",8,100,30,2,key=f"bins_{name}")
    minbin=st.sidebar.number_input("Minimum entries per profile bin",4,1000,20,step=2,key=f"minbin_{name}")
    col=st.sidebar.color_picker("Plot color",COLOR,key=f"color_{name}")
    mono=st.sidebar.checkbox("Monochrome density",False,key=f"mono_{name}")
    lines=st.sidebar.checkbox("Join profile points",True,key=f"lines_{name}")
    errors=st.sidebar.checkbox("Show muon-bootstrap error bars",True,key=f"errors_{name}")
    marker_size=st.sidebar.slider("Marker size",4,15,7,key=f"marker_{name}")
    lw=st.sidebar.slider("Line width",1.,5.,2.,.5,key=f"lw_{name}")
    normalization=st.sidebar.selectbox("Histogram normalization",["Counts","Probability","Density"],key=f"norm_{name}")
    legend=st.sidebar.checkbox("Show legend",True,key=f"legend_{name}")
    show_stats=st.sidebar.checkbox("Show N, mean and RMS in legend",True,
                                   key=f"legendstats_{name}",
                                   help="Uses the original unbinned values (after filtering). For profiles this is y, not the plotted bin means. RMS is the spread around the mean.")
    legend_text=st.sidebar.text_input("Legend entry", "Entries" if spec.chart!="profile" else "Muon-bootstrap mean",
                                       key=f"legendname_{name}")
    legend_position=st.sidebar.selectbox("Legend position",
        ["Inside upper right","Inside upper left","Inside lower right",
         "Inside lower left","Outside right","Above plot"],
        key=f"legendpos_{name}")
    with st.sidebar.expander("Fonts and axis limits"):
        title_font=st.slider("Title font size",14,34,21,key=f"titlefont_{name}")
        axis_font=st.slider("Axes font size",10,22,15,key=f"axisfont_{name}")
        xmin=st.number_input("X min (0 = auto)",value=0.,format="%.6g",key=f"xmin_{name}")
        xmax=st.number_input("X max (0 = auto)",value=0.,format="%.6g",key=f"xmax_{name}")
        ymin=st.number_input("Y min (0 = auto)",value=0.,format="%.6g",key=f"ymin_{name}")
        ymax=st.number_input("Y max (0 = auto)",value=0.,format="%.6g",key=f"ymax_{name}")
    custom=dict(title=title,xlabel=xlabel,ylabel=ylabel,xscale=xscale,yscale=yscale,bins=bins,
        min_bin=minbin,color=col,mono=mono,lines=lines,errors=errors,marker_size=marker_size,
        line_width=lw,normalization=normalization,legend=legend,
        show_stats=show_stats,legend_position=legend_position,title_font=title_font,
        axis_font=axis_font,xmin=xmin if xmin!=0 or xmax!=0 else None,
        xmax=xmax if xmin!=0 or xmax!=0 else None,
        ymin=ymin if ymin!=0 or ymax!=0 else None,
        ymax=ymax if ymin!=0 or ymax!=0 else None)
    st.subheader(spec.title)
    fig,displayed,n,omitted=draw(data[spec.dataset],spec,custom)
    if fig is None:st.warning("No finite data for these filters and axis scales.");st.stop()
    if fig.data:
        details=format_plot_stats(data[spec.dataset],spec,custom) if show_stats else ""
        plot_legend=(legend_text or "Entries")+(("<br>"+details) if details else "")
        if spec.chart=="density":
            # Histogram2d uses a colorbar; add a separate legend-only marker
            # so its exact underlying x/y moments are visible and exportable.
            fig.data[0].showlegend=False
            if legend:
                fig.add_trace(go.Scatter(x=[None],y=[None],mode="markers",
                    marker=dict(color=col,size=8,symbol="square"),
                    name=plot_legend,showlegend=True,hoverinfo="skip"))
        else:
            fig.update_traces(name=plot_legend,showlegend=legend)
    left,right=st.columns([3.6,1.45])
    with left:
        st.plotly_chart(fig,use_container_width=True,config={"displaylogo":False,"scrollZoom":True})
        st.caption(f"Observations plotted: {n:,} / {len(data[spec.dataset]):,} available {spec.dataset} rows. "+
            (f"{omitted:,} nonpositive coordinate(s) omitted for logarithmic axes. " if omitted else "")+
            ("Profiles: 160 muon-block bootstrap resamples; bins need ≥4 distinct muons." if spec.chart=="profile" else ""))
    with right:
        st.markdown("##### What does this plot show?")
        st.write(spec.meaning)
        st.markdown("##### How is it calculated?")
        st.write(spec.calculation)
        st.markdown("##### Important limitations")
        st.write(spec.limitations)
        if spec.chart=="profile":st.caption("Error bars are muon-cluster bootstrap standard deviations, not fit errors.")
        st.caption("Numerical definitions: ion *births* at Garfield callback; electron endpoints at actual physical anodes.")
    if sec=="Individual anodes":
        with st.expander("Per-anode summary table — exact Garfield counters",expanded=True):
            st.caption("One row per (station, detector, global wire). Statistics across "
                       "muon–anode records are equal-record means; summed ion/electron "
                       "counts are exact. Mean delays are referenced to the corresponding "
                       "Heed primary, not to cathode collection.")
            adata=data["anode"]
            if len(adata):
                groupkeys=["station","det_id","global_wire"]
                table=(adata.groupby(groupkeys,dropna=False)
                       .agg(n_muon_anode_records=("muon_slot","nunique"),
                            ion_births_total=("ion_births","sum"),
                            collected_e_total=("collected_e","sum"),
                            mean_A_UD=("A_UD_exact","mean"),
                            mean_absA_UD=("absA_UD_exact","mean"),
                            mean_birth_delay_ns=("ion_birth_delay_ns","mean"),
                            mean_birth_sigma_ns=("sigma_t_birth_relative_ns","mean"))
                       .reset_index().sort_values(["station","det_id","global_wire"]))
                st.dataframe(table,use_container_width=True,hide_index=True)
                st.download_button("Download exact per-anode summary CSV",
                    table.to_csv(index=False).encode("utf-8"),
                    file_name=f"NA60_2100_anodes_{selected}.csv",mime="text/csv")
            else:
                st.info("No exact anode-census rows match these filters.")
    with st.expander("Population statistics and provenance",expanded=False):
        st.write({"selected_muons":count_mu,"Heed_primaries":count_pr,
                  "zero_ion_primaries":int((pr["ion_births_exact"]==0).sum()),
                  "zero_collected_muons":int((mu["collected_electrons"]==0).sum()),
                  "ratio_of_sums_G_ion":avg_ion,"ratio_of_sums_G_collected":avg_col,
                  "good_collecting_wire_census_muons":int(mu["wire_census_ok"].sum())})
        meta=folder/"census"/"v0_manifest.txt"
        if meta.exists():st.code(meta.read_text()[:12000])
    st.markdown("##### Export plot and underlying displayed data")
    a,b,c=st.columns(3)
    with a:st.download_button("Interactive HTML",as_html(fig),file_name=f"NA60_2100_{name}.html",mime="text/html")
    with b:
        csv_bytes=displayed.to_csv(index=False).encode("utf-8")
        st.download_button("Plot values CSV",csv_bytes,file_name=f"NA60_2100_{name}.csv",mime="text/csv")
    with c:
        st.caption("PNG/PDF via 'Prepare image exports' below.")
    if st.button("Prepare image exports (PNG and vector PDF)",key=f"export_{name}"):
        try:
            png=fig.to_image(format="png",width=1200,height=850,scale=2)
            pdf=fig.to_image(format="pdf",width=1200,height=850)
            st.download_button("Download PNG",data=png,file_name=f"NA60_2100_{name}.png",mime="image/png")
            st.download_button("Download vector PDF",data=pdf,file_name=f"NA60_2100_{name}.pdf",mime="application/pdf")
        except Exception as exc:
            st.error(f"Static export requires compatible Kaleido/Chromium: {exc}. HTML and CSV remain available.")
    with st.expander("Unavailable directly from these compact CSV exports"):
        st.write("A full 997-muon ion-azimuth histogram, individual primary-relative delta-w distributions, and all-ion drift/current trajectories are not present in the compact census. The ROOT K256 ion-birth rows can support a separately validated weighted-angle export; that is not the same as an exact per-ion census. No top/bottom sectors or 1800 V comparisons are displayed here.")
if __name__=="__main__":
    main()

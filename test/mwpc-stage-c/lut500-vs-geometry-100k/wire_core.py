"""Pure, tested transport/charge-sharing primitives for NA60+ P3 noFW.

IMPORTANT: LUT wire offsets are relative to each electron's *own* source cell.
For long, oblique tracks, the global wire indices may span arbitrarily many
physical anodes even though a single LUT endpoint is in {-1,0,+1}.
"""
from __future__ import annotations
import hashlib
import math
import numpy as np

PITCH = 3.0
HALF_GAP = 3.0
# Native Geant4 chamber gas-u aperture widths in mm, from
# pi_six_geometry_native_timing.py::WIDTHS.
STATION_WIDTHS = (307.2, 512.0, 512.0, 512.0, 512.0, 512.0)
THRESHOLDS = (0.01, 0.05, 0.10)


def stable_hash(*parts: object, person: bytes = b'mwpc-multianode') -> int:
    payload = '\x1f'.join(map(str, parts)).encode('utf-8')
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8, person=person).digest(), 'little')


def rng_for(*parts: object) -> np.random.Generator:
    return np.random.default_rng(stable_hash(*parts))


def n_anodes(station: int, pitch: float = PITCH) -> int:
    if station not in range(6):
        raise ValueError(f'Bad station: {station}')
    return int(math.floor(STATION_WIDTHS[station]/pitch + 0.5))


def nearest_wire(u_mm: np.ndarray, station: int, pitch: float = PITCH):
    """Native finite-lattice indexing, including MS0 half-pitch phase.

    Finite anode indices 0..n-1. For p=3, MS0 has an even count
    (centres shifted by 1.5mm relative to MS1/3/5, with local u origin 0).
    """
    u = np.asarray(u_mm, dtype=float)
    n = n_anodes(station,pitch)
    pos = u/pitch + (n-1)/2
    k = np.floor(pos + 0.5).astype(np.int64)
    local = u - (k-(n-1)/2)*pitch
    if np.any(np.abs(local) > pitch/2 + 1e-8):
        raise AssertionError('Cell folding failed')
    return k, local


def cell_ids(local_u_mm: np.ndarray, v_mm: np.ndarray):
    u = np.asarray(local_u_mm, dtype=float)
    v = np.asarray(v_mm, dtype=float)
    if u.shape != v.shape:
        raise ValueError('u/v lengths differ')
    if np.any(~np.isfinite(u)) or np.any(~np.isfinite(v)):
        raise ValueError('nonfinite seed coordinates')
    if np.any(np.abs(u)>1.500001) or np.any(np.abs(v)>3.05):
        raise ValueError('seed outside LUT support')
    iu = np.clip(np.floor((u+1.5)/0.1).astype(np.int32),0,29)
    iv = np.clip(np.floor((v+3.)/0.2).astype(np.int32),0,29)
    return iu*30 + iv


def build_lut_kernel(iu, iv, collected, delta):
    """Make a conditional-on-collection P(delta|cell), 900x3."""
    iu = np.asarray(iu,dtype=int); iv=np.asarray(iv,dtype=int)
    ok = np.asarray(collected,dtype=bool)
    delta = np.asarray(delta,dtype=int)
    if len(iu)!=len(iv) or len(iv)!=len(ok) or len(ok)!=len(delta):
        raise ValueError('LUT length mismatch')
    if np.any((iu<0)|(iu>=30)|(iv<0)|(iv>=30)):
        raise ValueError('Unexpected LUT grid')
    if np.any(~np.isin(delta[ok],(-1,0,1))):
        raise ValueError('LUT endpoints outside supported offsets -1,0,+1')
    cell = iu[ok]*30+iv[ok]
    counts=np.zeros((900,3),dtype=np.int64)
    np.add.at(counts,(cell,delta[ok]+1),1)
    denom=counts.sum(axis=1)
    if np.any(denom==0):
        raise ValueError(f'LUT missing collected trials for {int((denom==0).sum())} cells')
    return counts/denom[:,None]


def draw_offsets(cells: np.ndarray, kernel: np.ndarray, rng):
    cells=np.asarray(cells,dtype=np.int32)
    if kernel.shape!=(900,3): raise ValueError('Expected 900x3 kernel')
    probabilities=kernel[cells]
    # Exactly -1,0,+1; never confuse these with physical wire IDs.
    x=rng.random(len(cells))
    return np.where(x<probabilities[:,0],-1,
           np.where(x<probabilities[:,0]+probabilities[:,1],0,1)).astype(np.int64)


def physical_wire_charges(wire_id,weights,station):
    wire=np.asarray(wire_id,dtype=np.int64)
    weights=np.asarray(weights,dtype=np.float64)
    if len(wire)!=len(weights) or np.any(weights<0) or not np.all(np.isfinite(weights)):
        raise ValueError('Bad wire/weights')
    n=n_anodes(station)
    valid=(wire>=0)&(wire<n)&(weights>0)
    dropped=float(weights[(~valid)&(weights>0)].sum())
    if not np.any(valid): return {},dropped
    ids,q=np.unique(wire[valid],return_inverse=True)
    total=np.bincount(q,weights=weights[valid],minlength=len(ids))
    return {int(k):float(v) for k,v in zip(ids,total)},dropped


def observables(wires: dict[int,float]) -> dict:
    q=np.asarray(sorted(wires.values(),reverse=True),dtype=float)
    total=float(q.sum())
    if total==0:
        return {'n_positive':0,'qtotal':0.,'q2_fraction':0.,'other_fraction':0.,
                **{f'n_gt{int(100*t)}pct':0 for t in THRESHOLDS}}
    # q2 = absolute charge fraction on SECOND strongest physical anode.
    second = float(q[1]/total) if len(q)>1 else 0.
    result={'n_positive':int(len(q)), 'qtotal':total,
            'q2_fraction':second,'other_fraction':float(1-q[0]/total)}
    for t in THRESHOLDS:
        result[f'n_gt{int(100*t)}pct']=int(np.count_nonzero(q>t*total))
    return result


def compare_one(u,v,weights,station,kernel,*,rng):
    """Both assignments use IDENTICAL electron starts and charge weights.

    Geometry: source-cell nearest-wire midline for each primary electron.
    LUT: one stochastic collected-wire offset for each same primary electron.
    The use of identical weights isolates wire assignment from gain fluctuations.
    """
    u=np.asarray(u,dtype=float); v=np.asarray(v,dtype=float)
    q=np.asarray(weights,dtype=float)
    if len(u)!=len(v) or len(u)!=len(q): raise ValueError('Muon seed arrays differ')
    source,local=nearest_wire(u,station)
    cells=cell_ids(local,v)
    offsets=draw_offsets(cells,kernel,rng)
    geo,geo_edge=physical_wire_charges(source,q,station)
    lut,lut_edge=physical_wire_charges(source+offsets,q,station)
    total=float(q.sum())
    if not np.isclose(sum(geo.values())+geo_edge,total,rtol=1.e-12,atol=1.e-9):
        raise AssertionError("Geometric physical-wire charge not conserved")
    if not np.isclose(sum(lut.values())+lut_edge,total,rtol=1.e-12,atol=1.e-9):
        raise AssertionError("LUT physical-wire charge not conserved")
    return geo,lut,geo_edge,lut_edge
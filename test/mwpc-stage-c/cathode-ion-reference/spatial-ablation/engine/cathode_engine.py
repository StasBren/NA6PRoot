"""Reusable cathode U/V response engine: birth records / cloud moments + an identified kernel -> unshaped signals.

Physics is exactly the validated K-model of s_kernel.py (frozen 2100 V transport paths, frozen spectral weighting LUT,
K0, Kc1, Ks1, Kc2, Ks2 and the hemisphere term). This module changes no physical model; it packages:

  KernelTable           a precomputed kernel directory, opened ONCE per process, with an explicit identity
                        (sha256 of K.npy, KPHI0.npy, kernel.json) checked against an expected id, plus the ids of the
                        transport paths and the weighting LUT it was built from (so a strip-geometry study can show
                        that only the weighting changed)
  Moments               per-cloud moments (from p3_birth_moments CSV, or from in-memory births)
  moments_from_births   numpy implementation of BirthMoments.hh (same definitions; tested against it)
  CathodeResponseEngine evaluate(moments) -> Response(ion, electron) on a 10 ns grid; implementation "vector"
                        (default) or "loop" (the validated per-cloud reference algorithm, kept for regression)
  Electronics           interchangeable shaper interface (CRRCn stand-in, impulse table); applied AFTER the engine
No development-muon index, run directory or S0 file is used in normal operation.
"""
from __future__ import annotations

import csv
import ctypes
import hashlib
import os
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

EFC = 1.602176634e-4          # fC per elementary charge
R_WIRE_UM = 15.0
MODELS = {"K1": (1, False), "K2": (2, False), "K1A": (1, True), "K2A": (2, True)}
MOMENT_COLS = ("primary", "wire", "t_ns", "N_fC", "r_eff_um", "D_U_mm", "D_V_mm", "C1", "S1", "C2", "S2", "A")


def _sha256(p: Path, chunk=1 << 24) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def s_perp(th_m, nmax):
    sm = np.sign(np.sin(th_m))
    sp = sm - sm.mean()
    for n in range(1, nmax + 1):
        for b in (np.cos(n * th_m), np.sin(n * th_m)):
            sp = sp - b * (sp @ b) / (b @ b)
    return sm, sp


class KernelTable:
    """A built kernel directory (s_kernel.py build). Opened once; read-only."""

    def __init__(self, path, expected_id: str | None = None, verify: bool = True, preload: bool = True):
        t0 = time.perf_counter()
        self.path = Path(path)
        self.meta = json.loads((self.path / "kernel.json").read_text())
        ax = self.meta["green_axes"]
        self.kernel_id = None
        if verify or expected_id:
            h = hashlib.sha256()
            for f in ("K.npy", "KPHI0.npy", "kernel.json"):
                h.update(_sha256(self.path / f).encode())
            self.kernel_id = h.hexdigest()
            if expected_id and expected_id != self.kernel_id:
                raise ValueError(f"kernel id mismatch: expected {expected_id}, found {self.kernel_id}")
        self.verify_s = time.perf_counter() - t0
        self.transport_id = ax.get("grid_vertices_sha256")      # trajectories (gas, field, wire geometry)
        self.weighting_id = ax.get("lut_sha256")                # electrode weighting (strip geometry)
        self.names = self.meta["basis"]
        self.wires = np.array(ax["wires"])
        self.rn = np.array(ax["r_um"]); self.D = np.array(ax["D_mm"]); self.tau = np.array(ax["tau_ns"])
        self.M, self.dD, self.dsub = len(ax["theta_deg"]), ax["dD_mm"], ax["d_sub"]
        self.pitch_mm, self.stereo_tan = ax["pitch_mm"], ax["stereo_tan"]
        K = np.load(self.path / "K.npy", mmap_mode="r")
        self.K = np.ascontiguousarray(K) if preload else K       # (basis, side, wire, r, D, tau) float32
        self.KP = np.load(self.path / "KPHI0.npy")                # (basis, side, wire, r, D) float64
        self.load_s = time.perf_counter() - t0 - self.verify_s
        self.nbytes = int(self.K.nbytes if preload else 0) + int(self.KP.nbytes)

    def describe(self):
        return {"path": str(self.path), "kernel_id": self.kernel_id, "transport_id": self.transport_id,
                "weighting_id": self.weighting_id, "basis": self.names, "shape": list(self.K.shape),
                "verify_s": self.verify_s, "load_s": self.load_s, "resident_bytes": self.nbytes}


@dataclass
class Moments:
    cols: dict
    base_U: int
    base_V: int
    source: str = ""

    @property
    def n(self):
        return len(self.cols["N_fC"])


def read_moments_csv(path) -> Moments:
    rows = list(csv.DictReader(open(path, newline="")))
    cols = {c: np.array([float(r[c]) for r in rows]) for c in MOMENT_COLS}
    return Moments(cols, int(rows[0]["base_U"]), int(rows[0]["base_V"]), str(path))


def moments_from_births(primary, u_cm, v_cm, w_cm, t_ns, weight, patch_shift_cm, u_in_mm, w_in_mm) -> Moments:
    """numpy version of BirthMoments.hh (T1 per primary; cloud = (primary, wire); same frame formulas)."""
    primary = np.asarray(primary, int); w = np.asarray(weight, float)
    uIn, wIn = .1 * u_in_mm, .1 * w_in_mm
    llround = lambda x: int(np.sign(x) * np.floor(abs(x) + 0.5))      # C++ std::llround (half away from zero)
    bu, bv = llround((wIn - .1 * uIn) / .17), llround((wIn + .1 * uIn) / .17)
    up, inv = np.unique(primary, return_inverse=True)
    tm = (np.bincount(inv, w * t_ns) / np.bincount(inv, w))[inv]
    t = t_ns + (tm - t_ns)
    u, v = 10 * np.asarray(u_cm), 10 * np.asarray(v_cm)
    wu = 10. * (np.asarray(w_cm) + wIn - .1 * patch_shift_cm - .17 * bu)
    wv = 10. * (np.asarray(w_cm) + wIn + .1 * patch_shift_cm - .17 * bv)
    k = np.rint(u / 3).astype(int)
    du = u - 3 * k
    r = 1000 * np.hypot(du, v)
    if np.any(r <= R_WIRE_UM):
        raise ValueError("ion birth on or inside the anode")
    th = np.arctan2(v, du)
    q = w * EFC
    keys, cl = np.unique(np.stack([primary, k], 1), axis=0, return_inverse=True)
    cl = cl.ravel()
    Nq = np.bincount(cl, q)
    m = lambda x: np.bincount(cl, q * x) / Nq
    cols = {"primary": keys[:, 0].astype(float), "wire": keys[:, 1].astype(float), "t_ns": m(t), "N_fC": Nq,
            "r_eff_um": np.exp(m(np.log(r))), "D_U_mm": m(wu), "D_V_mm": m(wv), "C1": m(np.cos(th)),
            "S1": m(np.sin(th)), "C2": m(np.cos(2 * th)), "S2": m(np.sin(2 * th)), "A": m(np.sign(np.sin(th)))}
    return Moments(cols, bu, bv, "births")


@dataclass
class Response:
    ion: np.ndarray            # (bins, 2 planes, strips) current, fC/ns, unshaped, frozen sign convention
    electron: np.ndarray       # same layout: APPROXIMATE instantaneous avalanche-electron step (current in its bin)
    labels: list
    t_center_ns: np.ndarray
    dt_ns: float
    meta: dict = field(default_factory=dict)

    def write_replay_format(self, out: Path, name: str):
        out = Path(out)
        for sub, arr, fn in (("replay", self.ion, "waveforms.csv"), ("prompt", self.electron, "prompt_waveforms.csv")):
            (out / name / sub).mkdir(parents=True, exist_ok=False)
            np.savetxt(out / name / sub / fn, np.column_stack((self.t_center_ns, arr.reshape(len(arr), -1))),
                       delimiter=",", header="time_center_ns," + ",".join(self.labels), comments="")


class CathodeResponseEngine:
    def __init__(self, kernel: KernelTable, model: str = "K2A", half_strips: int = 5, window_ns: float = 7000.0,
                 dt_ns: float = 10.0, chunk_clouds: int = 64):
        if model not in MODELS:
            raise ValueError(f"model {model} not in {list(MODELS)}")
        self.k, self.model = kernel, model
        self.nmax, self.hemi = MODELS[model]
        self.offsets = np.arange(-half_strips, half_strips + 1)
        self.dt = dt_ns
        self.nb = int(round(window_ns / dt_ns))
        if kernel.tau[-1] < window_ns:
            raise ValueError("kernel tau range shorter than the window")
        self.edges = np.arange(self.nb + 1) * dt_ns
        self.chunk = chunk_clouds
        th = (np.arange(kernel.M) + 0.5) * 2 * np.pi / kernel.M
        self.sm, self.sp = s_perp(th, self.nmax)
        self.cm = {n: (np.cos(n * th) @ self.sm, np.sin(n * th) @ self.sm) for n in (1, 2)}
        basis = ["K0"] + [f"K{c}{n}" for n in range(1, self.nmax + 1) for c in ("c", "s")] + ([f"KA{self.nmax}"] if self.hemi else [])
        self.bidx = [kernel.names.index(b) for b in basis]
        self.basis = basis

    # --- per-cloud coefficients and interpolation indices (shared by both implementations) ---
    def _prepare(self, m: Moments):
        k, X = self.k, m.cols
        if not set(X["wire"].astype(int)) <= set(k.wires.tolist()):
            raise ValueError("cloud on a wire outside the kernel table")
        if X["r_eff_um"].min() < k.rn[0] or X["r_eff_um"].max() > k.rn[-1]:
            raise ValueError("r_eff outside the radial grid (no extrapolation)")
        n = m.n
        coef = [np.ones(n)]
        a0s = self.sm.sum() / k.M * np.ones(n)
        for nn in range(1, self.nmax + 1):
            coef += [X[f"C{nn}"], X[f"S{nn}"]]
            a0s = a0s + (2 / k.M) * (X[f"C{nn}"] * self.cm[nn][0] + X[f"S{nn}"] * self.cm[nn][1])
        if self.hemi:
            coef.append((X["A"] - a0s) / (self.sp @ self.sm))
        C = np.stack(coef, 1) * X["N_fC"][:, None]
        wk = np.searchsorted(k.wires, X["wire"].astype(int))
        lr, ln = np.log(X["r_eff_um"]), np.log(k.rn)
        jr = np.clip(np.searchsorted(ln, lr) - 1, 0, len(k.rn) - 2)
        fr = (lr - ln[jr]) / (ln[jr + 1] - ln[jr])
        D = {}
        for si, s in enumerate(("U", "V")):
            x = (X[f"D_{s}_mm"] - k.D[0]) / k.dD
            jd = np.floor(x).astype(int)
            if (jd - k.dsub * self.offsets.max()).min() < 0 or (jd + 1 - k.dsub * self.offsets.min()).max() >= len(k.D):
                raise ValueError("strip-relative D outside the kernel table")
            D[si] = (jd, x - jd)
        return C, wk, jr, fr, D

    def evaluate(self, m: Moments, impl: str = "vector") -> Response:
        t0 = time.perf_counter()
        C, wk, jr, fr, D = self._prepare(m)
        ion = np.zeros((self.nb, 2, len(self.offsets)))
        ele = np.zeros_like(ion)
        tc = m.cols["t_ns"]
        {"loop": self._loop, "vector": self._vector, "compiled": self._compiled}[impl](C, wk, jr, fr, D, tc, ion, ele)
        labels = [f"U_{m.base_U + o}" for o in self.offsets] + [f"V_{m.base_V + o}" for o in self.offsets]
        return Response(ion / self.dt, ele / self.dt, labels, self.edges[:-1] + self.dt / 2, self.dt,
                        {"model": self.model, "impl": impl, "clouds": m.n, "kernel_id": self.k.kernel_id,
                         "transport_id": self.k.transport_id, "weighting_id": self.k.weighting_id,
                         "evaluate_s": time.perf_counter() - t0})

    def _loop(self, C, wk, jr, fr, D, tc, ion, ele):
        """The validated s_kernel.py algorithm, unchanged in order of operations (regression reference)."""
        k, tau, bidx = self.k, self.k.tau, self.bidx
        for si in range(2):
            jd, fd = D[si]
            for c in range(len(tc)):
                Kc = k.K[bidx, si, wk[c]]
                b0 = int(np.floor(tc[c] / self.dt))
                for oi, o in enumerate(self.offsets):
                    d0 = jd[c] - k.dsub * o
                    Q = np.zeros(len(tau)); P = 0.0
                    for dr, wr_ in ((0, 1 - fr[c]), (1, fr[c])):
                        for dd, wd_ in ((0, 1 - fd[c]), (1, fd[c])):
                            wgt = wr_ * wd_
                            Q += wgt * (C[c] @ np.asarray(Kc[:, jr[c] + dr, d0 + dd], dtype=np.float64))
                            P += wgt * (C[c] @ k.KP[bidx, si, wk[c], jr[c] + dr, d0 + dd])
                    qe = np.interp(self.edges - tc[c], tau, Q, left=0.0, right=Q[-1])
                    ion[:, si, oi] += np.diff(qe)
                    if b0 < self.nb:
                        ele[b0, si, oi] += P

    _lib = None

    @classmethod
    def compiled_library(cls):
        """libkernel_eval.so (build: g++ -O3 -std=c++17 -shared -fPIC -o libkernel_eval.so kernel_eval.cc)."""
        if cls._lib is None:
            p = os.environ.get("NA60_KERNEL_EVAL_LIB", str(Path(__file__).resolve().parent / "libkernel_eval.so"))
            lib = ctypes.CDLL(p)
            lib.kernel_eval.restype = ctypes.c_int
            cls._lib = lib
        return cls._lib

    def _compiled(self, C, wk, jr, fr, D, tc, ion, ele):
        k = self.k
        if not (isinstance(k.K, np.ndarray) and k.K.flags.c_contiguous and k.K.dtype == np.float32):
            raise ValueError("compiled implementation needs a preloaded, contiguous float32 kernel (preload=True)")
        lib = self.compiled_library()
        keep = []                                         # keep temporaries alive during the call
        def arr(a, t):
            x = np.ascontiguousarray(a, dtype=t); keep.append(x); return x.ctypes.data_as(ctypes.c_void_p)
        nbA, _, nw, nr, nD, nt = k.K.shape
        out_i = np.zeros((self.nb, 2, len(self.offsets))); out_e = np.zeros_like(out_i)
        rc = lib.kernel_eval(arr(k.K, np.float32), arr(k.KP, np.float64), ctypes.c_int(nbA), ctypes.c_int(nw),
                             ctypes.c_int(nr), ctypes.c_int(nD), ctypes.c_int(nt),
                             arr(self.bidx, np.int32), ctypes.c_int(len(self.bidx)), arr(k.tau, np.float64),
                             ctypes.c_int(len(tc)), arr(C, np.float64), arr(wk, np.int32), arr(jr, np.int32),
                             arr(fr, np.float64), arr(D[0][0], np.int32), arr(D[0][1], np.float64),
                             arr(D[1][0], np.int32), arr(D[1][1], np.float64), arr(tc, np.float64),
                             ctypes.c_int(len(self.offsets)), arr(self.offsets, np.int32), ctypes.c_int(k.dsub),
                             ctypes.c_int(len(self.edges)), arr(self.edges, np.float64), ctypes.c_double(self.dt),
                             out_i.ctypes.data_as(ctypes.c_void_p), out_e.ctypes.data_as(ctypes.c_void_p))
        if rc != 0:
            raise ValueError("compiled kernel evaluation: index outside the table")
        ion += out_i; ele += out_e

    def _vector(self, C, wk, jr, fr, D, tc, ion, ele):
        """Same basis, corners and time interpolation, vectorized over clouds x strips x corners (chunked)."""
        k, tau = self.k, self.k.tau
        bI = np.array(self.bidx)[:, None, None, None]   # gather without copying the kernel table
        O = len(self.offsets)
        dr = np.array([0, 0, 1, 1]); dd = np.array([0, 1, 0, 1])
        for si in range(2):
            jd, fd = D[si]
            for s0 in range(0, len(tc), self.chunk):
                sl = slice(s0, min(len(tc), s0 + self.chunk))
                n = sl.stop - sl.start
                wr = np.stack([1 - fr[sl], fr[sl]], 1)[:, dr]                 # (n, 4)
                wd = np.stack([1 - fd[sl], fd[sl]], 1)[:, dd]
                wcorner = wr * wd                                              # (n, 4)
                R = (jr[sl][:, None, None] + dr[None, None, :])                # (n, 1, 4)
                Dd = (jd[sl][:, None, None] - k.dsub * self.offsets[None, :, None] + dd[None, None, :])  # (n, O, 4)
                W = np.broadcast_to(wk[sl][:, None, None], (n, O, 4))
                R = np.broadcast_to(R, (n, O, 4))
                G = k.K[bI, si, W, R, Dd, :].astype(np.float64)                # (basis, n, O, 4, tau)
                Q = np.einsum("bnoxt,nb,nx->not", G, C[sl], wcorner, optimize=True)
                P = np.einsum("bnox,nb,nx->no", k.KP[bI, si, W, R, Dd], C[sl], wcorner, optimize=True)
                # np.interp(edges - t, tau, Q, left=0, right=Q[-1]) for every (cloud, strip)
                x = self.edges[None, :] - tc[sl][:, None]                      # (n, E)
                j = np.clip(np.searchsorted(tau, x, side="right") - 1, 0, len(tau) - 2)
                f = (x - tau[j]) / (tau[j + 1] - tau[j])
                Qt = Q.transpose(0, 2, 1)                                      # (n, tau, O)
                ci = np.arange(n)[:, None]
                q0, q1 = Qt[ci, j], Qt[ci, j + 1]                              # (n, E, O)
                qe = q0 + f[:, :, None] * (q1 - q0)
                qe = np.where((x < tau[0])[:, :, None], 0.0, qe)
                qe = np.where((x > tau[-1])[:, :, None], Qt[:, -1:, :], qe)
                ion[:, si, :] += np.diff(qe, axis=1).sum(0)                    # bin per cloud (as the loop), then sum
                b0 = np.floor(tc[sl] / self.dt).astype(int)
                ok = b0 < self.nb
                np.add.at(ele[:, si, :], b0[ok], P[ok])


class Electronics:
    """Interchangeable electronics response, applied to Response arrays (unshaped current, fC/ns per dt bin)."""

    def __init__(self, impulse_t_ns, impulse_h, label):
        h = np.asarray(impulse_h, float)
        self.t, self.h, self.label = np.asarray(impulse_t_ns, float), h / h.max(), label

    @classmethod
    def crrc_standin(cls, tp_ns, order=3, dt=1.0, tmax=None):
        t = np.arange(0.0, tmax or max(12 * tp_ns, 2000.0), dt)
        tau = tp_ns / order
        return cls(t, (t / tau) ** order * np.exp(-t / tau), f"CR-RC^{order} stand-in Tp={tp_ns:g} ns (PROVISIONAL, not VMM3a)")

    @classmethod
    def from_table(cls, path):
        t, h = np.loadtxt(path, unpack=True)
        return cls(t, h, f"impulse table {path}")

    def shape(self, current, dt_in):
        dt = self.t[1] - self.t[0]
        kk = int(round(dt_in / dt))
        up = np.repeat(current.reshape(len(current), -1), kk, axis=0) * dt
        n = up.shape[0] + len(self.h)
        nf = 1 << int(np.ceil(np.log2(n)))
        out = np.fft.irfft(np.fft.rfft(up, nf, axis=0) * np.fft.rfft(self.h, nf)[:, None], nf, axis=0)[:up.shape[0]]
        return out.reshape((up.shape[0],) + current.shape[1:])

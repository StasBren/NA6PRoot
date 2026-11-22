// V0: archive-only statistics on IMMUTABLE 2100 V K256 ion archives. No weighting solver, no transport.
// ROOT-free analysis core (I/O lives in v0_archive_stats.cc).
//
// CONVENTIONS (checked against CathodeIonArchive.hh and the cathode_ion_reference Garfield source):
//  * Patch frame: anode k at u = k*pitch, v = 0; cathodes at v = +-0.3 cm. theta = atan2(v, u - k*pitch).
//  * "wire" of an ion = ions.local_wire_index / global_wire_index = llround(u_patch/pitch):
//    the anode NEAREST TO THE ION BIRTH. It is NOT the wire that collected the avalanche electrons.
//    The collecting wire is only known from microscopic_wire_seed.csv (collected_e per primary and wire).
//  * Hemisphere A = (W_up - W_down)/(W_up + W_down), up: 0 < theta < pi, down: -pi < theta < 0; births exactly on
//    theta = 0 or pi are excluded, matching A_ion_hemi in microscopic_wire_seed.csv.
//  * Ion yield N_hat of a (primary, wire) fragment = sum of inverse_inclusion_weight (Horvitz-Thompson).
//    It is NEVER derived from collected electrons. A primary with ion_births_exact = 0 has no ion rows but is kept.
//  * Reservoir: K ions sampled uniformly without replacement PER HEED PRIMARY (not per wire). A wire fragment of a
//    primary can therefore have exact births > 0 and zero sampled ions; such fragments are reported, not dropped.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>
namespace V0 {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = 3.14159265358979323846;

struct Ion { int prim, lwire, gwire; double u, v, w, t, wt; };                    // cm, cm, cm, ns
struct Prim { int idx = 0, cluster = -1, success = 1; double u = 0, v = 0, w = 0, t = 0; long long Nexact = 0, nSaved = 0, collected = 0; };
struct SeedRow { int prim = 0, gwire = 0; long long collected = 0, births = 0; double A = kNaN, meanT = kNaN, sigT = kNaN, meanWmm = kNaN, sigWmm = kNaN; };
struct Muon {
  std::string path, id, channel; int slot = -1, station = 0, det = 0, base = 0, K = 0; double pitch = 0.3, thetaU = 0, thetaW = 0, wInMm = 0;
  std::vector<Prim> prims; std::vector<Ion> ions; bool haveSeed = false; std::vector<SeedRow> seed;
  // synthetic only: law-level truth per primary (expected hemisphere asymmetry of the generating distribution)
  std::map<int, double> lawA;
};

// descriptors: 0 A (hemisphere), 1 <sin>, 2 <cos>, 3 <dy> um, 4 <dx> um, 5 <sin2>, 6 <cos2>
constexpr int kND = 7;
static const char* kDesc[kND] = {"A", "sin", "cos", "dy_um", "dx_um", "sin2", "cos2"};

struct Frag {   // one (muon, primary, nearest-birth global wire)
  int muon = 0, prim = 0, gwire = 0, lwire = 0, nres = 0;
  double Nhat = 0, NprimExact = 0, neffRes = 0, fpc = 0;
  double x[kND], se[kND];
  double rbar = kNaN, tbar = kNaN, sigt = kNaN, tbarSE = kNaN, wbar = kNaN, sigw = kNaN;   // t relative to t_primary [ns], w [cm]
  // exact census from microscopic_wire_seed.csv (NaN / -1 if absent)
  long long exBirths = -1, exCollected = -1; double exA = kNaN, exTbar = kNaN, exSigt = kNaN, exSigwMm = kNaN;
  double startDuMm = 0, startVmm = 0;
  Frag() { for (int k = 0; k < kND; ++k) { x[k] = kNaN; se[k] = kNaN; } }
};
struct PrimRec {   // one (muon, primary), all wire fragments aggregated
  int muon = 0, prim = 0, nFrag = 0, nFragExact = 0; long long Nexact = 0, nSaved = 0, collected = 0;
  double Nhat = 0, x[kND], se[kND], tbar = kNaN, exA = kNaN, startDuMm = 0, startVmm = 0, tPrim = 0;
  int birthWireArgmax = -999, collWireArgmax = -999;
  PrimRec() { for (int k = 0; k < kND; ++k) { x[k] = kNaN; se[k] = kNaN; } }
};

inline void features(double dx, double dy, double f[kND]) {
  const double a = std::atan2(dy, dx);
  f[0] = (dy > 0) ? 1. : (dy < 0 ? -1. : kNaN);   // exact theta=0/pi excluded from A
  f[1] = std::sin(a); f[2] = std::cos(a); f[3] = dy * 1e4; f[4] = dx * 1e4; f[5] = std::sin(2 * a); f[6] = std::cos(2 * a);
}
inline double wvarUnbiased(double Sw, double Sw2, double SS) { const double d = Sw - Sw2 / Sw; return d > 1e-12 ? SS / d : 0.; }

inline void buildFragments(const std::vector<Muon>& mus, std::vector<Frag>& F, std::vector<PrimRec>& P) {
  F.clear(); P.clear();
  for (size_t im = 0; im < mus.size(); ++im) {
    const Muon& m = mus[im];
    std::map<int, const Prim*> pr; for (auto& p : m.prims) pr[p.idx] = &p;
    std::map<std::pair<int, int>, std::vector<const Ion*>> G; for (auto& i : m.ions) G[{i.prim, i.gwire}].push_back(&i);
    std::map<std::pair<int, int>, const SeedRow*> S; if (m.haveSeed) for (auto& s : m.seed) S[{s.prim, s.gwire}] = &s;
    // every (primary, wire) known from ions OR from the exact census
    std::map<std::pair<int, int>, int> keys; for (auto& kv : G) keys[kv.first] = 1; for (auto& kv : S) if (kv.second->births > 0) keys[kv.first] = 1;
    std::map<int, std::vector<size_t>> fragsOfPrim;
    for (auto& kk : keys) {
      Frag g; g.muon = (int)im; g.prim = kk.first.first; g.gwire = kk.first.second; g.lwire = g.gwire - m.base;
      auto ip = pr.find(g.prim); const double tP = ip != pr.end() ? ip->second->t : 0.;
      if (ip != pr.end()) { g.NprimExact = (double)ip->second->Nexact; g.startDuMm = 10. * (ip->second->u - std::llround(ip->second->u / m.pitch) * m.pitch); g.startVmm = 10. * ip->second->v; }
      auto is = S.find(kk.first);
      if (is != S.end()) { g.exBirths = is->second->births; g.exCollected = is->second->collected; g.exA = is->second->A; g.exTbar = is->second->meanT - tP; g.exSigt = is->second->sigT; g.exSigwMm = is->second->sigWmm; }
      auto ig = G.find(kk.first);
      if (ig != G.end()) {
        const auto& v = ig->second; g.nres = (int)v.size(); const double uw = g.lwire * m.pitch;
        double Sw = 0, Sw2 = 0, S1[kND] = {0}, Wk[kND] = {0}, sr = 0, sw = 0, st = 0;
        for (auto* p : v) { double f[kND]; features(p->u - uw, p->v, f); Sw += p->wt; Sw2 += p->wt * p->wt;
          for (int k = 0; k < kND; ++k) if (!std::isnan(f[k])) { S1[k] += p->wt * f[k]; Wk[k] += p->wt; }
          sr += p->wt * std::hypot(p->u - uw, p->v); sw += p->wt * p->w; st += p->wt * (p->t - tP); }
        g.Nhat = Sw; g.neffRes = Sw * Sw / Sw2;
        const double Nw = g.exBirths > 0 ? (double)g.exBirths : g.Nhat;   // fragment population size (exact if known)
        g.fpc = Nw > 1 ? std::max(0., (Nw - g.nres) / (Nw - 1.)) : 0.;
        for (int k = 0; k < kND; ++k) g.x[k] = Wk[k] > 0 ? S1[k] / Wk[k] : kNaN;
        g.rbar = 1e4 * sr / Sw; g.wbar = sw / Sw; g.tbar = st / Sw;
        double SS[kND] = {0}, SSt = 0, SSw = 0;
        for (auto* p : v) { double f[kND]; features(p->u - uw, p->v, f); for (int k = 0; k < kND; ++k) if (!std::isnan(f[k])) SS[k] += p->wt * (f[k] - g.x[k]) * (f[k] - g.x[k]);
          SSt += p->wt * (p->t - tP - g.tbar) * (p->t - tP - g.tbar); SSw += p->wt * (p->w - g.wbar) * (p->w - g.wbar); }
        g.sigt = std::sqrt(wvarUnbiased(Sw, Sw2, SSt)); g.sigw = std::sqrt(wvarUnbiased(Sw, Sw2, SSw));
        g.tbarSE = std::sqrt(g.sigt * g.sigt / std::max(1., g.neffRes) * g.fpc);
        for (int k = 0; k < kND; ++k) g.se[k] = std::sqrt(wvarUnbiased(Sw, Sw2, SS[k]) / std::max(1., g.neffRes) * g.fpc);
        if (!std::isnan(g.x[0])) g.se[0] = std::sqrt(std::max(0., 1. - g.x[0] * g.x[0]) / std::max(1., g.neffRes) * g.fpc);
      }
      fragsOfPrim[g.prim].push_back(F.size()); F.push_back(g);
    }
    // primary records: ALL primaries of the primaries tree, including zero-ion and zero-collected ones
    std::map<int, std::pair<long long, int>> collBest, birthBest;
    if (m.haveSeed) for (auto& s : m.seed) { auto& c = collBest[s.prim]; if (s.collected > c.first) c = {s.collected, s.gwire}; auto& b = birthBest[s.prim]; if (s.births > b.first) b = {s.births, s.gwire}; }
    for (auto& p : m.prims) {
      PrimRec r; r.muon = (int)im; r.prim = p.idx; r.Nexact = p.Nexact; r.nSaved = p.nSaved; r.collected = p.collected; r.tPrim = p.t;
      r.startDuMm = 10. * (p.u - std::llround(p.u / m.pitch) * m.pitch); r.startVmm = 10. * p.v;
      if (collBest.count(p.idx)) r.collWireArgmax = collBest[p.idx].second;
      if (birthBest.count(p.idx)) r.birthWireArgmax = birthBest[p.idx].second;
      double S[kND] = {0}, W[kND] = {0}, V[kND] = {0}, sT = 0, wT = 0, sxA = 0, wxA = 0;
      for (size_t fi : fragsOfPrim[p.idx]) { const Frag& g = F[fi]; r.nFrag++; if (g.exBirths >= 0) r.nFragExact++; r.Nhat += g.Nhat;
        for (int k = 0; k < kND; ++k) if (!std::isnan(g.x[k])) { S[k] += g.Nhat * g.x[k]; W[k] += g.Nhat; V[k] += g.Nhat * g.Nhat * g.se[k] * g.se[k]; }
        if (!std::isnan(g.tbar)) { sT += g.Nhat * g.tbar; wT += g.Nhat; }
        if (g.exBirths > 0 && !std::isnan(g.exA)) { sxA += g.exBirths * g.exA; wxA += g.exBirths; } }
      for (int k = 0; k < kND; ++k) if (W[k] > 0) { r.x[k] = S[k] / W[k]; r.se[k] = std::sqrt(V[k]) / W[k]; }
      if (wT > 0) r.tbar = sT / wT;
      if (wxA > 0) r.exA = sxA / wxA;
      P.push_back(r);
    }
  }
}

// ---------------- small statistics helpers ----------------
inline std::vector<double> ranks(const std::vector<double>& v) {
  std::vector<size_t> id(v.size()); std::iota(id.begin(), id.end(), 0); std::sort(id.begin(), id.end(), [&](size_t a, size_t b) { return v[a] < v[b]; });
  std::vector<double> r(v.size()); for (size_t i = 0; i < id.size();) { size_t j = i; while (j + 1 < id.size() && v[id[j + 1]] == v[id[i]]) ++j; for (size_t k = i; k <= j; ++k) r[id[k]] = 0.5 * (i + j); i = j + 1; } return r;
}
inline double pearson(const std::vector<double>& a, const std::vector<double>& b) {
  const double n = (double)a.size(); if (n < 3) return kNaN; double ma = 0, mb = 0; for (size_t i = 0; i < a.size(); ++i) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n;
  double sab = 0, saa = 0, sbb = 0; for (size_t i = 0; i < a.size(); ++i) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
  return saa > 0 && sbb > 0 ? sab / std::sqrt(saa * sbb) : kNaN;
}
inline double spearman(const std::vector<double>& a, const std::vector<double>& b) { return pearson(ranks(a), ranks(b)); }
inline std::vector<double> quantileEdges(std::vector<double> v, int nb) { std::vector<double> e; if (v.empty()) return e; std::sort(v.begin(), v.end()); for (int i = 1; i < nb; ++i) e.push_back(v[std::min(v.size() - 1, (size_t)(v.size() * (double)i / nb))]); return e; }
inline int binOf(const std::vector<double>& e, double x) { return (int)(std::upper_bound(e.begin(), e.end(), x) - e.begin()); }
inline double quantile(std::vector<double> v, double q) { if (v.empty()) return kNaN; std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, (size_t)(q * (v.size() - 1) + 0.5))]; }

// ---------------- independence / effective-number test ----------------
// Units (a muon, or a (muon,wire)) contain members i with weight N_i, value x_i, sampling error se_i and a conditioning class.
// Model: x_i = mu(class, N_i) + e_i, Var(e_i) = s2(class, N_i), Corr(e_i, e_j) = rho within a unit.
//  X = sum N x / sum N ; E[X] = sum N mu_c / S ; Var(X) = (1-rho) p + rho q,  p = sum N^2 s2 / S^2, q = (sum N s)^2 / S^2.
//  Observable y = (X - E[X])^2 - sum N^2 se^2 / S^2 (sampling noise removed). Least squares y ~ a p + b q, rho = b/(a+b).
struct Member { double N, x, se; int cls; double zN; double mu = 0, s2 = 0; };   // N = weight in the unit average, zN = conditioning size
struct Unit { int key; std::vector<Member> m; };
struct IndepResult { int units = 0; double obs = 0, pred = 0, a = kNaN, b = kNaN, rho = kNaN, rhoLo = kNaN, rhoHi = kNaN; std::vector<double> obsByNeff, predByNeff; std::vector<int> nByNeff; };
// Conditional mean and variance of x given (class, log zN), zN = the PRIMARY's ion yield (a wire fragment inherits it): quantile bins in log N within each class (>= ~60 members
// per bin, at most 20 bins), noise-corrected bin variance, then piecewise-linear interpolation in log N between bin
// medians. Interpolation avoids the bias of a coarse bin whose variance changes strongly with N.
inline void assignMoments(std::vector<Unit>& U) {
  std::map<int, std::vector<Member*>> byCls; for (auto& u : U) for (auto& m : u.m) byCls[m.cls].push_back(&m);
  for (auto& kv : byCls) { auto& v = kv.second; std::sort(v.begin(), v.end(), [](const Member* a, const Member* b) { return a->zN < b->zN; });
    const int nb = std::max(1, std::min(20, (int)v.size() / 60)); std::vector<double> lx, mu, s2;
    for (int b = 0; b < nb; ++b) { const size_t i0 = v.size() * b / nb, i1 = v.size() * (b + 1) / nb; double s = 0, ss = 0, se2 = 0; const double n = double(i1 - i0); if (n < 2) continue;
      for (size_t i = i0; i < i1; ++i) { s += v[i]->x; ss += v[i]->x * v[i]->x; se2 += v[i]->se * v[i]->se; }
      const double m = s / n; lx.push_back(std::log(std::max(1e-9, v[(i0 + i1) / 2]->zN))); mu.push_back(m); s2.push_back(std::max(0., (ss - n * m * m) / (n - 1) - se2 / n)); }
    for (auto* m : v) { const double l = std::log(std::max(1e-9, m->zN)); size_t j = std::upper_bound(lx.begin(), lx.end(), l) - lx.begin();
      if (lx.empty()) continue; if (j == 0) { m->mu = mu[0]; m->s2 = s2[0]; } else if (j >= lx.size()) { m->mu = mu.back(); m->s2 = s2.back(); }
      else { const double f = (l - lx[j - 1]) / std::max(1e-12, lx[j] - lx[j - 1]); m->mu = mu[j - 1] + f * (mu[j] - mu[j - 1]); m->s2 = s2[j - 1] + f * (s2[j] - s2[j - 1]); } } }
}
inline IndepResult indepTest(std::vector<Unit> U, const std::vector<double>& neffEdges, unsigned long long seed = 88172645463325252ULL) {
  assignMoments(U);
  struct Row { double y, p, q, neff; int key; }; std::vector<Row> R;
  for (auto& u : U) { double S = 0, SX = 0, SE = 0, P = 0, Q = 0, noise = 0, S2 = 0; for (auto& m : u.m) { S += m.N; SX += m.N * m.x; SE += m.N * m.mu; P += m.N * m.N * m.s2; Q += m.N * std::sqrt(m.s2); noise += m.N * m.N * m.se * m.se; S2 += m.N * m.N; }
    if (S <= 0) continue; const double d = SX / S - SE / S; R.push_back({d * d - noise / (S * S), P / (S * S), (Q / S) * (Q / S), S * S / S2, u.key}); }
  IndepResult r; r.units = (int)R.size(); if (R.empty()) return r;
  auto fit = [&](const std::vector<size_t>& idx, double& a, double& b) { double spp = 0, spq = 0, sqq = 0, syp = 0, syq = 0; for (size_t i : idx) { spp += R[i].p * R[i].p; spq += R[i].p * R[i].q; sqq += R[i].q * R[i].q; syp += R[i].y * R[i].p; syq += R[i].y * R[i].q; }
    const double det = spp * sqq - spq * spq; if (std::fabs(det) < 1e-300) { a = b = kNaN; return; } a = (syp * sqq - syq * spq) / det; b = (spp * syq - spq * syp) / det; };
  std::vector<size_t> all(R.size()); std::iota(all.begin(), all.end(), 0); fit(all, r.a, r.b); r.rho = r.b / (r.a + r.b);
  for (auto& x : R) { r.obs += x.y; r.pred += x.p; } r.obs /= R.size(); r.pred /= R.size();
  std::map<int, std::vector<size_t>> byKey; for (size_t i = 0; i < R.size(); ++i) byKey[R[i].key].push_back(i); std::vector<int> keys; for (auto& kv : byKey) keys.push_back(kv.first);
  std::vector<double> rs; unsigned long long s = seed;
  for (int bs = 0; bs < 400; ++bs) { std::vector<size_t> idx; for (size_t k = 0; k < keys.size(); ++k) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; auto& v = byKey[keys[s % keys.size()]]; idx.insert(idx.end(), v.begin(), v.end()); }
    double a, b; fit(idx, a, b); if (std::isfinite(a + b) && (a + b) != 0) rs.push_back(b / (a + b)); }
  r.rhoLo = quantile(rs, 0.16); r.rhoHi = quantile(rs, 0.84);
  const int nb = (int)neffEdges.size() + 1; r.obsByNeff.assign(nb, 0); r.predByNeff.assign(nb, 0); r.nByNeff.assign(nb, 0);
  for (auto& x : R) { int b = binOf(neffEdges, x.neff); r.obsByNeff[b] += x.y; r.predByNeff[b] += x.p; r.nByNeff[b]++; }
  for (int b = 0; b < nb; ++b) if (r.nByNeff[b]) { r.obsByNeff[b] /= r.nByNeff[b]; r.predByNeff[b] /= r.nByNeff[b]; }
  return r;
}
}  // namespace V0

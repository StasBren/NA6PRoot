// Operation-cost benchmark for candidate production digitizer stages (no physics content).
// Table CONTENTS are synthetic; table SIZES, loop counts and memory access patterns follow the real chain:
//   * 52.8 Heed primaries per crossing on average (2100 V V0: 51,550 primaries / 977 muons, user-reported)
//   * Fast4 tables: LUT500 900 cells x 500 trials, J/psi donor 539,422 seeds (slides 14, 17)
//   * joint 2100 V start-cell donor: 51,550 records in 72 cells (|du| 0.25 mm x v 0.5 mm, V0 section [8])
//   * cathode kernel: 2 planes x 3 angular harmonics x 601 strip-offset bins x 400 time bins (2.5 ns), one shaper
// What is measured: ns per crossing of each stage, single thread, after warm-up, median of repetitions.
// Not measured: Heed itself (needs Garfield++), Geant4, I/O, real kernel accuracy.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {
struct Rng {   // xoshiro256**
  uint64_t s[4];
  explicit Rng(uint64_t seed) { for (auto& x : s) { seed += 0x9e3779b97f4a7c15ULL; uint64_t z = seed; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; x = z ^ (z >> 31); } }
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  uint64_t next() { const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
  double uni() { return (next() >> 11) * 0x1.0p-53; }
  uint32_t below(uint32_t n) { return (uint32_t)(((next() >> 32) * (uint64_t)n) >> 32); }
  int poisson(double mu) { if (mu < 30) { double L = std::exp(-mu), p = 1; int k = 0; do { ++k; p *= uni(); } while (p > L); return k - 1; }
    double u1 = uni(), u2 = uni(); return std::max(0, (int)std::lround(mu + std::sqrt(mu) * std::sqrt(-2 * std::log(u1 + 1e-300)) * std::cos(6.283185307179586 * u2))); }
};

struct Primary { float du, v, w, t; };                    // start offset from nearest anode [mm], v [mm], w [mm], t [ns]
struct LutRec { int8_t dk; float t, dw, phi; };            // Fast4 transport record (1800 V LUT500 layout)
struct DonorRec { float N, sinA, cosA, sigw; };             // Fast4 avalanche donor record
struct JointRec { int8_t dk; float N, sinA, cosA, delay, dw; };   // joint 2100 V transport+avalanche record (archive primary)
struct WireSum { int wire; double N, Ns, Nc, Nt, Ntt, Nw; };
struct WireTuple { int8_t dk; float N, Ns, Nc, tm, ts, w; };     // Strategy B: per-wire record

constexpr int kNcellJ = 72, kStrips = 13, kHarm = 3, kSbin = 601, kTbin = 400, kSamples = 16;
constexpr double kSstep = 0.02 /*mm*/, kTstep = 2.5 /*ns*/;

struct Tables {
  std::vector<LutRec> lut;            // 900*500
  std::vector<DonorRec> donor;        // 539422
  std::vector<JointRec> joint;        // 51550 grouped by cell
  std::vector<uint32_t> jStart, jCount;
  std::vector<WireTuple> wtab;        // Strategy B: 72 track cells x 2000 records x 3 wires
  std::vector<float> K;               // [plane][harm][sbin][tbin]
  std::vector<float> K2;              // second time derivative of harmonic 0, [plane][sbin][tbin]
};

Tables makeTables() {
  Tables T; Rng r(7);
  T.lut.resize(900 * 500); for (auto& x : T.lut) x = {(int8_t)(r.uni() < 0.04 ? (r.uni() < .5 ? -1 : 1) : 0), (float)(30 * r.uni()), (float)(0.1 * r.uni()), (float)(6.28 * r.uni())};
  T.donor.resize(539422); for (auto& x : T.donor) x = {(float)(3000 * r.uni()), (float)(2 * r.uni() - 1), (float)(2 * r.uni() - 1), (float)(0.05 * r.uni())};
  T.joint.resize(51550); T.jStart.resize(kNcellJ); T.jCount.resize(kNcellJ);
  for (int c = 0; c < kNcellJ; ++c) { T.jStart[c] = (uint32_t)(51550ull * c / kNcellJ); T.jCount[c] = (uint32_t)(51550ull * (c + 1) / kNcellJ) - T.jStart[c]; }
  for (auto& x : T.joint) x = {(int8_t)(r.uni() < 0.01 ? 1 : 0), (float)(3000 * r.uni()), (float)(2 * r.uni() - 1), (float)(2 * r.uni() - 1), (float)(40 * r.uni()), (float)(0.1 * r.uni())};
  T.wtab.resize(72 * 2000 * 3); for (auto& x : T.wtab) x = {(int8_t)(r.below(3) - 1), (float)(150000 * r.uni()), (float)r.uni(), (float)r.uni(), (float)(25 + 10 * r.uni()), (float)(10 * r.uni()), (float)r.uni()};
  T.K.resize((size_t)2 * kHarm * kSbin * kTbin); for (auto& x : T.K) x = (float)(r.uni() - 0.3);
  T.K2.resize((size_t)2 * kSbin * kTbin); for (auto& x : T.K2) x = (float)(1e-3 * (r.uni() - 0.5));
  return T;
}

inline int jcell(float du, float v) { const int a = std::min(5, (int)(std::fabs(du) / 0.25f)); const int b = std::max(0, std::min(11, (int)std::floor((v + 3.f) / 0.5f))); return a * 12 + b; }

// ---- stage: Heed surrogate (clusters along a 6 mm chord, cluster size from a small table) ----
int heedSurrogate(Rng& r, double uEntry, double slope, std::vector<Primary>& P) {
  static const float sizeCdf[8] = {.65f, .80f, .87f, .91f, .935f, .95f, .96f, 1.f};
  P.clear(); const int ncl = r.poisson(52.8 / 2.305);   // mean cluster size of the table below = 2.305 electrons
  for (int c = 0; c < ncl; ++c) { const float v = (float)(-3 + 6 * r.uni()); const float u = (float)(uEntry + slope * v); const float t = (float)(0.02 * r.uni());
    const float x = (float)r.uni(); int sz = 1; while (sz < 8 && x > sizeCdf[sz - 1]) ++sz; if (sz == 8) sz = 8 + (int)r.below(20);
    const float du = u - 3.f * std::round(u / 3.f), w = (float)(0.1 * r.uni());
    for (int k = 0; k < sz; ++k) P.push_back({du, v, w, t}); }
  return (int)P.size();
}
// ---- Strategy A (Fast4 two-table, 1800 V layout): LUT500 draw + J/psi donor draw per primary ----
void fast4Draws(Rng& r, const Tables& T, const std::vector<Primary>& P, int baseWire, std::vector<WireSum>& W) {
  W.clear();
  for (const auto& p : P) { const int iu = std::min(29, (int)((p.du + 1.5f) / 0.1f)), iv = std::min(29, std::max(0, (int)((p.v + 3.f) / 0.2f)));
    const LutRec& L = T.lut[(size_t)(iu * 30 + iv) * 500 + r.below(500)]; const DonorRec& D = T.donor[r.below((uint32_t)T.donor.size())];
    const int wire = baseWire + L.dk; WireSum* s = nullptr; for (auto& x : W) if (x.wire == wire) { s = &x; break; }
    if (!s) { W.push_back({wire, 0, 0, 0, 0, 0, 0}); s = &W.back(); }
    const double t = p.t + L.t; s->N += D.N; s->Ns += D.N * D.sinA; s->Nc += D.N * D.cosA; s->Nt += D.N * t; s->Ntt += D.N * t * t; s->Nw += D.N * (p.w + L.dw); }
}
// ---- Strategy A (2100 V joint start-cell donor): one draw per primary ----
void jointDraws(Rng& r, const Tables& T, const std::vector<Primary>& P, int baseWire, std::vector<WireSum>& W) {
  W.clear();
  for (const auto& p : P) { const int c = jcell(p.du, p.v); const JointRec& J = T.joint[T.jStart[c] + r.below(T.jCount[c])];
    const int wire = baseWire + J.dk; WireSum* s = nullptr; for (auto& x : W) if (x.wire == wire) { s = &x; break; }
    if (!s) { W.push_back({wire, 0, 0, 0, 0, 0, 0}); s = &W.back(); }
    const double t = p.t + J.delay; s->N += J.N; s->Ns += J.N * J.sinA; s->Nc += J.N * J.cosA; s->Nt += J.N * t; s->Ntt += J.N * t * t; s->Nw += J.N * (p.w + J.dw); }
}
// ---- Strategy B: one conditional per-crossing draw of up to 3 correlated per-wire tuples ----
void directWire(Rng& r, const Tables& T, double uEntry, int baseWire, std::vector<WireSum>& W) {
  W.clear(); const int cell = (int)std::min(71., std::fabs(uEntry) / 1.5 * 72); const size_t rec = ((size_t)cell * 2000 + r.below(2000)) * 3;
  for (int k = 0; k < 3; ++k) { const WireTuple& x = T.wtab[rec + k]; if (x.N <= 0) continue; const double N = x.N;
    W.push_back({baseWire + x.dk, N, x.Ns * N, x.Nc * N, N * x.tm, N * (x.tm * x.tm + x.ts * x.ts), N * x.w}); }
}
// ---- cathode response: shaped-kernel evaluation into strip sample buffers ----
// buffer[plane][strip][sample]; samples on a fixed grid around the expected peak (tStart + j*dtS)
inline float kAt(const Tables& T, int plane, int h, int sb, double tt) {
  if (tt <= 0) return 0.f;
  const double x = tt / kTstep; int i = (int)x; if (i >= kTbin - 1) i = kTbin - 2; const float f = (float)(x - i);
  const float* row = &T.K[(((size_t)plane * kHarm + h) * kSbin + sb) * kTbin]; return row[i] + f * (row[i + 1] - row[i]);
}
inline float k2At(const Tables& T, int plane, int sb, double tt) {
  if (tt <= 0) return 0.f;
  int i = std::min(kTbin - 1, (int)(tt / kTstep)); return T.K2[((size_t)plane * kSbin + sb) * kTbin + i];
}
void evalWire(const Tables& T, double N, double Ns, double Nc, double tMean, double tVar, double wMm, int wireIdx, bool secondOrder,
              float buf[2][kStrips][kSamples], double tStart, double dtS) {
  for (int plane = 0; plane < 2; ++plane) {
    const double s0 = (plane == 0 ? wMm - 0.3 * wireIdx : wMm + 0.3 * wireIdx);   // stereo strip coordinate of the source [mm]
    const int centre = (int)std::lround(s0 / 1.7);
    for (int k = 0; k < kStrips; ++k) {
      const double ds = (centre + k - kStrips / 2) * 1.7 - s0; int sb = (int)std::lround(ds / kSstep) + kSbin / 2; sb = std::max(0, std::min(kSbin - 1, sb));
      for (int j = 0; j < kSamples; ++j) { const double tt = tStart + j * dtS - tMean;
        float y = (float)(N * kAt(T, plane, 0, sb, tt) + Ns * kAt(T, plane, 1, sb, tt) + Nc * kAt(T, plane, 2, sb, tt));
        if (secondOrder) y += (float)(0.5 * tVar * N * k2At(T, plane, sb, tt));
        buf[plane][k][j] += y; } } }
}
struct Hit { int plane, strip; float amp, tpk, tthr; };
int digitize(float buf[2][kStrips][kSamples], double tStart, double dtS, float thr, Hit* out) {
  int n = 0;
  for (int p = 0; p < 2; ++p) for (int k = 0; k < kStrips; ++k) { const float* y = buf[p][k]; int jm = 0; for (int j = 1; j < kSamples; ++j) if (y[j] > y[jm]) jm = j;
      if (y[jm] <= thr) continue; float a = y[jm], tp = (float)(tStart + jm * dtS);
      if (jm > 0 && jm < kSamples - 1) { const float d = y[jm - 1] - 2 * y[jm] + y[jm + 1]; if (d < 0) { const float o = 0.5f * (y[jm - 1] - y[jm + 1]) / d; a = y[jm] - 0.25f * (y[jm - 1] - y[jm + 1]) * o; tp += o * (float)dtS; } }
      int jt = 0; while (jt < jm && y[jt] <= thr) ++jt; float tt = (float)(tStart + jt * dtS); if (jt > 0) tt -= (float)dtS * (y[jt] - thr) / (y[jt] - y[jt - 1]);
      out[n++] = {p, k, a, tp, tt}; }
  return n;
}

double nowNs() { return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

struct Result { std::string name; double ns; };
}  // namespace

int main(int argc, char** argv) {
  const int nCross = argc > 1 ? std::atoi(argv[1]) : 20000, reps = argc > 2 ? std::atoi(argv[2]) : 7;
  const double tp = argc > 3 ? std::atof(argv[3]) : 100.;   // peaking time sets the sample grid only
  const Tables T = makeTables();
  std::printf("tables: LUT500 %zu, J/psi donor %zu, joint-2100 %zu, B-table %zu, kernel %.1f MB\n", T.lut.size(), T.donor.size(), T.joint.size(), T.wtab.size(), (T.K.size() + T.K2.size()) * 4 / 1e6);
  // fixed crossing inputs (entry u, slope) so every stage sees identical work
  std::vector<double> uE(nCross), sl(nCross); { Rng r(11); for (int i = 0; i < nCross; ++i) { uE[i] = -1.5 + 3 * r.uni(); sl[i] = 0.4 * (r.uni() - 0.5); } }
  std::vector<std::vector<Primary>> prim(nCross); { Rng r(13); for (int i = 0; i < nCross; ++i) heedSurrogate(r, uE[i], sl[i], prim[i]); }
  double meanP = 0; for (auto& p : prim) meanP += p.size(); meanP /= nCross;
  std::vector<WireSum> W; W.reserve(16); volatile double sink = 0;
  const double dtS = tp / 8., tStart = 30. + tp - 8 * dtS;   // 16 samples spanning the expected peak region
  auto bench = [&](const char* name, auto&& fn) { std::vector<double> v; for (int k = 0; k < reps + 1; ++k) { const double a = nowNs(); fn(); const double b = nowNs(); if (k) v.push_back((b - a) / nCross); }
    std::sort(v.begin(), v.end()); return Result{name, v[v.size() / 2]}; };
  std::vector<Result> R;
  R.push_back(bench("heed surrogate (Poisson clusters + size table)", [&] { Rng r(1); std::vector<Primary> P; P.reserve(256); double s = 0; for (int i = 0; i < nCross; ++i) s += heedSurrogate(r, uE[i], sl[i], P); sink = sink + s; }));
  R.push_back(bench("A: Fast4 draws (LUT500 + J/psi donor) + wire sums", [&] { Rng r(2); double s = 0; for (int i = 0; i < nCross; ++i) { fast4Draws(r, T, prim[i], 100, W); s += W[0].N; } sink = sink + s; }));
  R.push_back(bench("A: joint 2100 V start-cell donor + wire sums", [&] { Rng r(3); double s = 0; for (int i = 0; i < nCross; ++i) { jointDraws(r, T, prim[i], 100, W); s += W[0].N; } sink = sink + s; }));
  R.push_back(bench("B: direct per-wire correlated draw", [&] { Rng r(4); double s = 0; for (int i = 0; i < nCross; ++i) { directWire(r, T, uE[i], 100, W); s += W.empty() ? 0 : W[0].N; } sink = sink + s; }));
  // pre-generate per-wire sums for the kernel stages so only the kernel is timed
  std::vector<std::vector<WireSum>> wires(nCross); { Rng r(5); for (int i = 0; i < nCross; ++i) { jointDraws(r, T, prim[i], 100, W); wires[i] = W; } }
  double meanW = 0; for (auto& w : wires) meanW += w.size(); meanW /= nCross;
  float buf[2][kStrips][kSamples]; Hit hits[2 * kStrips];
  R.push_back(bench("kernel T4: per wire, mean time only", [&] { double s = 0; for (int i = 0; i < nCross; ++i) { std::memset(buf, 0, sizeof buf);
      for (auto& w : wires[i]) evalWire(T, w.N, w.Ns, w.Nc, w.Nt / w.N, 0, w.Nw / w.N, w.wire - 100, false, buf, tStart, dtS); s += buf[0][6][8]; } sink = sink + s; }));
  R.push_back(bench("kernel T3: per wire, mean + variance term", [&] { double s = 0; for (int i = 0; i < nCross; ++i) { std::memset(buf, 0, sizeof buf);
      for (auto& w : wires[i]) { const double m = w.Nt / w.N; evalWire(T, w.N, w.Ns, w.Nc, m, w.Ntt / w.N - m * m, w.Nw / w.N, w.wire - 100, true, buf, tStart, dtS); } s += buf[0][6][8]; } sink = sink + s; }));
  // T2: one kernel evaluation per primary (each with its own time); use joint draws to give each primary a record
  std::vector<std::vector<JointRec>> pj(nCross); { Rng r(6); for (int i = 0; i < nCross; ++i) for (auto& p : prim[i]) { const int c = jcell(p.du, p.v); pj[i].push_back(T.joint[T.jStart[c] + r.below(T.jCount[c])]); } }
  R.push_back(bench("kernel T2: per primary, own time", [&] { double s = 0; for (int i = 0; i < nCross; ++i) { std::memset(buf, 0, sizeof buf);
      for (size_t k = 0; k < prim[i].size(); ++k) { const auto& J = pj[i][k]; const auto& p = prim[i][k]; evalWire(T, J.N, J.N * J.sinA, J.N * J.cosA, p.t + J.delay, 0, p.w + J.dw, J.dk, false, buf, tStart, dtS); } s += buf[0][6][8]; } sink = sink + s; }));
  float proto[2][kStrips][kSamples]; for (int p = 0; p < 2; ++p) for (int k = 0; k < kStrips; ++k) for (int j = 0; j < kSamples; ++j) proto[p][k][j] = (float)std::exp(-0.5 * (j - 8 - 0.3 * k) * (j - 8 - 0.3 * k) / 9.) * (k + 1);
  R.push_back(bench("digitize: peak, parabolic refine, threshold time", [&] { double s = 0; for (int i = 0; i < nCross; ++i) { std::memcpy(buf, proto, sizeof buf); buf[i & 1][i % kStrips][i % kSamples] += 0.5f;
      s += digitize(buf, tStart, dtS, 2.f, hits); } sink = sink + s; }));
  std::printf("crossings %d, primaries/crossing %.2f, wires with ions/crossing %.2f, strips %d per plane, %d samples, reps %d (median)\n", nCross, meanP, meanW, kStrips, kSamples, reps);
  for (auto& r : R) std::printf("  %-55s %9.1f ns/crossing\n", r.name.c_str(), r.ns);
  std::printf("checksum %g\n", (double)sink);
  return 0;
}

// Compiled K2A cathode response on ABSOLUTE strips, from in-memory full-birth cloud moments (milestone 2 adapter).
//
//   auto k = K2AStrip::KernelTable::Open(KERNEL_DIR, EXPECTED_KERNEL_ID);   // once per process, identity checked
//   K2AStrip::Evaluator eval(*k, {"K2A", 5, 7000., 10.});                     // model, half strips, window, dt
//   K2AStrip::Result r = eval.Evaluate(event);                                // event: K2AInline::Event (milestone 1)
//   if (r.usable) ... r.U.strips, r.U.ion (strip-major, fC/ns), r.U.electron ...
//
// Physics and numerics are the validated offline path, unchanged:
//   * kernel evaluation: spatial-ablation/engine/kernel_eval.cc, compiled in unmodified (same basis sum, 4 (r, D)
//     corners, np.interp time interpolation, electron step in bin floor(t/dt)); kernel = the s_kernel.py directory
//     (K.npy, KPHI0.npy, kernel.json), identity = cathode_engine.KernelTable.kernel_id;
//   * per-cloud coefficients and interpolation indices: CathodeResponseEngine._prepare, operation for operation;
//   * recentering: the E2a/E2b convention of k2a_recentered_k2a.py. Per Heed primary p,
//       m_U(p) = floor(Dbar_U(p) / P + 1/2), Dbar_U(p) = sum_i (N_i / sum_j N_j) D_U,i over the primary's clouds
//       (numpy pairwise-sum order), likewise V; P = kernel pitch_mm. Clouds of primaries with the same (m_U, m_V) form
//       one group, evaluated with D' = D - P m and base + m, groups in ascending (m_U, m_V) order; each group's
//       response (fC per bin / dt) is added to its absolute strip labels. Exact in the frozen model (w-invariant
//       field, identical strips, linear superposition); whole primaries keep their census.
//   * an explicit shift table (e.g. the E2b groups.json) can be imposed for reproduction tests.
// Output: per plane the sorted absolute strip indices that any group window touches, with ion current and the
// APPROXIMATE avalanche-electron step on the 10 ns grid (bins x strips), and full provenance: every cloud's
// (primary, wire) -> (m_U, m_V) group -> strip windows.
//
// Statuses (per cloud and OR-ed per event). Nothing is clamped and no cloud is dropped: if ANY cloud carries a
// blocking status, the event is not usable and NO waveform is produced (the statuses and the blocked charge are).
//   INPUT_BLOCKED          the moment event itself is blocked (K2AInline blocking status)       BLOCKING
//   WIRE_NOT_IN_KERNEL     cloud wire not tabulated                                            BLOCKING
//   R_EFF_BELOW_GRID / R_EFF_ABOVE_GRID   r_eff outside the radial grid (no extrapolation)     BLOCKING
//   D_U_OUTSIDE_TABLE / D_V_OUTSIDE_TABLE strip-relative D', all +-half_strips strips, outside the D table  BLOCKING
//   NONFINITE_MOMENT       non-finite moment                                                   BLOCKING
//   NONPOSITIVE_CHARGE     N_fC <= 0                                                           BLOCKING
//   T_NEGATIVE             cloud time < 0 (before the evaluator clock)                         BLOCKING
//   ENGINE_INDEX_ERROR     kernel_eval refused an index (cannot happen after the checks)       BLOCKING
//   SHIFT_UNDEFINED        the primary's Dbar is not finite (zero/non-finite charge or D)      BLOCKING
//   T_AFTER_WINDOW         cloud time >= window: its response lies after the readout window     informational
//   SHIFT_NEAR_TIE         Dbar/P + 1/2 within 1e-9 of an integer (either m is exact physics)   informational
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "K2AKernelIO.hh"
#include "MomentAccumulator.hh"

extern "C" int kernel_eval(const float* K, const double* KP, int nbasisAll, int nwire, int nr, int nD, int ntau,
                           const int* bidx, int nb, const double* tau, int ncl, const double* C, const int* wk,
                           const int* jr, const double* fr, const int* jdU, const double* fdU, const int* jdV,
                           const double* fdV, const double* tcl, int nO, const int* offs, int dsub, int nE,
                           const double* edges, double dt, double* ion, double* ele);

namespace K2AStrip {
constexpr const char* kVersion = "k2a-strip-evaluator v1";

enum Status : std::uint32_t {
  OK = 0,
  INPUT_BLOCKED = 1u << 0,
  WIRE_NOT_IN_KERNEL = 1u << 1,
  R_EFF_BELOW_GRID = 1u << 2,
  R_EFF_ABOVE_GRID = 1u << 3,
  D_U_OUTSIDE_TABLE = 1u << 4,
  D_V_OUTSIDE_TABLE = 1u << 5,
  NONFINITE_MOMENT = 1u << 6,
  NONPOSITIVE_CHARGE = 1u << 7,
  T_NEGATIVE = 1u << 8,
  ENGINE_INDEX_ERROR = 1u << 9,
  T_AFTER_WINDOW = 1u << 10,
  SHIFT_NEAR_TIE = 1u << 11,
  SHIFT_UNDEFINED = 1u << 12,
};
constexpr std::uint32_t kInformational = T_AFTER_WINDOW | SHIFT_NEAR_TIE;
constexpr std::uint32_t kBlocking = ~kInformational & ((1u << 13) - 1);

inline std::string StatusNames(std::uint32_t s) {
  if (s == 0) return "OK";
  static const std::pair<std::uint32_t, const char*> n[] = {
      {INPUT_BLOCKED, "INPUT_BLOCKED"}, {WIRE_NOT_IN_KERNEL, "WIRE_NOT_IN_KERNEL"}, {R_EFF_BELOW_GRID, "R_EFF_BELOW_GRID"},
      {R_EFF_ABOVE_GRID, "R_EFF_ABOVE_GRID"}, {D_U_OUTSIDE_TABLE, "D_U_OUTSIDE_TABLE"}, {D_V_OUTSIDE_TABLE, "D_V_OUTSIDE_TABLE"},
      {NONFINITE_MOMENT, "NONFINITE_MOMENT"}, {NONPOSITIVE_CHARGE, "NONPOSITIVE_CHARGE"}, {T_NEGATIVE, "T_NEGATIVE"},
      {ENGINE_INDEX_ERROR, "ENGINE_INDEX_ERROR"}, {T_AFTER_WINDOW, "T_AFTER_WINDOW"}, {SHIFT_NEAR_TIE, "SHIFT_NEAR_TIE"},
      {SHIFT_UNDEFINED, "SHIFT_UNDEFINED"}};
  std::string o;
  for (const auto& [b, name] : n) if (s & b) o += (o.empty() ? "" : "|") + std::string(name);
  return o;
}

// numpy's pairwise summation (umath pairwise_sum, PW_BLOCKSIZE 128, 8-way unrolled), so that Dbar and therefore the
// integer shifts are reproduced exactly from the same moments.
inline double NumpySum(const double* a, std::size_t n) {
  if (n < 8) { double r = 0.; for (std::size_t i = 0; i < n; ++i) r += a[i]; return r; }
  if (n <= 128) {
    double r[8];
    for (int j = 0; j < 8; ++j) r[j] = a[j];
    std::size_t i = 8;
    for (; i < n - (n % 8); i += 8) for (int j = 0; j < 8; ++j) r[j] += a[i + j];
    double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; ++i) res += a[i];
    return res;
  }
  std::size_t n2 = n / 2; n2 -= n2 % 8;
  return NumpySum(a, n2) + NumpySum(a + n2, n - n2);
}

// ------------------------------------------------------------------------------------------------------ kernel table
class KernelTable {
 public:
  static std::unique_ptr<KernelTable> Open(const std::string& dir, const std::string& expected_id, bool verify = true,
                                           bool preload = false) {
    using clk = std::chrono::steady_clock;
    auto k = std::unique_ptr<KernelTable>(new KernelTable);
    const auto t0 = clk::now();
    k->dir = dir;
    if (verify || !expected_id.empty()) {
      std::string cat;
      for (const char* f : {"K.npy", "KPHI0.npy", "kernel.json"}) cat += K2AKernelIO::Sha256::OfFile(dir + "/" + f);
      k->kernel_id = K2AKernelIO::Sha256::OfString(cat);
      if (!expected_id.empty() && expected_id != k->kernel_id)
        throw std::runtime_error("kernel id mismatch: expected " + expected_id + ", found " + k->kernel_id);
    }
    const auto t1 = clk::now();
    k->meta = K2AKernelIO::Json::ParseFile(dir + "/kernel.json");
    const auto& ax = k->meta["green_axes"];
    for (const auto& b : k->meta["basis"].arr) k->basis.push_back(b.Str());
    for (double w : ax["wires"].Nums()) k->wires.push_back(int(w));
    k->rn = ax["r_um"].Nums(); k->D = ax["D_mm"].Nums(); k->tau = ax["tau_ns"].Nums();
    k->M = int(ax["theta_deg"].arr.size()); k->dD = ax["dD_mm"].Num(); k->dsub = int(ax["d_sub"].Num());
    k->pitch_mm = ax["pitch_mm"].Num(); k->stereo_tan = ax["stereo_tan"].Num();
    if (ax.Has("grid_vertices_sha256")) k->transport_id = ax["grid_vertices_sha256"].Str();
    if (ax.Has("lut_sha256")) k->weighting_id = ax["lut_sha256"].Str();
    k->K = K2AKernelIO::NpyArray::Open(dir + "/K.npy", preload);
    k->KP = K2AKernelIO::NpyArray::Open(dir + "/KPHI0.npy", true);
    const auto& s = k->K.shape;
    if (s.size() != 6 || k->K.descr != "<f4" || k->KP.descr != "<f8" || k->KP.shape.size() != 5)
      throw std::runtime_error("kernel arrays: unexpected rank or dtype");
    if (s[0] != (long long)k->basis.size() || s[1] != 2 || s[2] != (long long)k->wires.size() || s[3] != (long long)k->rn.size() ||
        s[4] != (long long)k->D.size() || s[5] != (long long)k->tau.size())
      throw std::runtime_error("kernel K.npy shape does not match kernel.json axes");
    for (int i = 0; i < 5; ++i) if (k->KP.shape[i] != s[i]) throw std::runtime_error("KPHI0.npy shape does not match K.npy");
    for (std::size_t i = 1; i < k->rn.size(); ++i) if (!(k->rn[i] > k->rn[i - 1])) throw std::runtime_error("r grid not increasing");
    for (std::size_t i = 1; i < k->tau.size(); ++i) if (!(k->tau[i] > k->tau[i - 1])) throw std::runtime_error("tau grid not increasing");
    const auto t2 = clk::now();
    k->verify_s = std::chrono::duration<double>(t1 - t0).count();
    k->load_s = std::chrono::duration<double>(t2 - t1).count();
    return k;
  }
  int BasisIndex(const std::string& n) const {
    for (std::size_t i = 0; i < basis.size(); ++i) if (basis[i] == n) return int(i);
    throw std::runtime_error("kernel has no basis function " + n);
  }
  std::size_t ResidentBytes() const { return (K.Mapped() ? 0 : K.Bytes()) + KP.Bytes(); }

  std::string dir, kernel_id, transport_id, weighting_id;
  K2AKernelIO::Json meta;
  std::vector<std::string> basis;
  std::vector<int> wires;
  std::vector<double> rn, D, tau;
  int M = 0, dsub = 0;
  double dD = 0, pitch_mm = 0, stereo_tan = 0;
  K2AKernelIO::NpyArray K, KP;
  double verify_s = 0, load_s = 0;

 private:
  KernelTable() = default;
};

struct Config {
  std::string model = "K2A";
  int half_strips = 5;
  double window_ns = 7000.0;
  double dt_ns = 10.0;
};

struct CloudRecord {
  int primary = 0, wire = 0, m_U = 0, m_V = 0;
  double N_fC = 0, t_ns = 0, r_eff_um = 0, D_U_mm = 0, D_V_mm = 0, D_U_local_mm = 0, D_V_local_mm = 0;
  std::uint32_t status = OK;
};
struct ShiftRecord {
  int primary = 0, m_U = 0, m_V = 0;
  double Dbar_U_mm = 0, Dbar_V_mm = 0, margin_U = 0, margin_V = 0;   // margin: distance of Dbar/P + 1/2 to an integer
  bool imposed = false;                                                // from an explicit shift table
  bool undefined = false;                                              // Dbar not finite: SHIFT_UNDEFINED
  int computed_m_U = 0, computed_m_V = 0;
};
struct GroupRecord {
  int m_U = 0, m_V = 0, clouds = 0;
  std::vector<int> primaries;
  double N_fC = 0;
  int first_U = 0, first_V = 0, n_strips = 0;                          // window: first + [0, n_strips)
};
struct Plane {
  std::vector<int> strips;            // ascending absolute strip indices
  std::vector<double> ion, electron;  // strip-major: [s * nbins + bin], fC/ns
};
struct Result {
  std::uint32_t status = OK;
  bool usable = false;
  int base_U = 0, base_V = 0, nbins = 0;
  double dt_ns = 0, window_ns = 0;
  std::string model, kernel_id;
  std::vector<CloudRecord> clouds;
  std::vector<ShiftRecord> shifts;
  std::vector<GroupRecord> groups;
  Plane U, V;
  double charge_fC = 0, blocked_charge_fC = 0;
  int shifts_differing_from_computed = 0;
  double evaluate_s = 0;
};

// ---------------------------------------------------------------------------------------------------------- evaluator
class Evaluator {
 public:
  Evaluator(const KernelTable& k, const Config& c) : k_(k), c_(c) {
    if (c.model == "K1") { nmax_ = 1; hemi_ = false; } else if (c.model == "K2") { nmax_ = 2; hemi_ = false; }
    else if (c.model == "K1A") { nmax_ = 1; hemi_ = true; } else if (c.model == "K2A") { nmax_ = 2; hemi_ = true; }
    else throw std::runtime_error("model " + c.model + " not in K1, K2, K1A, K2A");
    if (c.half_strips < 0 || !(c.dt_ns > 0) || !(c.window_ns > 0)) throw std::runtime_error("bad evaluator configuration");
    for (int o = -c.half_strips; o <= c.half_strips; ++o) offs_.push_back(o);
    nb_ = int(std::nearbyint(c.window_ns / c.dt_ns));                 // Python round(): half to even
    if (k.tau.back() < c.window_ns) throw std::runtime_error("kernel tau range shorter than the window");
    for (int i = 0; i <= nb_; ++i) edges_.push_back(i * c.dt_ns);
    // s_perp (cathode_engine.s_perp), theta = (i + 0.5) 2 pi / M
    const int M = k.M;
    std::vector<double> th(M), sm(M), sp(M);
    for (int i = 0; i < M; ++i) th[i] = (i + 0.5) * 2 * M_PI / M;
    double smsum = 0;
    for (int i = 0; i < M; ++i) { const double s = std::sin(th[i]); sm[i] = double((s > 0) - (s < 0)); }
    smsum = NumpySum(sm.data(), sm.size());
    for (int i = 0; i < M; ++i) sp[i] = sm[i] - smsum / M;
    for (int n = 1; n <= nmax_; ++n)
      for (int cs = 0; cs < 2; ++cs) {
        std::vector<double> b(M);
        for (int i = 0; i < M; ++i) b[i] = cs ? std::sin(n * th[i]) : std::cos(n * th[i]);
        double spb = 0, bb = 0;
        for (int i = 0; i < M; ++i) { spb += sp[i] * b[i]; bb += b[i] * b[i]; }
        for (int i = 0; i < M; ++i) sp[i] = sp[i] - b[i] * spb / bb;
      }
    spsm_ = 0; for (int i = 0; i < M; ++i) spsm_ += sp[i] * sm[i];
    a0base_ = smsum / M;
    for (int n = 1; n <= 2; ++n) {
      double cc = 0, ss = 0;
      for (int i = 0; i < M; ++i) { cc += std::cos(n * th[i]) * sm[i]; ss += std::sin(n * th[i]) * sm[i]; }
      cm_[n][0] = cc; cm_[n][1] = ss;
    }
    std::vector<std::string> names = {"K0"};
    for (int n = 1; n <= nmax_; ++n) { names.push_back("Kc" + std::to_string(n)); names.push_back("Ks" + std::to_string(n)); }
    if (hemi_) names.push_back("KA" + std::to_string(nmax_));
    for (const auto& n : names) bidx_.push_back(k.BasisIndex(n));
    lnr_.resize(k.rn.size());
    for (std::size_t i = 0; i < k.rn.size(); ++i) lnr_[i] = std::log(k.rn[i]);
  }

  // shifts_override: primary -> (m_U, m_V) (reproduction of a given grouping); nullptr = the production rule.
  Result Evaluate(const K2AInline::Event& ev, const std::map<int, std::pair<int, int>>* shifts_override = nullptr) const {
    const auto t0 = std::chrono::steady_clock::now();
    Result r;
    r.base_U = ev.base_U; r.base_V = ev.base_V; r.nbins = nb_; r.dt_ns = c_.dt_ns; r.window_ns = c_.window_ns;
    r.model = c_.model; r.kernel_id = k_.kernel_id;
    const double P = k_.pitch_mm;
    // --- per-primary integer shifts (clouds are ordered by (primary, wire), as in the moments) ---
    std::map<int, std::pair<int, int>> shift;
    for (std::size_t i = 0; i < ev.clouds.size();) {
      std::size_t j = i;
      while (j < ev.clouds.size() && ev.clouds[j].primary == ev.clouds[i].primary) ++j;
      std::vector<double> N, w, wu, wv;
      for (std::size_t q = i; q < j; ++q) N.push_back(ev.clouds[q].N_fC);
      const double Ns = NumpySum(N.data(), N.size());
      for (std::size_t q = i; q < j; ++q) { const double x = ev.clouds[q].N_fC / Ns; wu.push_back(x * ev.clouds[q].D_U_mm); wv.push_back(x * ev.clouds[q].D_V_mm); }
      ShiftRecord s; s.primary = ev.clouds[i].primary;
      s.Dbar_U_mm = NumpySum(wu.data(), wu.size()); s.Dbar_V_mm = NumpySum(wv.data(), wv.size());
      const double xu = s.Dbar_U_mm / P + .5, xv = s.Dbar_V_mm / P + .5;
      s.undefined = !(Ns > 0) || !(std::abs(xu) < 1e6) || !(std::abs(xv) < 1e6);   // also catches NaN; no UB in the cast
      s.computed_m_U = s.undefined ? 0 : int(std::floor(xu)); s.computed_m_V = s.undefined ? 0 : int(std::floor(xv));
      s.margin_U = std::abs(xu - std::nearbyint(xu)); s.margin_V = std::abs(xv - std::nearbyint(xv));
      s.m_U = s.computed_m_U; s.m_V = s.computed_m_V;
      if (shifts_override) {
        auto it = shifts_override->find(s.primary);
        if (it == shifts_override->end()) throw std::runtime_error("shift table lacks primary " + std::to_string(s.primary));
        s.imposed = true; s.m_U = it->second.first; s.m_V = it->second.second;
        if (s.m_U != s.computed_m_U || s.m_V != s.computed_m_V) ++r.shifts_differing_from_computed;
      }
      shift[s.primary] = {s.m_U, s.m_V};
      r.shifts.push_back(s);
      i = j;
    }
    // --- per-cloud domain checks (no clamping) ---
    const bool inputOk = ev.UsableForCathodeSignal();
    const int nD = int(k_.D.size()), hs = c_.half_strips;
    for (std::size_t q = 0; q < ev.clouds.size(); ++q) {
      const auto& c = ev.clouds[q];
      CloudRecord cr; cr.primary = c.primary; cr.wire = c.wire; cr.N_fC = c.N_fC; cr.t_ns = c.t_ns; cr.r_eff_um = c.r_eff_um;
      cr.D_U_mm = c.D_U_mm; cr.D_V_mm = c.D_V_mm;
      cr.m_U = shift[c.primary].first; cr.m_V = shift[c.primary].second;
      cr.D_U_local_mm = c.D_U_mm - P * cr.m_U; cr.D_V_local_mm = c.D_V_mm - P * cr.m_V;
      std::uint32_t st = inputOk ? OK : INPUT_BLOCKED;
      for (double x : {c.N_fC, c.t_ns, c.r_eff_um, c.D_U_mm, c.D_V_mm, c.C1, c.S1, c.C2, c.S2, c.A})
        if (!std::isfinite(x)) st |= NONFINITE_MOMENT;
      if (!(c.N_fC > 0)) st |= NONPOSITIVE_CHARGE;
      if (std::find(k_.wires.begin(), k_.wires.end(), c.wire) == k_.wires.end()) st |= WIRE_NOT_IN_KERNEL;
      if (c.r_eff_um < k_.rn.front()) st |= R_EFF_BELOW_GRID;
      if (c.r_eff_um > k_.rn.back()) st |= R_EFF_ABOVE_GRID;
      for (int si = 0; si < 2; ++si) {
        const double x = ((si ? cr.D_V_local_mm : cr.D_U_local_mm) - k_.D[0]) / k_.dD;
        const double jd = std::floor(x);
        if (!std::isfinite(x) || jd - k_.dsub * hs < 0 || jd + 1 + k_.dsub * hs >= nD) st |= si ? D_V_OUTSIDE_TABLE : D_U_OUTSIDE_TABLE;
      }
      if (c.t_ns < 0) st |= T_NEGATIVE;
      if (c.t_ns >= c_.window_ns) st |= T_AFTER_WINDOW;
      for (const auto& s : r.shifts) {
        if (s.primary != c.primary) continue;
        if (s.undefined && !s.imposed) st |= SHIFT_UNDEFINED;
        else if (!s.imposed && (s.margin_U < 1e-9 || s.margin_V < 1e-9)) st |= SHIFT_NEAR_TIE;
      }
      cr.status = st;
      r.status |= st;
      r.charge_fC += c.N_fC;
      if (st & kBlocking) r.blocked_charge_fC += c.N_fC;
      r.clouds.push_back(cr);
    }
    // --- groups (ascending (m_U, m_V)); provenance recorded even for a blocked event ---
    std::map<std::pair<int, int>, std::vector<std::size_t>> grp;
    for (std::size_t q = 0; q < r.clouds.size(); ++q) grp[{r.clouds[q].m_U, r.clouds[q].m_V}].push_back(q);
    for (const auto& [m, idx] : grp) {
      GroupRecord g; g.m_U = m.first; g.m_V = m.second; g.clouds = int(idx.size());
      std::set<int> ps;
      for (auto q : idx) { ps.insert(r.clouds[q].primary); g.N_fC += r.clouds[q].N_fC; }
      g.primaries.assign(ps.begin(), ps.end());
      g.first_U = ev.base_U + m.first - hs; g.first_V = ev.base_V + m.second - hs; g.n_strips = 2 * hs + 1;
      r.groups.push_back(g);
    }
    r.usable = (r.status & kBlocking) == 0;
    if (!r.usable) { r.evaluate_s = Seconds(t0); return r; }
    // --- evaluation per group with the unchanged compiled kernel evaluator ---
    std::map<int, std::vector<double>> acc[2][2];   // [plane][ion/electron] absolute strip -> bins
    const int nO = int(offs_.size()), nbas = int(bidx_.size());
    const auto& Ks = k_.K.shape;
    std::vector<double> ion, ele;
    for (const auto& [m, idx] : grp) {
      const int n = int(idx.size());
      std::vector<double> C(std::size_t(n) * nbas), fr(n), fdU(n), fdV(n), tcl(n);
      std::vector<int> wk(n), jr(n), jdU(n), jdV(n);
      for (int a = 0; a < n; ++a) {
        const auto& c = ev.clouds[idx[a]];
        const auto& cr = r.clouds[idx[a]];
        // coefficients: CathodeResponseEngine._prepare
        std::vector<double> coef = {1.0};
        double a0s = a0base_ * 1.0;
        const double Cn[3][2] = {{0, 0}, {c.C1, c.S1}, {c.C2, c.S2}};
        for (int nn = 1; nn <= nmax_; ++nn) {
          coef.push_back(Cn[nn][0]); coef.push_back(Cn[nn][1]);
          a0s = a0s + (2.0 / k_.M) * (Cn[nn][0] * cm_[nn][0] + Cn[nn][1] * cm_[nn][1]);
        }
        if (hemi_) coef.push_back((c.A - a0s) / spsm_);
        for (int b = 0; b < nbas; ++b) C[std::size_t(a) * nbas + b] = coef[b] * c.N_fC;
        wk[a] = int(std::find(k_.wires.begin(), k_.wires.end(), c.wire) - k_.wires.begin());
        const double lr = std::log(c.r_eff_um);
        int j = int(std::lower_bound(lnr_.begin(), lnr_.end(), lr) - lnr_.begin()) - 1;   // searchsorted(left) - 1
        j = std::clamp(j, 0, int(lnr_.size()) - 2);                                       // np.clip (index only)
        jr[a] = j; fr[a] = (lr - lnr_[j]) / (lnr_[j + 1] - lnr_[j]);
        const double xu = (cr.D_U_local_mm - k_.D[0]) / k_.dD, xv = (cr.D_V_local_mm - k_.D[0]) / k_.dD;
        jdU[a] = int(std::floor(xu)); fdU[a] = xu - jdU[a];
        jdV[a] = int(std::floor(xv)); fdV[a] = xv - jdV[a];
        tcl[a] = c.t_ns;
      }
      ion.assign(std::size_t(nb_) * 2 * nO, 0.); ele.assign(ion.size(), 0.);
      const int rc = kernel_eval(k_.K.Data<float>(), k_.KP.Data<double>(), int(Ks[0]), int(Ks[2]), int(Ks[3]), int(Ks[4]), int(Ks[5]),
                                 bidx_.data(), nbas, k_.tau.data(), n, C.data(), wk.data(), jr.data(), fr.data(), jdU.data(),
                                 fdU.data(), jdV.data(), fdV.data(), tcl.data(), nO, offs_.data(), k_.dsub, nb_ + 1, edges_.data(),
                                 c_.dt_ns, ion.data(), ele.data());
      if (rc != 0) {
        for (auto q : idx) r.clouds[q].status |= ENGINE_INDEX_ERROR;
        r.status |= ENGINE_INDEX_ERROR; r.usable = false; r.blocked_charge_fC = r.charge_fC;
        r.U = Plane(); r.V = Plane(); r.evaluate_s = Seconds(t0); return r;
      }
      for (int si = 0; si < 2; ++si) {
        const int base = (si ? ev.base_V + m.second : ev.base_U + m.first);
        for (int oi = 0; oi < nO; ++oi) {
          for (int e2 = 0; e2 < 2; ++e2) {
            auto& v = acc[si][e2][base + offs_[oi]];
            if (v.empty()) v.assign(nb_, 0.);
            const auto& src = e2 ? ele : ion;
            for (int b = 0; b < nb_; ++b) v[b] = v[b] + src[(std::size_t(b) * 2 + si) * nO + oi] / c_.dt_ns;
          }
        }
      }
    }
    for (int si = 0; si < 2; ++si) {
      Plane& p = si ? r.V : r.U;
      for (const auto& [s, v] : acc[si][0]) {
        p.strips.push_back(s);
        p.ion.insert(p.ion.end(), v.begin(), v.end());
        const auto& e = acc[si][1].at(s);
        p.electron.insert(p.electron.end(), e.begin(), e.end());
      }
    }
    r.evaluate_s = Seconds(t0);
    return r;
  }
  int Bins() const { return nb_; }
  const std::vector<double>& Edges() const { return edges_; }

 private:
  static double Seconds(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  const KernelTable& k_;
  Config c_;
  int nmax_ = 2, nb_ = 0;
  bool hemi_ = true;
  std::vector<int> offs_, bidx_;
  std::vector<double> edges_, lnr_;
  double spsm_ = 0, a0base_ = 0, cm_[3][2] = {{0, 0}, {0, 0}, {0, 0}};
};

// --------------------------------------------------------------------------------------------- moments CSV (tests)
// Reads a BirthMoments / K2AInline / s_kernel.py moments CSV (by column name) into a K2AInline::Event, for tests and
// for reproduction of offline results. The production path passes the in-memory Event directly.
inline K2AInline::Event ReadMomentsCsv(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot read " + path);
  std::string line;
  std::getline(f, line);
  std::vector<std::string> h;
  { std::stringstream ss(line); std::string x; while (std::getline(ss, x, ',')) h.push_back(x); }
  auto col = [&](const std::string& n) {
    for (std::size_t i = 0; i < h.size(); ++i) if (h[i] == n) return int(i);
    throw std::runtime_error(path + ": no column " + n);
  };
  const int ip = col("primary"), iw = col("wire"), it = col("t_ns"), iN = col("N_fC"), ir = col("r_eff_um"), iu = col("D_U_mm"),
            iv = col("D_V_mm"), c1 = col("C1"), s1 = col("S1"), c2 = col("C2"), s2 = col("S2"), iA = col("A"), bu = col("base_U"),
            bv = col("base_V");
  K2AInline::Event e;
  bool first = true;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::vector<std::string> v;
    { std::stringstream ss(line); std::string x; while (std::getline(ss, x, ',')) v.push_back(x); }
    if (v.size() < h.size()) throw std::runtime_error(path + ": short row");
    auto d = [&](int i) { return std::strtod(v[i].c_str(), nullptr); };
    K2AInline::Cloud c;
    c.primary = int(std::lround(d(ip))); c.wire = int(std::lround(d(iw)));
    c.t_ns = d(it); c.N_fC = d(iN); c.r_eff_um = d(ir); c.D_U_mm = d(iu); c.D_V_mm = d(iv);
    c.C1 = d(c1); c.S1 = d(s1); c.C2 = d(c2); c.S2 = d(s2); c.A = d(iA);
    const int b1 = int(std::lround(d(bu))), b2 = int(std::lround(d(bv)));
    if (first) { e.base_U = b1; e.base_V = b2; first = false; }
    else if (b1 != e.base_U || b2 != e.base_V) throw std::runtime_error(path + ": base_U/base_V differ between rows");
    e.clouds.push_back(c);
  }
  for (std::size_t i = 1; i < e.clouds.size(); ++i)
    if (std::make_pair(e.clouds[i - 1].primary, e.clouds[i - 1].wire) >= std::make_pair(e.clouds[i].primary, e.clouds[i].wire))
      throw std::runtime_error(path + ": clouds not strictly ordered by (primary, wire)");
  return e;
}

// replay format (k2a_recentered_k2a.py / frozen replay): replay/waveforms.csv and prompt/prompt_waveforms.csv,
// union labels U then V ascending, time_center_ns = edges + dt/2; plus clouds.csv, groups.csv, shifts.csv.
inline void WriteResult(const std::string& dir, const Result& r, const std::vector<double>& edges) {
  auto open = [](const std::string& p) { std::ofstream f(p); if (!f) throw std::runtime_error("cannot write " + p); f << std::setprecision(17); return f; };
  {
    auto f = open(dir + "/clouds.csv");
    f << "primary,wire,m_U,m_V,N_fC,t_ns,r_eff_um,D_U_mm,D_V_mm,D_U_local_mm,D_V_local_mm,status,status_names\n";
    for (const auto& c : r.clouds)
      f << c.primary << ',' << c.wire << ',' << c.m_U << ',' << c.m_V << ',' << c.N_fC << ',' << c.t_ns << ',' << c.r_eff_um << ','
        << c.D_U_mm << ',' << c.D_V_mm << ',' << c.D_U_local_mm << ',' << c.D_V_local_mm << ',' << c.status << ',' << StatusNames(c.status) << '\n';
  }
  {
    auto f = open(dir + "/shifts.csv");
    f << "primary,m_U,m_V,computed_m_U,computed_m_V,imposed,Dbar_U_mm,Dbar_V_mm,margin_U,margin_V\n";
    for (const auto& s : r.shifts)
      f << s.primary << ',' << s.m_U << ',' << s.m_V << ',' << s.computed_m_U << ',' << s.computed_m_V << ',' << int(s.imposed) << ','
        << s.Dbar_U_mm << ',' << s.Dbar_V_mm << ',' << s.margin_U << ',' << s.margin_V << '\n';
  }
  {
    auto f = open(dir + "/groups.csv");
    f << "group,m_U,m_V,clouds,N_fC,first_U,first_V,n_strips,primaries\n";
    for (const auto& g : r.groups) {
      f << "g_" << g.m_U << '_' << g.m_V << ',' << g.m_U << ',' << g.m_V << ',' << g.clouds << ',' << g.N_fC << ',' << g.first_U << ','
        << g.first_V << ',' << g.n_strips << ',';
      for (std::size_t i = 0; i < g.primaries.size(); ++i) f << (i ? " " : "") << g.primaries[i];
      f << '\n';
    }
  }
  if (!r.usable) return;
  for (int e2 = 0; e2 < 2; ++e2) {
    const std::string sub = dir + (e2 ? "/prompt" : "/replay");
    std::filesystem::create_directories(sub);
    auto f = open(sub + (e2 ? "/prompt_waveforms.csv" : "/waveforms.csv"));
    f << "time_center_ns";
    for (int si = 0; si < 2; ++si) for (int s : (si ? r.V : r.U).strips) f << ',' << (si ? "V_" : "U_") << s;
    f << '\n';
    for (int b = 0; b < r.nbins; ++b) {
      f << edges[b] + r.dt_ns / 2;
      for (int si = 0; si < 2; ++si) {
        const Plane& p = si ? r.V : r.U;
        const auto& a = e2 ? p.electron : p.ion;
        for (std::size_t s = 0; s < p.strips.size(); ++s) f << ',' << a[s * r.nbins + b];
      }
      f << '\n';
    }
  }
}
}  // namespace K2AStrip

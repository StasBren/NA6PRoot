// T0 (SYNTHETIC, no ROOT, no Garfield): incremental K2AInline::MomentAccumulator vs the validated batch reference
// BirthMoments::Compute on identical birth streams, plus statuses, census and half-cell boundary rounding.
//
//   moment_accumulator_t0 WORKDIR          (WORKDIR must not exist; receives the writer outputs and t0_report.json)
//
// Agreement classes (reported separately):
//   BITWISE    primary, wire, rows, base_U/V, N_fC, r_eff_um, D_U_mm, D_V_mm, C1, S1, C2, S2, A, min_clearance_um,
//              t_raw_rms_ns: the same expressions in the same order -> identical bits are REQUIRED.
//   T1 TIME    inline t_ns = per-primary sum(w t)/sum(w) (bitwise vs an independent recomputation in the same order);
//              batch t_ns = mean over ions of t + (mean - t). |inline - batch| <= 4 (n + 2) u T, T = max |t| of the cloud
//              and its primary, u = 2^-53 (two roundings per ion, recursive summation of n terms and one division:
//              standard gamma_n bound, Higham 2002 eq. 4.4, taken with a factor-2 margin).
//   T1 SPREAD  inline t_rms_ns = 0 by definition; batch t_rms_ns is the rounding noise of sqrt(<t^2> - <t>^2) of
//              collapsed times: <= sqrt(8 (n + 2) u) T.
//   RAW TIME   t_raw_mean_ns and the per-primary t_raw_rms_ns vs long-double two-pass references: mean within
//              (4 n + 16) u mean|t|; RMS within min(sqrt(dV), dV / rms_ref) with dV = 2 (4 n + 16) u <t^2>.
//   CENSUS     integer births per primary/cloud exact; N_fC vs (sum of weights) e within (n + 1) u relative.
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>

#include "../../spatial-ablation/BirthMoments.hh"
#include "../MomentAccumulator.hh"

namespace {
constexpr double kU = 1.1102230246251565e-16;   // 2^-53
int gFail = 0, gChecks = 0;
std::map<std::string, long> gBitwise;            // field -> bitwise-equal count
std::map<std::string, long> gCount;
std::map<std::string, double> gWorst;            // field -> worst |diff| / tolerance

void Check(bool ok, const std::string& what) {
  ++gChecks;
  if (!ok) { ++gFail; std::cout << "FAIL: " << what << "\n"; }
}
bool Same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
void Bit(const std::string& f, double a, double b, const std::string& ctx) {
  ++gCount[f];
  if (Same(a, b)) ++gBitwise[f];
  else { ++gFail; std::cout << "FAIL bitwise " << f << " " << ctx << ": " << std::setprecision(17) << a << " vs " << b << "\n"; }
  ++gChecks;
}
void Tol(const std::string& f, double a, double b, double tol, const std::string& ctx) {
  ++gCount[f]; ++gChecks;
  if (Same(a, b)) ++gBitwise[f];
  const double d = std::abs(a - b);
  gWorst[f] = std::max(gWorst[f], tol > 0 ? d / tol : (d > 0 ? INFINITY : 0.));
  if (!(d <= tol)) { ++gFail; std::cout << "FAIL tol " << f << " " << ctx << ": " << std::setprecision(17) << a << " vs " << b
                                        << " |d|=" << d << " tol=" << tol << "\n"; }
}

struct Synth { std::vector<BirthMoments::Birth> births; std::map<int, long long> perPrimary; std::vector<int> primaries; };

// u_cm whose 10 u_cm is exactly target_mm (walks the neighbouring doubles of target_mm / 10)
double ExactUcm(double target_mm) {
  double x = target_mm / 10;
  for (int i = 0; i < 64 && 10 * x != target_mm; ++i) x = std::nextafter(x, 10 * x < target_mm ? INFINITY : -INFINITY);
  if (10 * x != target_mm) throw std::runtime_error("cannot place an exact half-cell tie");
  return x;
}

Synth Make(std::uint64_t seed, bool k256weights, int nPrim) {
  Synth s;
  std::mt19937_64 g(seed);
  std::uniform_real_distribution<double> U(0, 1);
  const long long sizes[] = {0, 1, 2, 3, 17, 255, 256, 257, 1000, 4096, 30000};
  int pid = 3;
  for (int ip = 0; ip < nPrim; ++ip) {
    pid += 1 + int(g() % 3);                           // original Heed indices: increasing, with gaps
    const long long n = sizes[g() % (sizeof sizes / sizeof sizes[0])];
    s.primaries.push_back(pid); s.perPrimary[pid] = n;
    const int k0 = int(g() % 7) - 3;
    const double w0 = -0.4 + 0.8 * U(g), t0 = 5 + 60 * U(g), wt = k256weights ? std::max(1., (n * (1 + 3 * U(g))) / 256.) : 1.;
    for (long long j = 0; j < n; ++j) {
      const int k = std::clamp(k0 + int(g() % 3) - 1, -4, 4);
      const double r_mm = 0.0161 * std::pow(1.4 / 0.0161, U(g));     // 16.1 um .. 1.4 mm, log-uniform
      const double ph = 2 * M_PI * U(g);
      double u_mm = 3 * k + r_mm * std::cos(ph), v_mm = r_mm * std::sin(ph);
      v_mm = std::clamp(v_mm, -2.99, 2.99);
      const double t = t0 + (U(g) < .5 ? 0.3 * U(g) : -2 * std::log(1 - U(g)));   // prompt core + exponential tail
      s.births.push_back({pid, u_mm / 10, v_mm / 10, w0 + 0.02 * (U(g) - .5), t, wt});
    }
  }
  // exact half-cell ties (inner and patch-edge boundaries), appended to the last two non-empty primaries
  int tieP = -1;
  for (auto it = s.primaries.rbegin(); it != s.primaries.rend(); ++it) if (s.perPrimary[*it] > 0) { tieP = *it; break; }
  for (double tmm : {1.5, -1.5, 4.5, -4.5, 7.5, -7.5, 10.5, -10.5}) {
    s.births.push_back({tieP, ExactUcm(tmm), 0.05, 0.01, 40.0, s.births.back().weight});
    ++s.perPrimary[tieP];
  }
  // reorder so each primary's births are contiguous and in creation order (the generator's callback order)
  std::stable_sort(s.births.begin(), s.births.end(), [](const auto& a, const auto& b) { return a.primary < b.primary; });
  return s;
}

K2AInline::Event Inline(const Synth& s, const K2AInline::MuonFrame& f) {
  K2AInline::MomentAccumulator acc(f);
  std::size_t i = 0;
  for (int p : s.primaries) {
    acc.BeginPrimary(p);
    std::uint64_t ord = 0;
    for (; i < s.births.size() && s.births[i].primary == p; ++i, ++ord) {
      const auto& b = s.births[i];
      if (b.weight == 1.0) acc.AddCallbackBirth(b.u_cm, b.v_cm, b.w_cm, b.t_ns, ord);
      else acc.AddBirth(b.u_cm, b.v_cm, b.w_cm, b.t_ns, b.weight);
    }
    acc.EndPrimary(std::uint64_t(s.perPrimary.at(p)));
  }
  return acc.Finish();
}

void Compare(const std::string& name, const Synth& s, const K2AInline::MuonFrame& f) {
  const auto ref = BirthMoments::Compute(s.births, {f.patch_shift_cm, f.u_in_mm, f.w_in_mm});
  const auto ev = Inline(s, f);
  Check(ev.UsableForCathodeSignal() && (ev.status & ~K2AInline::ZERO_ION_PRIMARY) == 0,
        name + ": event usable, no flag other than the informational ZERO_ION_PRIMARY (got " + K2AInline::StatusNames(ev.status) + ")");
  Check(ev.base_U == ref.base_U && ev.base_V == ref.base_V, name + ": base_U/base_V");
  Check(ev.clouds.size() == ref.clouds.size(), name + ": cloud count");
  Check(ev.births == std::uint64_t(ref.rows), name + ": total births vs batch rows");
  // independent per-primary references (same order sums for the T1 mean; long double two-pass for raw RMS)
  std::map<int, std::pair<double, double>> pm;
  std::map<int, std::vector<std::pair<double, double>>> pt;   // primary -> (t, w)
  for (const auto& b : s.births) { auto& a = pm[b.primary]; a.first += b.weight * b.t_ns; a.second += b.weight; pt[b.primary].push_back({b.t_ns, b.weight}); }
  for (const auto& p : ev.primaries) {
    Check(p.births_seen == std::uint64_t(s.perPrimary.at(p.primary)) && p.callback_births == p.births_seen,
          name + ": primary census " + std::to_string(p.primary));
    if (p.births_seen == 0) { Check(p.status == K2AInline::ZERO_ION_PRIMARY && p.n_clouds == 0, name + ": zero-ion primary flagged"); continue; }
    const auto& a = pm.at(p.primary);
    Bit("primary_T1_mean", p.t_mean_ns, a.first / a.second, name + " p" + std::to_string(p.primary));
    long double sw = 0, swt = 0, sv = 0, st2 = 0, tmax = 0;
    for (auto [t, w] : pt.at(p.primary)) { sw += w; swt += (long double)w * t; tmax = std::max<long double>(tmax, std::abs(t)); }
    const long double m = swt / sw;
    for (auto [t, w] : pt.at(p.primary)) { sv += w * (t - m) * (t - m); st2 += (long double)w * t * t; }
    const double rmsRef = double(std::sqrt(sv / sw)), n = double(p.births_seen);
    const double dV = 2 * (4 * n + 16) * kU * double(st2 / sw);
    Tol("primary_t_raw_rms_ns", p.t_raw_rms_ns, rmsRef, std::min(std::sqrt(dV), rmsRef > 0 ? dV / rmsRef : INFINITY), name);
  }
  // per-cloud references for raw time and the error bounds
  std::map<std::pair<int, int>, std::vector<const BirthMoments::Birth*>> members;
  for (const auto& b : s.births) members[{b.primary, int(std::nearbyint(10 * b.u_cm / 3))}].push_back(&b);
  for (std::size_t j = 0; j < std::min(ev.clouds.size(), ref.clouds.size()); ++j) {
    const auto& c = ev.clouds[j]; const auto& r = ref.clouds[j];
    const std::string ctx = name + " (" + std::to_string(r.primary) + "," + std::to_string(r.wire) + ")";
    Check(c.primary == r.primary && c.wire == r.wire, ctx + ": cloud identity/order");
    Check(c.births == std::uint64_t(r.rows), ctx + ": cloud births");
    Check(c.status == K2AInline::OK, ctx + ": cloud status");
    Bit("N_fC", c.N_fC, r.N_fC, ctx); Bit("r_eff_um", c.r_eff_um, r.r_eff_um, ctx);
    Bit("D_U_mm", c.D_U_mm, r.D_U_mm, ctx); Bit("D_V_mm", c.D_V_mm, r.D_V_mm, ctx);
    Bit("C1", c.C1, r.C1, ctx); Bit("S1", c.S1, r.S1, ctx); Bit("C2", c.C2, r.C2, ctx); Bit("S2", c.S2, r.S2, ctx);
    Bit("A", c.A, r.A, ctx); Bit("min_clearance_um", c.min_clearance_um, r.min_clearance_um, ctx);
    Bit("t_raw_rms_ns", c.t_raw_rms_ns, r.t_raw_rms_ns, ctx);
    const auto& mem = members.at({r.primary, r.wire});
    const double n = double(mem.size());
    Check(mem.size() == c.births, ctx + ": independent cloud membership count");
    double T = std::abs(pm.at(r.primary).first / pm.at(r.primary).second), sumW = 0;
    long double sw = 0, swt = 0, sat = 0;
    for (const auto* b : mem) { T = std::max(T, std::abs(b->t_ns)); sumW += b->weight; sw += b->weight; swt += (long double)b->weight * b->t_ns; sat += (long double)b->weight * std::abs(b->t_ns); }
    Bit("cloud_t_ns_equals_primary_T1", c.t_ns, pm.at(r.primary).first / pm.at(r.primary).second, ctx);
    Tol("t_ns_vs_batch", c.t_ns, r.t_ns, 4 * (n + 2) * kU * T, ctx);
    Check(Same(c.t_rms_ns, 0.), ctx + ": collapsed t_rms_ns is exactly 0");
    Tol("batch_t_rms_ns_rounding", r.t_rms_ns, 0., std::sqrt(8 * (n + 2) * kU) * T, ctx);
    Tol("t_raw_mean_ns", c.t_raw_mean_ns, double(swt / sw), (4 * n + 16) * kU * double(sat / sw), ctx);
    Tol("N_fC_vs_census", c.N_fC, sumW * BirthMoments::kEfC, (n + 1) * kU * sumW * BirthMoments::kEfC, ctx);
  }
}

// half-cell boundary rounding: K2A rule (rint, half to even) in BOTH implementations; the generator census rule differs
void Ties() {
  struct T { double mm; int k2a; } tab[] = {{1.5, 0}, {-1.5, 0}, {4.5, 2}, {-4.5, -2}, {7.5, 2}, {-7.5, -2}, {10.5, 4}, {-10.5, -4}};
  for (auto [mm, k] : tab) {
    const double u = ExactUcm(mm);
    Check(10 * u / 3 == std::trunc(10 * u / 3) + std::copysign(.5, mm), "tie construction is an exact half-cell boundary");
    K2AInline::MomentAccumulator a({0, 0, 0});
    a.BeginPrimary(0); a.AddBirth(u, 0.05, 0, 1); a.EndPrimary(1);
    const auto e = a.Finish();
    const auto r = BirthMoments::Compute({{0, u, 0.05, 0, 1, 1}}, {0, 0, 0});
    Check(e.clouds.size() == 1 && e.clouds[0].wire == k, "inline tie u=" + std::to_string(mm) + " mm -> wire " + std::to_string(k));
    Check(r.clouds.size() == 1 && r.clouds[0].wire == k, "batch tie u=" + std::to_string(mm) + " mm -> wire " + std::to_string(k));
    const long long gen = std::llround(u / 0.3);   // generator census / archive local_wire_index rule
    std::cout << "  tie u_mm=" << std::setw(5) << mm << ": K2A wire " << std::setw(2) << k << " (inline " << e.clouds[0].wire
              << ", batch " << r.clouds[0].wire << "); generator llround(u/0.3) = " << gen << (gen != k ? "  <- differs" : "") << "\n";
  }
}

bool Throws(const std::vector<BirthMoments::Birth>& b) {
  try { BirthMoments::Compute(b, {0, 0, 0}); } catch (const std::exception&) { return true; }
  return false;
}

void Statuses(const std::filesystem::path& W) {
  using namespace K2AInline;
  auto one = [](double u, double v, double t, double w, std::uint64_t census) {
    MomentAccumulator a({0, 0, 0}); a.BeginPrimary(0); a.AddBirth(0.3, 0.05, 0, 1); a.AddBirth(u, v, 0, t, w); a.EndPrimary(census); return a.Finish();
  };
  // inside anode: r = 10 um (batch throws)
  auto e = one(0.3 + 0.0006, 0.0008, 1, 1, 2);
  Check((e.status & INSIDE_ANODE) && !e.UsableForCathodeSignal(), "INSIDE_ANODE flagged and blocking");
  Check(Throws({{0, 0.3 + 0.0006, 0.0008, 0, 1, 1}}), "batch refuses the inside-anode birth");
  Check(!WriteOutputs((W / "blocked_inside_anode").string(), e) && !std::filesystem::exists(W / "blocked_inside_anode.moments.csv") &&
        std::filesystem::exists(W / "blocked_inside_anode.clouds.csv"), "blocked event: diagnostics written, NO moments file");
  for (auto [lab, t, w] : {std::tuple{"NaN time", NAN, 1.}, {"inf time", INFINITY, 1.}, {"zero weight", 1., 0.}, {"negative weight", 1., -1.}}) {
    auto x = one(0.31, 0.05, t, w, 2);
    Check((x.status & NONFINITE) && !x.UsableForCathodeSignal() && x.births == 2 && x.primaries[0].callback_births == 2 &&
          !(x.status & CENSUS_MISMATCH), std::string("NONFINITE flagged, blocking, still counted in the census: ") + lab);
    Check(Throws({{0, 0.31, 0.05, 0, t, w}}), std::string("batch refuses: ") + lab);
  }
  auto o = one(1.5, 0.05, 1, 1, 2);                      // u = 15 mm -> K2A wire 5
  Check((o.status & OUTSIDE_PATCH) && !o.UsableForCathodeSignal(), "OUTSIDE_PATCH flagged and blocking (batch has no patch check)");
  auto c = one(0.31, 0.05, 1, 1, 3);
  Check((c.status & CENSUS_MISMATCH) && !c.UsableForCathodeSignal(), "CENSUS_MISMATCH flagged and blocking");
  { MomentAccumulator a({0, 0, 0}); a.AddBirth(0.3, .05, 0, 1); a.BeginPrimary(0); a.AddBirth(0.31, .05, 0, 1); a.EndPrimary(1);
    auto x = a.Finish(); Check((x.status & SEQUENCE_ERROR) && x.orphan_births == 1 && !x.UsableForCathodeSignal(), "birth before BeginPrimary -> SEQUENCE_ERROR"); }
  { MomentAccumulator a({0, 0, 0}); a.BeginPrimary(2); a.EndPrimary(0); a.BeginPrimary(2); a.EndPrimary(0);
    auto x = a.Finish(); Check((x.status & SEQUENCE_ERROR) && !x.UsableForCathodeSignal(), "repeated primary index -> SEQUENCE_ERROR"); }
  { MomentAccumulator a({0, 0, 0}); a.BeginPrimary(0); a.AddCallbackBirth(0.3, .05, 0, 1, 0); a.AddCallbackBirth(0.3, .06, 0, 1, 2); a.EndPrimary(2);
    auto x = a.Finish(); Check((x.status & SEQUENCE_ERROR) && !x.UsableForCathodeSignal(), "callback ordinal gap -> SEQUENCE_ERROR"); }
  { MomentAccumulator a({0, 0, 0}); a.BeginPrimary(0); a.AddBirth(0.3, .05, 0, 1);
    auto x = a.Finish(); Check((x.status & SEQUENCE_ERROR) && x.primaries.size() == 1 && !x.UsableForCathodeSignal(), "Finish with an open primary -> SEQUENCE_ERROR"); }
  { MomentAccumulator a({0, 0, 0}); a.BeginPrimary(0); a.EndPrimary(0); a.BeginPrimary(1); a.AddBirth(0.3, .05, 0, 1); a.EndPrimary(1);
    auto x = a.Finish();
    Check(x.status == ZERO_ION_PRIMARY && x.UsableForCathodeSignal() && x.clouds.size() == 1, "ZERO_ION_PRIMARY is informational only");
    Check(WriteOutputs((W / "usable_zero_ion").string(), x) && std::filesystem::exists(W / "usable_zero_ion.moments.csv"), "usable event writes the moments file"); }
}

void WriterRoundTrip(const std::filesystem::path& W, const Synth& s, const K2AInline::MuonFrame& f) {
  const auto ev = Inline(s, f);
  const auto ref = BirthMoments::Compute(s.births, {f.patch_shift_cm, f.u_in_mm, f.w_in_mm});
  K2AInline::WriteOutputs((W / "unit").string(), ev);
  std::ostringstream b; BirthMoments::Write(b, ref);
  std::ifstream in(W / "unit.moments.csv"); std::stringstream a; a << in.rdbuf();
  // same schema and the same text in every column except t_ns, t_rms_ns (T1 representation, see header)
  std::istringstream A(a.str()), B(b.str()); std::string la, lb; long rows = 0, diffCols = 0;
  while (std::getline(A, la) && std::getline(B, lb)) {
    std::vector<std::string> ca, cb; std::string x;
    for (std::istringstream s1(la); std::getline(s1, x, ',');) ca.push_back(x);
    for (std::istringstream s2(lb); std::getline(s2, x, ',');) cb.push_back(x);
    Check(ca.size() == cb.size(), "moments.csv column count");
    for (std::size_t j = 0; j < ca.size(); ++j) if (j != 2 && j != 3) { if (ca[j] != cb[j]) { ++diffCols; std::cout << "  col " << j << ": " << ca[j] << " vs " << cb[j] << "\n"; } }
    ++rows;
  }
  Check(rows == long(ref.clouds.size()) + 1 && diffCols == 0, "engine-facing moments.csv: BirthMoments schema; identical text outside t_ns/t_rms_ns");
}
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::runtime_error("usage: moment_accumulator_t0 WORKDIR");
    const std::filesystem::path W = argv[1];
    if (std::filesystem::exists(W)) throw std::runtime_error("WORKDIR exists");
    std::filesystem::create_directories(W);
    const double uIn = 123.4567, phase = -0.8765;
    const K2AInline::MuonFrame f{.1 * uIn - .1 * phase, uIn, -45.678};
    std::cout << "T0 tie rounding (u exactly on half-cell boundaries):\n";
    Ties();
    const auto unit = Make(20261125, false, 60), k256 = Make(77, true, 40);
    std::cout << "T0 unit-weight event: " << unit.births.size() << " births, " << unit.primaries.size() << " primaries\n";
    Compare("unit", unit, f);
    std::cout << "T0 K256-weight event: " << k256.births.size() << " births, " << k256.primaries.size() << " primaries\n";
    Compare("k256", k256, f);
    Compare("frame0", unit, {0, 0, 0});
    Statuses(W);
    WriterRoundTrip(W, unit, f);
    const auto ev = Inline(unit, f);
    std::ofstream j(W / "t0_report.json");
    j << std::setprecision(6) << "{\"test\": \"T0 synthetic\", \"checks\": " << gChecks << ", \"failures\": " << gFail
      << ", \"accumulator_state_bytes_peak_unit_event\": " << ev.state_bytes_peak << ", \"clouds_unit_event\": " << ev.clouds.size()
      << ", \"births_unit_event\": " << ev.births << ", \"bitwise\": {";
    bool first = true;
    for (auto& [k, n] : gCount) { j << (first ? "" : ", ") << "\"" << k << "\": [" << gBitwise[k] << ", " << n << "]"; first = false; }
    j << "}, \"worst_diff_over_tolerance\": {";
    first = true;
    for (auto& [k, v] : gWorst) { j << (first ? "" : ", ") << "\"" << k << "\": " << v; first = false; }
    j << "}}\n";
    std::cout << "agreement (bitwise-equal / compared; worst |diff|/tolerance where a tolerance applies):\n";
    for (auto& [k, n] : gCount)
      std::cout << "  " << std::left << std::setw(32) << k << std::right << std::setw(6) << gBitwise[k] << " / " << std::setw(6) << n
                << (gWorst.count(k) ? "   worst " + [&] { std::ostringstream o; o << std::setprecision(3) << gWorst[k]; return o.str(); }() : std::string("   (bitwise required)")) << "\n";
    std::cout << "accumulator state (unit event, " << ev.clouds.size() << " clouds, " << ev.births << " births): "
              << ev.state_bytes_peak << " bytes peak\n";
    std::cout << (gFail ? "FAIL" : "PASS") << " T0: " << gChecks - gFail << "/" << gChecks << " checks\n";
    return gFail ? 1 : 0;
  } catch (const std::exception& e) { std::cerr << "T0 ERROR: " << e.what() << "\n"; return 2; }
}

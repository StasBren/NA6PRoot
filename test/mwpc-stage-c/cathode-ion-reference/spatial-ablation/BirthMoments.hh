// Birth-record -> per-cloud K2A moments, without ion transport.
//
// Conventions reproduced from the FROZEN exporter p3_export_full_spectral_paths.cc (wire branch), whose first path
// vertex is exactly the archived birth (it sets points.front() = {b.u, b.v, b.w, b.t}):
//   u_mm = 10 u_patch_cm, v_mm = 10 v_absolute_cm, t = t_birth_ns, weight = inverse_inclusion_weight,
//   uCentralIndex = llround((w_in_cm - 0.1 u_in_cm) / 0.17), vCentralIndex = llround((w_in_cm + 0.1 u_in_cm) / 0.17),
//   wU_mm = 10 (w + w_in_cm - 0.1 patch_shift_cm - 0.17 uCentralIndex),
//   wV_mm = 10 (w + w_in_cm + 0.1 patch_shift_cm - 0.17 vCentralIndex),
// and from the T-ladder / S-ladder tools:
//   T1 (t_time_collapse_export.py --mode primary): every ion of a Heed primary at the primary's weight-mean birth time,
//   cloud = (Heed primary, wire k = rint(u_mm / 3)), q = weight x e (fC), and the moment definitions of
//   s_kernel.py moments (charge-weighted means; r_eff = exp<log r>; A = <sign(sin theta)>).
// Pure C++17, no ROOT, so the same code is unit-tested synthetically and used by the ROOT front end.
#pragma once
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace BirthMoments {
constexpr double kEfC = 1.602176634e-4;
constexpr double kWireRadiusUm = 15.0;

struct Birth { int primary; double u_cm, v_cm, w_cm, t_ns, weight; };
struct Frame { double patch_shift_cm, u_in_mm, w_in_mm; };
struct Cloud {
  int primary = 0, wire = 0;
  double t_ns = 0, t_rms_ns = 0, t_raw_rms_ns = 0, N_fC = 0, r_eff_um = 0, D_U_mm = 0, D_V_mm = 0;
  double C1 = 0, S1 = 0, C2 = 0, S2 = 0, A = 0, min_clearance_um = 0;
  long long rows = 0;
};
struct Result { std::vector<Cloud> clouds; int base_U = 0, base_V = 0; long long rows = 0; double weight_sum = 0; };

inline void Require(bool ok, const std::string& m) { if (!ok) throw std::runtime_error(m); }

inline Result Compute(const std::vector<Birth>& births, const Frame& f) {
  Result res;
  const double uIn = .1 * f.u_in_mm, wIn = .1 * f.w_in_mm;
  res.base_U = int(std::llround((wIn - .1 * uIn) / .17));
  res.base_V = int(std::llround((wIn + .1 * uIn) / .17));
  // T1: ion-weight mean birth time per primary, then the same per-ion shift as the collapse tool (t + (mean - t))
  std::map<int, std::pair<double, double>> pt;   // primary -> (sum w t, sum w)
  for (const auto& b : births) { auto& a = pt[b.primary]; a.first += b.weight * b.t_ns; a.second += b.weight; }
  struct Acc { double q = 0, t = 0, t2 = 0, traw = 0, traw2 = 0, lr = 0, du = 0, dv = 0, c1 = 0, s1 = 0, c2 = 0, s2 = 0,
               a = 0, rmin = std::numeric_limits<double>::infinity(); long long n = 0; };
  std::map<std::pair<int, int>, Acc> acc;        // ordered like numpy unique over (primary, wire)
  for (const auto& b : births) {
    Require(std::isfinite(b.u_cm) && std::isfinite(b.v_cm) && std::isfinite(b.w_cm) && std::isfinite(b.t_ns) &&
            b.weight > 0, "invalid birth record");
    const double u = 10 * b.u_cm, v = 10 * b.v_cm;
    const double wu = 10. * (b.w_cm + wIn - .1 * f.patch_shift_cm - .17 * res.base_U);
    const double wv = 10. * (b.w_cm + wIn + .1 * f.patch_shift_cm - .17 * res.base_V);
    const auto& p = pt[b.primary];
    const double tm = p.first / p.second;
    const double t = b.t_ns + (tm - b.t_ns);
    const int k = int(std::nearbyint(u / 3));      // numpy rint: round half to even
    const double du = u - 3 * k;
    const double r = 1000 * std::hypot(du, v);
    Require(r > kWireRadiusUm, "ion birth on or inside the anode");
    const double th = std::atan2(v, du);
    const double q = b.weight * kEfC;
    auto& a = acc[{b.primary, k}];
    a.q += q; a.t += q * t; a.t2 += q * t * t; a.traw += q * b.t_ns; a.traw2 += q * b.t_ns * b.t_ns;
    a.lr += q * std::log(r); a.du += q * wu; a.dv += q * wv;
    a.c1 += q * std::cos(th); a.s1 += q * std::sin(th); a.c2 += q * std::cos(2 * th); a.s2 += q * std::sin(2 * th);
    const double s = std::sin(th);
    a.a += q * double((s > 0) - (s < 0));
    a.rmin = std::min(a.rmin, r); ++a.n;
    ++res.rows; res.weight_sum += b.weight;
  }
  for (const auto& [key, a] : acc) {
    Cloud c; c.primary = key.first; c.wire = key.second; c.N_fC = a.q; c.rows = a.n;
    c.t_ns = a.t / a.q; c.t_rms_ns = std::sqrt(std::max(0., a.t2 / a.q - c.t_ns * c.t_ns));
    const double tr = a.traw / a.q; c.t_raw_rms_ns = std::sqrt(std::max(0., a.traw2 / a.q - tr * tr));
    c.r_eff_um = std::exp(a.lr / a.q); c.D_U_mm = a.du / a.q; c.D_V_mm = a.dv / a.q;
    c.C1 = a.c1 / a.q; c.S1 = a.s1 / a.q; c.C2 = a.c2 / a.q; c.S2 = a.s2 / a.q; c.A = a.a / a.q;
    c.min_clearance_um = a.rmin - kWireRadiusUm;
    res.clouds.push_back(c);
  }
  return res;
}

inline void Write(std::ostream& o, const Result& r) {
  o << std::setprecision(17)
    << "primary,wire,t_ns,t_rms_ns,N_fC,r_eff_um,D_U_mm,D_V_mm,C1,S1,C2,S2,A,min_clearance_um,base_U,base_V,"
       "t_raw_rms_ns,rows\n";
  for (const auto& c : r.clouds)
    o << c.primary << ',' << c.wire << ',' << c.t_ns << ',' << c.t_rms_ns << ',' << c.N_fC << ',' << c.r_eff_um << ','
      << c.D_U_mm << ',' << c.D_V_mm << ',' << c.C1 << ',' << c.S1 << ',' << c.C2 << ',' << c.S2 << ',' << c.A << ','
      << c.min_clearance_um << ',' << r.base_U << ',' << r.base_V << ',' << c.t_raw_rms_ns << ',' << c.rows << '\n';
}
}  // namespace BirthMoments

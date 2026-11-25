// In-memory FULL-BIRTH K2A cloud moments, accumulated birth by birth inside the Garfield++ avalanche path.
//
// Production purpose: the generator's ionisation callback hands every positive-ion birth to this accumulator BEFORE any
// reservoir sampling, so the cloud moments of the complete birth population are formed without ion_births.root,
// trajectory exports or CSV. Pure C++17, no ROOT, no Garfield (unit-tested synthetically, T0).
//
// Definitions are those of the validated batch path spatial-ablation/BirthMoments.hh (frozen 2100 V strip frame,
// tan(stereo) = 0.1, strip pitch 0.17 cm), copied expression by expression so that, for births arriving in the same
// order, every spatial moment is the SAME floating-point computation:
//   u_mm = 10 u_patch_cm, v_mm = 10 v_absolute_cm, q = weight x e (fC)
//   base_U = llround((w_in_cm - 0.1 u_in_cm) / 0.17), base_V = llround((w_in_cm + 0.1 u_in_cm) / 0.17)   (muon ENTRY)
//   wU_mm = 10 (w + w_in_cm - 0.1 patch_shift_cm - 0.17 base_U), wV_mm = 10 (w + w_in_cm + 0.1 patch_shift_cm - 0.17 base_V)
//   cloud = (Heed primary, wire k = rint(u_mm / 3)), rint = round half to EVEN (numpy rint of the validated pipeline;
//           NOT the generator's llround(u / 0.3) census rule, which rounds half away from zero: the two differ only for a
//           birth exactly on a half-cell boundary, and K2A clouds use the rint rule)
//   r = 1000 hypot(u_mm - 3k, v_mm) um, theta = atan2(v_mm, u_mm - 3k)
//   charge-weighted means of log r (-> r_eff = exp), wU, wV, cos/sin theta, cos/sin 2 theta, sign(sin theta) (-> A),
//   min r - 15 um (min clearance), raw birth-time mean and RMS.
//
// Timing (T1): the K2A time of every cloud of a Heed primary is the primary's weight-mean birth time
//   t_ns = sum(w t) / sum(w) over ALL births of the primary (same summation order as the batch code -> bitwise equal to
//   its per-primary mean). The batch code reaches the same quantity per ion as t + (mean - t) and averages it again, so
//   the batch cloud t_ns agrees only to rounding. The collapsed per-cloud spread t_rms_ns is 0 BY DEFINITION here
//   (the batch code shows rounding noise). The RAW timing spread is kept, not collapsed:
//     Cloud::t_raw_mean_ns / t_raw_rms_ns   raw birth-time mean / RMS of the cloud's own births (t_raw_rms_ns: same formula
//                                           and order as BirthMoments -> bitwise equal)
//     PrimaryRecord::t_mean_ns / t_raw_rms_ns   T1 time and raw birth-time RMS of the whole Heed primary.
//
// Electron step: no extra data. The validated K2A electron term is the instantaneous step at the T1 time evaluated on
// these same cloud moments (electron-step kernel); it is not recomputed here.
//
// Statuses (bit flags; nothing is clamped, nothing is silently dropped):
//   INSIDE_ANODE     birth with r <= 15 um (BirthMoments throws)                          BLOCKING
//   NONFINITE        non-finite u, v, w or t, or weight not finite and > 0                BLOCKING
//   OUTSIDE_PATCH    K2A wire |k| > 4 (outside the nine-wire analytic patch)              BLOCKING
//                    (the generator itself still ABORTS on such a birth, exactly as before; this flag only guarantees
//                     that the accumulator can never hand out moments for one)
//   CENSUS_MISMATCH  callback birth count reported at EndPrimary != births received        BLOCKING
//   SEQUENCE_ERROR   birth outside Begin/EndPrimary, non-increasing primary index, callback ordinal gap,
//                    Finish with an open primary, use after Finish                          BLOCKING
//   ZERO_ION_PRIMARY primary with no ion birth (no cloud)                                  informational
// Event::UsableForCathodeSignal() is false when ANY blocking flag is set on the event; the writer then does NOT write the
// engine-facing moments file (only diagnostics), so invalid moments cannot reach a cathode-signal calculation by
// accident. Engine-side domain statuses (wire not in kernel, r_eff / D outside the Green table, time window) belong to
// the evaluator and are not decided here.
//
// Memory: O(clouds). One small accumulator per wire of the CURRENT primary (<= 9 inside the patch), finalised into a
// Cloud record at EndPrimary; nothing per birth is retained. StateBytes() accounts the container capacities exactly.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace K2AInline {
constexpr double kEfC = 1.602176634e-4;     // elementary charge in fC (as BirthMoments)
constexpr double kWireRadiusUm = 15.0;      // anode radius
constexpr double kStereoTan = .1;           // frozen 2100 V strip frame (BirthMoments::StripGeometry defaults)
constexpr double kStripPitchCm = .17;
constexpr int kPatchHalfWires = 4;          // nine-wire analytic patch, local wires -4..4
constexpr const char* kVersion = "k2a-inline-moments v1";

enum Status : std::uint32_t {
  OK = 0,
  INSIDE_ANODE = 1u << 0,
  NONFINITE = 1u << 1,
  OUTSIDE_PATCH = 1u << 2,
  CENSUS_MISMATCH = 1u << 3,
  SEQUENCE_ERROR = 1u << 4,
  ZERO_ION_PRIMARY = 1u << 5,
};
constexpr std::uint32_t kBlocking = INSIDE_ANODE | NONFINITE | OUTSIDE_PATCH | CENSUS_MISMATCH | SEQUENCE_ERROR;

inline std::string StatusNames(std::uint32_t s) {
  if (s == 0) return "OK";
  std::string o;
  const std::pair<std::uint32_t, const char*> n[] = {{INSIDE_ANODE, "INSIDE_ANODE"}, {NONFINITE, "NONFINITE"},
      {OUTSIDE_PATCH, "OUTSIDE_PATCH"}, {CENSUS_MISMATCH, "CENSUS_MISMATCH"}, {SEQUENCE_ERROR, "SEQUENCE_ERROR"},
      {ZERO_ION_PRIMARY, "ZERO_ION_PRIMARY"}};
  for (const auto& [b, name] : n)
    if (s & b) o += (o.empty() ? "" : "|") + std::string(name);
  if (s & ~(kBlocking | ZERO_ION_PRIMARY)) o += (o.empty() ? "" : "|") + std::string("UNKNOWN");
  return o;
}

struct MuonFrame { double patch_shift_cm = 0, u_in_mm = 0, w_in_mm = 0; };   // as CathodeRef::Archive::Muon writes them

struct Cloud {
  int primary = 0, wire = 0;
  double N_fC = 0;
  double t_ns = 0;            // T1 time = primary weight-mean birth time
  double t_rms_ns = 0;        // collapsed spread: 0 by definition
  double t_raw_mean_ns = 0, t_raw_rms_ns = 0;
  double r_eff_um = 0, D_U_mm = 0, D_V_mm = 0, C1 = 0, S1 = 0, C2 = 0, S2 = 0, A = 0, min_clearance_um = 0;
  std::uint64_t births = 0;
  double weight_sum = 0;
  std::uint32_t status = OK;
};

struct PrimaryRecord {
  int primary = 0;
  std::uint64_t callback_births = 0, births_seen = 0;
  double weight_sum = 0, t_mean_ns = 0, t_raw_rms_ns = 0;
  int n_clouds = 0;
  std::uint32_t status = OK;
};

struct Event {
  MuonFrame frame;
  int base_U = 0, base_V = 0;
  std::vector<Cloud> clouds;            // ordered by (primary, wire), like the batch code
  std::vector<PrimaryRecord> primaries; // in BeginPrimary order (strictly increasing original Heed index)
  std::uint64_t births = 0, orphan_births = 0;
  double weight_sum = 0;
  std::uint32_t status = OK;            // OR of all cloud / primary / sequence flags
  std::size_t state_bytes_peak = 0;
  bool UsableForCathodeSignal() const { return (status & kBlocking) == 0; }
};

class MomentAccumulator {
  struct Acc {
    int wire = 0;
    double q = 0, traw = 0, traw2 = 0, lr = 0, du = 0, dv = 0, c1 = 0, s1 = 0, c2 = 0, s2 = 0, a = 0,
           rmin = std::numeric_limits<double>::infinity(), wsum = 0;
    std::uint64_t n = 0;
    std::uint32_t status = OK;
  };
  Event ev_;
  double wIn_ = 0;
  bool open_ = false, finished_ = false;
  int last_primary_ = std::numeric_limits<int>::min();
  PrimaryRecord cur_;
  double pwt_ = 0, pw_ = 0, pwt2_ = 0;   // current primary: sum w t, sum w, sum w t^2
  std::vector<Acc> act_;                 // current primary's per-wire accumulators (sorted by wire)

  void Flag(std::uint32_t s) { ev_.status |= s; }
  void Track() {
    const std::size_t b = sizeof(*this) + ev_.clouds.capacity() * sizeof(Cloud) +
                          ev_.primaries.capacity() * sizeof(PrimaryRecord) + act_.capacity() * sizeof(Acc);
    ev_.state_bytes_peak = std::max(ev_.state_bytes_peak, b);
  }
  Acc& Find(int k) {
    auto it = std::lower_bound(act_.begin(), act_.end(), k, [](const Acc& x, int w) { return x.wire < w; });
    if (it == act_.end() || it->wire != k) { Acc a; a.wire = k; it = act_.insert(it, a); Track(); }
    return *it;
  }

 public:
  explicit MomentAccumulator(const MuonFrame& f) {
    ev_.frame = f;
    const double uIn = .1 * f.u_in_mm, wIn = .1 * f.w_in_mm;
    wIn_ = wIn;
    ev_.base_U = int(std::llround((wIn - kStereoTan * uIn) / kStripPitchCm));
    ev_.base_V = int(std::llround((wIn + kStereoTan * uIn) / kStripPitchCm));
    act_.reserve(2 * kPatchHalfWires + 1);
    Track();
  }

  // heed_primary_index: the ORIGINAL index of the Heed primary (archive primary_index); strictly increasing.
  void BeginPrimary(int heed_primary_index) {
    if (finished_ || open_ || heed_primary_index <= last_primary_) Flag(SEQUENCE_ERROR);
    if (open_) EndPrimary(cur_.births_seen);   // close it so nothing is lost; the event is already blocked
    open_ = true; last_primary_ = std::max(last_primary_, heed_primary_index);
    cur_ = PrimaryRecord{}; cur_.primary = heed_primary_index;
    pwt_ = pw_ = pwt2_ = 0; act_.clear();
  }

  void AddBirth(double u_cm, double v_cm, double w_cm, double t_ns, double weight = 1.0) {
    if (!open_ || finished_) { Flag(SEQUENCE_ERROR); ++ev_.orphan_births; ++ev_.births; return; }
    ++cur_.births_seen; ++ev_.births;
    if (!(std::isfinite(u_cm) && std::isfinite(v_cm) && std::isfinite(w_cm) && std::isfinite(t_ns) &&
          std::isfinite(weight) && weight > 0)) {
      cur_.status |= NONFINITE; Flag(NONFINITE); return;
    }
    pwt_ += weight * t_ns; pw_ += weight; pwt2_ += weight * t_ns * t_ns;
    // --- per-ion expressions exactly as BirthMoments::Compute ---
    const double u = 10 * u_cm, v = 10 * v_cm;
    const double wu = 10. * (w_cm + wIn_ - kStereoTan * ev_.frame.patch_shift_cm - kStripPitchCm * ev_.base_U);
    const double wv = 10. * (w_cm + wIn_ + kStereoTan * ev_.frame.patch_shift_cm - kStripPitchCm * ev_.base_V);
    const int k = int(std::nearbyint(u / 3));      // round half to even (default FE_TONEAREST)
    const double du = u - 3 * k;
    const double r = 1000 * std::hypot(du, v);
    const double th = std::atan2(v, du);
    const double q = weight * kEfC;
    Acc& a = Find(k);
    if (!(r > kWireRadiusUm)) { a.status |= INSIDE_ANODE; Flag(INSIDE_ANODE); }
    if (k < -kPatchHalfWires || k > kPatchHalfWires) { a.status |= OUTSIDE_PATCH; Flag(OUTSIDE_PATCH); }
    a.q += q; a.traw += q * t_ns; a.traw2 += q * t_ns * t_ns;
    a.lr += q * std::log(r); a.du += q * wu; a.dv += q * wv;
    a.c1 += q * std::cos(th); a.s1 += q * std::sin(th); a.c2 += q * std::cos(2 * th); a.s2 += q * std::sin(2 * th);
    const double s = std::sin(th);
    a.a += q * double((s > 0) - (s < 0));
    a.rmin = std::min(a.rmin, r); ++a.n; a.wsum += weight;
  }

  // Generator form: ordinal is the BirthAccumulator callback ordinal within the primary (0, 1, 2, ...).
  void AddCallbackBirth(double u_cm, double v_cm, double w_cm, double t_ns, std::uint64_t ordinal) {
    if (open_ && ordinal != cur_.births_seen) { cur_.status |= SEQUENCE_ERROR; Flag(SEQUENCE_ERROR); }
    AddBirth(u_cm, v_cm, w_cm, t_ns, 1.0);
  }

  // callback_births: the generator's own exact birth census of this primary (BirthAccumulator::total).
  void EndPrimary(std::uint64_t callback_births) {
    if (!open_ || finished_) { Flag(SEQUENCE_ERROR); return; }
    open_ = false;
    cur_.callback_births = callback_births;
    if (callback_births != cur_.births_seen) cur_.status |= CENSUS_MISMATCH;
    cur_.weight_sum = pw_;
    if (cur_.births_seen == 0) cur_.status |= ZERO_ION_PRIMARY;
    const double tm = pw_ > 0 ? pwt_ / pw_ : std::numeric_limits<double>::quiet_NaN();
    cur_.t_mean_ns = tm;
    cur_.t_raw_rms_ns = pw_ > 0 ? std::sqrt(std::max(0., pwt2_ / pw_ - tm * tm)) : std::numeric_limits<double>::quiet_NaN();
    for (const auto& a : act_) {
      Cloud c; c.primary = cur_.primary; c.wire = a.wire; c.N_fC = a.q; c.births = a.n; c.weight_sum = a.wsum;
      c.status = a.status | (cur_.status & (NONFINITE | SEQUENCE_ERROR | CENSUS_MISMATCH));
      c.t_ns = tm; c.t_rms_ns = 0.;
      const double tr = a.traw / a.q;
      c.t_raw_mean_ns = tr; c.t_raw_rms_ns = std::sqrt(std::max(0., a.traw2 / a.q - tr * tr));
      c.r_eff_um = std::exp(a.lr / a.q); c.D_U_mm = a.du / a.q; c.D_V_mm = a.dv / a.q;
      c.C1 = a.c1 / a.q; c.S1 = a.s1 / a.q; c.C2 = a.c2 / a.q; c.S2 = a.s2 / a.q; c.A = a.a / a.q;
      c.min_clearance_um = a.rmin - kWireRadiusUm;
      cur_.status |= c.status & ~ZERO_ION_PRIMARY;
      ev_.clouds.push_back(c);
      ++cur_.n_clouds;
    }
    Flag(cur_.status);
    ev_.weight_sum += pw_;
    ev_.primaries.push_back(cur_);
    act_.clear();
    Track();
  }

  Event Finish() {
    if (finished_) { Flag(SEQUENCE_ERROR); return ev_; }
    if (open_) { Flag(SEQUENCE_ERROR); EndPrimary(cur_.births_seen); }
    finished_ = true;
    std::stable_sort(ev_.clouds.begin(), ev_.clouds.end(),
                     [](const Cloud& x, const Cloud& y) { return std::make_pair(x.primary, x.wire) < std::make_pair(y.primary, y.wire); });
    Track();
    return ev_;
  }
  std::size_t StateBytes() const { return ev_.state_bytes_peak; }
};

// ---------------------------------------------------------------------------------------------------------------------
// Output. PREFIX.inline.json always; PREFIX.clouds.csv and PREFIX.primaries.csv (all records, with statuses) always;
// PREFIX.moments.csv (engine-facing, the BirthMoments::Write schema, so existing K2A tools read it unchanged) and
// PREFIX.moments.csv.frame.json ONLY when the event is usable for a cathode signal.
inline std::string Hex(std::uint32_t s) { std::ostringstream o; o << "0x" << std::hex << s; return o.str(); }

inline void WriteMomentsCsv(std::ostream& o, const Event& e) {
  o << std::setprecision(17)
    << "primary,wire,t_ns,t_rms_ns,N_fC,r_eff_um,D_U_mm,D_V_mm,C1,S1,C2,S2,A,min_clearance_um,base_U,base_V,"
       "t_raw_rms_ns,rows\n";
  for (const auto& c : e.clouds)
    o << c.primary << ',' << c.wire << ',' << c.t_ns << ',' << c.t_rms_ns << ',' << c.N_fC << ',' << c.r_eff_um << ','
      << c.D_U_mm << ',' << c.D_V_mm << ',' << c.C1 << ',' << c.S1 << ',' << c.C2 << ',' << c.S2 << ',' << c.A << ','
      << c.min_clearance_um << ',' << e.base_U << ',' << e.base_V << ',' << c.t_raw_rms_ns << ',' << c.births << '\n';
}

inline void WriteCloudsCsv(std::ostream& o, const Event& e) {
  o << std::setprecision(17)
    << "primary,wire,births,weight_sum,N_fC,t_ns,t_rms_ns,t_raw_mean_ns,t_raw_rms_ns,r_eff_um,D_U_mm,D_V_mm,C1,S1,C2,S2,A,"
       "min_clearance_um,status,status_names\n";
  for (const auto& c : e.clouds)
    o << c.primary << ',' << c.wire << ',' << c.births << ',' << c.weight_sum << ',' << c.N_fC << ',' << c.t_ns << ','
      << c.t_rms_ns << ',' << c.t_raw_mean_ns << ',' << c.t_raw_rms_ns << ',' << c.r_eff_um << ',' << c.D_U_mm << ','
      << c.D_V_mm << ',' << c.C1 << ',' << c.S1 << ',' << c.C2 << ',' << c.S2 << ',' << c.A << ',' << c.min_clearance_um
      << ',' << c.status << ',' << StatusNames(c.status) << '\n';
}

inline void WritePrimariesCsv(std::ostream& o, const Event& e) {
  o << std::setprecision(17) << "primary,callback_births,births_seen,weight_sum,t_mean_ns,t_raw_rms_ns,n_clouds,status,status_names\n";
  for (const auto& p : e.primaries)
    o << p.primary << ',' << p.callback_births << ',' << p.births_seen << ',' << p.weight_sum << ',' << p.t_mean_ns << ','
      << p.t_raw_rms_ns << ',' << p.n_clouds << ',' << p.status << ',' << StatusNames(p.status) << '\n';
}

// extra_json: already-formatted '"key": value' pairs (may be empty), appended to the summary object.
inline bool WriteOutputs(const std::string& prefix, const Event& e, const std::string& extra_json = "") {
  auto must = [](std::ofstream& f, const std::string& p) { if (!f) throw std::runtime_error("K2A inline: cannot write " + p); };
  { std::ofstream f(prefix + ".clouds.csv"); must(f, prefix + ".clouds.csv"); WriteCloudsCsv(f, e); f.close(); must(f, prefix + ".clouds.csv"); }
  { std::ofstream f(prefix + ".primaries.csv"); must(f, prefix + ".primaries.csv"); WritePrimariesCsv(f, e); f.close(); must(f, prefix + ".primaries.csv"); }
  const bool usable = e.UsableForCathodeSignal();
  if (usable) {
    std::ofstream f(prefix + ".moments.csv"); must(f, prefix + ".moments.csv"); WriteMomentsCsv(f, e); f.close(); must(f, prefix + ".moments.csv");
    std::ofstream g(prefix + ".moments.csv.frame.json"); must(g, "frame.json");
    g << std::setprecision(17) << "{\"patch_shift_cm\": " << e.frame.patch_shift_cm << ", \"u_in_mm\": " << e.frame.u_in_mm
      << ", \"w_in_mm\": " << e.frame.w_in_mm << ", \"base_U\": " << e.base_U << ", \"base_V\": " << e.base_V
      << ", \"stereo_tan\": " << kStereoTan << ", \"strip_pitch_mm\": " << 1.7 << "}\n";
    g.close(); must(g, "frame.json");
  }
  std::uint64_t zero = 0;
  for (const auto& p : e.primaries) zero += (p.status & ZERO_ION_PRIMARY) != 0;
  std::ofstream j(prefix + ".inline.json"); must(j, prefix + ".inline.json");
  j << std::setprecision(17) << "{\"tool\": \"" << kVersion << "\", \"usable_for_cathode_signal\": " << (usable ? "true" : "false")
    << ", \"status\": " << e.status << ", \"status_names\": \"" << StatusNames(e.status) << "\", \"status_hex\": \"" << Hex(e.status)
    << "\", \"n_primaries\": " << e.primaries.size() << ", \"n_primaries_zero_ions\": " << zero << ", \"n_clouds\": " << e.clouds.size()
    << ", \"births\": " << e.births << ", \"orphan_births\": " << e.orphan_births << ", \"weight_sum\": " << e.weight_sum
    << ", \"base_U\": " << e.base_U << ", \"base_V\": " << e.base_V << ", \"patch_shift_cm\": " << e.frame.patch_shift_cm
    << ", \"u_in_mm\": " << e.frame.u_in_mm << ", \"w_in_mm\": " << e.frame.w_in_mm
    << ", \"accumulator_state_bytes_peak\": " << e.state_bytes_peak
    << ", \"moments_file\": " << (usable ? "\"" + prefix + ".moments.csv\"" : std::string("null"))
    << (extra_json.empty() ? "" : ", " + extra_json) << "}\n";
  j.close(); must(j, prefix + ".inline.json");
  return usable;
}
}  // namespace K2AInline

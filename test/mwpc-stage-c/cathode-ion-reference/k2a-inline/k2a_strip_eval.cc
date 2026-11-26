// Standalone compiled K2A absolute-strip evaluator (no Garfield, no ROOT, no Python): kernel opened once, moments
// evaluated, timing and peak memory measured independently of the avalanche generator.
//
//   k2a_strip_eval --kernel DIR (--kernel-id SHA | --no-verify) --moments FILE.csv --out NEWDIR
//                  [--shifts groups.json] [--model K2A] [--half-strips 5] [--window-ns 7000] [--repeat N] [--preload]
//                  [--displace-event M_U M_V]
//
// FILE.csv      any moments CSV with the BirthMoments / K2AInline / s_kernel.py columns (tests and reproduction only;
//               production passes the in-memory K2AInline::Event to K2AStrip::Evaluator directly)
// --shifts      impose the per-primary (m_U, m_V) of a groups.json ("primary_shift"), e.g. the E2b grouping; the
//               computed production shifts are still reported and compared
// --repeat N    evaluate N times (cold = first, warm = the rest); every repetition must be bit-identical
// --preload     copy K.npy into memory instead of mapping it read-only
// --displace-event M_U M_V  TEST ONLY: add M_U P to every D_U and M_V P to every D_V (the same event displaced by
//               M strips); the response must equal the undisplaced one on labels shifted by M
// Writes NEWDIR/replay/waveforms.csv, NEWDIR/prompt/prompt_waveforms.csv (fC/ns, union labels U then V ascending),
// clouds.csv, shifts.csv, groups.csv, result.json. Exit 0 usable, 4 blocked by a domain/status flag, 2 error.
#include <sys/resource.h>

#include <cstring>
#include <filesystem>
#include <iostream>

#include "K2AStripEvaluator.hh"

namespace {
long PeakRssKiB() { rusage ru{}; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss; }
bool SameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}
std::string Js(const std::string& s) { std::string o = "\""; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; } return o + "\""; }
}  // namespace

int main(int argc, char** argv) {
  try {
    const long rss0 = PeakRssKiB();
    std::string kdir, kid, mom, out, shifts;
    K2AStrip::Config cfg;
    int repeat = 1, dU = 0, dV = 0;
    bool verify = true, preload = false, displace = false;
    for (int i = 1; i < argc; ++i) {
      const std::string a = argv[i];
      auto next = [&]() -> std::string { if (i + 1 >= argc) throw std::runtime_error(a + " needs a value"); return argv[++i]; };
      if (a == "--kernel") kdir = next(); else if (a == "--kernel-id") kid = next(); else if (a == "--no-verify") verify = false;
      else if (a == "--moments") mom = next(); else if (a == "--out") out = next(); else if (a == "--shifts") shifts = next();
      else if (a == "--model") cfg.model = next(); else if (a == "--half-strips") cfg.half_strips = std::stoi(next());
      else if (a == "--window-ns") cfg.window_ns = std::stod(next()); else if (a == "--repeat") repeat = std::stoi(next());
      else if (a == "--preload") preload = true;
      else if (a == "--displace-event") { dU = std::stoi(next()); dV = std::stoi(next()); displace = true; }
      else throw std::runtime_error("unknown option " + a);
    }
    if (kdir.empty() || mom.empty() || out.empty() || repeat < 1) throw std::runtime_error("need --kernel, --moments, --out (see header)");
    if (verify && kid.empty()) throw std::runtime_error("give --kernel-id SHA (or --no-verify for a timing-only run)");
    if (std::filesystem::exists(out)) throw std::runtime_error("refusing: " + out + " exists");
    const auto k = K2AStrip::KernelTable::Open(kdir, verify ? kid : "", verify, preload);
    const long rssKernel = PeakRssKiB();
    K2AStrip::Evaluator ev(*k, cfg);
    auto event = K2AStrip::ReadMomentsCsv(mom);
    if (displace)
      for (auto& c : event.clouds) { c.D_U_mm += k->pitch_mm * dU; c.D_V_mm += k->pitch_mm * dV; }
    std::map<int, std::pair<int, int>> table;
    if (!shifts.empty()) {
      const auto j = K2AKernelIO::Json::ParseFile(shifts)["primary_shift"];
      for (const auto& [p, v] : j.obj) table[std::stoi(p)] = {int(v.arr.at(0).Num()), int(v.arr.at(1).Num())};
    }
    std::vector<double> times;
    K2AStrip::Result first;
    bool identical = true;
    for (int r = 0; r < repeat; ++r) {
      auto res = ev.Evaluate(event, shifts.empty() ? nullptr : &table);
      times.push_back(res.evaluate_s);
      if (r == 0) first = std::move(res);
      else identical = identical && SameBits(res.U.ion, first.U.ion) && SameBits(res.V.ion, first.V.ion) &&
                       SameBits(res.U.electron, first.U.electron) && SameBits(res.V.electron, first.V.electron) &&
                       res.U.strips == first.U.strips && res.V.strips == first.V.strips && res.status == first.status;
    }
    const long rssEnd = PeakRssKiB();
    std::filesystem::create_directories(out);
    K2AStrip::WriteResult(out, first, ev.Edges());
    double warm = 0; std::vector<double> w(times.begin() + 1, times.end());
    std::sort(w.begin(), w.end());
    if (!w.empty()) warm = w[w.size() / 2];
    double qU = 0, qV = 0;
    for (double x : first.U.ion) qU += x * first.dt_ns;
    for (double x : first.V.ion) qV += x * first.dt_ns;
    std::ofstream j(out + "/result.json");
    j << std::setprecision(10) << "{\"tool\": \"" << K2AStrip::kVersion << "\", \"moments\": " << Js(mom)
      << ", \"moments_sha256\": \"" << K2AKernelIO::Sha256::OfFile(mom) << "\", \"kernel\": " << Js(kdir)
      << ", \"kernel_id\": \"" << k->kernel_id << "\", \"kernel_verified\": " << (verify ? "true" : "false")
      << ", \"transport_id\": \"" << k->transport_id << "\", \"weighting_id\": \"" << k->weighting_id << "\", \"model\": \"" << cfg.model
      << "\", \"half_strips\": " << cfg.half_strips << ", \"window_ns\": " << cfg.window_ns << ", \"dt_ns\": " << cfg.dt_ns
      << ", \"shifts_from\": " << (shifts.empty() ? std::string("\"production rule floor(Dbar/P + 1/2)\"") : Js(shifts))
      << ", \"shifts_differing_from_computed\": " << first.shifts_differing_from_computed
      << ", \"displace_event_TEST_ONLY\": " << (displace ? "[" + std::to_string(dU) + ", " + std::to_string(dV) + "]" : std::string("null"))
      << ", \"status\": " << first.status << ", \"status_names\": \"" << K2AStrip::StatusNames(first.status)
      << "\", \"usable\": " << (first.usable ? "true" : "false") << ", \"clouds\": " << first.clouds.size()
      << ", \"primaries\": " << first.shifts.size() << ", \"groups\": " << first.groups.size()
      << ", \"strips_U\": " << first.U.strips.size() << ", \"strips_V\": " << first.V.strips.size() << ", \"bins\": " << first.nbins
      << std::setprecision(17) << ", \"cloud_charge_fC\": " << first.charge_fC << ", \"blocked_charge_fC\": " << first.blocked_charge_fC
      << ", \"window_ion_charge_U_fC\": " << qU << ", \"window_ion_charge_V_fC\": " << qV << std::setprecision(6)
      << ", \"timing\": {\"kernel_verify_s\": " << k->verify_s << ", \"kernel_load_s\": " << k->load_s
      << ", \"evaluate_cold_s\": " << times[0] << ", \"evaluate_warm_median_s\": " << (w.empty() ? std::string("null") : [&] { std::ostringstream o; o << std::setprecision(6) << warm; return o.str(); }())
      << ", \"repeat\": " << repeat
      << ", \"repeats_bit_identical\": " << (identical ? "true" : "false") << "}, \"memory\": {\"peak_rss_kib_start\": " << rss0
      << ", \"peak_rss_kib_after_kernel\": " << rssKernel << ", \"peak_rss_kib_end\": " << rssEnd
      << ", \"kernel_K_bytes\": " << k->K.Bytes() << ", \"kernel_K_mapped\": " << (k->K.Mapped() ? "true" : "false")
      << ", \"kernel_resident_copy_bytes\": " << k->ResidentBytes() << "}}\n";
    std::cout << (first.usable ? "PASS" : "BLOCKED") << " k2a_strip_eval " << cfg.model << ": " << first.clouds.size() << " clouds, "
              << first.groups.size() << " groups -> U " << first.U.strips.size() << " / V " << first.V.strips.size()
              << " absolute strips; status " << K2AStrip::StatusNames(first.status) << "; evaluate cold " << times[0] * 1e3
              << " ms, warm " << (w.empty() ? std::string("n/a") : std::to_string(warm * 1e3)) << " ms; peak RSS " << rssEnd << " KiB" << (identical ? "" : "; REPEATS DIFFER") << "\n";
    if (!identical) return 2;
    return first.usable ? 0 : 4;
  } catch (const std::exception& e) { std::cerr << "k2a_strip_eval: " << e.what() << "\n"; return 2; }
}

// Read-only birth-record -> K2A cloud moments for one archived 2100 V muon, WITHOUT ion transport.
//
//   p3_birth_moments ion_births.root OUT_MOMENTS.csv [--repeat N]
//   options (single and batch): --stereo-tan T --strip-pitch-mm P   readout strip geometry of the strip frame
//                               (default 0.1 / 1.7 = frozen reference); each output also gets <out>.frame.json
//                               (patch_shift_cm, u_in_mm, w_in_mm, base_U/V, geometry) for re-striping exports
//   p3_birth_moments --batch LIST.txt OUTDIR            persistent process over many events: LIST lines are
//                                                       "EVENT_ID PATH/ion_births.root"; writes OUTDIR/EVENT_ID.moments.csv
//                                                       and OUTDIR/batch_timing.json (process start-up, per-event cold/warm)
//
// Reading and validation use the FROZEN archive reader Plane2100Archive.hh unchanged (weights = exact/saved per
// primary, complete census, primary association). The physical-frame metadata (muons: patch_shift_cm, u_in_mm,
// w_in_mm) and the local/physical consistency check (ions: u_chamber_cm, w_chamber_cm) are read exactly as the frozen
// exporter's ChamberFrame does. Moments: BirthMoments.hh. Writes OUT_MOMENTS.csv and OUT_MOMENTS.csv.json with the
// timing split (open+read+validate, moment construction), rows, peak RSS and the archive size.
#include <TFile.h>
#include <TTree.h>
#include <sys/resource.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "BirthMoments.hh"
#include "Plane2100Archive.hh"

namespace {
struct OneResult { BirthMoments::Result res; double readS = 0, momS = 0; long long nTrue = 0; BirthMoments::Frame frame{}; };

BirthMoments::StripGeometry gGeom;
double gPitchMm = 1.7;   // as given on the command line (10 * pitch_cm is not exactly 1.7 in double)

void WriteFrame(const std::string& path, const OneResult& r) {
  std::ofstream f(path);
  f << std::setprecision(17) << "{\"patch_shift_cm\": " << r.frame.patch_shift_cm << ", \"u_in_mm\": " << r.frame.u_in_mm
    << ", \"w_in_mm\": " << r.frame.w_in_mm << ", \"base_U\": " << r.res.base_U << ", \"base_V\": " << r.res.base_V
    << ", \"stereo_tan\": " << gGeom.stereo_tan << ", \"strip_pitch_mm\": " << gPitchMm << "}\n";
  BirthMoments::Require(bool(f), "frame write failed");
}

OneResult ProcessOne(const std::string& in) {
  using clk = std::chrono::steady_clock;
  OneResult o;
  const auto t0 = clk::now();
  NA60Plane::Archive archive(in);
  BirthMoments::Require(archive.mode == "weighted" && archive.Entries() > 0, "need a complete weighted archive");
  std::unique_ptr<TFile> f(TFile::Open(in.c_str(), "READ"));
  BirthMoments::Require(f && !f->IsZombie(), "cannot reopen archive for frame metadata");
  auto* mu = dynamic_cast<TTree*>(f->Get("muons"));
  auto* ions = dynamic_cast<TTree*>(f->Get("ions"));
  BirthMoments::Require(mu && ions && mu->GetEntries() == 1, "missing muons/ions tree");
  for (const char* s : {"patch_shift_cm", "u_in_mm", "w_in_mm"})
    BirthMoments::Require(mu->GetBranch(s) != nullptr, std::string("missing muon branch ") + s);
  for (const char* s : {"u_chamber_cm", "w_chamber_cm"})
    BirthMoments::Require(ions->GetBranch(s) != nullptr, std::string("missing ion branch ") + s);
  BirthMoments::Frame fr{};
  mu->SetBranchAddress("patch_shift_cm", &fr.patch_shift_cm);
  mu->SetBranchAddress("u_in_mm", &fr.u_in_mm);
  mu->SetBranchAddress("w_in_mm", &fr.w_in_mm);
  mu->GetEntry(0);
  double uCh = 0, wCh = 0;
  ions->SetBranchAddress("u_chamber_cm", &uCh);
  ions->SetBranchAddress("w_chamber_cm", &wCh);
  std::vector<BirthMoments::Birth> births;
  births.reserve(archive.Entries());
  for (long long j = 0; j < archive.Entries(); ++j) {
    const auto b = archive.Get(j);
    BirthMoments::Require(ions->GetEntry(j) > 0, "physical ion location row is absent");
    BirthMoments::Require(std::abs(uCh - (b.u + fr.patch_shift_cm)) < 1e-8 &&
                          std::abs(wCh - (b.w + .1 * fr.w_in_mm)) < 1e-8, "local/physical ion-coordinate mapping failed");
    births.push_back({b.primary, b.u, b.v, b.w, b.t, b.weight});
  }
  archive.CheckComplete();
  const auto t1 = clk::now();
  o.res = BirthMoments::Compute(births, fr, gGeom);
  o.frame = fr;
  const auto t2 = clk::now();
  o.readS = std::chrono::duration<double>(t1 - t0).count();
  o.momS = std::chrono::duration<double>(t2 - t1).count();
  o.nTrue = archive.nTrue;
  BirthMoments::Require(std::abs(o.res.weight_sum - double(archive.nTrue)) < 1e-6 * std::max(1., double(archive.nTrue)),
                        "weighted census not conserved");
  return o;
}

long RssKiB() { rusage ru{}; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss; }
}  // namespace

int main(int argc, char** argv) {
  const auto tProcess = std::chrono::steady_clock::now();
  try {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
      const std::string x = argv[i];
      if (x == "--stereo-tan" && i + 1 < argc) gGeom.stereo_tan = std::stod(argv[++i]);
      else if (x == "--strip-pitch-mm" && i + 1 < argc) { gPitchMm = std::stod(argv[++i]); gGeom.pitch_cm = .1 * gPitchMm; }
      else args.push_back(x);
    }
    BirthMoments::Require(gGeom.stereo_tan > 0 && gGeom.stereo_tan < 1 && gGeom.pitch_cm > 0, "bad strip geometry");
    argc = int(args.size()) + 1;
    std::vector<char*> av{argv[0]};
    for (auto& x : args) av.push_back(x.data());
    argv = av.data();
    if (argc >= 4 && std::string(argv[1]) == "--batch") {
      const std::string list = argv[2]; const std::filesystem::path outdir = argv[3];
      BirthMoments::Require(!std::filesystem::exists(outdir), "output directory exists");
      std::filesystem::create_directories(outdir);
      std::ifstream L(list); BirthMoments::Require(bool(L), "cannot read " + list);
      std::string id, path;
      std::ofstream j(outdir / "batch_timing.json");
      j << std::setprecision(10) << "{\"tool\": \"p3_birth_moments --batch\", \"events\": [\n";
      int n = 0; double sumRead = 0, sumMom = 0;
      while (L >> id >> path) {
        const auto r = ProcessOne(path);
        const auto out = outdir / (id + ".moments.csv");
        BirthMoments::Require(!std::filesystem::exists(out), "duplicate event id " + id);
        std::ofstream o(out); BirthMoments::Write(o, r.res); o.close(); BirthMoments::Require(bool(o), "write failed");
        WriteFrame(out.string() + ".frame.json", r);
        j << (n ? ",\n" : "") << "  {\"id\": \"" << id << "\", \"archive\": \"" << path << "\", \"archive_bytes\": "
          << std::filesystem::file_size(path) << ", \"rows\": " << r.res.rows << ", \"clouds\": " << r.res.clouds.size()
          << ", \"read_validate_s\": " << r.readS << ", \"moments_s\": " << r.momS << ", \"warm\": " << (n ? "true" : "false")
          << ", \"peak_rss_kib_so_far\": " << RssKiB() << "}";
        sumRead += r.readS; sumMom += r.momS; ++n;
        std::cout << "event " << id << ": rows=" << r.res.rows << " clouds=" << r.res.clouds.size() << " read+validate "
                  << r.readS * 1e3 << " ms, moments " << r.momS * 1e3 << " ms" << (n == 1 ? " (cold: first event)" : "") << "\n";
      }
      const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - tProcess).count();
      j << "\n], \"n_events\": " << n << ", \"sum_read_validate_s\": " << sumRead << ", \"sum_moments_s\": " << sumMom
        << ", \"process_wall_s\": " << wall << ", \"peak_rss_kib\": " << RssKiB() << "}\n";
      std::cout << "PASS birth moments batch events=" << n << " process wall " << wall << " s\n";
      return 0;
    }
    BirthMoments::Require(argc == 3 || (argc == 5 && std::string(argv[3]) == "--repeat"),
                          "usage: p3_birth_moments ion_births.root OUT_MOMENTS.csv [--repeat N] | --batch LIST OUTDIR");
    const std::string in = argv[1], out = argv[2];
    const int repeat = argc == 5 ? std::stoi(argv[4]) : 1;
    BirthMoments::Require(!std::filesystem::exists(out), "output exists");
    double readS = 1e300, momS = 1e300;
    BirthMoments::Result res;
    OneResult last;
    for (int rep = 0; rep < repeat; ++rep) {
      auto r = ProcessOne(in);
      readS = std::min(readS, r.readS); momS = std::min(momS, r.momS); last = r; res = std::move(r.res);
    }
    WriteFrame(out + ".frame.json", last);
    std::ofstream o(out);
    BirthMoments::Write(o, res);
    o.close();
    BirthMoments::Require(bool(o), "write failed");
    std::ofstream j(out + ".json");
    j << std::setprecision(10) << "{\"tool\": \"p3_birth_moments\", \"archive\": \"" << in << "\", \"archive_bytes\": "
      << std::filesystem::file_size(in) << ", \"rows\": " << res.rows << ", \"clouds\": " << res.clouds.size()
      << ", \"repeat\": " << repeat << ", \"read_validate_s_min\": " << readS << ", \"moments_s_min\": " << momS
      << ", \"peak_rss_kib\": " << RssKiB() << ", \"moments_bytes\": " << std::filesystem::file_size(out) << "}\n";
    std::cout << "PASS birth moments rows=" << res.rows << " clouds=" << res.clouds.size() << " read+validate "
              << readS * 1e3 << " ms, moments " << momS * 1e3 << " ms (min over " << repeat << ")\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
}

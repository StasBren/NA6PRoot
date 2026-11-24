// Read-only birth-record -> K2A cloud moments for one archived 2100 V muon, WITHOUT ion transport.
//
//   p3_birth_moments ion_births.root OUT_MOMENTS.csv [--repeat N]
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

int main(int argc, char** argv) {
  try {
    BirthMoments::Require(argc == 3 || (argc == 5 && std::string(argv[3]) == "--repeat"),
                          "usage: p3_birth_moments ion_births.root OUT_MOMENTS.csv [--repeat N]");
    const std::string in = argv[1], out = argv[2];
    const int repeat = argc == 5 ? std::stoi(argv[4]) : 1;
    BirthMoments::Require(!std::filesystem::exists(out), "output exists");
    using clk = std::chrono::steady_clock;
    double readS = 1e300, momS = 1e300;
    BirthMoments::Result res;
    for (int rep = 0; rep < repeat; ++rep) {
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
      res = BirthMoments::Compute(births, fr);
      const auto t2 = clk::now();
      readS = std::min(readS, std::chrono::duration<double>(t1 - t0).count());
      momS = std::min(momS, std::chrono::duration<double>(t2 - t1).count());
      BirthMoments::Require(std::abs(res.weight_sum - double(archive.nTrue)) < 1e-6 * std::max(1., double(archive.nTrue)),
                            "weighted census not conserved");
    }
    std::ofstream o(out);
    BirthMoments::Write(o, res);
    o.close();
    BirthMoments::Require(bool(o), "write failed");
    rusage ru{};
    getrusage(RUSAGE_SELF, &ru);
    std::ofstream j(out + ".json");
    j << std::setprecision(10) << "{\"tool\": \"p3_birth_moments\", \"archive\": \"" << in << "\", \"archive_bytes\": "
      << std::filesystem::file_size(in) << ", \"rows\": " << res.rows << ", \"clouds\": " << res.clouds.size()
      << ", \"repeat\": " << repeat << ", \"read_validate_s_min\": " << readS << ", \"moments_s_min\": " << momS
      << ", \"peak_rss_kib\": " << ru.ru_maxrss << ", \"moments_bytes\": " << std::filesystem::file_size(out) << "}\n";
    std::cout << "PASS birth moments rows=" << res.rows << " clouds=" << res.clouds.size() << " read+validate "
              << readS * 1e3 << " ms, moments " << momS * 1e3 << " ms (min over " << repeat << ")\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 2; }
}

// Read-only dump of a 2100 V cathode-ion archive (full or weighted) to text, for the INDEPENDENT T1 references.
//
//   k2a_archive_dump ion_births.root OUTDIR
//
// Reads the ROOT trees directly with TTree branch addresses (deliberately NOT through Plane2100Archive.hh or
// BirthMoments.hh, so the numpy/fsum reference built on this dump shares no code with the batch path). Writes
//   OUTDIR/births.csv     primary_index,callback_ordinal,u_patch_cm,v_absolute_cm,w_relative_cm,t_birth_ns,
//                         inverse_inclusion_weight,local_wire_index,garfield_seed     (archive row order)
//   OUTDIR/primaries.csv  primary_index,cluster_index,garfield_seed,successful,ion_births_exact,ion_births_saved,
//                         inverse_inclusion_probability
//   OUTDIR/muon.json      storage_mode, reservoir, heed_seed, crossing_id, station, patch_shift_cm, u_in_mm, w_in_mm,
//                         base_global_wire, rows
// Doubles are written with 17 significant digits (exact round trip). The archive is opened READ-only.
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
void Req(bool ok, const std::string& m) { if (!ok) throw std::runtime_error(m); }
template <class T>
void Addr(TTree* t, const char* n, T* x) {
  Req(t->GetBranch(n) != nullptr, std::string("missing branch ") + n);
  Req(t->SetBranchAddress(n, x) >= 0, std::string("bad branch type ") + n);
}
std::string Js(const std::string& s) {
  std::string o = "\"";
  for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; }
  return o + "\"";
}
}  // namespace

int main(int argc, char** argv) {
  try {
    Req(argc == 3, "usage: k2a_archive_dump ion_births.root OUTDIR");
    const std::filesystem::path out = argv[2];
    Req(!std::filesystem::exists(out), "OUTDIR exists");
    std::unique_ptr<TFile> f(TFile::Open(argv[1], "READ"));
    Req(f && !f->IsZombie(), std::string("cannot open ") + argv[1]);
    auto* mu = dynamic_cast<TTree*>(f->Get("muons"));
    auto* pr = dynamic_cast<TTree*>(f->Get("primaries"));
    auto* io = dynamic_cast<TTree*>(f->Get("ions"));
    Req(mu && pr && io && mu->GetEntries() == 1, "missing muons/primaries/ions tree");
    std::filesystem::create_directories(out);
    std::string *mode = nullptr, *id = nullptr;
    int K = 0, heed = 0, station = 0, base = 0;
    double ps = 0, ui = 0, wi = 0;
    Addr(mu, "storage_mode", &mode); Addr(mu, "crossing_id", &id); Addr(mu, "reservoir_per_primary", &K);
    Addr(mu, "heed_seed", &heed); Addr(mu, "station", &station); Addr(mu, "base_global_wire", &base);
    Addr(mu, "patch_shift_cm", &ps); Addr(mu, "u_in_mm", &ui); Addr(mu, "w_in_mm", &wi);
    mu->GetEntry(0);
    int p = 0, cl = 0, gs = 0, ok = 0;
    long long ne = 0, ns = 0;
    double ipi = 0;
    Addr(pr, "primary_index", &p); Addr(pr, "cluster_index", &cl); Addr(pr, "garfield_seed", &gs); Addr(pr, "successful", &ok);
    Addr(pr, "ion_births_exact", &ne); Addr(pr, "ion_births_saved", &ns); Addr(pr, "inverse_inclusion_probability", &ipi);
    {
      std::ofstream o(out / "primaries.csv");
      o << std::setprecision(17) << "primary_index,cluster_index,garfield_seed,successful,ion_births_exact,ion_births_saved,inverse_inclusion_probability\n";
      for (Long64_t j = 0; j < pr->GetEntries(); ++j) {
        Req(pr->GetEntry(j) > 0, "primaries row unreadable");
        o << p << ',' << cl << ',' << gs << ',' << ok << ',' << ne << ',' << ns << ',' << ipi << '\n';
      }
      Req(bool(o), "primaries.csv write failed");
    }
    int ip = 0, lw = 0, igs = 0;
    long long ord = 0;
    double u = 0, v = 0, w = 0, t = 0, wt = 0;
    Addr(io, "primary_index", &ip); Addr(io, "callback_ordinal", &ord); Addr(io, "u_patch_cm", &u); Addr(io, "v_absolute_cm", &v);
    Addr(io, "w_relative_cm", &w); Addr(io, "t_birth_ns", &t); Addr(io, "inverse_inclusion_weight", &wt);
    Addr(io, "local_wire_index", &lw); Addr(io, "garfield_seed", &igs);
    {
      std::ofstream o(out / "births.csv");
      o << std::setprecision(17) << "primary_index,callback_ordinal,u_patch_cm,v_absolute_cm,w_relative_cm,t_birth_ns,inverse_inclusion_weight,local_wire_index,garfield_seed\n";
      for (Long64_t j = 0; j < io->GetEntries(); ++j) {
        Req(io->GetEntry(j) > 0, "ions row unreadable");
        o << ip << ',' << ord << ',' << u << ',' << v << ',' << w << ',' << t << ',' << wt << ',' << lw << ',' << igs << '\n';
      }
      Req(bool(o), "births.csv write failed");
    }
    std::ofstream j(out / "muon.json");
    j << std::setprecision(17) << "{\"archive\": " << Js(argv[1]) << ", \"storage_mode\": " << Js(*mode) << ", \"reservoir\": " << K
      << ", \"heed_seed\": " << heed << ", \"crossing_id\": " << Js(*id) << ", \"station\": " << station
      << ", \"base_global_wire\": " << base << ", \"patch_shift_cm\": " << ps << ", \"u_in_mm\": " << ui << ", \"w_in_mm\": " << wi
      << ", \"n_primaries\": " << pr->GetEntries() << ", \"rows\": " << io->GetEntries() << "}\n";
    Req(bool(j), "muon.json write failed");
    std::cout << "PASS archive dump mode=" << *mode << " primaries=" << pr->GetEntries() << " rows=" << io->GetEntries() << "\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << "k2a_archive_dump: " << e.what() << "\n"; return 2; }
}

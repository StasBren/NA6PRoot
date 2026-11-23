// V0 archive statistics (READ-ONLY on ion_births.root and microscopic_wire_seed.csv). Conventions: see V0Core.hh.
//
// Real archives (needs -DV0_ROOT):
//   v0_archive_stats --list FILES.txt --exclude-index-csv sealed_validation20.csv [--out-dir DIR] [--max-muons N]
//   FILES.txt: one muon_NNNN/ion_births.root path per line; microscopic_wire_seed.csv is read from the same directory.
//   Holdout guard: --exclude-index-csv (column "index") and/or --exclude-slots 0101,0435 are REQUIRED, unless
//   --no-holdout-exclusion is given explicitly (then every output is labelled NOT-FOR-CALIBRATION).
//   optional: --ion-dump-muons N   raw reservoir ions of the first N muons -> v0_ion_examples.csv
//   optional: --lut500 garfield_electrons.csv[,more.csv]   1800 V LUT500 raw drift trials, binned like the 2100 V map, never merged
//   v0_archive_stats --dump-schema muon_NNNN/ion_births.root   (prints every branch type and checks the required ones)
// Synthetic estimator validation (no ROOT needed; NOT NA60+ data):
//   v0_archive_stats --synthetic N_MUONS [--rho-side B] [--seed S] [--out-dir DIR]
#include "V0Core.hh"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#ifdef V0_ROOT
#include <TBranch.h>
#include <TBranchElement.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>
#include <memory>
#endif
using namespace V0;
using Args = std::map<std::string, std::string>;
static Args parseArgs(int c, char** v) { Args m; for (int i = 1; i < c; ++i) { std::string a = v[i]; if (a.rfind("--", 0)) continue; a = a.substr(2); if (i + 1 < c && std::string(v[i + 1]).rfind("--", 0)) m[a] = v[++i]; else m[a] = "1"; } return m; }
static std::vector<std::string> splitCsv(const std::string& s) { std::vector<std::string> o; std::string t; std::stringstream ss(s); while (std::getline(ss, t, ',')) o.push_back(t); if (!s.empty() && s.back() == ',') o.push_back(""); return o; }
static double numOrNaN(const std::string& s) { if (s.empty()) return kNaN; try { return std::stod(s); } catch (...) { return kNaN; } }
static int slotOf(const std::string& p) { auto k = p.rfind("muon_"); if (k == std::string::npos) return -1; size_t j = k + 5; std::string d; while (j < p.size() && isdigit((unsigned char)p[j])) d += p[j++]; return d.empty() ? -1 : std::stoi(d); }
static std::string dirOf(const std::string& p) { auto k = p.find_last_of('/'); return k == std::string::npos ? std::string(".") : p.substr(0, k); }

// ---------------- microscopic_wire_seed.csv (exact census, written by the frozen Garfield source) ----------------
static bool readSeedCsv(const std::string& path, Muon& m, std::string& why) {
  std::ifstream f(path); if (!f) { why = "absent"; return false; }
  std::string line; if (!std::getline(f, line)) { why = "empty"; return false; }
  auto h = splitCsv(line); std::map<std::string, size_t> H; for (size_t i = 0; i < h.size(); ++i) H[h[i]] = i;
  for (const char* c : {"crossing_id", "primary_index", "wire_index_global", "collected_e", "ion_births", "A_ion_hemi", "mean_t_birth_ns", "sigma_t_birth_ns", "sigma_w_birth_mm", "mean_w_birth_global_mm"})
    if (!H.count(c)) { why = std::string("missing column ") + c; return false; }
  while (std::getline(f, line)) { if (line.empty()) continue; auto r = splitCsv(line); if (r.size() < h.size()) { why = "short row"; return false; }
    if (!m.id.empty() && r[H["crossing_id"]] != m.id) { why = "crossing_id differs from ROOT muons tree"; return false; }
    SeedRow s; s.prim = std::stoi(r[H["primary_index"]]); s.gwire = std::stoi(r[H["wire_index_global"]]); s.collected = std::stoll(r[H["collected_e"]]); s.births = std::stoll(r[H["ion_births"]]);
    s.A = numOrNaN(r[H["A_ion_hemi"]]); s.meanT = numOrNaN(r[H["mean_t_birth_ns"]]); s.sigT = numOrNaN(r[H["sigma_t_birth_ns"]]); s.sigWmm = numOrNaN(r[H["sigma_w_birth_mm"]]); s.meanWmm = numOrNaN(r[H["mean_w_birth_global_mm"]]);
    if (H.count("mean_t_collected_ns")) { s.meanTcol = numOrNaN(r[H["mean_t_collected_ns"]]); s.sigTcol = numOrNaN(r[H["sigma_t_collected_ns"]]); }
    m.seed.push_back(s); }
  long long a = 0, b = 0; for (auto& s : m.seed) a += s.births; for (auto& p : m.prims) b += p.Nexact;
  if (a != b) { why = "sum of exact per-wire births != sum of primaries.ion_births_exact"; m.seed.clear(); return false; }
  long long ca = 0, cb = 0; for (auto& s : m.seed) ca += s.collected; for (auto& p : m.prims) cb += p.collected;
  if (ca != cb) { why = "sum of per-wire collected_e != sum of primaries.collected_electrons"; m.seed.clear(); return false; }
  return true;
}

// ---------------- microscopic_seed.csv (per Heed primary, written by the same frozen source) ----------------
// columns used: primary_index, ion_births, collected_e, attached, other. attached/other are avalanche ELECTRON losses.
static bool readPrimaryCsv(const std::string& path, Muon& m, std::string& why) {
  std::ifstream f(path); if (!f) { why = "absent"; return false; }
  std::string line; if (!std::getline(f, line)) { why = "empty"; return false; }
  auto h = splitCsv(line); std::map<std::string, size_t> H; for (size_t i = 0; i < h.size(); ++i) H[h[i]] = i;
  for (const char* c : {"crossing_id", "primary_index", "ion_births", "collected_e", "attached", "other"}) if (!H.count(c)) { why = std::string("missing column ") + c; return false; }
  size_t n = 0;
  while (std::getline(f, line)) { if (line.empty()) continue; auto r = splitCsv(line); if (r.size() < h.size()) { why = "short row"; return false; }
    if (r[H["crossing_id"]] != m.id) { why = "crossing_id differs"; return false; }
    const int ip = std::stoi(r[H["primary_index"]]); if (ip < 0 || ip >= (int)m.prims.size()) { why = "primary_index out of range"; return false; }
    Prim& p = m.prims[ip]; if (std::stoll(r[H["ion_births"]]) != p.Nexact || std::stoll(r[H["collected_e"]]) != p.collected) { why = "ion_births/collected_e differ from ROOT primaries"; return false; }
    p.attached = std::stoll(r[H["attached"]]); p.otherLost = std::stoll(r[H["other"]]); ++n; }
  if (n != m.prims.size()) { why = "row count differs from ROOT primaries"; return false; }
  return true;
}

#ifdef V0_ROOT
// Branch types are those written by CathodeIonArchive.hh (identical sha256 on the overnight, independent-muon and
// v0-pilot branches): /I -> Int_t, /L -> Long64_t, /D -> Double_t, std::string objects for crossing_id/source_channel/storage_mode.
struct BranchSpec { const char* tree; const char* name; const char* type; };
static const BranchSpec kRequired[] = {
  {"muons", "crossing_id", "string"}, {"muons", "source_channel", "string"}, {"muons", "storage_mode", "string"},
  {"muons", "station", "Int_t"}, {"muons", "det_id", "Int_t"}, {"muons", "base_global_wire", "Int_t"}, {"muons", "reservoir_per_primary", "Int_t"},
  {"muons", "pitch_cm", "Double_t"}, {"muons", "theta_u_deg", "Double_t"}, {"muons", "theta_w_deg", "Double_t"}, {"muons", "w_in_mm", "Double_t"},
  {"primaries", "primary_index", "Int_t"}, {"primaries", "cluster_index", "Int_t"}, {"primaries", "successful", "Int_t"},
  {"primaries", "u_patch_cm", "Double_t"}, {"primaries", "v_absolute_cm", "Double_t"}, {"primaries", "w_relative_cm", "Double_t"}, {"primaries", "t_primary_ns", "Double_t"},
  {"primaries", "ion_births_exact", "Long64_t"}, {"primaries", "ion_births_saved", "Long64_t"}, {"primaries", "collected_electrons", "Long64_t"},
  {"ions", "primary_index", "Int_t"}, {"ions", "local_wire_index", "Int_t"}, {"ions", "global_wire_index", "Int_t"},
  {"ions", "u_patch_cm", "Double_t"}, {"ions", "v_absolute_cm", "Double_t"}, {"ions", "w_relative_cm", "Double_t"}, {"ions", "t_birth_ns", "Double_t"},
  {"ions", "inverse_inclusion_weight", "Double_t"}};
static std::string leafType(TTree* t, const char* name) {
  TBranch* b = t->GetBranch(name); if (!b) return "MISSING";
  if (auto* be = dynamic_cast<TBranchElement*>(b)) return be->GetClassName() == std::string("string") ? "string" : be->GetClassName();
  auto* l = dynamic_cast<TLeaf*>(b->GetListOfLeaves()->At(0)); return l ? l->GetTypeName() : "UNKNOWN";
}
static void checkSchema(TFile& f, const std::string& path) {
  for (const auto& r : kRequired) { auto* t = dynamic_cast<TTree*>(f.Get(r.tree)); if (!t) throw std::runtime_error(std::string("missing tree ") + r.tree + " in " + path);
    const std::string got = leafType(t, r.name); if (got != r.type) throw std::runtime_error(std::string("schema mismatch ") + r.tree + "." + r.name + ": expected " + r.type + ", found " + got + " in " + path); }
}
static void dumpSchema(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ")); if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  for (const char* tn : {"muons", "clusters", "primaries", "ions", "electron_endpoints", "electron_control"}) {
    auto* t = dynamic_cast<TTree*>(f->Get(tn)); if (!t) { std::cout << tn << ": ABSENT\n"; continue; }
    std::cout << tn << "  entries=" << t->GetEntries() << "\n"; for (auto* o : *t->GetListOfBranches()) std::cout << "   " << o->GetName() << " : " << leafType(t, o->GetName()) << "\n"; }
  bool ok = true; for (const auto& r : kRequired) { auto* t = dynamic_cast<TTree*>(f->Get(r.tree)); const std::string got = t ? leafType(t, r.name) : "NO TREE";
    if (got != r.type) { ok = false; std::cout << "MISMATCH " << r.tree << "." << r.name << " expected " << r.type << " found " << got << "\n"; } }
  std::cout << (ok ? "PASS" : "FAIL") << " V0 required-branch schema check (" << sizeof(kRequired) / sizeof(kRequired[0]) << " branches)\n";
}
static void bindAddr(TTree* t, const char* n, void* addr) { if (t->SetBranchAddress(n, addr) < 0) throw std::runtime_error(std::string("SetBranchAddress failed: ") + t->GetName() + "." + n); }
template <class T> static void bindT(TTree* t, const char* n, T* addr) { if (t->SetBranchAddress(n, addr) < 0) throw std::runtime_error(std::string("SetBranchAddress failed: ") + t->GetName() + "." + n); }
static Muon readMuon(const std::string& path) {
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ")); if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
  checkSchema(*f, path);
  auto *mu = dynamic_cast<TTree*>(f->Get("muons")), *pr = dynamic_cast<TTree*>(f->Get("primaries")), *io = dynamic_cast<TTree*>(f->Get("ions"));
  if (mu->GetEntries() != 1) throw std::runtime_error("muons tree must have exactly one entry in " + path);
  Muon m; m.path = path; m.slot = slotOf(path);
  std::string *id = nullptr, *ch = nullptr, *mode = nullptr; bindT(mu, "crossing_id", &id); bindT(mu, "source_channel", &ch); bindT(mu, "storage_mode", &mode);
  Int_t st = 0, det = 0, base = 0, K = 0; Double_t pitch = 0, thU = 0, thW = 0, wIn = 0;
  bindAddr(mu, "station", &st); bindAddr(mu, "det_id", &det); bindAddr(mu, "base_global_wire", &base); bindAddr(mu, "reservoir_per_primary", &K);
  bindAddr(mu, "pitch_cm", &pitch); bindAddr(mu, "theta_u_deg", &thU); bindAddr(mu, "theta_w_deg", &thW); bindAddr(mu, "w_in_mm", &wIn);
  mu->GetEntry(0); m.id = id ? *id : ""; m.channel = ch ? *ch : ""; m.station = st; m.det = det; m.base = base; m.K = K; m.pitch = pitch; m.thetaU = thU; m.thetaW = thW; m.wInMm = wIn;
  if (!mode || *mode != "weighted") throw std::runtime_error("storage_mode is not 'weighted' in " + path);
  if (std::fabs(m.pitch - 0.3) > 1e-12) throw std::runtime_error("pitch_cm != 0.3 (not the frozen 2100 V reference) in " + path);
  Int_t pi = 0, pc = 0, ok = 0; Double_t pu = 0, pv = 0, pw = 0, pt = 0; Long64_t N = 0, n = 0, col = 0;
  bindAddr(pr, "primary_index", &pi); bindAddr(pr, "cluster_index", &pc); bindAddr(pr, "successful", &ok);
  bindAddr(pr, "u_patch_cm", &pu); bindAddr(pr, "v_absolute_cm", &pv); bindAddr(pr, "w_relative_cm", &pw); bindAddr(pr, "t_primary_ns", &pt);
  bindAddr(pr, "ion_births_exact", &N); bindAddr(pr, "ion_births_saved", &n); bindAddr(pr, "collected_electrons", &col);
  for (Long64_t j = 0; j < pr->GetEntries(); ++j) { pr->GetEntry(j); if (pi != (Int_t)j) throw std::runtime_error("primary_index != row number in " + path);
    Prim p; p.idx = pi; p.cluster = pc; p.success = ok; p.u = pu; p.v = pv; p.w = pw; p.t = pt; p.Nexact = N; p.nSaved = n; p.collected = col; m.prims.push_back(p); }
  Int_t ip = 0, lw = 0, gw = 0; Double_t u = 0, v = 0, w = 0, t = 0, wt = 0;
  bindAddr(io, "primary_index", &ip); bindAddr(io, "local_wire_index", &lw); bindAddr(io, "global_wire_index", &gw);
  bindAddr(io, "u_patch_cm", &u); bindAddr(io, "v_absolute_cm", &v); bindAddr(io, "w_relative_cm", &w); bindAddr(io, "t_birth_ns", &t); bindAddr(io, "inverse_inclusion_weight", &wt);
  m.ions.reserve(io->GetEntries()); std::vector<long long> rows(m.prims.size(), 0);
  for (Long64_t j = 0; j < io->GetEntries(); ++j) { io->GetEntry(j);
    if (ip < 0 || ip >= (Int_t)m.prims.size()) throw std::runtime_error("ion row with invalid primary_index in " + path);
    if (lw != (Int_t)std::llround(u / m.pitch) || gw != m.base + lw) throw std::runtime_error("ion wire index is not the nearest anode to the birth in " + path);
    const Prim& p = m.prims[ip]; const double expect = double(p.Nexact) / double(p.nSaved);
    if (!(wt > 0) || std::fabs(wt - expect) > 1e-8 * std::max(1., expect)) throw std::runtime_error("inverse_inclusion_weight != ion_births_exact/ion_births_saved in " + path);
    ++rows[ip]; m.ions.push_back({ip, lw, gw, u, v, w, t, wt}); }
  for (size_t i = 0; i < rows.size(); ++i) if (rows[i] != m.prims[i].nSaved) throw std::runtime_error("ion rows per primary != ion_births_saved in " + path);
  return m;
}
#endif

// ---------------- synthetic archives with KNOWN realized and law-level truth ----------------
// Each primary's N ions are REALIZED in full (finite-ion physical fluctuation), the exact per-(primary,wire) census is
// computed from the realized cloud (as Garfield writes microscopic_wire_seed.csv), and only then a K-reservoir is drawn.
static std::vector<Muon> synthesize(int nMu, double sideBias, unsigned long long seed, int K) {
  std::mt19937_64 g(seed), gRes(seed ^ 0x5bd1e995ULL); std::uniform_real_distribution<double> U(0, 1); std::normal_distribution<double> Nn(0, 1);
  std::vector<Muon> out; const double pitch = 0.3;
  for (int im = 0; im < nMu; ++im) {
    Muon m; m.path = "synthetic"; m.slot = 9000 + im; m.thetaU = 40 * (U(g) - 0.5); m.thetaW = 40 * (U(g) - 0.5); m.id = "synthetic|" + std::to_string(im); m.channel = "synthetic"; m.pitch = pitch; m.base = 100; m.K = K; m.haveSeed = true;
    const double z = U(g) < 0.5 ? -1 : 1; const int M = 8 + (int)(30 * U(g));
    for (int ip = 0; ip < M; ++ip) {
      Prim p; p.idx = ip; p.cluster = ip; p.v = (2 * U(g) - 1) * 0.29; p.u = (2 * U(g) - 1) * 0.15; p.t = 3 * U(g);
      const bool noIon = U(g) < 0.04; long long N = noIon ? 0 : std::max<long long>(1, (long long)std::gamma_distribution<double>(0.8, 3000. / 0.8)(g));
      p.Nexact = N; p.collected = (N == 0 || U(g) < 0.15) ? 0 : std::max<long long>(1, N / 12);
      // arrival hemisphere follows the start side (geometry), plus a muon-level latent shift controlled by sideBias
      const bool fromAbove = (p.v > 0) ? (U(g) < 0.85) : (U(g) < 0.15); const bool up = U(g) < (fromAbove ? 0.8 : 0.2) + sideBias * z;
      const double mu = (up ? 1 : -1) * kPi / 2 + 0.5 * Nn(g), kap = 3.5 / (1. + N / 3000.);
      const double D = std::gamma_distribution<double>(3.5, 7.5)(g) + 4. * (up ? 1 : -1), sigIn = 1.5 + 0.5 * U(g);
      const bool split = U(g) < 0.15; const double fSplit = split ? 0.1 + 0.4 * U(g) : 0.; const int other = U(g) < 0.5 ? -1 : 1;
      double Alaw = 0, Z = 0; for (int k = 0; k < 2000; ++k) { const double th = (k + 0.5) * 2 * kPi / 2000 - kPi, w = std::exp(kap * std::cos(th - mu)); Z += w; Alaw += (std::sin(th) > 0 ? 1 : -1) * w; }
      m.lawA[ip] = Alaw / Z;
      // realize the full cloud
      struct R { double u, v, w, t; int wire; }; std::vector<R> cloud; cloud.reserve(N);
      const double wC = 0.02 * Nn(g), wS = 0.004 * (1 + std::log1p(N / 1000.));
      for (long long k = 0; k < N; ++k) { double th; do { th = (2 * U(g) - 1) * kPi; } while (U(g) > std::exp(kap * (std::cos(th - mu) - 1)));
        const int wire = (split && U(g) < fSplit) ? other : 0; const double r = 0.0018 * std::exp(0.45 * Nn(g));
        cloud.push_back({wire * pitch + r * std::cos(th), r * std::sin(th), wC + wS * Nn(g), p.t + D + sigIn * Nn(g), wire}); }
      // exact census per wire (what microscopic_wire_seed.csv holds)
      std::map<int, std::array<double, 7>> ex;   // births, up, down, St, St2, Sw, Sw2
      for (auto& c : cloud) { auto& e = ex[c.wire]; e[0] += 1; if (c.v > 0) e[1] += 1; else if (c.v < 0) e[2] += 1; e[3] += c.t; e[4] += c.t * c.t; e[5] += c.w; e[6] += c.w * c.w; }
      p.attached = (long long)(N * 0.02 * U(g)); p.otherLost = 0;
      long long colLeft = p.collected;
      for (auto& kv : ex) { auto& e = kv.second; SeedRow s; s.prim = ip; s.gwire = m.base + kv.first; s.births = (long long)e[0];
        s.collected = (kv.first == 0) ? colLeft : 0; s.A = (e[1] + e[2]) > 0 ? (e[1] - e[2]) / (e[1] + e[2]) : kNaN;
        s.meanT = e[3] / e[0]; s.sigT = std::sqrt(std::max(0., e[4] / e[0] - s.meanT * s.meanT)); if (s.collected > 0) { s.meanTcol = p.t + D - 0.4; s.sigTcol = 0.6; } const double mw = e[5] / e[0]; s.meanWmm = 10 * mw; s.sigWmm = 10 * std::sqrt(std::max(0., e[6] / e[0] - mw * mw)); m.seed.push_back(s); }
      if (N == 0 && p.collected == 0) { /* no seed row, as in the Garfield writer */ }
      // K-reservoir per primary (algorithm R), independent RNG
      std::vector<size_t> res; for (size_t k = 0; k < cloud.size(); ++k) { if ((int)res.size() < K) res.push_back(k); else { std::uniform_int_distribution<size_t> d(0, k); size_t j = d(gRes); if ((int)j < K) res[j] = k; } }
      p.nSaved = (long long)res.size(); const double wt = res.empty() ? 0 : double(N) / res.size();
      for (size_t k : res) { auto& c = cloud[k]; m.ions.push_back({ip, c.wire, m.base + c.wire, c.u, c.v, c.w, c.t, wt}); }
      m.prims.push_back(p);
    }
    out.push_back(m);
  }
  return out;
}

static double mean(const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return v.empty() ? kNaN : s / v.size(); }
static double rms(const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return v.empty() ? kNaN : std::sqrt(s / v.size()); }

int main(int argc, char** argv) {
  try {
  Args A = parseArgs(argc, argv); const std::string od = A.count("out-dir") ? A["out-dir"] : "v0_out";
  const bool synthetic = A.count("synthetic");
  std::vector<Muon> mus; std::set<int> excl; std::vector<std::string> skipped, seedProblems; bool calibrationOK = true;
  if (synthetic) mus = synthesize(std::stoi(A["synthetic"]), A.count("rho-side") ? std::stod(A["rho-side"]) : 0., A.count("seed") ? std::stoull(A["seed"]) : 12345ULL, 256);
  else {
#ifdef V0_ROOT
    if (A.count("dump-schema")) { dumpSchema(A["dump-schema"]); return 0; }
    if (A.count("exclude-index-csv")) { std::ifstream f(A["exclude-index-csv"]); if (!f) throw std::runtime_error("cannot read --exclude-index-csv"); std::string line; std::getline(f, line); auto h = splitCsv(line); int col = -1; for (size_t i = 0; i < h.size(); ++i) if (h[i] == "index") col = (int)i; if (col < 0) throw std::runtime_error("no 'index' column in exclusion CSV");
      while (std::getline(f, line)) { if (line.empty()) continue; auto r = splitCsv(line); if ((int)r.size() > col) excl.insert(std::stoi(r[col])); } }
    if (A.count("exclude-slots")) for (auto& s : splitCsv(A["exclude-slots"])) if (!s.empty()) excl.insert(std::stoi(s));
    if (excl.empty()) { if (!A.count("no-holdout-exclusion")) throw std::runtime_error("refusing to run: give --exclude-index-csv sealed_validation20.csv (or --exclude-slots), or explicitly --no-holdout-exclusion"); calibrationOK = false; }
    std::ifstream fl(A.at("list")); if (!fl) throw std::runtime_error("cannot read --list"); std::string line; const size_t mx = A.count("max-muons") ? std::stoul(A["max-muons"]) : (size_t)1 << 30;
    while (std::getline(fl, line) && mus.size() < mx) { if (line.empty()) continue; const int s = slotOf(line); if (excl.count(s)) { skipped.push_back(line); continue; }
      Muon m = readMuon(line); std::string why; m.haveSeed = readSeedCsv(dirOf(line) + "/microscopic_wire_seed.csv", m, why); if (!m.haveSeed) seedProblems.push_back(line + " : wire census: " + why);
      std::string why2; if (!readPrimaryCsv(dirOf(line) + "/microscopic_seed.csv", m, why2)) seedProblems.push_back(line + " : primary census: " + why2); mus.push_back(std::move(m)); }
#else
    std::cerr << "built without ROOT: use --synthetic, or rebuild with -DV0_ROOT\n"; return 2;
#endif
  }
  if (mus.empty()) throw std::runtime_error("no muons");
  if (std::system(("mkdir -p '" + od + "'").c_str())) return 2;
  if (synthetic) {   // round-trip test of the census-CSV reader using the exact column layout written by the Garfield source
    const std::string p = od + "/selftest_microscopic_wire_seed.csv"; std::ofstream f(p); f << std::setprecision(15);
    f << "crossing_id,station,theta_u_deg,theta_w_deg,primary_index,cluster_index,wire_index_global,wire_u_global_mm,collected_e,ion_births,ion_up,ion_down,ion_boundary,ion_top,ion_bottom,ion_side,A_ion_hemi,A_ion_TB,mean_w_birth_global_mm,sigma_w_birth_mm,mean_t_birth_ns,sigma_t_birth_ns,mean_t_collected_ns,sigma_t_collected_ns\n";
    for (auto& s : mus[0].seed) { f << mus[0].id << ",0,0,0," << s.prim << "," << s.prim << "," << s.gwire << ",0," << s.collected << "," << s.births << ",0,0,0,0,0,0,";
      if (std::isnan(s.A)) f << ","; else f << s.A << ","; f << "," << s.meanWmm << "," << s.sigWmm << "," << s.meanT << "," << s.sigT << ",,\n"; }
    f.close(); Muon t; t.id = mus[0].id; t.prims = mus[0].prims; std::string why; const bool ok = readSeedCsv(p, t, why); bool same = ok && t.seed.size() == mus[0].seed.size();
    for (size_t i = 0; same && i < t.seed.size(); ++i) { const auto &a = t.seed[i], &b = mus[0].seed[i]; same = a.prim == b.prim && a.gwire == b.gwire && a.births == b.births && a.collected == b.collected && std::fabs(a.A - b.A) < 1e-12 && std::fabs(a.sigT - b.sigT) < 1e-9; }
    std::cout << (same ? "PASS" : "FAIL") << " census-CSV reader round trip (" << t.seed.size() << " rows" << (ok ? "" : ", " + why) << ")\n"; if (!same) return 3; }
  std::vector<Frag> F; std::vector<PrimRec> P; buildFragments(mus, F, P);
  std::ofstream S(od + "/v0_summary.txt"); auto out = [&](const std::ostringstream& s) { S << s.str(); std::cout << s.str(); };
  const std::string tag = synthetic ? "   [SYNTHETIC ARCHIVES: estimator validation only, NOT NA60+ data]" : (calibrationOK ? "" : "   [NO HOLDOUT EXCLUSION: NOT FOR CALIBRATION]");
  int nSeed = 0; for (auto& m : mus) nSeed += m.haveSeed;
  { std::ostringstream s; s << std::setprecision(5);
    s << "V0 archive statistics: " << mus.size() << " muons, " << P.size() << " Heed primaries, " << F.size() << " (primary, nearest-birth wire) fragments" << tag << "\n";
    s << "  exact per-wire census (microscopic_wire_seed.csv) available for " << nSeed << "/" << mus.size() << " muons; holdout slots excluded: " << skipped.size() << "\n";
    s << "  NOTE 'wire' below = anode nearest to the ION BIRTH (ions.local_wire_index); collecting wire is taken only from collected_e in the census CSV.\n";
    double sv = 0, sw = 0; for (auto& m : mus) for (auto& i : m.ions) { sv += i.wt * i.v; sw += i.wt; } s << "  ion-weighted mean v_absolute = " << sv / sw << " cm (anode plane at v=0)\n"; out(s); }

  // ===== [1] primary census: ion births vs collected electrons =====
  { std::ostringstream s; s << std::setprecision(5); long long n = 0, bad = 0, z00 = 0, z01 = 0, z10 = 0, z11 = 0, resBad = 0, htBad = 0; double yAll = 0, yZeroCol = 0;
    std::map<std::pair<int, int>, double> wsum; for (size_t im = 0; im < mus.size(); ++im) for (auto& i : mus[im].ions) wsum[{(int)im, i.prim}] += i.wt;
    for (size_t im = 0; im < mus.size(); ++im) for (auto& p : mus[im].prims) { ++n; if (p.success != 1) ++bad; const bool ion = p.Nexact > 0, col = p.collected > 0;
      (ion ? (col ? z11 : z10) : (col ? z01 : z00))++; yAll += p.Nexact; if (!col) yZeroCol += p.Nexact;
      if (p.nSaved != std::min<long long>(p.Nexact, mus[im].K)) ++resBad; const double w = wsum[{(int)im, p.idx}]; if (std::fabs(w - p.Nexact) > 1e-6 * std::max(1., (double)p.Nexact)) ++htBad; }
    s << "\n[1] Primary census (primaries tree, ALL Heed primaries incl. zero-ion ones)\n  primaries " << n << "  unsuccessful " << bad << "\n";
    s << "  ion_births_exact>0 & collected>0: " << z11 << "   ion>0 & collected=0: " << z10 << "   ion=0 & collected>0: " << z01 << "   ion=0 & collected=0: " << z00 << "\n";
    s << "  fraction collected=0: " << double(z10 + z00) / n << "   fraction ion_births_exact=0: " << double(z01 + z00) / n << "   share of total ion yield in collected=0 primaries: " << yZeroCol / yAll << "\n";
    s << "  reservoir checks: n_saved != min(N,K): " << resBad << "   sum(weights) != ion_births_exact: " << htBad << "\n"; out(s); }

  // ===== [2] nearest-birth wire vs collecting wire; reservoir coverage of wire fragments =====
  { std::ostringstream s; s << std::setprecision(5); long long np = 0, diff = 0, fr = 0, frUnsampled = 0; double yEx = 0, yUns = 0, yNoCol = 0; std::vector<double> rel;
    for (auto& r : P) if (r.collWireArgmax != -999 && r.birthWireArgmax != -999) { ++np; if (r.collWireArgmax != r.birthWireArgmax) ++diff; }
    for (auto& g : F) if (g.exBirths > 0) { ++fr; yEx += g.exBirths; if (g.nres == 0) { ++frUnsampled; yUns += g.exBirths; } if (g.exCollected == 0) yNoCol += g.exBirths; if (g.nres > 0) rel.push_back(g.Nhat / g.exBirths - 1); }
    std::vector<double> rw; double wTot = 0, wFar = 0; for (auto& m : mus) for (auto& i : m.ions) { const double r = std::hypot(i.u - i.lwire * m.pitch, i.v); wTot += i.wt; if (r > 0.05) wFar += i.wt; }
    s << "\n[2] Wires (needs census CSV for exact parts)\n  primaries with both argmax wires known: " << np << "; argmax(collected_e) != argmax(ion births): " << diff << " (" << (np ? double(diff) / np : kNaN) << ")\n";
    s << "  exact fragments with births>0: " << fr << "; with ZERO sampled reservoir ions: " << frUnsampled << " carrying " << (yEx > 0 ? yUns / yEx : kNaN) << " of exact yield\n";
    s << "  share of exact ion yield on a (primary, wire) whose collected_e = 0: " << (yEx > 0 ? yNoCol / yEx : kNaN) << "\n";
    s << "  per-fragment yield estimate N_hat/N_exact - 1 (sampled fragments): mean " << mean(rel) << "  RMS " << rms(rel) << "  (reservoir is per primary, so wire shares are estimated)\n";
    s << "  ion weight born > 0.5 mm from its nearest anode: " << wFar / wTot << "\n"; out(s); }

  // ===== [3] timing: nested decomposition of t_birth - t_primary =====
  // t_ij = T_i + dt_ij. Levels (ion-yield weighted): ions within a (primary, wire) fragment | fragments within a primary |
  // primaries within a muon | muons. Also: spread around a per-(muon, wire) mean, i.e. what per-wire time moments keep.
  { std::ostringstream s; s << std::setprecision(5);
    s << "\n[3] Birth time relative to the Heed primary time, x = t_birth - t_primary [ns]; ion-yield weighted RMS per level\n";
    for (int ex = 1; ex >= 0; --ex) {
      struct FT { int muon, prim, wire; double N, m, s2, n2; }; std::vector<FT> v;
      for (auto& g : F) { const double N = ex ? (double)g.exBirths : g.Nhat, m = ex ? g.exTbar : g.tbar, sg = ex ? g.exSigt : g.sigt;
        if (!(N > 0) || std::isnan(m) || std::isnan(sg)) continue; v.push_back({g.muon, g.prim, g.gwire, N, m, sg * sg, ex ? 0. : g.tbarSE * g.tbarSE}); }
      if (v.empty()) continue;
      double W = 0, S = 0; std::map<std::pair<int, int>, std::array<double, 2>> pm, wm; std::map<int, std::array<double, 2>> mm;
      for (auto& f : v) { W += f.N; S += f.N * f.m; auto& a1 = pm[{f.muon, f.prim}]; a1[0] += f.N; a1[1] += f.N * f.m; auto& a2 = wm[{f.muon, f.wire}]; a2[0] += f.N; a2[1] += f.N * f.m; auto& a3 = mm[f.muon]; a3[0] += f.N; a3[1] += f.N * f.m; }
      const double M = S / W; double wf = 0, bfp = 0, bpm = 0, bm = 0, aroundWire = 0, noise = 0;
      for (auto& f : v) { const auto& P = pm[{f.muon, f.prim}]; const auto& Wr = wm[{f.muon, f.wire}]; const double mp = P[1] / P[0], mw = Wr[1] / Wr[0];
        wf += f.N * f.s2; bfp += f.N * (f.m - mp) * (f.m - mp); aroundWire += f.N * (f.s2 + (f.m - mw) * (f.m - mw)); noise += f.N * f.n2; }
      for (auto& kv : pm) { const double mp = kv.second[1] / kv.second[0], mu = mm[kv.first.first][1] / mm[kv.first.first][0]; bpm += kv.second[0] * (mp - mu) * (mp - mu); }
      for (auto& kv : mm) { const double mu = kv.second[1] / kv.second[0]; bm += kv.second[0] * (mu - M) * (mu - M); }
      wf /= W; bfp /= W; bpm /= W; bm /= W; aroundWire /= W; noise /= W; const double tot = wf + bfp + bpm + bm;
      s << (ex ? "  EXACT census (" : "  K256 reservoir (") << v.size() << " fragments): mean " << M << "  total RMS " << std::sqrt(tot) << "\n"
        << "    within fragment (one avalanche, one wire) " << std::sqrt(wf) << "  | between fragments of one primary " << std::sqrt(bfp)
        << "  | between primaries of one muon " << std::sqrt(bpm) << "  | between muons " << std::sqrt(bm) << "\n"
        << "    variance shares: " << wf / tot << " / " << bfp / tot << " / " << bpm / tot << " / " << bm / tot
        << "   residual RMS with one time per primary " << std::sqrt(wf + bfp) << "   residual RMS around per-(muon,wire) mean " << std::sqrt(aroundWire) << "\n";
      if (!ex) s << "    (reservoir: between-level terms include sampling noise of fragment means, RMS " << std::sqrt(noise) << " ns)\n"; }
    std::vector<double> Nv; for (auto& g : F) if (g.exBirths > 0 && !std::isnan(g.exSigt)) Nv.push_back((double)g.exBirths);
    if (!Nv.empty()) { auto e = quantileEdges(Nv, 5); s << "  exact within-fragment sigma_t by exact-births quintile (median [p10, p90] ns):";
      for (int b = 0; b < 5; ++b) { std::vector<double> x; for (auto& g : F) if (g.exBirths > 0 && !std::isnan(g.exSigt) && binOf(e, (double)g.exBirths) == b) x.push_back(g.exSigt);
        s << "  " << quantile(x, .5) << " [" << quantile(x, .1) << ", " << quantile(x, .9) << "]"; } s << "\n"; }
    out(s); }

  // ===== [4] reservoir uncertainty measured against the exact census =====
  { std::ostringstream s; s << std::setprecision(5); std::vector<double> dA, zA, seA, dT, zT, dS;
    for (auto& g : F) { if (g.nres == 0 || g.exBirths <= 0) continue; if (!std::isnan(g.exA) && !std::isnan(g.x[0])) { dA.push_back(g.x[0] - g.exA); seA.push_back(g.se[0]); if (g.se[0] > 0) zA.push_back((g.x[0] - g.exA) / g.se[0]); }
      if (!std::isnan(g.exTbar)) { dT.push_back(g.tbar - g.exTbar); if (g.tbarSE > 0) zT.push_back((g.tbar - g.exTbar) / g.tbarSE); } if (!std::isnan(g.exSigt) && g.nres > 1) dS.push_back(g.sigt - g.exSigt); }
    s << "\n[4] Reservoir (K256) error measured against the EXACT census of the same avalanches (sampling error only; no physics)\n";
    if (dA.empty()) s << "  no census CSV: reservoir error cannot be measured directly\n";
    else { s << "  hemisphere A: fragments " << dA.size() << "  mean(est-exact) " << mean(dA) << "  RMS " << rms(dA) << "  predicted RMS from SE " << rms(seA) << "  RMS of z=(est-exact)/SE (SE>0) " << rms(zA) << " (n=" << zA.size() << ")\n";
      s << "  mean birth time: mean(est-exact) " << mean(dT) << " ns  RMS " << rms(dT) << " ns  RMS z " << rms(zT) << ";  sigma_t est-exact: mean " << mean(dS) << " RMS " << rms(dS) << " ns\n"; }
    if (synthetic) { std::vector<double> dPhys, pred; for (auto& r : P) { if (r.Nexact <= 0 || std::isnan(r.exA)) continue; const double L = mus[r.muon].lawA.at(r.prim); dPhys.push_back(r.exA - L); pred.push_back(std::sqrt(std::max(0., 1 - L * L) / r.Nexact)); }
      s << "  [synthetic only] PHYSICAL finite-ion fluctuation, realized cloud A minus generating-law A per primary: RMS " << rms(dPhys) << "  (iid-binomial expectation " << rms(pred) << ")\n"
        << "  -> the SE of the estimator describes reservoir sampling around the REALIZED avalanche; it does not and should not include the realized-vs-law term.\n"; }
    out(s); }

  // ===== [5] conditional descriptor distributions =====
  std::vector<double> NvF; for (auto& g : F) if (g.nres > 0) NvF.push_back(g.Nhat); const auto NeF = quantileEdges(NvF, 10);
  std::vector<double> du; for (auto& g : F) du.push_back(std::fabs(g.startDuMm)); const auto duE = quantileEdges(du, 3);
  std::vector<double> dl; for (auto& g : F) if (g.nres > 0) dl.push_back(g.tbar); const auto dlE = quantileEdges(dl, 5);
  { std::ostringstream s; s << std::setprecision(4);
    s << "\n[5] Conditional per-fragment descriptor distributions (K256 estimates; RMS has reservoir noise subtracted)\n      columns: n, then [mean, RMS] for A, <sin>, <cos>, <dy> um, <dx> um, then fraction |A|>0.9\n";
    auto table = [&](const std::string& name, int nb, auto key) { s << "  by " << name << "\n"; for (int b = 0; b < nb; ++b) { double sm[5] = {0}, sq[5] = {0}, ns[5] = {0}; int n = 0, sat = 0;
        for (auto& g : F) { if (g.nres == 0 || key(g) != b) continue; ++n; sat += std::fabs(g.x[0]) > 0.9; for (int k = 0; k < 5; ++k) { if (std::isnan(g.x[k])) continue; sm[k] += g.x[k]; sq[k] += g.x[k] * g.x[k]; ns[k] += g.se[k] * g.se[k]; } }
        if (!n) continue; s << "    " << std::setw(3) << b << ": " << std::setw(6) << n; for (int k = 0; k < 5; ++k) { const double mm = sm[k] / n; s << "  [" << mm << ", " << std::sqrt(std::max(0., (sq[k] - ns[k]) / n - mm * mm)) << "]"; } s << "  " << double(sat) / n << "\n"; } };
    table("fragment N_hat decile", 10, [&](const Frag& g) { return binOf(NeF, g.Nhat); });
    table("primary start: 0..2 = v<0 & |du| tercile, 3..5 = v>0 & |du| tercile", 6, [&](const Frag& g) { return (g.startVmm > 0 ? 3 : 0) + binOf(duE, std::fabs(g.startDuMm)); });
    table("delay (t_bar - t_primary) quintile", 5, [&](const Frag& g) { return binOf(dlE, g.tbar); });
    std::map<std::pair<int, int>, long long> colOf; for (auto& r : P) colOf[{r.muon, r.prim}] = r.collected;
    std::vector<double> cv; for (auto& r : P) if (r.collected > 0) cv.push_back((double)r.collected); const auto cE = quantileEdges(cv, 4);
    table("primary collected electrons: 0 = none collected, 1..4 = quartile of collected > 0", 5, [&](const Frag& g) { const long long c = colOf[{g.muon, g.prim}]; return c <= 0 ? 0 : 1 + binOf(cE, (double)c); });
    table("collecting-wire relation (census): 0 = this wire collected electrons of the primary, 1 = it collected none, 2 = census absent", 3,
          [&](const Frag& g) { return g.exCollected < 0 ? 2 : (g.exCollected > 0 ? 0 : 1); });
    std::vector<double> th; for (auto& m : mus) th.push_back(std::fabs(m.thetaU)); const auto thE = quantileEdges(th, 3);
    table("muon |theta_u| tercile (across-wire incidence)", 3, [&](const Frag& g) { return binOf(thE, std::fabs(mus[g.muon].thetaU)); });
    s << "  note: arrival azimuth of the primary electron is not stored in K256 archives (electron endpoints exist only for the first primary of FULL archives); start geometry is used instead.\n";
    out(s); }

  // ===== [6] rank correlations =====
  { std::vector<std::string> nm = {"A", "sin", "cos", "dy", "dx", "|A|", "logN", "logCol", "delay", "sigt", "rbar", "sigw", "startV", "|du|", "thU", "thW"}; std::vector<std::vector<double>> c(nm.size()); int n = 0;
    for (auto& g : F) { if (g.nres < 30 || std::isnan(g.x[0])) continue; ++n; const Muon& m = mus[g.muon];
      double v[] = {g.x[0], g.x[1], g.x[2], g.x[3], g.x[4], std::fabs(g.x[0]), std::log10(g.Nhat), std::log10(1. + (double)m.prims[g.prim].collected), g.tbar, g.sigt, g.rbar, g.sigw, g.startVmm, std::fabs(g.startDuMm), m.thetaU, m.thetaW}; for (size_t k = 0; k < nm.size(); ++k) c[k].push_back(v[k]); }
    std::ostringstream s; s << "\n[6] Spearman rank correlations over fragments with >=30 reservoir ions (n=" << n << ")\n       "; for (auto& x : nm) s << std::setw(7) << x; s << "\n";
    for (size_t i = 0; i < nm.size(); ++i) { s << std::setw(7) << nm[i]; for (size_t j = 0; j < nm.size(); ++j) s << std::setw(7) << std::fixed << std::setprecision(2) << spearman(c[i], c[j]); s << "\n"; } s << std::defaultfloat; out(s); }

  // ===== [7] effective-number / independence tests =====
  { std::ostringstream s; s << std::setprecision(4); const std::vector<double> neffE = {2, 3, 5, 8, 12};
        std::vector<double> duP; for (auto& r : P) duP.push_back(std::fabs(r.startDuMm)); const auto duEP = quantileEdges(duP, 3);
        auto cls = [&](double v, double d) { return (v > 0 ? 3 : 0) + binOf(duEP, std::fabs(d)); };
    s << "\n[7] Independent-avalanche hypothesis. obs/pred = observed Var / prediction for independent members; rho = fitted equicorrelation (68% bootstrap over muons)\n";
    s << "    muon level: members are PRIMARIES (all wire fragments of a primary aggregated first). wire level: members are (primary, wire) fragments of one muon and wire.\n";
    auto line = [&](const std::string& lab, const IndepResult& r, bool showNeff) { s << "  " << std::left << std::setw(46) << lab << std::right << " units " << std::setw(5) << r.units << "  obs/pred " << std::setw(7) << (r.pred > 0 ? r.obs / r.pred : kNaN)
      << "  rho " << std::setw(8) << r.rho << " [" << r.rhoLo << ", " << r.rhoHi << "]\n"; if (showNeff) { s << "      obs/pred by N_eff [<2,2-3,3-5,5-8,8-12,>12]:"; for (size_t b = 0; b < r.obsByNeff.size(); ++b) s << ' ' << (r.predByNeff[b] > 0 ? r.obsByNeff[b] / r.predByNeff[b] : kNaN) << '(' << r.nByNeff[b] << ')'; s << "\n"; } };
    for (int k = 0; k < 5; ++k) {
      for (int cond = 0; cond < 2; ++cond) { std::map<int, Unit> U; for (auto& r : P) { if (!(r.Nhat > 0) || std::isnan(r.x[k])) continue; auto& u = U[r.muon]; u.key = r.muon;
          u.m.push_back({r.Nhat, r.x[k], r.se[k], cond ? cls(r.startVmm, r.startDuMm) : 0, r.Nhat}); }
        std::vector<Unit> V; for (auto& kv : U) if (kv.second.m.size() >= 3) V.push_back(kv.second);
        line(std::string(kDesc[k]) + (cond ? " muon | N, start(v-side,|du|)" : " muon | N"), indepTest(V, neffE), k == 0); }
      std::map<std::pair<int, int>, Unit> W; for (auto& g : F) { if (g.nres == 0 || std::isnan(g.x[k])) continue; auto& u = W[{g.muon, g.gwire}]; u.key = g.muon; u.m.push_back({g.Nhat, g.x[k], g.se[k], 0, g.NprimExact > 0 ? g.NprimExact : g.Nhat}); }
      std::vector<Unit> V; for (auto& kv : W) if (kv.second.m.size() >= 2) V.push_back(kv.second); line(std::string(kDesc[k]) + " per wire | N_primary", indepTest(V, neffE), false); }
    { std::map<int, Unit> U; for (auto& r : P) { if (r.Nexact <= 0 || std::isnan(r.exA)) continue; auto& u = U[r.muon]; u.key = r.muon; u.m.push_back({(double)r.Nexact, r.exA, 0., 0, (double)r.Nexact}); }
      std::vector<Unit> V; for (auto& kv : U) if (kv.second.m.size() >= 3) V.push_back(kv.second); if (!V.empty()) line("A EXACT census, muon | N", indepTest(V, neffE), true); }
    // slide-style summary
    std::vector<double> aP, aM; std::map<int, std::pair<double, double>> am; for (auto& r : P) if (r.Nhat > 0 && !std::isnan(r.x[0])) { aP.push_back(r.x[0]); am[r.muon].first += r.Nhat * r.x[0]; am[r.muon].second += r.Nhat; }
    for (auto& kv : am) aM.push_back(kv.second.first / kv.second.second);
    s << "  RMS of hemisphere A: per primary " << rms(aP) << "  per muon (yield-weighted) " << rms(aM) << "   (1800 V slide: 0.87 / 0.38; different HV and campaign, not comparable as calibration)\n"; out(s); }

  // ===== [8] 2100 V transport and timing map by primary start position (voltage-consistency audit) =====
  // Cells: |du| (start offset from nearest anode) in 6 bins of 0.25 mm x v_start in 12 bins of 0.5 mm. 2100 V quantities
  // come only from the K256 archives and their exact census; the optional 1800 V LUT500 is binned identically and
  // reported in separate columns. Observables differ: 2100 V "onset" = mean endpoint time of the avalanche electrons
  // collected (AvalancheElectron), 1800 V "drift" = endpoint time of one DriftElectron. They are never merged.
  {
    const int NU = 6, NV = 12; auto cellOf = [&](double duMm, double vMm) { const int a = std::min(NU - 1, (int)(std::fabs(duMm) / 0.25)); const int b = std::max(0, std::min(NV - 1, (int)std::floor((vMm + 3.) / 0.5))); return a * NV + b; };
    struct C2 { double n = 0, zc = 0, zi = 0, inc = 0, multi = 0, same = 0, nOn = 0, on = 0, on2 = 0, nIo = 0, io = 0, io2 = 0, nOff = 0, off = 0, off2 = 0, nSig = 0, sig = 0, nDw = 0, dw = 0, dw2 = 0, nAtt = 0, att = 0, nDir = 0, sn = 0, cs = 0; };
    struct C1 { double n = 0, nc = 0, nT = 0, t = 0, t2 = 0, nDw = 0, dw = 0, dw2 = 0, same = 0, sn = 0, cs = 0; };
    std::vector<C2> c2(NU * NV); std::vector<C1> c1(NU * NV); bool haveLut = false; std::string lutNote;
    struct PT { double on, io, sig, dw; int cell, thBin, nBin, mu; }; std::vector<PT> pts;
    std::map<std::pair<int, int>, std::vector<const Frag*>> fr; for (auto& g : F) fr[{g.muon, g.prim}].push_back(&g);
    std::vector<double> thv; for (auto& m : mus) thv.push_back(std::fabs(m.thetaU)); const auto thE = quantileEdges(thv, 3);
    std::vector<double> nv; for (auto& r : P) if (r.Nexact > 0) nv.push_back((double)r.Nexact); const auto nE = quantileEdges(nv, 5);
    for (auto& r : P) {
      const Muon& m = mus[r.muon]; const Prim& p = m.prims[r.prim]; const int cell = cellOf(r.startDuMm, r.startVmm); C2& c = c2[cell]; c.n += 1;
      c.zc += p.collected == 0; c.zi += p.Nexact == 0; c.inc += (p.Nexact > 0 && p.collected == 0);
      if (p.attached >= 0) { const double tot = double(p.attached + p.otherLost + p.collected); if (tot > 0) { c.nAtt += 1; c.att += p.attached / tot; } }
      if (!m.haveSeed) continue;
      double B = 0, Bt = 0, Bt2 = 0, Bw = 0, Cc = 0, Ct = 0, best = 0, second = 0; int bestW = -999;
      for (auto* g : fr[{r.muon, r.prim}]) { if (g->exBirths > 0 && !std::isnan(g->exTbar) && !std::isnan(g->exSigt)) { B += g->exBirths; Bt += g->exBirths * g->exTbar; Bt2 += g->exBirths * (g->exSigt * g->exSigt + g->exTbar * g->exTbar); if (!std::isnan(g->exMeanWmm)) Bw += g->exBirths * g->exMeanWmm; }
        if (g->exCollected > 0) { Cc += g->exCollected; if (!std::isnan(g->exTcol)) Ct += g->exCollected * g->exTcol; if (g->exCollected > best) { second = best; best = g->exCollected; bestW = g->gwire; } else if (g->exCollected > second) second = g->exCollected; } }
      const int startWire = m.base + (int)std::llround(p.u / m.pitch);
      if (Cc > 0) { const double on = Ct / Cc; c.nOn += 1; c.on += on; c.on2 += on * on; c.same += bestW == startWire; c.multi += second > 0.05 * Cc; }
      double io = kNaN, sig = kNaN, dw = kNaN;
      if (B > 0) { io = Bt / B; sig = std::sqrt(std::max(0., Bt2 / B - io * io)); dw = Bw / B - (m.wInMm + 10. * p.w); c.nIo += 1; c.io += io; c.io2 += io * io; c.nSig += 1; c.sig += sig; c.nDw += 1; c.dw += dw; c.dw2 += dw * dw;
        if (Cc > 0) { const double off = io - Ct / Cc; c.nOff += 1; c.off += off; c.off2 += off * off; } }
      if (!std::isnan(r.x[1]) && !std::isnan(r.x[2])) { c.nDir += 1; c.sn += r.x[1]; c.cs += r.x[2]; }
      pts.push_back({Cc > 0 ? Ct / Cc : kNaN, io, sig, dw, cell, binOf(thE, std::fabs(m.thetaU)), p.Nexact > 0 ? binOf(nE, (double)p.Nexact) : -1, r.muon});
    }
    if (A.count("lut500")) {   // 1800 V LUT500 raw garfield_electrons.csv (DriftElectron, attachment on). Never merged with 2100 V numbers.
      for (auto& path : splitCsv(A["lut500"])) { std::ifstream f(path); if (!f) throw std::runtime_error("cannot read LUT500 file " + path);
        std::string line; std::getline(f, line); auto h = splitCsv(line); std::map<std::string, size_t> H; for (size_t i = 0; i < h.size(); ++i) H[h[i]] = i;
        for (const char* col : {"u_start_mm", "v_start_mm", "w_start_mm", "collected", "wire_index", "w_end_mm", "drift_time_ns", "phi_hit_deg"}) if (!H.count(col)) throw std::runtime_error(std::string("LUT500 file lacks column ") + col);
        long long rows = 0;
        while (std::getline(f, line)) { if (line.empty()) continue; auto r = splitCsv(line); if (r.size() < h.size()) continue; ++rows;
          const double u = std::stod(r[H["u_start_mm"]]), v = std::stod(r[H["v_start_mm"]]); const double du = u - 3. * std::round(u / 3.); C1& c = c1[cellOf(du, v)]; c.n += 1;
          const int col = (int)numOrNaN(r[H["collected"]]); if (col != 1) continue; c.nc += 1; const double t = numOrNaN(r[H["drift_time_ns"]]); if (!std::isnan(t)) { c.nT += 1; c.t += t; c.t2 += t * t; }
          const double dw = numOrNaN(r[H["w_end_mm"]]) - numOrNaN(r[H["w_start_mm"]]); if (!std::isnan(dw)) { c.nDw += 1; c.dw += dw; c.dw2 += dw * dw; }
          c.same += (int)numOrNaN(r[H["wire_index"]]) == (int)std::llround(u / 3.); const double ph = numOrNaN(r[H["phi_hit_deg"]]) * kPi / 180.; if (!std::isnan(ph)) { c.sn += std::sin(ph); c.cs += std::cos(ph); } }
        lutNote += path + " rows=" + std::to_string(rows) + "; "; }
      haveLut = true; }
    auto mr = [](double s, double s2, double n, double& m, double& r) { m = n > 0 ? s / n : kNaN; r = n > 1 ? std::sqrt(std::max(0., s2 / n - m * m)) : kNaN; };
    std::ostringstream s; s << std::setprecision(4);
    s << "\n[8] 2100 V transport/timing by primary START position (voltage-consistency audit). Delays relative to t_primary [ns].\n";
    s << "    onset = mean endpoint time of collected avalanche electrons (exact census); ion = mean ion-birth time; sigma = within-avalanche ion-birth RMS; dw = mean ion w minus start w [mm]\n";
    { C2 t; for (auto& c : c2) { t.n += c.n; t.zc += c.zc; t.zi += c.zi; t.inc += c.inc; t.multi += c.multi; t.same += c.same; t.nOn += c.nOn; t.on += c.on; t.on2 += c.on2; t.nIo += c.nIo; t.io += c.io; t.io2 += c.io2; t.nOff += c.nOff; t.off += c.off; t.off2 += c.off2; t.nSig += c.nSig; t.sig += c.sig; t.nDw += c.nDw; t.dw += c.dw; t.dw2 += c.dw2; t.nAtt += c.nAtt; t.att += c.att; }
      double om, orr, im, ir, fm, fr2, dm, dr; mr(t.on, t.on2, t.nOn, om, orr); mr(t.io, t.io2, t.nIo, im, ir); mr(t.off, t.off2, t.nOff, fm, fr2); mr(t.dw, t.dw2, t.nDw, dm, dr);
      std::vector<double> onv, sgv; for (auto& q : pts) { if (!std::isnan(q.on)) onv.push_back(q.on); if (!std::isnan(q.sig)) sgv.push_back(q.sig); }
      s << "  all primaries " << t.n << ": P(collected=0) " << t.zc / t.n << "  P(ion=0) " << t.zi / t.n << "  P(ion>0 & collected=0) " << t.inc / t.n
        << "  mean attached share of avalanche electrons " << (t.nAtt > 0 ? t.att / t.nAtt : kNaN) << "\n";
      s << "  onset delay: mean " << om << "  RMS " << orr << "  p10/p50/p90 " << quantile(onv, .1) << "/" << quantile(onv, .5) << "/" << quantile(onv, .9)
        << "   ion-birth delay: mean " << im << " RMS " << ir << "   ion minus onset: mean " << fm << " RMS " << fr2 << "\n";
      s << "  within-avalanche sigma_t: p10/p50/p90 " << quantile(sgv, .1) << "/" << quantile(sgv, .5) << "/" << quantile(sgv, .9) << "   dw: mean " << dm << " RMS " << dr
        << "   main collecting wire = nearest-to-start wire: " << (t.nOn > 0 ? t.same / t.nOn : kNaN) << "   >5% on a second wire: " << (t.nOn > 0 ? t.multi / t.nOn : kNaN) << "\n"; }
    auto condTable = [&](const char* name, int nb, auto key) { s << "  onset/ion/sigma/dw by " << name << " (n, onset mean RMS, ion mean, sigma median, dw RMS)\n";
      for (int b = 0; b < nb; ++b) { std::vector<double> on, io, sg, dw; for (auto& q : pts) { if (key(q) != b) continue; if (!std::isnan(q.on)) on.push_back(q.on); if (!std::isnan(q.io)) io.push_back(q.io); if (!std::isnan(q.sig)) sg.push_back(q.sig); if (!std::isnan(q.dw)) dw.push_back(q.dw); }
        if (on.empty() && io.empty()) continue; double m1 = mean(on), r1 = 0; for (double x : on) r1 += (x - m1) * (x - m1); r1 = on.size() > 1 ? std::sqrt(r1 / on.size()) : kNaN;
        double md = mean(dw), rd = 0; for (double x : dw) rd += (x - md) * (x - md); rd = dw.size() > 1 ? std::sqrt(rd / dw.size()) : kNaN;
        s << "    " << std::setw(2) << b << ": " << std::setw(6) << on.size() << "  " << m1 << " " << r1 << "  " << mean(io) << "  " << quantile(sg, .5) << "  " << rd << "\n"; } };
    condTable("|theta_u| tercile", 3, [](const PT& q) { return q.thBin; });
    condTable("ion-yield quintile", 5, [](const PT& q) { return q.nBin; });
    s << "  onset delay mean [ns] by start cell: rows v_start (-3..3 mm, 0.5 mm), columns |du| (0..1.5 mm, 0.25 mm)" << (haveLut ? "; second block = 1800 V LUT500 drift time, third = 2100 minus 1800" : "") << "\n";
    for (int blk = 0; blk < (haveLut ? 3 : 1); ++blk) { s << (blk == 0 ? "   2100 V onset\n" : blk == 1 ? "   1800 V LUT500 single-electron drift time\n" : "   difference (2100 onset - 1800 drift); different observables, see header\n");
      for (int b = 0; b < NV; ++b) { s << "    v " << std::setw(5) << -3. + 0.5 * b << ":"; for (int a = 0; a < NU; ++a) { const C2& x = c2[a * NV + b]; const C1& y = c1[a * NV + b]; const double m2 = x.nOn > 0 ? x.on / x.nOn : kNaN, m1 = y.nT > 0 ? y.t / y.nT : kNaN;
          s << std::setw(9) << (blk == 0 ? m2 : blk == 1 ? m1 : m2 - m1); } s << "\n"; } }
    if (haveLut) { double w = 0, d18 = 0, d21 = 0, p18 = 0, p21 = 0, dw18 = 0, dw21 = 0, s18 = 0, s21 = 0;
      for (int k = 0; k < NU * NV; ++k) { const C2& x = c2[k]; const C1& y = c1[k]; if (x.nOn < 1 || y.nT < 2) continue; double m, r; w += x.nOn; d21 += x.on; mr(y.t, y.t2, y.nT, m, r); d18 += x.nOn * m;
        p21 += x.zc / x.n * x.nOn; p18 += (1. - y.nc / y.n) * x.nOn; mr(x.dw, x.dw2, x.nDw, m, r); if (std::isfinite(r)) dw21 += x.nOn * r * r; mr(y.dw, y.dw2, y.nDw, m, r); if (std::isfinite(r)) dw18 += x.nOn * r * r;
        s21 += x.same; s18 += y.nc > 0 ? x.nOn * y.same / y.nc : 0; }
      s << "  occupancy-weighted (2100 V start positions) comparison over cells present in both:\n    mean delay 2100 onset " << d21 / w << " vs 1800 LUT500 drift " << d18 / w << " ns;  P(no collection) 2100 " << p21 / w << " vs 1800 " << p18 / w
        << ";  within-cell RMS dw 2100 " << std::sqrt(dw21 / w) << " vs 1800 " << std::sqrt(dw18 / w) << " mm;  collected on start-nearest wire 2100 " << s21 / w << " vs 1800 " << s18 / w << "\n  LUT500 inputs: " << lutNote << "\n"; }
    else s << "  no --lut500 file given: no 1800 V comparison made\n";
    s << "  NOT available at 2100 V in K256 archives: individual electron drift lines and electron arrival azimuth at the wire (ion-cloud <sin>,<cos> are a different observable)\n";
    out(s);
    std::ofstream f(od + "/v0_transport_map.csv"); f << std::setprecision(8) << "cell,du_lo_mm,du_hi_mm,v_lo_mm,v_hi_mm,v2100_n_primaries,v2100_p_collected0,v2100_p_ion0,v2100_p_ion_pos_collected0,v2100_n_onset,v2100_onset_mean_ns,v2100_onset_rms_ns,v2100_ion_delay_mean_ns,v2100_ion_minus_onset_mean_ns,v2100_sigma_t_mean_ns,v2100_dw_mean_mm,v2100_dw_rms_mm,v2100_same_wire_frac,v2100_multiwire_frac,v2100_cloud_sin_mean,v2100_cloud_cos_mean,v2100_attached_share,lut1800_n,lut1800_p_not_collected,lut1800_drift_mean_ns,lut1800_drift_rms_ns,lut1800_dw_mean_mm,lut1800_dw_rms_mm,lut1800_same_wire_frac,lut1800_hit_sin_mean,lut1800_hit_cos_mean\n";
    for (int a = 0; a < NU; ++a) for (int b = 0; b < NV; ++b) { const int k = a * NV + b; const C2& x = c2[k]; const C1& y = c1[k]; double om, orr, dm, dr, tm, tr, ym, yr; mr(x.on, x.on2, x.nOn, om, orr); mr(x.dw, x.dw2, x.nDw, dm, dr); mr(y.t, y.t2, y.nT, tm, tr); mr(y.dw, y.dw2, y.nDw, ym, yr);
      f << k << ',' << 0.25 * a << ',' << 0.25 * (a + 1) << ',' << -3 + 0.5 * b << ',' << -2.5 + 0.5 * b << ',' << x.n << ',' << (x.n ? x.zc / x.n : kNaN) << ',' << (x.n ? x.zi / x.n : kNaN) << ',' << (x.n ? x.inc / x.n : kNaN) << ',' << x.nOn << ',' << om << ',' << orr << ','
        << (x.nIo ? x.io / x.nIo : kNaN) << ',' << (x.nOff ? x.off / x.nOff : kNaN) << ',' << (x.nSig ? x.sig / x.nSig : kNaN) << ',' << dm << ',' << dr << ',' << (x.nOn ? x.same / x.nOn : kNaN) << ',' << (x.nOn ? x.multi / x.nOn : kNaN) << ','
        << (x.nDir ? x.sn / x.nDir : kNaN) << ',' << (x.nDir ? x.cs / x.nDir : kNaN) << ',' << (x.nAtt ? x.att / x.nAtt : kNaN) << ',';
      if (haveLut) f << y.n << ',' << (y.n ? 1 - y.nc / y.n : kNaN) << ',' << tm << ',' << tr << ',' << ym << ',' << yr << ',' << (y.nc ? y.same / y.nc : kNaN) << ',' << (y.nc ? y.sn / y.nc : kNaN) << ',' << (y.nc ? y.cs / y.nc : kNaN) << '\n';
      else f << ",,,,,,,,\n"; }
  }

  // ===== exports =====
  { std::ofstream f(od + "/v0_fragments.csv"); f << std::setprecision(10) << "muon_slot,primary,global_wire_nearest_birth,n_reservoir,N_hat,neff_res,A,se_A,sin,se_sin,cos,se_cos,dy_um,se_dy,dx_um,se_dx,sin2,cos2,rbar_um,delay_ns,se_delay,sigt_ns,wbar_cm,sigw_cm,exact_births,exact_collected_e,exact_A,exact_delay_ns,exact_sigt_ns,exact_sigw_mm,start_v_mm,start_du_mm,primary_collected_e,primary_ion_births_exact,exact_onset_ns,exact_sigt_collected_ns\n";
    for (auto& g : F) { f << mus[g.muon].slot << ',' << g.prim << ',' << g.gwire << ',' << g.nres << ',' << g.Nhat << ',' << g.neffRes; for (int k = 0; k < 5; ++k) f << ',' << g.x[k] << ',' << g.se[k];
      f << ',' << g.x[5] << ',' << g.x[6] << ',' << g.rbar << ',' << g.tbar << ',' << g.tbarSE << ',' << g.sigt << ',' << g.wbar << ',' << g.sigw << ',' << g.exBirths << ',' << g.exCollected << ',' << g.exA << ',' << g.exTbar << ',' << g.exSigt << ',' << g.exSigwMm << ',' << g.startVmm << ',' << g.startDuMm << ',' << mus[g.muon].prims[g.prim].collected << ',' << g.NprimExact << ',' << g.exTcol << ',' << g.exSigTcol << '\n'; } }
  { std::ofstream f(od + "/v0_primaries.csv"); f << std::setprecision(10) << "muon_slot,primary,ion_births_exact,n_saved,collected_electrons,n_fragments,N_hat,A,se_A,sin,cos,dy_um,dx_um,delay_ns,exact_A,start_v_mm,start_du_mm,t_primary_ns,argmax_collected_wire,argmax_birth_wire,synthetic_law_A,start_nearest_wire,exact_delay_ns,exact_sigt_prim_ns,exact_onset_ns,collected_frac_max_wire\n";
    for (auto& r : P) f << mus[r.muon].slot << ',' << r.prim << ',' << r.Nexact << ',' << r.nSaved << ',' << r.collected << ',' << r.nFrag << ',' << r.Nhat << ',' << r.x[0] << ',' << r.se[0] << ',' << r.x[1] << ',' << r.x[2] << ',' << r.x[3] << ',' << r.x[4] << ',' << r.tbar << ',' << r.exA << ',' << r.startVmm << ',' << r.startDuMm << ',' << r.tPrim << ',' << r.collWireArgmax << ',' << r.birthWireArgmax << ',' << (mus[r.muon].lawA.count(r.prim) ? mus[r.muon].lawA.at(r.prim) : kNaN) << ',' << r.startWire << ',' << r.exDelay << ',' << r.exSigtPrim << ',' << r.exOnset << ',' << r.collFracMax << '\n'; }
  if (A.count("ion-dump-muons")) {   // raw reservoir ions of the first N muons, for within-avalanche time plots
    const size_t nd = std::min(mus.size(), (size_t)std::stoul(A["ion-dump-muons"])); std::ofstream f(od + "/v0_ion_examples.csv"); f << std::setprecision(10) << "muon_slot,primary,global_wire_nearest_birth,t_minus_tprim_ns,weight,dx_um,dy_um\n";
    for (size_t im = 0; im < nd; ++im) { const Muon& m = mus[im]; std::map<int, double> tp; for (auto& p : m.prims) tp[p.idx] = p.t;
      for (auto& i : m.ions) f << m.slot << ',' << i.prim << ',' << i.gwire << ',' << i.t - tp[i.prim] << ',' << i.wt << ',' << 1e4 * (i.u - i.lwire * m.pitch) << ',' << 1e4 * i.v << '\n'; } }
  { std::ofstream f(od + "/v0_wires.csv"); f << std::setprecision(10) << "muon_slot,global_wire,n_fragments,N_hat,N_exact,collected_e_exact,Neff,sumN_A,sumN_sin,sumN_cos,sumN_dy_um,sumN_dx_um,sumN_delay,sumN_delay2,sumN_w,sumN_w2,sumN_A_delay,sumNexact_Aexact\n";
    std::map<std::pair<int, int>, std::vector<const Frag*>> U; for (auto& g : F) U[{g.muon, g.gwire}].push_back(&g);
    for (auto& kv : U) { double N = 0, N2 = 0, Nx = 0, col = 0, s[5] = {0}, t = 0, tt = 0, w = 0, ww = 0, at = 0, ax = 0; bool ex = true;
      for (auto* g : kv.second) { if (g->exBirths < 0) ex = false; else { Nx += g->exBirths; col += g->exCollected; if (!std::isnan(g->exA)) ax += g->exBirths * g->exA; } if (g->nres == 0) continue; N += g->Nhat; N2 += g->Nhat * g->Nhat;
        for (int k = 0; k < 5; ++k) if (!std::isnan(g->x[k])) s[k] += g->Nhat * g->x[k]; t += g->Nhat * g->tbar; tt += g->Nhat * (g->tbar * g->tbar + g->sigt * g->sigt); w += g->Nhat * g->wbar; ww += g->Nhat * (g->wbar * g->wbar + g->sigw * g->sigw); at += g->Nhat * g->x[0] * g->tbar; }
      f << mus[kv.first.first].slot << ',' << kv.first.second << ',' << kv.second.size() << ',' << N << ',' << (ex ? Nx : kNaN) << ',' << (ex ? col : kNaN) << ',' << (N2 > 0 ? N * N / N2 : 0.); for (int k = 0; k < 5; ++k) f << ',' << s[k];
      f << ',' << t << ',' << tt << ',' << w << ',' << ww << ',' << at << ',' << (ex ? ax : kNaN) << '\n'; } }
  { std::ofstream f(od + "/v0_muons.csv"); f << std::setprecision(10) << "muon_slot,crossing_id,channel,station,det_id,theta_u_deg,theta_w_deg,n_primaries,n_zero_ion,n_zero_collected,ion_births_exact,collected_electrons,A_hat,A_exact,census_csv\n";
    for (size_t im = 0; im < mus.size(); ++im) { const Muon& m = mus[im]; long long z0 = 0, zc = 0, N = 0, C = 0; for (auto& p : m.prims) { z0 += p.Nexact == 0; zc += p.collected == 0; N += p.Nexact; C += p.collected; }
      double a = 0, w = 0, ax = 0, wx = 0; for (auto& r : P) if ((size_t)r.muon == im) { if (r.Nhat > 0 && !std::isnan(r.x[0])) { a += r.Nhat * r.x[0]; w += r.Nhat; } if (r.Nexact > 0 && !std::isnan(r.exA)) { ax += r.Nexact * r.exA; wx += r.Nexact; } }
      f << m.slot << ',' << m.id << ',' << m.channel << ',' << m.station << ',' << m.det << ',' << m.thetaU << ',' << m.thetaW << ',' << m.prims.size() << ',' << z0 << ',' << zc << ',' << N << ',' << C << ',' << (w > 0 ? a / w : kNaN) << ',' << (wx > 0 ? ax / wx : kNaN) << ',' << m.haveSeed << '\n'; } }
  { std::ofstream f(od + "/v0_manifest.txt"); f << "tool: v0_archive_stats (V0, archive-only; no weighting field, no transport)\nmode: " << (synthetic ? "SYNTHETIC (estimator validation, not NA60+ data)" : "archives") << "\ncalibration_use_allowed: " << (synthetic ? "no" : (calibrationOK ? "yes (holdout excluded)" : "NO (no holdout exclusion)")) << "\nmuons_used: " << mus.size() << "\n";
    f << "excluded_holdout_slots:"; for (int e : excl) f << ' ' << std::setw(4) << std::setfill('0') << e; f << std::setfill(' ') << "\nskipped_files:\n"; for (auto& s : skipped) f << "  " << s << "\n";
    f << "census_csv_problems:\n"; for (auto& s : seedProblems) f << "  " << s << "\n";
    f << "debug_only_slots_note: 0301,0408,0446,0435,0511,0101 were inspected during development; they are allowed in population statistics, never in sealed validation\n";
    f << "conventions: wire = nearest anode to ion birth; A excludes theta=0,pi births; t relative to primaries.t_primary_ns; N_hat = sum of inverse inclusion weights\n"; }
  if (!seedProblems.empty()) std::cerr << "WARNING: census CSV unusable for " << seedProblems.size() << " muons (see v0_manifest.txt)\n";
  return 0;
  } catch (const std::exception& e) { std::cerr << "V0 FAILED: " << e.what() << "\n"; return 2; }
}

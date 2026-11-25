// SYNTHETIC STAND-IN for the few ROOT TFile/TTree calls used by the archive writer and readers - NOT ROOT.
// A "ROOT file" is a text file: per tree a header and rows; doubles are written with 17 significant digits, so every
// value reads back bit for bit (strtod of a %.17g string is exact). Tests only.
#pragma once
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Long64_t = long long;
using Int_t = int;
class TFile;
class TBranch {};
class TObject { public: virtual ~TObject() = default; };

class TTree : public TObject {
 public:
  struct Col { std::string name; char type; void* addr = nullptr; std::vector<std::string> v; };
  TTree(const char* name, const char* title) : name_(name), title_(title) {}
  void SetDirectory(TFile* f);
  void SetAutoFlush(long long) {}
  TBranch* Branch(const char* n, int* a, const char*) { return Add(n, 'I', a); }
  TBranch* Branch(const char* n, long long* a, const char*) { return Add(n, 'L', a); }
  TBranch* Branch(const char* n, double* a, const char*) { return Add(n, 'D', a); }
  TBranch* Branch(const char* n, std::string* a) { return Add(n, 'S', a); }
  int Fill() {
    for (auto& c : cols_) {
      char b[64];
      switch (c.type) {
        case 'I': std::snprintf(b, sizeof b, "%d", *static_cast<int*>(c.addr)); c.v.push_back(b); break;
        case 'L': std::snprintf(b, sizeof b, "%lld", *static_cast<long long*>(c.addr)); c.v.push_back(b); break;
        case 'D': std::snprintf(b, sizeof b, "%.17g", *static_cast<double*>(c.addr)); c.v.push_back(b); break;
        default: c.v.push_back(Esc(*static_cast<std::string*>(c.addr)));
      }
    }
    return int(++n_);
  }
  int Write() { return 1; }
  Long64_t GetEntries() const { return n_; }
  TBranch* GetBranch(const char* n) { static TBranch b; return Find(n) ? &b : nullptr; }
  int SetBranchAddress(const char* n, int* a) { return Set(n, 'I', a); }
  int SetBranchAddress(const char* n, long long* a) { return Set(n, 'L', a); }
  int SetBranchAddress(const char* n, double* a) { return Set(n, 'D', a); }
  int SetBranchAddress(const char* n, std::string** a) {
    if (!*a) { owned_.push_back(std::make_unique<std::string>()); *a = owned_.back().get(); }
    return Set(n, 'S', *a);
  }
  int GetEntry(Long64_t i) {
    if (i < 0 || i >= n_) return 0;
    for (auto& c : cols_) {
      if (!c.addr) continue;
      const std::string& s = c.v[std::size_t(i)];
      switch (c.type) {
        case 'I': *static_cast<int*>(c.addr) = std::atoi(s.c_str()); break;
        case 'L': *static_cast<long long*>(c.addr) = std::atoll(s.c_str()); break;
        case 'D': *static_cast<double*>(c.addr) = std::strtod(s.c_str(), nullptr); break;
        default: *static_cast<std::string*>(c.addr) = Unesc(s);
      }
    }
    return 8;
  }
  // serialisation (used by TFile)
  void Save(std::ostream& o) const {
    o << "TREE\t" << name_ << '\t' << cols_.size() << '\t' << n_ << '\n';
    for (const auto& c : cols_) o << "COL\t" << c.name << '\t' << c.type << '\n';
    for (Long64_t i = 0; i < n_; ++i) {
      for (std::size_t j = 0; j < cols_.size(); ++j) o << (j ? "\t" : "") << cols_[j].v[std::size_t(i)];
      o << '\n';
    }
  }
  static TTree* Load(std::istream& in, const std::string& head) {
    std::istringstream h(head); std::string tag, name; std::size_t nc = 0; long long nr = 0;
    h >> tag >> name >> nc >> nr;
    auto* t = new TTree(name.c_str(), "");
    std::string l;
    for (std::size_t j = 0; j < nc; ++j) {
      std::getline(in, l); std::istringstream c(l); std::string ct, cn; char ty = 0; c >> ct >> cn >> ty;
      t->cols_.push_back({cn, ty, nullptr, {}});
    }
    for (long long i = 0; i < nr; ++i) {
      std::getline(in, l);
      std::size_t p = 0;
      for (std::size_t j = 0; j < nc; ++j) {
        const std::size_t q = l.find('\t', p);
        t->cols_[j].v.push_back(l.substr(p, q == std::string::npos ? std::string::npos : q - p));
        p = q == std::string::npos ? l.size() : q + 1;
      }
    }
    t->n_ = nr;
    return t;
  }
  const std::string& Name() const { return name_; }

 private:
  Col* Find(const std::string& n) { for (auto& c : cols_) if (c.name == n) return &c; return nullptr; }
  TBranch* Add(const char* n, char ty, void* a) { cols_.push_back({n, ty, a, {}}); static TBranch b; return &b; }
  int Set(const char* n, char ty, void* a) {
    Col* c = Find(n);
    if (!c) return -5;
    if (c->type != ty) return -2;
    c->addr = a; return 0;
  }
  static std::string Esc(const std::string& s) { std::string o; for (char ch : s) o += (ch == '\t' ? ' ' : ch == '\n' ? ' ' : ch); return o.empty() ? std::string("\x01") : o; }
  static std::string Unesc(const std::string& s) { return s == "\x01" ? std::string() : s; }
  std::string name_, title_;
  std::vector<Col> cols_;
  Long64_t n_ = 0;
  std::vector<std::unique_ptr<std::string>> owned_;
};

class TFile : public TObject {
 public:
  static TFile* Open(const char* path, const char* opt, const char* = "", int = 0) {
    auto* f = new TFile; f->path_ = path; const std::string o = opt;
    if (o == "RECREATE") { std::ofstream t(path); f->zombie_ = !t; f->write_ = true; return f; }
    std::ifstream in(path);
    std::string l;
    if (!in || !std::getline(in, l) || l != "STUBROOT1") { f->zombie_ = true; return f; }
    while (std::getline(in, l)) if (l.rfind("TREE\t", 0) == 0) f->trees_.push_back(TTree::Load(in, l));
    return f;
  }
  ~TFile() override { for (auto* t : trees_) delete t; }
  bool IsZombie() const { return zombie_; }
  void SetCompressionLevel(int) {}
  void cd() {}
  void Register(TTree* t) { trees_.push_back(t); }
  int Write() {
    if (!write_) return 0;
    std::ofstream o(path_); o << "STUBROOT1\n";
    for (auto* t : trees_) t->Save(o);
    if (!o) throw std::runtime_error("stub ROOT write failed: " + path_);
    return 1;
  }
  void Close() {}
  TObject* Get(const char* name) { for (auto* t : trees_) if (t->Name() == name) return t; return nullptr; }
 private:
  std::string path_;
  bool zombie_ = false, write_ = false;
  std::vector<TTree*> trees_;
};
inline void TTree::SetDirectory(TFile* f) { if (f) f->Register(this); }

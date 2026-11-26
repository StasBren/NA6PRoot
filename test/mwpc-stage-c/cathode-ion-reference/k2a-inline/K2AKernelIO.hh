// Kernel-directory I/O for the compiled K2A strip evaluator: minimal JSON reader, read-only .npy mapping, SHA-256.
// No ROOT, no Python, no third-party library. Used once per process to open an s_kernel.py kernel directory
// (K.npy float32 (basis, side, wire, r, D, tau), KPHI0.npy float64 (basis, side, wire, r, D), kernel.json).
// kernel_id = sha256 over the concatenated lowercase hex sha256 of K.npy, KPHI0.npy, kernel.json (exactly
// cathode_engine.KernelTable.kernel_id).
#pragma once
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace K2AKernelIO {

// ------------------------------------------------------------------------------------------------ SHA-256 (FIPS 180-4)
class Sha256 {
 public:
  Sha256() { Reset(); }
  void Reset() {
    static const std::uint32_t h0[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::memcpy(h_, h0, sizeof h_); len_ = 0; n_ = 0;
  }
  void Update(const void* p, std::size_t n) {
    const auto* b = static_cast<const unsigned char*>(p);
    len_ += n;
    while (n) {
      const std::size_t k = std::min<std::size_t>(64 - n_, n);
      std::memcpy(buf_ + n_, b, k); n_ += k; b += k; n -= k;
      if (n_ == 64) { Block(buf_); n_ = 0; }
    }
  }
  std::string HexDigest() {
    const std::uint64_t bits = len_ * 8;
    const unsigned char pad = 0x80, z = 0;
    Update(&pad, 1);
    while (n_ != 56) Update(&z, 1);
    unsigned char L[8];
    for (int i = 0; i < 8; ++i) L[i] = (unsigned char)(bits >> (56 - 8 * i));
    Update(L, 8);
    char out[65];
    for (int i = 0; i < 8; ++i) std::snprintf(out + 8 * i, 9, "%08x", h_[i]);
    return std::string(out, 64);
  }
  static std::string OfFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + path);
    Sha256 h; std::vector<char> b(1 << 22);
    while (f) { f.read(b.data(), std::streamsize(b.size())); h.Update(b.data(), std::size_t(f.gcount())); }
    return h.HexDigest();
  }
  static std::string OfString(const std::string& s) { Sha256 h; h.Update(s.data(), s.size()); return h.HexDigest(); }

 private:
  static std::uint32_t R(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
  void Block(const unsigned char* p) {
    static const std::uint32_t K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01,
        0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
        0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08,
        0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
    std::uint32_t w[64];
    for (int i = 0; i < 16; ++i) w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16) | (std::uint32_t(p[4 * i + 2]) << 8) | p[4 * i + 3];
    for (int i = 16; i < 64; ++i) {
      const std::uint32_t s0 = R(w[i - 15], 7) ^ R(w[i - 15], 18) ^ (w[i - 15] >> 3), s1 = R(w[i - 2], 17) ^ R(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3], e = h_[4], f = h_[5], g = h_[6], h = h_[7];
    for (int i = 0; i < 64; ++i) {
      const std::uint32_t t1 = h + (R(e, 6) ^ R(e, 11) ^ R(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i];
      const std::uint32_t t2 = (R(a, 2) ^ R(a, 13) ^ R(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
      h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h_[0] += a; h_[1] += b; h_[2] += c; h_[3] += d; h_[4] += e; h_[5] += f; h_[6] += g; h_[7] += h;
  }
  std::uint32_t h_[8];
  unsigned char buf_[64];
  std::size_t n_;
  std::uint64_t len_;
};

// --------------------------------------------------------------------------------------------------- minimal JSON
struct Json {
  enum Type { Null, Bool, Number, String, Array, Object } type = Null;
  bool b = false;
  double num = 0;
  std::string str;
  std::vector<Json> arr;
  std::map<std::string, Json> obj;

  const Json& operator[](const std::string& k) const {
    auto it = obj.find(k);
    if (type != Object || it == obj.end()) throw std::runtime_error("JSON: missing key '" + k + "'");
    return it->second;
  }
  bool Has(const std::string& k) const { return type == Object && obj.count(k); }
  double Num() const { if (type != Number) throw std::runtime_error("JSON: not a number"); return num; }
  const std::string& Str() const { if (type != String) throw std::runtime_error("JSON: not a string"); return str; }
  std::vector<double> Nums() const {
    if (type != Array) throw std::runtime_error("JSON: not an array");
    std::vector<double> v; for (const auto& x : arr) v.push_back(x.Num()); return v;
  }

  static Json Parse(const std::string& s) { std::size_t i = 0; Json j = Value(s, i); Ws(s, i); if (i != s.size()) throw std::runtime_error("JSON: trailing data"); return j; }
  static Json ParseFile(const std::string& p) {
    std::ifstream f(p); if (!f) throw std::runtime_error("cannot read " + p);
    std::stringstream ss; ss << f.rdbuf(); return Parse(ss.str());
  }

 private:
  static void Ws(const std::string& s, std::size_t& i) { while (i < s.size() && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r' || s[i] == '\t')) ++i; }
  static void Expect(const std::string& s, std::size_t& i, const char* lit) {
    const std::size_t n = std::strlen(lit);
    if (s.compare(i, n, lit) != 0) throw std::runtime_error(std::string("JSON: expected ") + lit);
    i += n;
  }
  static std::string Str(const std::string& s, std::size_t& i) {
    if (s[i] != '"') throw std::runtime_error("JSON: expected string");
    ++i; std::string o;
    while (i < s.size() && s[i] != '"') {
      char c = s[i++];
      if (c != '\\') { o += c; continue; }
      if (i >= s.size()) break;
      const char e = s[i++];
      switch (e) {
        case 'n': o += '\n'; break; case 't': o += '\t'; break; case 'r': o += '\r'; break; case 'b': o += '\b'; break;
        case 'f': o += '\f'; break; case 'u': {
          unsigned cp = std::stoul(s.substr(i, 4), nullptr, 16); i += 4;
          if (cp < 0x80) o += char(cp); else if (cp < 0x800) { o += char(0xc0 | (cp >> 6)); o += char(0x80 | (cp & 0x3f)); }
          else { o += char(0xe0 | (cp >> 12)); o += char(0x80 | ((cp >> 6) & 0x3f)); o += char(0x80 | (cp & 0x3f)); }
          break; }
        default: o += e;
      }
    }
    if (i >= s.size()) throw std::runtime_error("JSON: unterminated string");
    ++i; return o;
  }
  static Json Value(const std::string& s, std::size_t& i) {
    Ws(s, i);
    if (i >= s.size()) throw std::runtime_error("JSON: unexpected end");
    Json j;
    const char c = s[i];
    if (c == '{') {
      j.type = Object; ++i; Ws(s, i);
      if (s[i] == '}') { ++i; return j; }
      for (;;) {
        Ws(s, i); std::string k = Str(s, i); Ws(s, i); Expect(s, i, ":");
        j.obj[k] = Value(s, i); Ws(s, i);
        if (s[i] == ',') { ++i; continue; }
        Expect(s, i, "}"); return j;
      }
    }
    if (c == '[') {
      j.type = Array; ++i; Ws(s, i);
      if (s[i] == ']') { ++i; return j; }
      for (;;) {
        j.arr.push_back(Value(s, i)); Ws(s, i);
        if (s[i] == ',') { ++i; continue; }
        Expect(s, i, "]"); return j;
      }
    }
    if (c == '"') { j.type = String; j.str = Str(s, i); return j; }
    if (s.compare(i, 4, "true") == 0) { i += 4; j.type = Bool; j.b = true; return j; }
    if (s.compare(i, 5, "false") == 0) { i += 5; j.type = Bool; return j; }
    if (s.compare(i, 4, "null") == 0) { i += 4; return j; }
    if (s.compare(i, 3, "NaN") == 0) { i += 3; j.type = Number; j.num = NAN; return j; }
    if (s.compare(i, 8, "Infinity") == 0) { i += 8; j.type = Number; j.num = INFINITY; return j; }
    if (s.compare(i, 9, "-Infinity") == 0) { i += 9; j.type = Number; j.num = -INFINITY; return j; }
    char* end = nullptr;
    j.num = std::strtod(s.c_str() + i, &end);     // correctly rounded: the same double as Python's float(text)
    if (end == s.c_str() + i) throw std::runtime_error("JSON: bad value");
    i = std::size_t(end - s.c_str()); j.type = Number;
    return j;
  }
};

// ----------------------------------------------------------------------------------------- read-only .npy mapping
// Supports npy format 1.x/2.x/3.x, little-endian '<f4' / '<f8', C order. The data are mmap'ed read-only (shared page
// cache, nothing copied) unless preload is requested; the mapping is released with the object.
class NpyArray {
 public:
  NpyArray() = default;
  NpyArray(const NpyArray&) = delete;
  NpyArray& operator=(const NpyArray&) = delete;
  NpyArray(NpyArray&& o) noexcept { *this = std::move(o); }
  NpyArray& operator=(NpyArray&& o) noexcept {
    std::swap(map_, o.map_); std::swap(len_, o.len_); std::swap(data_, o.data_); std::swap(shape, o.shape);
    std::swap(descr, o.descr); std::swap(copy_, o.copy_); return *this;
  }
  ~NpyArray() { if (map_) munmap(map_, len_); }

  static NpyArray Open(const std::string& path, bool preload) {
    NpyArray a;
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("cannot open " + path);
    struct stat st{};
    if (fstat(fd, &st) != 0) { ::close(fd); throw std::runtime_error("cannot stat " + path); }
    a.len_ = std::size_t(st.st_size);
    a.map_ = mmap(nullptr, a.len_, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd);
    if (a.map_ == MAP_FAILED) { a.map_ = nullptr; throw std::runtime_error("cannot map " + path); }
    const auto* p = static_cast<const unsigned char*>(a.map_);
    if (a.len_ < 10 || std::memcmp(p, "\x93NUMPY", 6) != 0) throw std::runtime_error(path + ": not a .npy file");
    const int major = p[6];
    std::size_t hl, off;
    if (major == 1) { hl = p[8] | (std::size_t(p[9]) << 8); off = 10; }
    else { hl = p[8] | (std::size_t(p[9]) << 8) | (std::size_t(p[10]) << 16) | (std::size_t(p[11]) << 24); off = 12; }
    const std::string h(reinterpret_cast<const char*>(p + off), hl);
    auto field = [&](const std::string& k) {
      const std::size_t q = h.find("'" + k + "'");
      if (q == std::string::npos) throw std::runtime_error(path + ": npy header lacks " + k);
      return q + k.size() + 2;
    };
    std::size_t q = h.find('\'', h.find(':', field("descr")));
    a.descr = h.substr(q + 1, h.find('\'', q + 1) - q - 1);
    q = h.find(':', field("fortran_order"));
    if (h.compare(h.find_first_not_of(' ', q + 1), 5, "False") != 0) throw std::runtime_error(path + ": Fortran order not supported");
    q = h.find('(', field("shape"));
    const std::string sh = h.substr(q + 1, h.find(')', q) - q - 1);
    std::stringstream ss(sh); std::string tok;
    while (std::getline(ss, tok, ',')) if (tok.find_first_not_of(" ") != std::string::npos) a.shape.push_back(std::stoll(tok));
    if (a.descr != "<f4" && a.descr != "<f8") throw std::runtime_error(path + ": dtype " + a.descr + " not supported");
    std::size_t n = 1; for (auto d : a.shape) n *= std::size_t(d);
    const std::size_t es = a.descr == "<f4" ? 4 : 8;
    if (off + hl + n * es != a.len_) throw std::runtime_error(path + ": size does not match its header");
    a.data_ = p + off + hl;
    if (reinterpret_cast<std::uintptr_t>(a.data_) % es) throw std::runtime_error(path + ": misaligned data");
    if (preload) {
      a.copy_.assign(a.data_, a.data_ + n * es);
      munmap(a.map_, a.len_); a.map_ = nullptr; a.len_ = 0;
      a.data_ = a.copy_.data();
    }
    return a;
  }
  template <class T> const T* Data() const {
    if ((sizeof(T) == 4 && descr != "<f4") || (sizeof(T) == 8 && descr != "<f8")) throw std::runtime_error("npy dtype mismatch");
    return reinterpret_cast<const T*>(data_);
  }
  std::size_t Bytes() const { std::size_t n = descr == "<f4" ? 4 : 8; for (auto d : shape) n *= std::size_t(d); return n; }
  bool Mapped() const { return map_ != nullptr; }
  std::vector<long long> shape;
  std::string descr;

 private:
  void* map_ = nullptr;
  std::size_t len_ = 0;
  const unsigned char* data_ = nullptr;
  std::vector<unsigned char> copy_;
};
}  // namespace K2AKernelIO

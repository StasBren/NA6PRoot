// SYNTHETIC STAND-IN for the Garfield++ classes used by the 2100 V generator - NOT Garfield, NOT physics.
// Purpose: compile the frozen and the K2A-inline generator sources in an environment without Garfield++/ROOT and run
// them end to end on a deterministic TOY avalanche, to test the plumbing (option parsing, sink composition, census,
// archive writing, statuses, exit codes, harness and comparator). Every number produced with these stubs is synthetic.
//   toy Heed:      clusters every 0.02 cm along the muon direction, 1..3 electrons each, seeded by the Heed seed
//   toy avalanche: ion births around the nearest wire of the start point, count and positions from the Garfield seed
//                  (Random::SetEngine), log-uniform radius 16 um .. 0.5 mm, exponential time tail; one collected
//                  electron endpoint per ion + 1
//   K2A_STUB_INSIDE_ANODE=1   primary 1 additionally gets one birth at r = 10 um (accumulator must block)
//   K2A_STUB_OUTSIDE_PATCH=1  primary 1 additionally gets one birth on wire 5 (the generator must abort)
//   K2A_STUB_ZERO_PRIMARY=N   primary N gets no birth
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace Garfield {
class RandomEngineRoot {
 public:
  explicit RandomEngineRoot(unsigned int s) : seed(s) {}
  unsigned int seed;
};
namespace Random {
inline RandomEngineRoot*& Engine() { static RandomEngineRoot* e = nullptr; return e; }
inline void SetEngine(RandomEngineRoot& e) { Engine() = &e; }
inline std::uint64_t Seed() { return Engine() ? Engine()->seed : 0u; }
}  // namespace Random

class Medium {
 public:
  virtual ~Medium() = default;
  bool IsIonisable() const { return true; }
};
class MediumMagboltz : public Medium {
 public:
  template <class... A> void SetComposition(A&&...) {}
  void SetTemperature(double) {}
  void SetPressure(double) {}
  void SetMaxElectronEnergy(double) {}
  bool Initialise(bool = false) { return true; }
  template <class... A> void EnablePenningTransfer(A&&...) {}
  template <class... A> bool LoadIonMobility(A&&...) { return true; }
};
class ComponentAnalyticField {
 public:
  template <class... A> double WeightingPotential(A&&...) { return 0.; }
  void SetMedium(Medium*) {}
  template <class... A> void AddWire(A&&...) {}
  template <class... A> void AddPlaneY(A&&...) {}
  template <class... A> void SetMagneticField(A&&...) {}
  template <class... A> void AddReadout(A&&...) {}
};
class ComponentConstant {
 public:
  template <class... A> void SetArea(A&&...) {}
  void SetMedium(Medium* m) { medium = m; }
  template <class... A> void SetElectricField(A&&...) {}
  Medium* medium = nullptr;
};
class Sensor {
 public:
  void AddComponent(ComponentConstant* c) { medium_ = c->medium; }
  void AddComponent(ComponentAnalyticField*) {}
  template <class... A> void SetArea(A&&...) {}
  template <class... A> void AddElectrode(A&&...) {}
  template <class... A> void SetTimeWindow(A&&...) {}
  void ClearSignal() {}
  double GetElectronSignal(const std::string&, int) { return 0.; }
  double GetIonSignal(const std::string&, int) { return 0.; }
  Medium* GetMedium(double, double, double) { static MediumMagboltz gas; return medium_ ? medium_ : &gas; }
 private:
  Medium* medium_ = nullptr;
};

class TrackHeed {
 public:
  void SetSensor(Sensor*) {}
  void SetParticle(const std::string&) {}
  void SetMomentum(double) {}
  void EnableDeltaElectronTransport() {}
  bool NewTrack(double x, double y, double z, double t, double dx, double dy, double dz) {
    x0_ = x; y0_ = y; z0_ = z; t0_ = t; d_[0] = dx; d_[1] = dy; d_[2] = dz; k_ = 0;
    rng_.seed(Random::Seed() * 0x9e3779b97f4a7c15ULL + 11);
    return true;
  }
  bool GetCluster(double& x, double& y, double& z, double& t, int& n, double& e, double& extra) {
    if (k_ >= 70) return false;
    const double s = 0.02 * (k_ + 0.5);
    std::uniform_real_distribution<double> U(0, 1);
    x = x0_ + s * d_[0]; y = y0_ + s * d_[1]; z = z0_ + s * d_[2]; t = t0_ + s / 30.;
    n = 1 + int(rng_() % 3); e = 30.; extra = 0.; ++k_;
    els_.clear();
    for (int i = 0; i < n; ++i) els_.push_back({x + 0.001 * (U(rng_) - .5), y + 0.001 * (U(rng_) - .5), z + 0.001 * (U(rng_) - .5), t});
    return true;
  }
  bool GetElectron(int i, double& x, double& y, double& z, double& t, double& e, double& dx, double& dy, double& dz) {
    if (i < 0 || i >= int(els_.size())) return false;
    x = els_[i][0]; y = els_[i][1]; z = els_[i][2]; t = els_[i][3]; e = 0.1; dx = dy = dz = 0.;
    return true;
  }
 private:
  double x0_ = 0, y0_ = 0, z0_ = 0, t0_ = 0, d_[3]{};
  int k_ = 0;
  std::mt19937_64 rng_;
  std::vector<std::vector<double>> els_;
};

class AvalancheMicroscopic {
 public:
  using Handler = void (*)(double, double, double, double, int, int, Medium*);
  explicit AvalancheMicroscopic(Sensor*) {}
  void DisableAvalancheSizeLimit() {}
  void EnableSignalCalculation(bool = true) {}
  void SetUserHandleIonisation(Handler h) { h_ = h; }
  bool AvalancheElectron(double x, double y, double z, double t, double, double, double, double) {
    std::mt19937_64 g(Random::Seed() * 0xbf58476d1ce4e5b9ULL + 7);
    std::uniform_real_distribution<double> U(0, 1);
    const int k = int(std::lround(x / 0.3));
    static int primaryCounter = -1; ++primaryCounter;   // counts avalanches in this process (stub-only switches)
    const char* zp = std::getenv("K2A_STUB_ZERO_PRIMARY");
    int n = 40 + int(g() % 900);
    if (zp && std::atoi(zp) == primaryCounter) n = 0;
    ends_.clear(); ni_ = 0;
    for (int j = 0; j < n; ++j) {
      const double r = 0.0016 * std::pow(0.05 / 0.0016, U(g)), ph = 6.283185307179586 * U(g);
      Emit(k * 0.3 + r * std::cos(ph), std::clamp(r * std::sin(ph), -0.29, 0.29), z + 0.002 * (U(g) - .5),
           t + 0.2 + (U(g) < .5 ? 0.3 * U(g) : -1.5 * std::log(1 - U(g))));
    }
    if (primaryCounter == 1 && std::getenv("K2A_STUB_INSIDE_ANODE")) Emit(k * 0.3 + 0.0006, 0.0008, z, t + 1);
    if (primaryCounter == 1 && std::getenv("K2A_STUB_OUTSIDE_PATCH")) Emit(5 * 0.3, 0.05, z, t + 1);
    for (int j = 0; j <= ni_; ++j) {
      const double ph = 6.283185307179586 * U(g);
      ends_.push_back({x, y, z, t, 0.1, k * 0.3 + 0.0015 * std::cos(ph), 0.0015 * std::sin(ph), z, t + 1 + U(g), 1.});
    }
    (void)y;
    return true;
  }
  void GetAvalancheSize(int& ne, int& ni) const { ne = int(ends_.size()); ni = ni_; }
  std::size_t GetNumberOfElectronEndpoints() const { return ends_.size(); }
  void GetElectronEndpoint(std::size_t j, double& x0, double& y0, double& z0, double& t0, double& e0, double& x1, double& y1,
                           double& z1, double& t1, double& e1, int& status) const {
    const auto& v = ends_[j];
    x0 = v[0]; y0 = v[1]; z0 = v[2]; t0 = v[3]; e0 = v[4]; x1 = v[5]; y1 = v[6]; z1 = v[7]; t1 = v[8]; e1 = v[9]; status = -5;
  }
 private:
  void Emit(double u, double v, double w, double t) { ++ni_; if (h_) h_(u, v, w, t, 0, 0, nullptr); }
  Handler h_ = nullptr;
  int ni_ = 0;
  std::vector<std::vector<double>> ends_;
};

class DriftLineRKF {
 public:
  DriftLineRKF() = default;
  explicit DriftLineRKF(Sensor*) {}
  template <class... A> void UseWeightingPotential(A&&...) {}
  template <class... A> void SetSensor(A&&...) {}
  template <class... A> void SetMaximumStepSize(A&&...) {}
  template <class... A> void EnableSignalCalculation(A&&...) {}
  template <class... A> void SetSignalAveragingOrder(A&&...) {}
  template <class... A> bool DriftIon(A&&...) { return false; }
  template <class... A> void GetEndPoint(A&...) {}
};
}  // namespace Garfield

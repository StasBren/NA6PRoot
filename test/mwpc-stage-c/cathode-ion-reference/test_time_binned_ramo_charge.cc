// Pure C++17 unit test: no Garfield, ROOT, meshes or muon input.
#include "TimeBinnedRamoCharge.hh"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void Same(double value, double expected, double tol = 1.e-11) {
  if (std::abs(value - expected) > tol)
    throw std::runtime_error("Ramo bin numerical check failed");
}
void Check(const NA60Stereo::TimeBinnedRamoCharge& p,
           const std::vector<double>& expected, double total, double tail) {
  const auto got = p.BinChargeFc();
  if (got.size() != expected.size()) throw std::runtime_error("bin length");
  for (std::size_t i = 0; i < got.size(); ++i) Same(got[i], expected[i]);
  Same(p.FullEndpointFc(), total);
  Same(p.TailFc(), tail);
  if (!p.CheckClosure(1.e-10)) throw std::runtime_error("Ramo closure failed");
}
}

int main() {
  try {
    {
      NA60Stereo::TimeBinnedRamoCharge p(10., 4);
      p.Add(3., 23., 20.);
      Check(p, {7., 10., 3., 0.}, 20., 0.);
    }
    {
      NA60Stereo::TimeBinnedRamoCharge p(10., 4);
      p.Add(20., 70., 50.);
      Check(p, {0., 0., 10., 10.}, 50., 30.);
    }
    {
      NA60Stereo::TimeBinnedRamoCharge p(10., 4);
      p.Add(0., 10., +10.);
      p.Add(10., 20., -4.);
      p.Add(45., 55., 20.);
      Check(p, {10., -4., 0., 0.}, 26., 20.);
    }
    {
      // Several RKF vertices, with a single track crossing many time bins.
      NA60Stereo::TimeBinnedRamoCharge p(10., 4);
      p.Add(0., 5., 4.);
      p.Add(5., 15., 10.);
      p.Add(15., 35., 20.);
      Check(p, {9., 10., 10., 5.}, 34., 0.);
    }
    std::cout << "PASS time-binned Ramo: partial bins, signed currents, "
                 "late arrivals, endpoint closure\n";
  } catch (const std::exception& error) {
    std::cerr << "FAIL time-binned Ramo: " << error.what() << "\n";
    return 2;
  }
  return 0;
}

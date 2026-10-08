#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Medium.hh"
#include "Garfield/Sensor.hh"

namespace na6p::mwpc {

constexpr double ElementaryChargeFc = 1.602176634e-4;

struct IonBirth {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
};

inline std::string StripLabel(const int i) {
  if (i < 0) return "strip_m" + std::to_string(-i);
  if (i > 0) return "strip_p" + std::to_string(i);
  return "strip_0";
}

// Endpoint form of Shockley-Ramo used throughout the Stage-B studies.
// qFc carries the charge sign; weighting potential is dimensionless.
inline double EndpointSignal(
    Garfield::ComponentAnalyticField& weighting,
    const std::string& label,
    const double qFc,
    const double u0, const double v0, const double w0,
    const double u1, const double v1, const double w1) {
  const double phi0 =
      weighting.WeightingPotential(u0, v0, w0, label);
  const double phi1 =
      weighting.WeightingPotential(u1, v1, w1, label);
  return qFc * (phi1 - phi0);
}

inline bool IonVelocity(
    Garfield::Sensor& sensor,
    const std::array<double, 3>& x,
    std::array<double, 3>& velocity) {
  double eu = 0., ev = 0., ew = 0., potential = 0.;
  Garfield::Medium* medium = nullptr;
  int status = 0;

  sensor.ElectricField(
      x[0], x[1], x[2],
      eu, ev, ew, potential, medium, status);
  if (status != 0 || !medium) return false;

  double bu = 0., bv = 0., bw = 0.;
  sensor.MagneticField(
      x[0], x[1], x[2],
      bu, bv, bw, status);
  if (status != 0) bu = bv = bw = 0.;

  return medium->IonVelocity(
      eu, ev, ew, bu, bv, bw,
      velocity[0], velocity[1], velocity[2]);
}

// Advance an avalanche ion only until the chosen observation time.
// This is the same local RK4 construction used by the late Stage-B
// finite-observation-time studies. Coordinates are (u,v,w), time is ns.
inline bool AdvanceIonRK4(
    Garfield::Sensor& sensor,
    std::array<double, 3>& x,
    double& t,
    const double targetTimeNs,
    const double maxDtNs) {
  if (targetTimeNs <= t) return true;

  while (t < targetTimeNs) {
    const double h = std::min(maxDtNs, targetTimeNs - t);

    std::array<double, 3> k1{}, k2{}, k3{}, k4{};
    if (!IonVelocity(sensor, x, k1)) return false;

    std::array<double, 3> x2 = x;
    for (int j = 0; j < 3; ++j) x2[j] += 0.5 * h * k1[j];
    if (!IonVelocity(sensor, x2, k2)) return false;

    std::array<double, 3> x3 = x;
    for (int j = 0; j < 3; ++j) x3[j] += 0.5 * h * k2[j];
    if (!IonVelocity(sensor, x3, k3)) return false;

    std::array<double, 3> x4 = x;
    for (int j = 0; j < 3; ++j) x4[j] += h * k3[j];
    if (!IonVelocity(sensor, x4, k4)) return false;

    for (int j = 0; j < 3; ++j) {
      x[j] += h * (k1[j] + 2. * k2[j] + 2. * k3[j] + k4[j]) / 6.;
    }
    t += h;

    if (!sensor.IsInside(x[0], x[1], x[2])) return false;
  }

  return true;
}

}  // namespace na6p::mwpc

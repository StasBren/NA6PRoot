// Numerical core for time-binned Shockley-Ramo *electrode* charge.
// Units: t ns; charge fC. No Garfield/ROOT dependency.
// The project follows Garfield Sensor's positive-ion convention Q=q*(phi1-phi0).
// Piecewise-linear phi(t) between supplied RKF vertices is an approximation:
// refine trajectory sampling to test early-time waveform stability.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace NA60Stereo {

class TimeBinnedRamoCharge {
 public:
  TimeBinnedRamoCharge(double dtNs, int bins)
      : dt_(dtNs), bin_(bins, 0.), fullDiff_(bins + 1, 0.) {
    if (!(std::isfinite(dt_) && dt_ > 0.) || bins < 1 || bins > 50000)
      throw std::invalid_argument("invalid Ramo time grid");
  }

  // Add the weighted induced charge of a single trajectory segment.
  // The charge rate is constant across this segment under linear phi(t).
  // Work scales with segment count, NOT with 10,000 bins per segment.
  void Add(double t0, double t1, double deltaChargeFc) {
    if (!std::isfinite(t0) || !std::isfinite(t1) ||
        !std::isfinite(deltaChargeFc) || !(t1 > t0))
      throw std::invalid_argument("invalid Ramo segment (time/charge)");
    endpointFc_ += deltaChargeFc;
    const double start = std::max(0., t0);
    const double end = std::min(WindowEndNs(), t1);
    if (!(end > start)) return;
    const double chargeRate = deltaChargeFc / (t1 - t0);
    const int first = std::clamp(int(std::floor(start / dt_)), 0, int(bin_.size()) - 1);
    // Use the last point strictly below end for exact bin-boundary semantics.
    const double lastT = std::nextafter(end, -std::numeric_limits<double>::infinity());
    const int last = std::clamp(int(std::floor(lastT / dt_)), 0, int(bin_.size()) - 1);
    if (first == last) {
      bin_[first] += chargeRate * (end - start);
    } else {
      bin_[first] += chargeRate * ((first + 1) * dt_ - start);
      bin_[last] += chargeRate * (end - last * dt_);
      if (last > first + 1) {
        fullDiff_[first + 1] += chargeRate * dt_;
        fullDiff_[last] -= chargeRate * dt_;
      }
    }
    windowFc_ += chargeRate * (end - start);
  }

  std::vector<double> BinChargeFc() const {
    auto out = bin_;
    double carry = 0.;
    for (std::size_t j = 0; j < out.size(); ++j) {
      carry += fullDiff_[j];
      out[j] += carry;
    }
    return out;
  }

  double FullEndpointFc() const { return endpointFc_; }
  double WindowChargeFc() const { return windowFc_; }
  double TailFc() const { return endpointFc_ - windowFc_; }
  double WindowEndNs() const { return dt_ * bin_.size(); }
  double DtNs() const { return dt_; }
  int Bins() const { return int(bin_.size()); }

  bool CheckClosure(double absTolFc, double relTol = 1.e-9) const {
    const auto charges = BinChargeFc();
    double sum = 0.;
    for (double q : charges) sum += q;
    const double tol = absTolFc +
        relTol * std::max(std::abs(windowFc_), std::abs(endpointFc_));
    return std::abs(sum - windowFc_) <= tol &&
           std::abs(windowFc_ + TailFc() - endpointFc_) <= tol;
  }

 private:
  double dt_;
  std::vector<double> bin_;
  std::vector<double> fullDiff_;
  double endpointFc_ = 0.;
  double windowFc_ = 0.;
};
}  // namespace NA60Stereo

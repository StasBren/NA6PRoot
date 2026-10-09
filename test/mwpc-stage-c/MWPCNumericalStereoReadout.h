#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentElmer.hh"

#include "../mwpc-avalanche/MWPCSignalUtils.h"

namespace na6p::mwpc {

// Numerical Stage-C2 stereo readout.
//
// A single FEM weighting map is solved for the central strip on each cathode.
// Because the ideal local chamber geometry (parallel cathodes + wires running
// along w) is translationally invariant along w, neighbouring parallel strips
// are obtained by translating the query point:
//
//   phi_k(u,v,w) = phi_0(u,v,w - k*p/cos(alpha)).
//
// This avoids one Elmer solve per strip while preserving the explicit anode
// wire conductor boundary present in the FEM map.  The finite FEM patch is
// still an approximation; domain misses are counted and invalidate the event.
class NumericalStereoReadout {
 public:
  struct Result {
    std::vector<double> qMinusFc;
    std::vector<double> qPlusFc;

    double qCathodeMinusFc = 0.;
    double qCathodePlusFc = 0.;
    double qMinusLocalSumFc = 0.;
    double qPlusLocalSumFc = 0.;

    // Kept for CSV compatibility with the old hybrid readout.
    // Numerical strip charges are direct Shockley-Ramo charges and are never
    // rescaled, so both are identically one.
    double minusScale = 1.;
    double plusScale = 1.;

    double minusCoverage = std::numeric_limits<double>::quiet_NaN();
    double plusCoverage = std::numeric_limits<double>::quiet_NaN();

    double xPlusAlphaCogMm = std::numeric_limits<double>::quiet_NaN();
    double xMinusAlphaCogMm = std::numeric_limits<double>::quiet_NaN();
    double uCogMm = std::numeric_limits<double>::quiet_NaN();
    double wCogMm = std::numeric_limits<double>::quiet_NaN();

    int segmentedMinusSegments = 0;
    int segmentedPlusSegments = 0;
    int skippedMinusSegments = 0;
    int skippedPlusSegments = 0;
    long long nonFiniteWeightingQueries = 0;

    std::string invalidReason;
    bool valid = false;
  };

  NumericalStereoReadout(
      Garfield::ComponentElmer& weighting,
      const double gapMinusCm,
      const double gapPlusCm,
      const double stripPitchCm,
      const double tanAlpha,
      const int halfStrips,
      const double mapUHalfSpanCm,
      const double mapWHalfSpanCm,
      const std::string& minusCentralLabel = "minus_strip_0",
      const std::string& plusCentralLabel = "plus_strip_0")
      : mWeighting(weighting),
        mGapMinusCm(gapMinusCm),
        mGapPlusCm(gapPlusCm),
        mStripPitchCm(stripPitchCm),
        mTanAlpha(tanAlpha),
        mHalfStrips(halfStrips),
        mMapUHalfSpanCm(mapUHalfSpanCm),
        mMapWHalfSpanCm(mapWHalfSpanCm),
        mMinusCentralLabel(minusCentralLabel),
        mPlusCentralLabel(plusCentralLabel),
        mAlpha(std::atan(tanAlpha)),
        mCosA(std::cos(mAlpha)),
        mSinA(std::sin(mAlpha)) {
    if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
        stripPitchCm <= 0. || tanAlpha <= 0. ||
        halfStrips < 1 || mapUHalfSpanCm <= 0. ||
        mapWHalfSpanCm <= 0. || std::abs(mCosA) < 1.e-12) {
      throw std::invalid_argument(
          "Invalid numerical stereo-readout geometry.");
    }

    mStripIds.reserve(2 * mHalfStrips + 1);
    for (int k = -mHalfStrips; k <= mHalfStrips; ++k) {
      mStripIds.push_back(k);
    }
  }

  int HalfStrips() const { return mHalfStrips; }
  double StripPitchCm() const { return mStripPitchCm; }
  double AlphaRad() const { return mAlpha; }
  double AlphaDeg() const { return 180. * mAlpha / Pi(); }
  const std::vector<int>& StripIds() const { return mStripIds; }

  const char* SharingDomainName() const {
    return "numerical-elmer-translated-central-strip";
  }

  Result EmptyResult() const {
    Result r;
    r.qMinusFc.assign(mStripIds.size(), 0.);
    r.qPlusFc.assign(mStripIds.size(), 0.);
    return r;
  }

  void AddSegment(
      Result& r,
      Garfield::ComponentAnalyticField& physicalField,
      const double qFc,
      const double u0, const double v0, const double w0,
      const double u1, const double v1, const double w1) {
    if (r.qMinusFc.size() != mStripIds.size() ||
        r.qPlusFc.size() != mStripIds.size()) {
      throw std::invalid_argument(
          "NumericalStereoReadout::Result was not created with EmptyResult().");
    }

    r.qCathodeMinusFc += EndpointSignal(
        physicalField, "cathode_minus", qFc,
        u0, v0, w0, u1, v1, w1);
    r.qCathodePlusFc += EndpointSignal(
        physicalField, "cathode_plus", qFc,
        u0, v0, w0, u1, v1, w1);

    for (std::size_t j = 0; j < mStripIds.size(); ++j) {
      const int k = mStripIds[j];
      const double wShiftCm =
          k * mStripPitchCm / mCosA;

      if (InsideMap(u0, v0, w0 - wShiftCm) &&
          InsideMap(u1, v1, w1 - wShiftCm)) {
        const double phi0 = mWeighting.WeightingPotential(
            u0, v0, w0 - wShiftCm, mMinusCentralLabel);
        const double phi1 = mWeighting.WeightingPotential(
            u1, v1, w1 - wShiftCm, mMinusCentralLabel);

        if (std::isfinite(phi0) && std::isfinite(phi1)) {
          r.qMinusFc[j] += qFc * (phi1 - phi0);
        } else {
          ++r.nonFiniteWeightingQueries;
        }
      } else {
        ++r.skippedMinusSegments;
      }

      if (InsideMap(u0, v0, w0 - wShiftCm) &&
          InsideMap(u1, v1, w1 - wShiftCm)) {
        const double phi0 = mWeighting.WeightingPotential(
            u0, v0, w0 - wShiftCm, mPlusCentralLabel);
        const double phi1 = mWeighting.WeightingPotential(
            u1, v1, w1 - wShiftCm, mPlusCentralLabel);

        if (std::isfinite(phi0) && std::isfinite(phi1)) {
          r.qPlusFc[j] += qFc * (phi1 - phi0);
        } else {
          ++r.nonFiniteWeightingQueries;
        }
      } else {
        ++r.skippedPlusSegments;
      }
    }

    ++r.segmentedMinusSegments;
    ++r.segmentedPlusSegments;
  }

  bool Finalize(Result& r) const {
    r.qMinusLocalSumFc = 0.;
    r.qPlusLocalSumFc = 0.;
    for (const double q : r.qMinusFc) r.qMinusLocalSumFc += q;
    for (const double q : r.qPlusFc) r.qPlusLocalSumFc += q;

    constexpr double MinDenom = 1.e-18;
    if (std::abs(r.qCathodeMinusFc) > MinDenom) {
      r.minusCoverage =
          r.qMinusLocalSumFc / r.qCathodeMinusFc;
    }
    if (std::abs(r.qCathodePlusFc) > MinDenom) {
      r.plusCoverage =
          r.qPlusLocalSumFc / r.qCathodePlusFc;
    }

    if (r.nonFiniteWeightingQueries > 0) {
      r.invalidReason = "non-finite numerical weighting-potential query";
      r.valid = false;
      return false;
    }

    if (r.skippedMinusSegments > 0 ||
        r.skippedPlusSegments > 0) {
      r.invalidReason = "translated strip query outside FEM map";
      r.valid = false;
      return false;
    }

    r.xPlusAlphaCogMm = CogMm(r.qMinusFc);
    r.xMinusAlphaCogMm = CogMm(r.qPlusFc);

    if (!std::isfinite(r.qMinusLocalSumFc) ||
        !std::isfinite(r.qPlusLocalSumFc) ||
        !std::isfinite(r.xPlusAlphaCogMm) ||
        !std::isfinite(r.xMinusAlphaCogMm) ||
        std::abs(mSinA) < 1.e-12 ||
        std::abs(mCosA) < 1.e-12) {
      r.invalidReason = "non-finite numerical strip response or CoG";
      r.valid = false;
      return false;
    }

    r.wCogMm =
        (r.xPlusAlphaCogMm + r.xMinusAlphaCogMm) / (2. * mCosA);
    r.uCogMm =
        (r.xMinusAlphaCogMm - r.xPlusAlphaCogMm) / (2. * mSinA);

    r.invalidReason.clear();
    r.valid = true;
    return true;
  }

 private:
  static constexpr double Pi() {
    return 3.14159265358979323846;
  }

  bool InsideMap(
      const double uCm,
      const double vCm,
      const double wCm) const {
    // Stay a little away from the artificial FEM side boundary.  The cathode
    // surfaces themselves are physical and may be approached by transported
    // ions, hence no analogous v guard is imposed.
    constexpr double SideGuardCm = 1.e-4;
    constexpr double VGuardCm = 5.e-7;

    return std::abs(uCm) <
               mMapUHalfSpanCm - SideGuardCm &&
           std::abs(wCm) <
               mMapWHalfSpanCm - SideGuardCm &&
           vCm >= -mGapMinusCm - VGuardCm &&
           vCm <= +mGapPlusCm + VGuardCm;
  }

  double CogMm(const std::vector<double>& q) const {
    double denom = 0.;
    double numerCm = 0.;

    for (std::size_t j = 0; j < q.size(); ++j) {
      const double amplitude = std::abs(q[j]);
      const double centerCm =
          mStripIds[j] * mStripPitchCm;
      denom += amplitude;
      numerCm += amplitude * centerCm;
    }

    if (denom <= 0.) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return 10. * numerCm / denom;
  }

  Garfield::ComponentElmer& mWeighting;

  double mGapMinusCm = 0.;
  double mGapPlusCm = 0.;
  double mStripPitchCm = 0.;
  double mTanAlpha = 0.;
  int mHalfStrips = 0;
  double mMapUHalfSpanCm = 0.;
  double mMapWHalfSpanCm = 0.;

  std::string mMinusCentralLabel;
  std::string mPlusCentralLabel;

  double mAlpha = 0.;
  double mCosA = 1.;
  double mSinA = 0.;

  std::vector<int> mStripIds;
};

}  // namespace na6p::mwpc

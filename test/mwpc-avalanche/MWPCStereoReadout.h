#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"

#include "MWPCSignalUtils.h"

namespace na6p::mwpc {

enum class StereoSharingDomain {
  // Stage-C default: the analytic segmented-strip weighting problem spans the
  // full cathode-to-cathode gas gap. This lets either cathode respond to
  // carriers anywhere in the chamber.
  FullCathodeGap,

  // Historical Stage-B construction: each segmented weighting problem spans
  // only wire plane -> selected cathode. Kept for regression/comparison.
  LegacyHalfGap
};

// Reusable form of the Stage-B3c.2 hybrid stereo model, with a Stage-C
// full-gap extension for complete muon tracks.
//
// Physical total coupling:
//   full analytic MWPC cell, including wires, via complete cathode labels.
//
// Segmented sharing:
//   ideal planar-strip weighting field on each cathode.
//
// The local segmented charges are rescaled side-by-side so their signed sum
// matches the corresponding wire-aware complete-cathode signal.
//
// Stereo convention:
//   minus cathode: +alpha family, x_plus  = w cos(alpha) - u sin(alpha)
//   plus  cathode: -alpha family, x_minus = w cos(alpha) + u sin(alpha)
class StereoReadout {
 public:
  struct Result {
    std::vector<double> qMinusFc;
    std::vector<double> qPlusFc;

    double qCathodeMinusFc = 0.;
    double qCathodePlusFc = 0.;
    double qMinusLocalSumFc = 0.;
    double qPlusLocalSumFc = 0.;
    double minusScale = 0.;
    double plusScale = 0.;

    double xPlusAlphaCogMm = std::numeric_limits<double>::quiet_NaN();
    double xMinusAlphaCogMm = std::numeric_limits<double>::quiet_NaN();
    double uCogMm = std::numeric_limits<double>::quiet_NaN();
    double wCogMm = std::numeric_limits<double>::quiet_NaN();

    int segmentedMinusSegments = 0;
    int segmentedPlusSegments = 0;
    int skippedMinusSegments = 0;
    int skippedPlusSegments = 0;

    std::string invalidReason;
    bool valid = false;
  };

  StereoReadout(
      const double gapMinusCm,
      const double gapPlusCm,
      const double stripPitchCm,
      const double stripWidthCm,
      const double tanAlpha,
      const int halfStrips,
      const StereoSharingDomain sharingDomain =
          StereoSharingDomain::FullCathodeGap)
      : mGapMinusCm(gapMinusCm),
        mGapPlusCm(gapPlusCm),
        mStripPitchCm(stripPitchCm),
        mStripWidthCm(stripWidthCm),
        mTanAlpha(tanAlpha),
        mHalfStrips(halfStrips),
        mSharingDomain(sharingDomain),
        mAlpha(std::atan(tanAlpha)),
        mCosA(std::cos(mAlpha)),
        mSinA(std::sin(mAlpha)) {
    if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
        stripPitchCm <= 0. || stripWidthCm <= 0. ||
        stripWidthCm > stripPitchCm + 1.e-12 ||
        tanAlpha <= 0. || halfStrips < 2) {
      throw std::invalid_argument("Invalid stereo-readout geometry.");
    }

    if (mSharingDomain == StereoSharingDomain::FullCathodeGap) {
      const double fullGapCm = mGapMinusCm + mGapPlusCm;

      // Stage-C extension of the same analytic strip model:
      // both real cathodes are the bounding planes, so the strip weighting
      // potential is defined throughout the complete gas volume.
      //
      // The auxiliary plane voltages only establish a valid analytic cell;
      // AddStripOnPlaneY constructs the actual unit weighting field.
      mWeightingMinus.AddPlaneY(
          -mGapMinusCm, 0., "minus_front");
      mWeightingMinus.AddPlaneY(
          +mGapPlusCm, 1., "minus_back");

      mWeightingPlus.AddPlaneY(
          -mGapMinusCm, 1., "plus_back");
      mWeightingPlus.AddPlaneY(
          +mGapPlusCm, 0., "plus_front");

      mStripIds.reserve(2 * mHalfStrips + 1);
      mLabelsMinus.reserve(2 * mHalfStrips + 1);
      mLabelsPlus.reserve(2 * mHalfStrips + 1);

      for (int k = -mHalfStrips; k <= mHalfStrips; ++k) {
        const double center = k * mStripPitchCm;
        const double sMin = center - 0.5 * mStripWidthCm;
        const double sMax = center + 0.5 * mStripWidthCm;

        const std::string lm = FamilyLabel("minus_strip", k);
        const std::string lp = FamilyLabel("plus_strip", k);

        mWeightingMinus.AddStripOnPlaneY(
            'x', -mGapMinusCm, sMin, sMax, lm, fullGapCm);
        mWeightingPlus.AddStripOnPlaneY(
            'x', +mGapPlusCm, sMin, sMax, lp, fullGapCm);

        mStripIds.push_back(k);
        mLabelsMinus.push_back(lm);
        mLabelsPlus.push_back(lp);
      }
    } else {
      // Exact historical Stage-B half-gap construction.
      mWeightingMinus.AddPlaneY(0., 1., "minus_back");
      mWeightingMinus.AddPlaneY(-mGapMinusCm, 0., "minus_front");

      // Mirror +v into local y with readout plane at -gap.
      mWeightingPlus.AddPlaneY(0., 1., "plus_back");
      mWeightingPlus.AddPlaneY(-mGapPlusCm, 0., "plus_front");

      mStripIds.reserve(2 * mHalfStrips + 1);
      mLabelsMinus.reserve(2 * mHalfStrips + 1);
      mLabelsPlus.reserve(2 * mHalfStrips + 1);

      for (int k = -mHalfStrips; k <= mHalfStrips; ++k) {
        const double center = k * mStripPitchCm;
        const double sMin = center - 0.5 * mStripWidthCm;
        const double sMax = center + 0.5 * mStripWidthCm;

        const std::string lm = FamilyLabel("minus_strip", k);
        const std::string lp = FamilyLabel("plus_strip", k);

        mWeightingMinus.AddStripOnPlaneY(
            'x', -mGapMinusCm, sMin, sMax, lm, mGapMinusCm);
        mWeightingPlus.AddStripOnPlaneY(
            'x', -mGapPlusCm, sMin, sMax, lp, mGapPlusCm);

        mStripIds.push_back(k);
        mLabelsMinus.push_back(lm);
        mLabelsPlus.push_back(lp);
      }
    }
  }

  int HalfStrips() const { return mHalfStrips; }
  double StripPitchCm() const { return mStripPitchCm; }
  double AlphaRad() const { return mAlpha; }
  double AlphaDeg() const { return 180. * mAlpha / Pi(); }
  const std::vector<int>& StripIds() const { return mStripIds; }

  const char* SharingDomainName() const {
    return mSharingDomain == StereoSharingDomain::FullCathodeGap
               ? "full-cathode-gap"
               : "legacy-half-gap";
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
          "StereoReadout::Result was not created with EmptyResult().");
    }

    r.qCathodeMinusFc += EndpointSignal(
        physicalField, "cathode_minus", qFc,
        u0, v0, w0, u1, v1, w1);
    r.qCathodePlusFc += EndpointSignal(
        physicalField, "cathode_plus", qFc,
        u0, v0, w0, u1, v1, w1);

    const double zMinus0 = w0 * mCosA - u0 * mSinA;
    const double zMinus1 = w1 * mCosA - u1 * mSinA;
    const double zPlus0 = w0 * mCosA + u0 * mSinA;
    const double zPlus1 = w1 * mCosA + u1 * mSinA;

    if (mSharingDomain == StereoSharingDomain::FullCathodeGap) {
      for (std::size_t j = 0; j < mStripIds.size(); ++j) {
        r.qMinusFc[j] += EndpointSignal(
            mWeightingMinus, mLabelsMinus[j], qFc,
            0., v0, zMinus0,
            0., v1, zMinus1);

        r.qPlusFc[j] += EndpointSignal(
            mWeightingPlus, mLabelsPlus[j], qFc,
            0., v0, zPlus0,
            0., v1, zPlus1);
      }
      ++r.segmentedMinusSegments;
      ++r.segmentedPlusSegments;
    } else {
      // Historical half-gap weighting problems are not defined beyond v=0.
      // Clip only in this legacy mode to avoid evaluating them outside their
      // analytic domain.
      Segment clipped;

      if (ClipToVRange(
              u0, v0, w0, u1, v1, w1,
              -mGapMinusCm, 0., clipped)) {
        const double zm0 =
            clipped.w0 * mCosA - clipped.u0 * mSinA;
        const double zm1 =
            clipped.w1 * mCosA - clipped.u1 * mSinA;

        for (std::size_t j = 0; j < mStripIds.size(); ++j) {
          r.qMinusFc[j] += EndpointSignal(
              mWeightingMinus, mLabelsMinus[j], qFc,
              0., clipped.v0, zm0,
              0., clipped.v1, zm1);
        }
        ++r.segmentedMinusSegments;
      } else {
        ++r.skippedMinusSegments;
      }

      if (ClipToVRange(
              u0, v0, w0, u1, v1, w1,
              0., mGapPlusCm, clipped)) {
        const double zp0 =
            clipped.w0 * mCosA + clipped.u0 * mSinA;
        const double zp1 =
            clipped.w1 * mCosA + clipped.u1 * mSinA;

        for (std::size_t j = 0; j < mStripIds.size(); ++j) {
          r.qPlusFc[j] += EndpointSignal(
              mWeightingPlus, mLabelsPlus[j], qFc,
              0., -clipped.v0, zp0,
              0., -clipped.v1, zp1);
        }
        ++r.segmentedPlusSegments;
      } else {
        ++r.skippedPlusSegments;
      }
    }
  }

  bool Finalize(Result& r) const {
    r.qMinusLocalSumFc = 0.;
    r.qPlusLocalSumFc = 0.;
    for (const double q : r.qMinusFc) r.qMinusLocalSumFc += q;
    for (const double q : r.qPlusFc) r.qPlusLocalSumFc += q;

    constexpr double MinDenom = 1.e-18;
    if (!std::isfinite(r.qCathodeMinusFc) ||
        !std::isfinite(r.qCathodePlusFc) ||
        !std::isfinite(r.qMinusLocalSumFc) ||
        !std::isfinite(r.qPlusLocalSumFc) ||
        std::abs(r.qMinusLocalSumFc) < MinDenom ||
        std::abs(r.qPlusLocalSumFc) < MinDenom) {
      r.invalidReason = "non-finite or vanishing local strip sum";
      r.valid = false;
      return false;
    }

    r.minusScale = r.qCathodeMinusFc / r.qMinusLocalSumFc;
    r.plusScale = r.qCathodePlusFc / r.qPlusLocalSumFc;

    for (double& q : r.qMinusFc) q *= r.minusScale;
    for (double& q : r.qPlusFc) q *= r.plusScale;

    r.xPlusAlphaCogMm = CogMm(r.qMinusFc);
    r.xMinusAlphaCogMm = CogMm(r.qPlusFc);

    if (!std::isfinite(r.xPlusAlphaCogMm) ||
        !std::isfinite(r.xMinusAlphaCogMm) ||
        std::abs(mSinA) < 1.e-12 ||
        std::abs(mCosA) < 1.e-12) {
      r.invalidReason = "non-finite projected CoG";
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
  struct Segment {
    double u0 = 0., v0 = 0., w0 = 0.;
    double u1 = 0., v1 = 0., w1 = 0.;
  };

  static bool ClipToVRange(
      const double u0, const double v0, const double w0,
      const double u1, const double v1, const double w1,
      const double vMin, const double vMax,
      Segment& out) {
    if (vMin > vMax) return false;

    constexpr double Eps = 1.e-14;
    const double dv = v1 - v0;

    double tLo = 0.;
    double tHi = 1.;

    if (std::abs(dv) < Eps) {
      if (v0 < vMin - Eps || v0 > vMax + Eps) return false;
    } else {
      double ta = (vMin - v0) / dv;
      double tb = (vMax - v0) / dv;
      if (ta > tb) std::swap(ta, tb);

      tLo = std::max(0., ta);
      tHi = std::min(1., tb);
      if (tHi < tLo - Eps) return false;
    }

    auto lerp = [](const double a, const double b, const double t) {
      return a + t * (b - a);
    };

    out.u0 = lerp(u0, u1, tLo);
    out.v0 = std::clamp(lerp(v0, v1, tLo), vMin, vMax);
    out.w0 = lerp(w0, w1, tLo);

    out.u1 = lerp(u0, u1, tHi);
    out.v1 = std::clamp(lerp(v0, v1, tHi), vMin, vMax);
    out.w1 = lerp(w0, w1, tHi);

    return true;
  }

  static constexpr double Pi() {
    return 3.14159265358979323846;
  }

  static std::string FamilyLabel(
      const std::string& prefix, const int k) {
    if (k < 0) return prefix + "_m" + std::to_string(-k);
    if (k > 0) return prefix + "_p" + std::to_string(k);
    return prefix + "_0";
  }

  double CogMm(const std::vector<double>& q) const {
    double denom = 0.;
    double numerCm = 0.;

    for (std::size_t j = 0; j < q.size(); ++j) {
      const double amplitude = std::abs(q[j]);
      const double centerCm = mStripIds[j] * mStripPitchCm;
      denom += amplitude;
      numerCm += amplitude * centerCm;
    }

    if (denom <= 0.) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return 10. * numerCm / denom;
  }

  double mGapMinusCm = 0.;
  double mGapPlusCm = 0.;
  double mStripPitchCm = 0.;
  double mStripWidthCm = 0.;
  double mTanAlpha = 0.;
  int mHalfStrips = 0;
  StereoSharingDomain mSharingDomain =
      StereoSharingDomain::FullCathodeGap;

  double mAlpha = 0.;
  double mCosA = 1.;
  double mSinA = 0.;

  Garfield::ComponentAnalyticField mWeightingMinus;
  Garfield::ComponentAnalyticField mWeightingPlus;

  std::vector<int> mStripIds;
  std::vector<std::string> mLabelsMinus;
  std::vector<std::string> mLabelsPlus;
};

}  // namespace na6p::mwpc

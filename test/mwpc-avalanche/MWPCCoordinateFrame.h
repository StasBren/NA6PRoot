#pragma once

#include <cmath>
#include <stdexcept>

namespace na6p::mwpc {

struct GlobalXYZ {
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

struct LocalUVW {
  double u = 0.;
  double v = 0.;
  double w = 0.;
};

inline GlobalXYZ operator+(const GlobalXYZ& a, const GlobalXYZ& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline GlobalXYZ operator-(const GlobalXYZ& a, const GlobalXYZ& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline GlobalXYZ operator*(const double s, const GlobalXYZ& a) {
  return {s * a.x, s * a.y, s * a.z};
}

inline double Dot(const GlobalXYZ& a, const GlobalXYZ& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline GlobalXYZ Cross(const GlobalXYZ& a, const GlobalXYZ& b) {
  return {
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x};
}

inline double Norm(const GlobalXYZ& a) {
  return std::sqrt(Dot(a, a));
}

// A rigid chamber frame.
//
// u: in the chamber plane, perpendicular to the anode wires
// v: chamber normal (wire-to-cathode / drift direction)
// w: in the chamber plane, along the anode wires
//
// The basis vectors are stored in GLOBAL coordinates.
// The constructor requires an orthonormal, right-handed frame:
//     u_hat x v_hat = w_hat.
//
// Garfield's analytic wire-cell coordinates will later be identified as
//     (x_G, y_G, z_G) = (u, v, w).
class ChamberFrame {
 public:
  ChamberFrame(const GlobalXYZ& originGlobal,
               const GlobalXYZ& uHatGlobal,
               const GlobalXYZ& vHatGlobal,
               const GlobalXYZ& wHatGlobal,
               const double tolerance = 1.e-10)
      : mOrigin(originGlobal),
        mU(uHatGlobal),
        mV(vHatGlobal),
        mW(wHatGlobal) {
    Validate(tolerance);
  }

  LocalUVW GlobalToLocalPoint(const GlobalXYZ& pGlobal) const {
    const auto d = pGlobal - mOrigin;
    return {Dot(d, mU), Dot(d, mV), Dot(d, mW)};
  }

  // Vectors (momentum, magnetic field, direction, ...) are rotated only.
  // They are not translated by the chamber origin.
  LocalUVW GlobalToLocalVector(const GlobalXYZ& aGlobal) const {
    return {Dot(aGlobal, mU), Dot(aGlobal, mV), Dot(aGlobal, mW)};
  }

  GlobalXYZ LocalToGlobalPoint(const LocalUVW& pLocal) const {
    return mOrigin + pLocal.u * mU + pLocal.v * mV + pLocal.w * mW;
  }

  GlobalXYZ LocalToGlobalVector(const LocalUVW& aLocal) const {
    return aLocal.u * mU + aLocal.v * mV + aLocal.w * mW;
  }

  const GlobalXYZ& OriginGlobal() const { return mOrigin; }
  const GlobalXYZ& UHatGlobal() const { return mU; }
  const GlobalXYZ& VHatGlobal() const { return mV; }
  const GlobalXYZ& WHatGlobal() const { return mW; }

 private:
  void Validate(const double tol) const {
    if (std::abs(Norm(mU) - 1.) > tol ||
        std::abs(Norm(mV) - 1.) > tol ||
        std::abs(Norm(mW) - 1.) > tol) {
      throw std::invalid_argument("ChamberFrame basis vectors must have unit length.");
    }

    if (std::abs(Dot(mU, mV)) > tol ||
        std::abs(Dot(mU, mW)) > tol ||
        std::abs(Dot(mV, mW)) > tol) {
      throw std::invalid_argument("ChamberFrame basis vectors must be orthogonal.");
    }

    const auto uv = Cross(mU, mV);
    const auto delta = uv - mW;
    if (Norm(delta) > tol) {
      throw std::invalid_argument(
          "ChamberFrame must be right-handed: u_hat x v_hat = w_hat.");
    }
  }

  GlobalXYZ mOrigin;
  GlobalXYZ mU;
  GlobalXYZ mV;
  GlobalXYZ mW;
};

// Nominal chamber orientation used by the current standalone sandbox:
//
// global X = horizontal = along wires       -> +w
// global Y = vertical   = across wires      -> +u
// global Z = beam       = chamber normal    -> +v
//
// Therefore:
//   u_hat = global +Y
//   v_hat = global +Z
//   w_hat = global +X
inline ChamberFrame MakeNominalNA60Frame(
    const GlobalXYZ& originGlobal = {}) {
  return ChamberFrame(
      originGlobal,
      {0., 1., 0.},  // +u = global +Y
      {0., 0., 1.},  // +v = global +Z
      {1., 0., 0.}); // +w = global +X
}

}  // namespace na6p::mwpc

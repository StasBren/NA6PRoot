#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "MWPCCoordinateFrame.h"

namespace {

bool Close(const double a, const double b, const double tol = 1.e-12) {
  return std::abs(a - b) < tol;
}

bool CheckLocal(const na6p::mwpc::LocalUVW& got,
                const na6p::mwpc::LocalUVW& expected,
                const std::string& label) {
  const bool ok = Close(got.u, expected.u) &&
                  Close(got.v, expected.v) &&
                  Close(got.w, expected.w);
  std::cout << (ok ? "PASS" : "FAIL") << "  " << label
            << "  got=(" << got.u << ", " << got.v << ", " << got.w << ")"
            << "  expected=(" << expected.u << ", " << expected.v << ", "
            << expected.w << ")\n";
  return ok;
}

bool CheckGlobal(const na6p::mwpc::GlobalXYZ& got,
                 const na6p::mwpc::GlobalXYZ& expected,
                 const std::string& label) {
  const bool ok = Close(got.x, expected.x) &&
                  Close(got.y, expected.y) &&
                  Close(got.z, expected.z);
  std::cout << (ok ? "PASS" : "FAIL") << "  " << label
            << "  got=(" << got.x << ", " << got.y << ", " << got.z << ")"
            << "  expected=(" << expected.x << ", " << expected.y << ", "
            << expected.z << ")\n";
  return ok;
}

}  // namespace

int main() {
  using namespace na6p::mwpc;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "\n=== MWPC COORDINATE-FRAME TEST ===\n";
  std::cout << "local u = across wires\n"
            << "local v = chamber normal / drift direction\n"
            << "local w = along wires\n"
            << "Garfield (x,y,z) will be identified with (u,v,w).\n\n";

  bool allOk = true;

  // Give the chamber a non-zero global origin so the test distinguishes
  // positions (which need translation) from vectors (which do not).
  const GlobalXYZ origin{10., 20., 30.};
  const auto frame = MakeNominalNA60Frame(origin);

  // This point was constructed as
  // origin + 1.2*u - 0.3*v + 2.5*w.
  const GlobalXYZ pointGlobal{12.5, 21.2, 29.7};
  const LocalUVW pointExpected{1.2, -0.3, 2.5};
  const auto pointLocal = frame.GlobalToLocalPoint(pointGlobal);
  allOk &= CheckLocal(pointLocal, pointExpected, "global point -> local");

  const auto pointRoundTrip = frame.LocalToGlobalPoint(pointLocal);
  allOk &= CheckGlobal(pointRoundTrip, pointGlobal,
                       "global point -> local -> global");

  // Nominal NA60 orientation:
  // global +Y (vertical B field) must be local +u.
  const GlobalXYZ bGlobal{0., 1., 0.};
  allOk &= CheckLocal(frame.GlobalToLocalVector(bGlobal),
                      {1., 0., 0.},
                      "global +Y magnetic field -> local +u");

  // Horizontal wire direction is global +X -> local +w.
  const GlobalXYZ wireGlobal{1., 0., 0.};
  allOk &= CheckLocal(frame.GlobalToLocalVector(wireGlobal),
                      {0., 0., 1.},
                      "global +X wire direction -> local +w");

  // Chamber normal / beam direction is global +Z -> local +v.
  const GlobalXYZ normalGlobal{0., 0., 1.};
  allOk &= CheckLocal(frame.GlobalToLocalVector(normalGlobal),
                      {0., 1., 0.},
                      "global +Z chamber normal -> local +v");

  // An arbitrary vector checks that vector transforms also round-trip.
  const GlobalXYZ vectorGlobal{0.3, -0.7, 1.1};
  const auto vectorLocal = frame.GlobalToLocalVector(vectorGlobal);
  const auto vectorRoundTrip = frame.LocalToGlobalVector(vectorLocal);
  allOk &= CheckGlobal(vectorRoundTrip, vectorGlobal,
                       "global vector -> local -> global");

  std::cout << "\n";
  if (!allOk) {
    std::cerr << "FAIL: coordinate-frame test failed.\n";
    return 2;
  }

  std::cout << "PASS: all coordinate-frame checks passed.\n";
  return 0;
}

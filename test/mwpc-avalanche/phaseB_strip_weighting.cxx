#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"

namespace {

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return std::atof(argv[i + 1]);
  }
  return defaultValue;
}

int ReadIntArg(const int argc, char** argv, const std::string& key,
               const int defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return std::atoi(argv[i + 1]);
  }
  return defaultValue;
}

std::string ReadStringArg(const int argc, char** argv, const std::string& key,
                          const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return argv[i + 1];
  }
  return defaultValue;
}

std::string StripLabel(const int i) {
  if (i < 0) return "strip_m" + std::to_string(-i);
  if (i > 0) return "strip_p" + std::to_string(i);
  return "strip_0";
}

}  // namespace

int main(int argc, char** argv) {
  // Chamber-local coordinates are mapped to Garfield as
  //   x = u : across wires
  //   y = v : chamber normal
  //   z = w : along wires
  //
  // The first Phase-B strip family is taken perpendicular to the anode wires,
  // i.e. the strips run along u (Garfield x) and are segmented in w
  // (Garfield z). This is deliberately the simplest idealised readout test.

  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double stripPitchCm =
      0.1 * ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthCm =
      0.1 * ReadArg(argc, argv, "--strip-width-mm", 1.7);
  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 8);
  const int scanSteps =
      ReadIntArg(argc, argv, "--scan-steps", 121);
  const double scanHalfRangeCm =
      0.1 * ReadArg(argc, argv, "--scan-half-range-mm", 8.0);
  const double uEvalCm =
      0.1 * ReadArg(argc, argv, "--u-mm", 0.0);
  const double vEvalCm =
      0.1 * ReadArg(argc, argv, "--v-mm", 0.0);
  const std::string side =
      ReadStringArg(argc, argv, "--side", "plus");
  const std::string outName =
      ReadStringArg(argc, argv, "--output", "phaseB_strip_weighting.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm ||
      halfStrips < 1 || scanSteps < 2 || scanHalfRangeCm <= 0.) {
    std::cerr << "Invalid Phase-B strip-weighting parameters.\n";
    return 2;
  }
  if (side != "plus" && side != "minus") {
    std::cerr << "--side must be 'plus' or 'minus'.\n";
    return 2;
  }

  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const double wireToCathodeGapCm =
      side == "plus" ? gapPlusCm : gapMinusCm;

  // This component is used only for weighting potentials, not for the
  // physical drift field. We intentionally use the analytic strip solution
  // with the wire-to-readout-cathode distance supplied explicitly.
  Garfield::ComponentAnalyticField weighting;
  weighting.AddPlaneY(planeY, 0., "");

  std::vector<std::string> labels;
  labels.reserve(2 * halfStrips + 1);

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double centerW = i * stripPitchCm;
    const double wMin = centerW - 0.5 * stripWidthCm;
    const double wMax = centerW + 0.5 * stripWidthCm;
    const std::string label = StripLabel(i);

    // direction='x': strip runs along Garfield x (= local u), while
    // smin/smax define its extent in Garfield z (= local w).
    weighting.AddStripOnPlaneY('x', planeY, wMin, wMax,
                               label, wireToCathodeGapCm);
    labels.push_back(label);
  }

  std::ofstream out(outName);
  out << "w_mm";
  for (int i = -halfStrips; i <= halfStrips; ++i) {
    out << "," << StripLabel(i);
  }
  out << ",sum_phi,central_fraction\n";

  std::cout << std::fixed << std::setprecision(5);
  std::cout << "\n=== PHASE B0 STRIP WEIGHTING-POTENTIAL PROFILE ===\n"
            << "readout side        : " << side << "\n"
            << "wire -> cathode gap : " << 10. * wireToCathodeGapCm << " mm\n"
            << "strip pitch         : " << 10. * stripPitchCm << " mm\n"
            << "strip width         : " << 10. * stripWidthCm << " mm\n"
            << "strip orientation   : along local u, segmented in local w\n"
            << "evaluation plane    : (u,v)=(" << 10. * uEvalCm << ", "
            << 10. * vEvalCm << ") mm\n"
            << "number of strips    : " << (2 * halfStrips + 1) << "\n\n"
            << "This is a weighting-field smoke test only; no avalanche or VMM"
               " response is included.\n\n"
            << "w[mm]    phi(strip0)    sum(phi)    central/sum\n";

  for (int iw = 0; iw < scanSteps; ++iw) {
    const double f = static_cast<double>(iw) / (scanSteps - 1);
    const double wCm = -scanHalfRangeCm +
                       2. * scanHalfRangeCm * f;

    double sumPhi = 0.;
    double phi0 = 0.;

    out << 10. * wCm;
    for (int i = -halfStrips; i <= halfStrips; ++i) {
      const std::string label = StripLabel(i);
      const double phi =
          weighting.WeightingPotential(uEvalCm, vEvalCm, wCm, label);
      out << "," << phi;
      sumPhi += phi;
      if (i == 0) phi0 = phi;
    }

    const double centralFraction = sumPhi > 0. ? phi0 / sumPhi : 0.;
    out << "," << sumPhi << "," << centralFraction << "\n";

    // Print a sparse subset so the terminal remains readable.
    if (iw == 0 || iw == scanSteps - 1 ||
        iw == scanSteps / 2 ||
        iw == scanSteps / 2 - scanSteps / 8 ||
        iw == scanSteps / 2 + scanSteps / 8) {
      std::cout << std::setw(7) << 10. * wCm << "    "
                << std::setw(11) << phi0 << "    "
                << std::setw(8) << sumPhi << "    "
                << std::setw(11) << centralFraction << "\n";
    }
  }

  // Also print the weighting-potential distribution over strips for an
  // avalanche centred at w = 0 on the wire plane.
  std::cout << "\nWeighting-potential sharing for a point at w = 0:\n"
            << "strip    center_w[mm]    phi_w    normalized\n";

  std::vector<double> phiAtZero;
  phiAtZero.reserve(labels.size());
  double totalAtZero = 0.;
  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double phi =
        weighting.WeightingPotential(uEvalCm, vEvalCm, 0., StripLabel(i));
    phiAtZero.push_back(phi);
    totalAtZero += phi;
  }

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double phi = phiAtZero[i + halfStrips];
    const double fraction = totalAtZero > 0. ? phi / totalAtZero : 0.;
    std::cout << std::setw(6) << i << "    "
              << std::setw(12) << 10. * i * stripPitchCm << "    "
              << std::setw(8) << phi << "    "
              << std::setw(10) << fraction << "\n";
  }

  std::cout << "\nWrote " << outName << "\n";
  return 0;
}

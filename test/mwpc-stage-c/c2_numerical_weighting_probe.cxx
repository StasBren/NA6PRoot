#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Garfield/ComponentElmer.hh"

namespace {

std::string ReadStringArg(const int argc, char** argv,
                          const std::string& key,
                          const std::string& defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return argv[i + 1];
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
  }
  return defaultValue;
}

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  return std::atof(
      ReadStringArg(argc, argv, key, std::to_string(defaultValue)).c_str());
}

bool Exists(const std::filesystem::path& p) {
  if (!std::filesystem::exists(p)) {
    std::cerr << "Missing required map file: " << p << "\n";
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const std::filesystem::path mapDir =
      ReadStringArg(argc, argv, "--map-dir", ".");
  const double gapMinusMm =
      ReadArg(argc, argv, "--gap-minus-mm", 2.5);
  const double gapPlusMm =
      ReadArg(argc, argv, "--gap-plus-mm", 2.5);
  const double wHalfMm =
      ReadArg(argc, argv, "--profile-half-w-mm", 4.0);
  const double wStepMm =
      ReadArg(argc, argv, "--profile-step-mm", 0.1);

  const auto meshDir = mapDir / "chamber";
  const auto header = meshDir / "mesh.header";
  const auto elements = meshDir / "mesh.elements";
  const auto nodes = meshDir / "mesh.nodes";
  const auto dielectrics = mapDir / "dielectrics.dat";
  const auto drift = meshDir / "drift.result";
  const auto wMinus = meshDir / "weight_minus_0.result";
  const auto wPlus = meshDir / "weight_plus_0.result";

  bool filesOk = true;
  for (const auto& p :
       {header, elements, nodes, dielectrics, drift, wMinus, wPlus}) {
    filesOk = Exists(p) && filesOk;
  }
  if (!filesOk) return 2;

  Garfield::ComponentElmer field(
      header.string(), elements.string(), nodes.string(),
      dielectrics.string(), drift.string(), "mm");

  if (!field.SetWeightingField(wMinus.string(), "minus_strip_0")) {
    std::cerr << "Failed to import minus-strip weighting map.\n";
    return 3;
  }
  if (!field.SetWeightingField(wPlus.string(), "plus_strip_0")) {
    std::cerr << "Failed to import plus-strip weighting map.\n";
    return 3;
  }

  // Garfield coordinates are cm even though the imported mesh is declared mm.
  const auto phi = [&](const double uMm, const double vMm,
                       const double wMm, const std::string& label) {
    return field.WeightingPotential(
        0.1 * uMm, 0.1 * vMm, 0.1 * wMm, label);
  };

  const double cathodeEpsMm = 0.02;
  const double wireProbeVMm = 0.05;

  const double phiMinusTarget =
      phi(0., -gapMinusMm + cathodeEpsMm, 0., "minus_strip_0");
  const double phiMinusOpposite =
      phi(0., +gapPlusMm - cathodeEpsMm, 0., "minus_strip_0");
  const double phiPlusTarget =
      phi(0., +gapPlusMm - cathodeEpsMm, 0., "plus_strip_0");
  const double phiPlusOpposite =
      phi(0., -gapMinusMm + cathodeEpsMm, 0., "plus_strip_0");
  const double phiMinusNearWire =
      phi(0., wireProbeVMm, 0., "minus_strip_0");
  const double phiPlusNearWire =
      phi(0., wireProbeVMm, 0., "plus_strip_0");

  std::cout << std::fixed << std::setprecision(8)
            << "\n=== STAGE C2a: NUMERICAL WEIGHTING-MAP PROBE ===\n"
            << "map directory                 : " << mapDir << "\n"
            << "minus target cathode phi      : " << phiMinusTarget << "\n"
            << "minus opposite cathode phi    : " << phiMinusOpposite << "\n"
            << "plus target cathode phi       : " << phiPlusTarget << "\n"
            << "plus opposite cathode phi     : " << phiPlusOpposite << "\n"
            << "minus phi near central wire   : " << phiMinusNearWire << "\n"
            << "plus phi near central wire    : " << phiPlusNearWire << "\n";

  const std::filesystem::path profileFile =
      mapDir / "central_strip_weighting_profile.csv";
  std::ofstream out(profileFile);
  if (!out) {
    std::cerr << "Could not create " << profileFile << "\n";
    return 4;
  }

  out << "w_mm,"
      << "phi_minus_vminus_quarter,phi_minus_v0,phi_minus_vplus_quarter,"
      << "phi_plus_vminus_quarter,phi_plus_v0,phi_plus_vplus_quarter\n";

  const double vMinusQuarter = -0.5 * gapMinusMm;
  const double vPlusQuarter = +0.5 * gapPlusMm;

  double minPhi = +1.e9;
  double maxPhi = -1.e9;
  double maxMirrorNearCathode = 0.;
  double maxMirrorFarCathode = 0.;
  bool allFinite = true;

  for (double w = -wHalfMm; w <= wHalfMm + 0.5 * wStepMm;
       w += wStepMm) {
    const double mmq =
        phi(0., vMinusQuarter, w, "minus_strip_0");
    const double mm0 =
        phi(0., 0.05, w, "minus_strip_0");
    const double mpq =
        phi(0., vPlusQuarter, w, "minus_strip_0");
    const double pmq =
        phi(0., vMinusQuarter, w, "plus_strip_0");
    const double pm0 =
        phi(0., 0.05, w, "plus_strip_0");
    const double ppq =
        phi(0., vPlusQuarter, w, "plus_strip_0");

    for (const double value : {mmq, mm0, mpq, pmq, pm0, ppq}) {
      allFinite = allFinite && std::isfinite(value);
      if (std::isfinite(value)) {
        minPhi = std::min(minPhi, value);
        maxPhi = std::max(maxPhi, value);
      }
    }

    // For the symmetric 2-sided reference geometry the reflected maps
    // should agree up to FEM discretisation error.
    maxMirrorNearCathode =
        std::max(maxMirrorNearCathode, std::abs(mmq - ppq));
    maxMirrorFarCathode =
        std::max(maxMirrorFarCathode, std::abs(mpq - pmq));

    out << w << ","
        << mmq << "," << mm0 << "," << mpq << ","
        << pmq << "," << pm0 << "," << ppq << "\n";
  }

  const bool boundarySanity =
      std::isfinite(phiMinusTarget) &&
      std::isfinite(phiPlusTarget) &&
      phiMinusTarget > 0.5 && phiMinusTarget < 1.05 &&
      phiPlusTarget > 0.5 && phiPlusTarget < 1.05 &&
      std::abs(phiMinusOpposite) < 0.1 &&
      std::abs(phiPlusOpposite) < 0.1 &&
      std::abs(phiMinusNearWire) < 0.2 &&
      std::abs(phiPlusNearWire) < 0.2;

  const bool profileSanity =
      allFinite && minPhi > -0.02 && maxPhi < 1.02 &&
      maxMirrorNearCathode < 0.02 &&
      maxMirrorFarCathode < 0.02;

  std::cout << "profile output                 : " << profileFile << "\n"
            << "profile phi range              : ["
            << minPhi << ", " << maxPhi << "]\n"
            << "max reflected diff (near)      : "
            << maxMirrorNearCathode << "\n"
            << "max reflected diff (far)       : "
            << maxMirrorFarCathode << "\n"
            << "boundary sanity                : "
            << (boundarySanity ? "PASS" : "CHECK") << "\n"
            << "profile sanity                 : "
            << (profileSanity ? "PASS" : "CHECK") << "\n"
            << "=========================================================\n";

  return (boundarySanity && profileSanity) ? 0 : 5;
}

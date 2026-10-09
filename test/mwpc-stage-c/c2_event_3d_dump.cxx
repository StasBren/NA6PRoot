#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

#include "../mwpc-avalanche/MWPCSignalUtils.h"

namespace {

constexpr double Pi = 3.14159265358979323846;

std::vector<na6p::mwpc::IonBirth> gIonBirths;

void RecordIonisation(const double x, const double y, const double z,
                      const double t, const int, const int,
                      Garfield::Medium*) {
  gIonBirths.push_back({x, y, z, t});
}

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return std::atof(argv[i + 1]);
    if (arg.rfind(prefix, 0) == 0) {
      return std::atof(arg.substr(prefix.size()).c_str());
    }
  }
  return defaultValue;
}

int ReadIntArg(const int argc, char** argv, const std::string& key,
               const int defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return std::atoi(argv[i + 1]);
    if (arg.rfind(prefix, 0) == 0) {
      return std::atoi(arg.substr(prefix.size()).c_str());
    }
  }
  return defaultValue;
}

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

double DegToRad(const double deg) {
  return deg * Pi / 180.;
}

bool FileExists(const std::filesystem::path& p) {
  if (!std::filesystem::exists(p)) {
    std::cerr << "Missing required file: " << p << "\n";
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.5);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 2.5);
  const double wirePitchCm =
      0.1 * ReadArg(argc, argv, "--wire-pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);
  const double bTesla = ReadArg(argc, argv, "--b-tesla", 0.0);
  const int halfWires =
      ReadIntArg(argc, argv, "--half-wires", 6);

  const std::string particle =
      ReadStringArg(argc, argv, "--particle", "mu-");
  const double momentumGeV =
      ReadArg(argc, argv, "--momentum-gev", 5.0);
  const double u0Cm =
      0.1 * ReadArg(argc, argv, "--u0-mm", 1.0);
  const double w0Cm =
      0.1 * ReadArg(argc, argv, "--w0-mm", 0.0);
  const double thetaUDeg =
      ReadArg(argc, argv, "--theta-u-deg", 10.0);
  const double thetaWDeg =
      ReadArg(argc, argv, "--theta-w-deg", 25.0);

  const int randomSeed =
      ReadIntArg(argc, argv, "--seed", 120001);
  const int detailedSeeds =
      ReadIntArg(argc, argv, "--detailed-seeds", 2);
  const int maxElectronPathsPerSeed =
      ReadIntArg(argc, argv, "--max-electron-paths-per-seed", 160);
  const int maxIonBirthPointsPerSeed =
      ReadIntArg(argc, argv, "--max-ion-births-per-seed", 5000);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::filesystem::path mapDir =
      ReadStringArg(argc, argv, "--map-dir", ".");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix", "c2_event3d");

  const double fieldUHalfMm =
      ReadArg(argc, argv, "--field-u-half-mm", 8.0);
  const int fieldNu =
      ReadIntArg(argc, argv, "--field-nu", 161);
  const int fieldNv =
      ReadIntArg(argc, argv, "--field-nv", 101);

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      halfWires < 2 || momentumGeV <= 0. ||
      detailedSeeds < 1 || maxElectronPathsPerSeed < 1 ||
      maxIonBirthPointsPerSeed < 1 || avalancheLimit < 1 ||
      fieldUHalfMm <= 0. || fieldNu < 3 || fieldNv < 3) {
    std::cerr << "Invalid C2 3D visualization parameters.\n";
    return 2;
  }

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  if (!gas.LoadIonMobility(mobilityFile)) {
    std::cerr << "Could not load ion mobility file: "
              << mobilityFile << "\n";
    return 3;
  }

  Garfield::ComponentAnalyticField physicalField;
  physicalField.SetMedium(&gas);
  for (int i = -halfWires; i <= halfWires; ++i) {
    physicalField.AddWire(
        i * wirePitchCm, 0., wireDiameterCm, hv, "anode");
  }
  physicalField.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  physicalField.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  physicalField.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&physicalField);
  const double uExtent =
      (halfWires + 0.5) * wirePitchCm;
  const double wExtent = 2.0;
  sensor.SetArea(
      -uExtent, -gapMinusCm, -wExtent,
      +uExtent, +gapPlusCm, +wExtent);

  // Numerical FEM map is sampled only for the visualization layers.
  const auto meshHeader = mapDir / "chamber" / "mesh.header";
  const auto meshElements = mapDir / "chamber" / "mesh.elements";
  const auto meshNodes = mapDir / "chamber" / "mesh.nodes";
  const auto dielectrics = mapDir / "dielectrics.dat";
  const auto driftResult = mapDir / "chamber" / "drift.result";
  const auto weightMinus = mapDir / "chamber" / "weight_minus_0.result";
  const auto weightPlus = mapDir / "chamber" / "weight_plus_0.result";

  bool filesOk = true;
  for (const auto& p :
       {meshHeader, meshElements, meshNodes, dielectrics,
        driftResult, weightMinus, weightPlus}) {
    filesOk = FileExists(p) && filesOk;
  }
  if (!filesOk) return 4;

  Garfield::ComponentElmer fem(
      meshHeader.string(), meshElements.string(), meshNodes.string(),
      dielectrics.string(), driftResult.string(), "mm");
  if (!fem.SetWeightingField(
          weightMinus.string(), "minus_strip_0") ||
      !fem.SetWeightingField(
          weightPlus.string(), "plus_strip_0")) {
    std::cerr << "Could not import weighting maps.\n";
    return 4;
  }

  const double tanU = std::tan(DegToRad(thetaUDeg));
  const double tanW = std::tan(DegToRad(thetaWDeg));
  const double marginCm =
      std::min(0.005, 0.02 * std::min(gapMinusCm, gapPlusCm));
  const double vStartCm = -gapMinusCm + marginCm;
  const double vExitCm = +gapPlusCm - marginCm;
  const double uStartCm = u0Cm + vStartCm * tanU;
  const double wStartCm = w0Cm + vStartCm * tanW;
  const double uExitCm = u0Cm + vExitCm * tanU;
  const double wExitCm = w0Cm + vExitCm * tanW;

  const double norm =
      std::sqrt(tanU * tanU + 1. + tanW * tanW);
  const double du = tanU / norm;
  const double dv = 1. / norm;
  const double dw = tanW / norm;

  std::ofstream muonOut(outputPrefix + "_muon_track.csv");
  std::ofstream clusterOut(outputPrefix + "_clusters.csv");
  std::ofstream primaryOut(outputPrefix + "_primary_electrons.csv");
  std::ofstream pathOut(outputPrefix + "_electron_paths.csv");
  std::ofstream ionBirthOut(outputPrefix + "_ion_births.csv");
  std::ofstream fieldOut(outputPrefix + "_field_slice_uv.csv");
  if (!muonOut || !clusterOut || !primaryOut ||
      !pathOut || !ionBirthOut || !fieldOut) {
    std::cerr << "Could not open C2 3D output files.\n";
    return 5;
  }

  muonOut
      << "point,u_mm,v_mm,w_mm\n";
  for (int i = 0; i <= 100; ++i) {
    const double f = i / 100.;
    muonOut
        << i << ","
        << 10. * (uStartCm + f * (uExitCm - uStartCm)) << ","
        << 10. * (vStartCm + f * (vExitCm - vStartCm)) << ","
        << 10. * (wStartCm + f * (wExitCm - wStartCm)) << "\n";
  }

  clusterOut
      << "cluster,u_mm,v_mm,w_mm,t_ns,n_primary,energy_eV\n";
  primaryOut
      << "cluster,seed,u_mm,v_mm,w_mm,t_ns,energy_eV\n";
  pathOut
      << "seed,path,point,u_mm,v_mm,w_mm,t_ns\n";
  ionBirthOut
      << "seed,birth,u_mm,v_mm,w_mm,t_ns\n";

  Garfield::RandomEngineRoot randomEngine(
      static_cast<unsigned int>(randomSeed));
  Garfield::Random::SetEngine(randomEngine);

  Garfield::TrackHeed heed;
  heed.SetSensor(&sensor);
  heed.SetParticle(particle);
  heed.SetMomentum(momentumGeV * 1.e9);
  heed.EnableDeltaElectronTransport();
  heed.NewTrack(
      uStartCm, vStartCm, wStartCm, 0.,
      du, dv, dw);

  int clusterIndex = 0;
  int seedIndex = 0;
  int detailedCount = 0;
  long long dumpedPaths = 0;
  long long dumpedPathPoints = 0;
  long long dumpedIonBirths = 0;

  double xc = 0., yc = 0., zc = 0., tc = 0.;
  double ec = 0., extra = 0.;
  int nc = 0;

  while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    clusterOut
        << clusterIndex << ","
        << 10. * xc << "," << 10. * yc << "," << 10. * zc << ","
        << tc << "," << nc << "," << ec << "\n";

    for (int ie = 0; ie < nc; ++ie) {
      double xe = 0., ye = 0., ze = 0., te = 0.;
      double ee = 0., dxe = 0., dye = 0., dze = 0.;
      if (!heed.GetElectron(
              ie, xe, ye, ze, te, ee, dxe, dye, dze)) {
        ++seedIndex;
        continue;
      }

      primaryOut
          << clusterIndex << "," << seedIndex << ","
          << 10. * xe << "," << 10. * ye << "," << 10. * ze << ","
          << te << "," << ee << "\n";

      if (detailedCount < detailedSeeds) {
        gIonBirths.clear();

        Garfield::AvalancheMicroscopic avalanche(&sensor);
        avalanche.EnableAvalancheSizeLimit(
            static_cast<unsigned int>(avalancheLimit));
        avalanche.EnableSignalCalculation(false);
        avalanche.EnableDriftLines();
        avalanche.SetUserHandleIonisation(RecordIonisation);

        avalanche.AvalancheElectron(
            xe, ye, ze, te, 0.1, 0., 0., 0.);

        const std::size_t nEndpoints =
            avalanche.GetNumberOfElectronEndpoints();
        const std::size_t pathStride =
            std::max<std::size_t>(
                1, (nEndpoints + maxElectronPathsPerSeed - 1) /
                       maxElectronPathsPerSeed);

        int localPath = 0;
        for (std::size_t iend = 0;
             iend < nEndpoints &&
             localPath < maxElectronPathsPerSeed;
             iend += pathStride, ++localPath) {
          const unsigned int np =
              avalanche.GetNumberOfElectronDriftLinePoints(iend);
          for (unsigned int ip = 0; ip < np; ++ip) {
            double x = 0., y = 0., z = 0., t = 0.;
            avalanche.GetElectronDriftLinePoint(
                x, y, z, t, ip, iend);
            pathOut
                << seedIndex << "," << localPath << "," << ip << ","
                << 10. * x << "," << 10. * y << "," << 10. * z << ","
                << t << "\n";
            ++dumpedPathPoints;
          }
          ++dumpedPaths;
        }

        const std::size_t birthStride =
            std::max<std::size_t>(
                1, (gIonBirths.size() + maxIonBirthPointsPerSeed - 1) /
                       maxIonBirthPointsPerSeed);
        int localBirth = 0;
        for (std::size_t ib = 0;
             ib < gIonBirths.size() &&
             localBirth < maxIonBirthPointsPerSeed;
             ib += birthStride, ++localBirth) {
          const auto& birth = gIonBirths[ib];
          ionBirthOut
              << seedIndex << "," << localBirth << ","
              << 10. * birth.u << ","
              << 10. * birth.v << ","
              << 10. * birth.w << ","
              << birth.t << "\n";
          ++dumpedIonBirths;
        }

        ++detailedCount;
      }

      ++seedIndex;
    }

    ++clusterIndex;
  }

  // Numerical FEM field slice at w=0.  Coordinates are exported in mm,
  // E components in V/cm, and potential in V.
  fieldOut
      << "iu,iv,u_mm,v_mm,w_mm,eu_Vcm,ev_Vcm,ew_Vcm,"
      << "E_mag_Vcm,potential_V,phi_minus_0,phi_plus_0,valid\n";

  const double wireRadiusMm = 5.e-3 * wireDiameterCm * 1.e4;
  for (int iu = 0; iu < fieldNu; ++iu) {
    const double uMm =
        -fieldUHalfMm +
        2. * fieldUHalfMm * iu / (fieldNu - 1.);
    for (int iv = 0; iv < fieldNv; ++iv) {
      const double vMm =
          -10. * gapMinusCm +
          10. * (gapMinusCm + gapPlusCm) * iv / (fieldNv - 1.);

      bool insideWire = false;
      const int nearest =
          static_cast<int>(
              std::lround(uMm / (10. * wirePitchCm)));
      if (std::abs(nearest) <= halfWires) {
        const double wireUMm =
            10. * nearest * wirePitchCm;
        insideWire =
            std::hypot(uMm - wireUMm, vMm) <= wireRadiusMm;
      }

      double eu = std::numeric_limits<double>::quiet_NaN();
      double ev = eu, ew = eu, emag = eu, potential = eu;
      double phiMinus = eu, phiPlus = eu;
      int valid = 0;

      if (!insideWire) {
        const double uCm = 0.1 * uMm;
        const double vCm = 0.1 * vMm;
        const auto e = fem.ElectricField(uCm, vCm, 0.);
        eu = e[0];
        ev = e[1];
        ew = e[2];
        potential = fem.ElectricPotential(uCm, vCm, 0.);
        phiMinus =
            fem.WeightingPotential(
                uCm, vCm, 0., "minus_strip_0");
        phiPlus =
            fem.WeightingPotential(
                uCm, vCm, 0., "plus_strip_0");
        emag = std::sqrt(eu * eu + ev * ev + ew * ew);
        valid =
            std::isfinite(emag) &&
            std::isfinite(potential) &&
            std::isfinite(phiMinus) &&
            std::isfinite(phiPlus);
      }

      fieldOut
          << iu << "," << iv << ","
          << uMm << "," << vMm << ",0,"
          << eu << "," << ev << "," << ew << ","
          << emag << "," << potential << ","
          << phiMinus << "," << phiPlus << ","
          << valid << "\n";
    }
  }

  std::cout << std::fixed << std::setprecision(4)
            << "\n=== STAGE C2 VISUALIZATION DUMP ===\n"
            << "event seed                    : " << randomSeed << "\n"
            << "muon p                         : " << momentumGeV << " GeV/c\n"
            << "theta_u / theta_w              : "
            << thetaUDeg << " / " << thetaWDeg << " deg\n"
            << "track start (u,v,w) mm         : "
            << 10. * uStartCm << ", "
            << 10. * vStartCm << ", "
            << 10. * wStartCm << "\n"
            << "track exit  (u,v,w) mm         : "
            << 10. * uExitCm << ", "
            << 10. * vExitCm << ", "
            << 10. * wExitCm << "\n"
            << "Heed clusters                  : " << clusterIndex << "\n"
            << "primary electrons              : " << seedIndex << "\n"
            << "detailed avalanche seeds       : " << detailedCount << "\n"
            << "electron drift paths exported  : " << dumpedPaths << "\n"
            << "electron path points exported  : " << dumpedPathPoints << "\n"
            << "ion-birth points exported      : " << dumpedIonBirths << "\n"
            << "field slice samples            : "
            << fieldNu * fieldNv << "\n"
            << "========================================\n";

  std::cout
      << "Wrote:\n"
      << "  " << outputPrefix << "_muon_track.csv\n"
      << "  " << outputPrefix << "_clusters.csv\n"
      << "  " << outputPrefix << "_primary_electrons.csv\n"
      << "  " << outputPrefix << "_electron_paths.csv\n"
      << "  " << outputPrefix << "_ion_births.csv\n"
      << "  " << outputPrefix << "_field_slice_uv.csv\n";

  return 0;
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

#include "../mwpc-avalanche/MWPCCoordinateFrame.h"

namespace {

constexpr double Pi = 3.14159265358979323846;

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

struct EventCounters {
  int clusters = 0;
  int primaryElectrons = 0;
  int seedsProcessed = 0;
  int zeroAvalanches = 0;
  long long avalancheElectrons = 0;
  long long avalancheIons = 0;
  long long collectedElectrons = 0;
  double energyLossEv = 0.;

  double sumGainWeightedSeedU = 0.;
  double sumGainWeightedSeedV = 0.;
  double sumGainWeightedSeedW = 0.;
  long long gainWeight = 0;

  std::set<int> activeWires;
  std::map<int, long long> wireCollectedElectrons;
};

}  // namespace

int main(int argc, char** argv) {
  // ------------------------------------------------------------------
  // Stage-C wire-only response.
  //
  // Physics chain retained:
  //   muon -> Heed clusters -> primary electrons -> microscopic avalanche
  //        -> electron endpoints collected on explicit analytic anode wires.
  //
  // Deliberately omitted:
  //   numerical cathode weighting, ion drift, strip charge sharing,
  //   electronics, thresholds and shaping.
  //
  // This isolates the coordinate measured by the anode-wire plane itself.
  // ------------------------------------------------------------------

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
  const int halfWires = ReadIntArg(argc, argv, "--half-wires", 6);
  const double wHalfCm =
      0.1 * ReadArg(argc, argv, "--w-half-mm", 20.0);

  const std::string particle =
      ReadStringArg(argc, argv, "--particle", "mu-");
  const double momentumGeV =
      ReadArg(argc, argv, "--momentum-gev", 5.0);

  // Truth coordinate is defined at the anode-wire plane v = 0.
  const double u0Cm =
      0.1 * ReadArg(argc, argv, "--u0-mm", 0.0);
  const double w0Cm =
      0.1 * ReadArg(argc, argv, "--w0-mm", 0.0);
  const double thetaUDeg =
      ReadArg(argc, argv, "--theta-u-deg", 0.0);
  const double thetaWDeg =
      ReadArg(argc, argv, "--theta-w-deg", 0.0);

  const int events =
      ReadIntArg(argc, argv, "--events", 1);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 160001);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);
  const int maxSeeds =
      ReadIntArg(argc, argv, "--max-seeds", 0);

  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "stageC_wire_response");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      halfWires < 2 || wHalfCm <= 0. ||
      momentumGeV <= 0. || events < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || maxSeeds < 0) {
    std::cerr << "Invalid Stage-C wire-response parameters.\n";
    return 2;
  }

  // Guard the coordinate convention used throughout the MWPC studies.
  const auto frame = na6p::mwpc::MakeNominalNA60Frame();
  const auto localFromGlobalX = frame.GlobalToLocalVector({1., 0., 0.});
  const auto localFromGlobalY = frame.GlobalToLocalVector({0., 1., 0.});
  const auto localFromGlobalZ = frame.GlobalToLocalVector({0., 0., 1.});
  constexpr double frameTol = 1.e-12;
  if (std::abs(localFromGlobalX.w - 1.) > frameTol ||
      std::abs(localFromGlobalY.u - 1.) > frameTol ||
      std::abs(localFromGlobalZ.v - 1.) > frameTol) {
    std::cerr << "Nominal ChamberFrame no longer matches Stage-C convention.\n";
    return 3;
  }

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);
  for (int i = -halfWires; i <= halfWires; ++i) {
    field.AddWire(i * wirePitchCm, 0., wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");

  // Nominal MNP33 mapping: global B_Y -> local B_u -> Garfield x.
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);

  const double uExtent = (halfWires + 0.5) * wirePitchCm;
  sensor.SetArea(
      -uExtent, -gapMinusCm, -wHalfCm,
      +uExtent, +gapPlusCm, +wHalfCm);

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

  if (std::abs(uStartCm) >= uExtent ||
      std::abs(uExitCm) >= uExtent ||
      std::abs(wStartCm) >= wHalfCm ||
      std::abs(wExitCm) >= wHalfCm) {
    std::cerr
        << "Requested track does not fit inside the wire-response sandbox.\n";
    return 4;
  }

  const double wireRadiusCm = 0.5 * wireDiameterCm;

  const std::string eventFile = outputPrefix + "_event_summary.csv";
  const std::string wireFile = outputPrefix + "_wire_summary.csv";

  std::ofstream eventOut(eventFile);
  std::ofstream wireOut(wireFile);
  if (!eventOut || !wireOut) {
    std::cerr << "Could not open Stage-C wire-response output files.\n";
    return 5;
  }

  eventOut
      << "event,random_seed,particle,momentum_GeV,"
      << "u0_mm,w0_mm,theta_u_deg,theta_w_deg,"
      << "u_start_mm,v_start_mm,w_start_mm,"
      << "u_exit_mm,v_exit_mm,w_exit_mm,"
      << "clusters,primary_electrons,seeds_processed,seed_truncated,"
      << "zero_avalanches,energy_loss_eV,"
      << "avalanche_electrons,avalanche_ions,"
      << "collected_electrons,active_wires,"
      << "gain_weighted_seed_u_mm,gain_weighted_seed_v_mm,"
      << "gain_weighted_seed_w_mm,"
      << "dominant_wire_index,dominant_wire_u_mm,dominant_wire_fraction,"
      << "wire_cog_u_mm,"
      << "binary_residual_mm,wire_cog_residual_mm,"
      << "dominant_tie\n";

  wireOut
      << "event,random_seed,u0_mm,wire_index,wire_u_mm,"
      << "collected_electrons,collected_fraction\n";

  std::cout << std::fixed << std::setprecision(4)
            << "\n=== STAGE C: WIRE-ONLY MUON RESPONSE ===\n"
            << "particle / momentum        : "
            << particle << " / " << momentumGeV << " GeV/c\n"
            << "truth at v=0 (u,w)         : "
            << 10. * u0Cm << ", " << 10. * w0Cm << " mm\n"
            << "angles (theta_u,theta_w)   : "
            << thetaUDeg << ", " << thetaWDeg << " deg\n"
            << "gap (-,+)                  : "
            << 10. * gapMinusCm << ", "
            << 10. * gapPlusCm << " mm\n"
            << "wire pitch / diameter      : "
            << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "B(local u)                 : "
            << bTesla << " T\n"
            << "events                     : "
            << events << "\n\n";

  for (int iev = 0; iev < events; ++iev) {
    const int randomSeed = baseSeed + iev;
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

    EventCounters c;
    int seedIndex = 0;
    bool truncated = false;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    double ec = 0., extra = 0.;
    int nc = 0;

    while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      ++c.clusters;
      c.primaryElectrons += nc;
      c.energyLossEv += ec;

      for (int ie = 0; ie < nc; ++ie) {
        if (maxSeeds > 0 && c.seedsProcessed >= maxSeeds) {
          truncated = true;
          break;
        }

        double xe = 0., ye = 0., ze = 0., te = 0.;
        double ee = 0., dxe = 0., dye = 0., dze = 0.;
        if (!heed.GetElectron(
                ie, xe, ye, ze, te, ee, dxe, dye, dze)) {
          ++seedIndex;
          continue;
        }

        ++c.seedsProcessed;

        Garfield::AvalancheMicroscopic avalanche(&sensor);
        avalanche.EnableAvalancheSizeLimit(
            static_cast<unsigned int>(avalancheLimit));
        avalanche.EnableSignalCalculation(false);

        const bool avalancheOk =
            avalanche.AvalancheElectron(
                xe, ye, ze, te, 0.1, 0., 0., 0.);

        int nAvalancheElectrons = 0;
        int nAvalancheIons = 0;
        avalanche.GetAvalancheSize(
            nAvalancheElectrons, nAvalancheIons);

        c.avalancheElectrons += nAvalancheElectrons;
        c.avalancheIons += nAvalancheIons;

        const long long seedGain =
            std::max(0, nAvalancheIons);
        c.sumGainWeightedSeedU += seedGain * xe;
        c.sumGainWeightedSeedV += seedGain * ye;
        c.sumGainWeightedSeedW += seedGain * ze;
        c.gainWeight += seedGain;

        if (!avalancheOk || nAvalancheIons <= 0) {
          ++c.zeroAvalanches;
        }

        const std::size_t nEndpoints =
            avalanche.GetNumberOfElectronEndpoints();

        for (std::size_t iend = 0; iend < nEndpoints; ++iend) {
          double uE0 = 0., vE0 = 0., wE0 = 0., tE0 = 0., eE0 = 0.;
          double uE1 = 0., vE1 = 0., wE1 = 0., tE1 = 0., eE1 = 0.;
          int status = 0;

          avalanche.GetElectronEndpoint(
              iend,
              uE0, vE0, wE0, tE0, eE0,
              uE1, vE1, wE1, tE1, eE1,
              status);

          const int nearest =
              static_cast<int>(std::lround(uE1 / wirePitchCm));
          if (std::abs(nearest) > halfWires) continue;

          const double wireU = nearest * wirePitchCm;
          const double rho = std::hypot(uE1 - wireU, vE1);
          if (rho < 1.5 * wireRadiusCm) {
            ++c.collectedElectrons;
            c.activeWires.insert(nearest);
            ++c.wireCollectedElectrons[nearest];
          }
        }

        ++seedIndex;
      }

      if (truncated) break;
    }

    const double gainWeightedSeedUMm =
        c.gainWeight > 0
            ? 10. * c.sumGainWeightedSeedU /
                  static_cast<double>(c.gainWeight)
            : std::numeric_limits<double>::quiet_NaN();
    const double gainWeightedSeedVMm =
        c.gainWeight > 0
            ? 10. * c.sumGainWeightedSeedV /
                  static_cast<double>(c.gainWeight)
            : std::numeric_limits<double>::quiet_NaN();
    const double gainWeightedSeedWMm =
        c.gainWeight > 0
            ? 10. * c.sumGainWeightedSeedW /
                  static_cast<double>(c.gainWeight)
            : std::numeric_limits<double>::quiet_NaN();

    int dominantWire = 0;
    long long dominantCount = -1;
    int dominantMultiplicity = 0;
    double wireCogUMm =
        std::numeric_limits<double>::quiet_NaN();

    if (c.collectedElectrons > 0) {
      double numerMm = 0.;
      for (const auto& [wireIndex, nCollected] :
           c.wireCollectedElectrons) {
        const double wireUMm = 10. * wireIndex * wirePitchCm;
        numerMm += wireUMm * static_cast<double>(nCollected);

        if (nCollected > dominantCount) {
          dominantCount = nCollected;
          dominantWire = wireIndex;
          dominantMultiplicity = 1;
        } else if (nCollected == dominantCount) {
          ++dominantMultiplicity;
        }
      }
      wireCogUMm =
          numerMm / static_cast<double>(c.collectedElectrons);
    }

    const double dominantWireUMm =
        10. * dominantWire * wirePitchCm;
    const double dominantFraction =
        c.collectedElectrons > 0 && dominantCount >= 0
            ? static_cast<double>(dominantCount) /
                  static_cast<double>(c.collectedElectrons)
            : std::numeric_limits<double>::quiet_NaN();

    const double truthUMm = 10. * u0Cm;
    const double binaryResidualMm =
        c.collectedElectrons > 0
            ? dominantWireUMm - truthUMm
            : std::numeric_limits<double>::quiet_NaN();
    const double wireCogResidualMm =
        std::isfinite(wireCogUMm)
            ? wireCogUMm - truthUMm
            : std::numeric_limits<double>::quiet_NaN();

    eventOut
        << iev << "," << randomSeed << ","
        << particle << "," << momentumGeV << ","
        << truthUMm << "," << 10. * w0Cm << ","
        << thetaUDeg << "," << thetaWDeg << ","
        << 10. * uStartCm << "," << 10. * vStartCm << ","
        << 10. * wStartCm << ","
        << 10. * uExitCm << "," << 10. * vExitCm << ","
        << 10. * wExitCm << ","
        << c.clusters << ","
        << c.primaryElectrons << ","
        << c.seedsProcessed << ","
        << (truncated ? 1 : 0) << ","
        << c.zeroAvalanches << ","
        << c.energyLossEv << ","
        << c.avalancheElectrons << ","
        << c.avalancheIons << ","
        << c.collectedElectrons << ","
        << c.activeWires.size() << ","
        << gainWeightedSeedUMm << ","
        << gainWeightedSeedVMm << ","
        << gainWeightedSeedWMm << ","
        << dominantWire << ","
        << dominantWireUMm << ","
        << dominantFraction << ","
        << wireCogUMm << ","
        << binaryResidualMm << ","
        << wireCogResidualMm << ","
        << (dominantMultiplicity > 1 ? 1 : 0)
        << "\n";

    for (const auto& [wireIndex, nCollected] :
         c.wireCollectedElectrons) {
      wireOut
          << iev << "," << randomSeed << ","
          << truthUMm << ","
          << wireIndex << ","
          << 10. * wireIndex * wirePitchCm << ","
          << nCollected << ","
          << (c.collectedElectrons > 0
                  ? static_cast<double>(nCollected) /
                        static_cast<double>(c.collectedElectrons)
                  : 0.)
          << "\n";
    }

    std::cout
        << "event " << iev
        << " seed=" << randomSeed
        << " primary_e=" << c.primaryElectrons
        << " avalanche_e=" << c.avalancheElectrons
        << " active_wires=" << c.activeWires.size()
        << " dominant=" << dominantWire
        << " f_dom=" << dominantFraction
        << " wireCoG=" << wireCogUMm << " mm"
        << "\n";
  }

  std::cout
      << "\nWrote:\n"
      << "  " << eventFile << "\n"
      << "  " << wireFile << "\n"
      << "\nMODEL BOUNDARY:\n"
      << "  - Heed muon ionisation: yes\n"
      << "  - microscopic electron avalanche: yes\n"
      << "  - anode-wire collection coordinate: yes\n"
      << "  - cathode strips / Elmer weighting: intentionally omitted\n"
      << "  - ion signal: intentionally omitted\n"
      << "  - electronics/noise/threshold/shaping: no\n";

  return 0;
}

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

#include "../mwpc-avalanche/MWPCCoordinateFrame.h"
#include "../mwpc-avalanche/MWPCSignalUtils.h"
#include "../mwpc-avalanche/MWPCStereoReadout.h"

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

struct EventCounters {
  int clusters = 0;
  int primaryElectrons = 0;
  int seedsProcessed = 0;
  int zeroAvalanches = 0;
  int lateElectronSegments = 0;
  long long avalancheElectrons = 0;
  long long avalancheIons = 0;
  long long recordedIonBirths = 0;
  long long ionTransportFailures = 0;
  long long collectedElectrons = 0;
  double energyLossEv = 0.;
  std::set<int> activeWires;
};

double DegToRad(const double deg) {
  return deg * Pi / 180.;
}

}  // namespace

int main(int argc, char** argv) {
  // ------------------------------------------------------------------
  // Stage C1a
  //
  // One complete local muon crossing:
  //   Heed clusters -> all conduction electrons -> microscopic avalanches
  //   -> electron + finite-time ion Shockley-Ramo segments
  //   -> Stage-B3c.2 hybrid stereo strip response -> event-level CoG.
  //
  // This is deliberately built from the coordinate and signal conventions
  // already established in Stages A/B.
  // ------------------------------------------------------------------

  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double wirePitchCm =
      0.1 * ReadArg(argc, argv, "--wire-pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);
  const double bTesla = ReadArg(argc, argv, "--b-tesla", 0.0);

  const double stripPitchCm =
      0.1 * ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthCm =
      0.1 * ReadArg(argc, argv, "--strip-width-mm", 1.7);
  const double tanAlpha =
      ReadArg(argc, argv, "--tan-alpha", 0.10);

  const int halfWires =
      ReadIntArg(argc, argv, "--half-wires", 6);
  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 12);

  const std::string particle =
      ReadStringArg(argc, argv, "--particle", "mu-");
  const double momentumGeV =
      ReadArg(argc, argv, "--momentum-gev", 5.0);

  // Truth coordinate is defined at the wire plane v=0.
  const double u0Cm =
      0.1 * ReadArg(argc, argv, "--u0-mm", 1.0);
  const double w0Cm =
      0.1 * ReadArg(argc, argv, "--w0-mm", 0.0);

  // Projected local angles measured by C0:
  // theta_u = atan2(p_u,p_v), theta_w = atan2(p_w,p_v).
  const double thetaUDeg =
      ReadArg(argc, argv, "--theta-u-deg", 0.0);
  const double thetaWDeg =
      ReadArg(argc, argv, "--theta-w-deg", 0.0);

  const int events =
      ReadIntArg(argc, argv, "--events", 1);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 120001);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);
  const int maxSeeds =
      ReadIntArg(argc, argv, "--max-seeds", 0);

  const double observationNs =
      ReadArg(argc, argv, "--observation-ns", 100.0);
  const double ionRkDtNs =
      ReadArg(argc, argv, "--ion-rk-dt-ns", 5.0);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "stageC1_full_muon");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm + 1.e-12 ||
      tanAlpha <= 0. || halfWires < 2 || halfStrips < 2 ||
      momentumGeV <= 0. || events < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || maxSeeds < 0 ||
      observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-C1 parameters.\n";
    return 2;
  }

  // ------------------------------------------------------------------
  // Existing Stage-A coordinate convention.
  // ------------------------------------------------------------------
  const auto nominalFrame = na6p::mwpc::MakeNominalNA60Frame();
  const auto localFromGlobalX =
      nominalFrame.GlobalToLocalVector({1., 0., 0.});
  const auto localFromGlobalY =
      nominalFrame.GlobalToLocalVector({0., 1., 0.});
  const auto localFromGlobalZ =
      nominalFrame.GlobalToLocalVector({0., 0., 1.});

  // Fail loudly if somebody changes the coordinate utility underneath C1.
  const double frameTol = 1.e-12;
  if (std::abs(localFromGlobalX.w - 1.) > frameTol ||
      std::abs(localFromGlobalY.u - 1.) > frameTol ||
      std::abs(localFromGlobalZ.v - 1.) > frameTol) {
    std::cerr << "Nominal ChamberFrame no longer matches Stage-C convention.\n";
    return 3;
  }

  // ------------------------------------------------------------------
  // Existing Stage-A physical cell.
  // ------------------------------------------------------------------
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  if (!gas.LoadIonMobility(mobilityFile)) {
    std::cerr << "Could not load ion mobility file: "
              << mobilityFile << "\n";
    return 4;
  }

  Garfield::ComponentAnalyticField physicalField;
  physicalField.SetMedium(&gas);

  for (int i = -halfWires; i <= halfWires; ++i) {
    physicalField.AddWire(
        i * wirePitchCm, 0., wireDiameterCm, hv, "anode");
  }
  physicalField.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  physicalField.AddPlaneY(+gapPlusCm, 0., "cathode_plus");

  // Nominal MNP33 mapping: global B_Y -> local B_u -> Garfield x.
  physicalField.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor physicalSensor;
  physicalSensor.AddComponent(&physicalField);

  const double uExtent =
      (halfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 4.0) * stripPitchCm;
  physicalSensor.SetArea(
      -uExtent, -gapMinusCm, -wExtent,
      +uExtent, +gapPlusCm, +wExtent);

  // ------------------------------------------------------------------
  // Reusable Stage-B3c.2 hybrid stereo kernel.
  // ------------------------------------------------------------------
  na6p::mwpc::StereoReadout readout(
      gapMinusCm, gapPlusCm,
      stripPitchCm, stripWidthCm,
      tanAlpha, halfStrips);

  const double alpha = readout.AlphaRad();
  const double cosA = std::cos(alpha);
  const double sinA = std::sin(alpha);

  // ------------------------------------------------------------------
  // Local track line.
  //
  // Truth (u0,w0) is specified at v=0.  Forward NA60+/DiCE transport is +v.
  // du/dv = tan(theta_u), dw/dv = tan(theta_w).
  // ------------------------------------------------------------------
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
      std::abs(wStartCm) >= wExtent ||
      std::abs(wExitCm) >= wExtent) {
    std::cerr
        << "Requested track does not fit inside the local sandbox.\n"
        << "Increase --half-wires/--half-strips or move (u0,w0).\n";
    return 5;
  }

  const double trueXPlusMm =
      10. * (w0Cm * cosA - u0Cm * sinA);
  const double trueXMinusMm =
      10. * (w0Cm * cosA + u0Cm * sinA);

  const double wireRadiusCm = 0.5 * wireDiameterCm;

  const std::string eventFile =
      outputPrefix + "_event_summary.csv";
  const std::string stripFile =
      outputPrefix + "_strip_summary.csv";
  const std::string seedFile =
      outputPrefix + "_seed_summary.csv";

  std::ofstream eventOut(eventFile);
  std::ofstream stripOut(stripFile);
  std::ofstream seedOut(seedFile);
  if (!eventOut || !stripOut || !seedOut) {
    std::cerr << "Could not open Stage-C1 output files.\n";
    return 6;
  }

  eventOut
      << "event,random_seed,particle,momentum_GeV,"
      << "u0_mm,w0_mm,theta_u_deg,theta_w_deg,"
      << "u_start_mm,v_start_mm,w_start_mm,"
      << "u_exit_mm,v_exit_mm,w_exit_mm,"
      << "clusters,primary_electrons,seeds_processed,seed_truncated,"
      << "zero_avalanches,energy_loss_eV,"
      << "avalanche_electrons,avalanche_ions,recorded_ion_births,"
      << "ion_transport_failures,late_electron_segments,"
      << "collected_electrons,active_wires,"
      << "q_cathode_minus_fC,q_cathode_plus_fC,"
      << "x_plus_true_mm,x_minus_true_mm,"
      << "x_plus_cog_mm,x_minus_cog_mm,"
      << "u_cog_mm,w_cog_mm,u_residual_mm,w_residual_mm,"
      << "readout_valid\n";

  stripOut
      << "event,side,strip,center_mm,Q_fC,A_fC\n";

  seedOut
      << "event,cluster,seed,seed_u_mm,seed_v_mm,seed_w_mm,seed_t_ns,"
      << "avalanche_electrons,avalanche_ions,ion_births,"
      << "electron_endpoints,collected_electrons\n";

  std::cout << std::fixed << std::setprecision(4)
            << "\n=== STAGE C1a: FULL MUON -> EVENT-LEVEL STRIP RESPONSE ===\n"
            << "particle / momentum        : "
            << particle << " / " << momentumGeV << " GeV/c\n"
            << "truth at v=0 (u,w)         : "
            << 10. * u0Cm << ", " << 10. * w0Cm << " mm\n"
            << "angles (theta_u,theta_w)   : "
            << thetaUDeg << ", " << thetaWDeg << " deg\n"
            << "track start (u,v,w)        : "
            << 10. * uStartCm << ", "
            << 10. * vStartCm << ", "
            << 10. * wStartCm << " mm\n"
            << "track exit  (u,v,w)        : "
            << 10. * uExitCm << ", "
            << 10. * vExitCm << ", "
            << 10. * wExitCm << " mm\n"
            << "B(local u)                 : "
            << bTesla << " T\n"
            << "observation time           : "
            << observationNs << " ns\n"
            << "stereo alpha               : "
            << readout.AlphaDeg() << " deg\n"
            << "events                      : "
            << events << "\n"
            << "max seeds (0=all)           : "
            << maxSeeds << "\n\n";

  for (int iev = 0; iev < events; ++iev) {
    const int randomSeed = baseSeed + iev;
    Garfield::RandomEngineRoot randomEngine(
        static_cast<unsigned int>(randomSeed));
    Garfield::Random::SetEngine(randomEngine);

    Garfield::TrackHeed heed;
    heed.SetSensor(&physicalSensor);
    heed.SetParticle(particle);
    heed.SetMomentum(momentumGeV * 1.e9);
    heed.EnableDeltaElectronTransport();

    heed.NewTrack(
        uStartCm, vStartCm, wStartCm, 0.,
        du, dv, dw);

    EventCounters counters;
    auto response = readout.EmptyResult();

    int clusterIndex = 0;
    int seedIndex = 0;
    bool truncated = false;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    double ec = 0., extra = 0.;
    int nc = 0;

    while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      ++counters.clusters;
      counters.primaryElectrons += nc;
      counters.energyLossEv += ec;

      for (int ie = 0; ie < nc; ++ie) {
        if (maxSeeds > 0 && counters.seedsProcessed >= maxSeeds) {
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

        ++counters.seedsProcessed;
        gIonBirths.clear();

        Garfield::AvalancheMicroscopic avalanche(&physicalSensor);
        avalanche.EnableAvalancheSizeLimit(
            static_cast<unsigned int>(avalancheLimit));
        avalanche.EnableSignalCalculation(false);
        avalanche.SetUserHandleIonisation(RecordIonisation);

        const bool avalancheOk =
            avalanche.AvalancheElectron(
                xe, ye, ze, te, 0.1, 0., 0., 0.);

        int nAvalancheElectrons = 0;
        int nAvalancheIons = 0;
        avalanche.GetAvalancheSize(
            nAvalancheElectrons, nAvalancheIons);

        counters.avalancheElectrons += nAvalancheElectrons;
        counters.avalancheIons += nAvalancheIons;
        counters.recordedIonBirths += gIonBirths.size();

        if (!avalancheOk || nAvalancheIons <= 0) {
          ++counters.zeroAvalanches;
        }

        int seedCollected = 0;
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

          if (tE0 >= observationNs) continue;
          if (tE1 > observationNs) {
            ++counters.lateElectronSegments;
            continue;
          }

          readout.AddSegment(
              response, physicalField,
              -na6p::mwpc::ElementaryChargeFc,
              uE0, vE0, wE0,
              uE1, vE1, wE1);

          const int nearest =
              static_cast<int>(std::lround(uE1 / wirePitchCm));
          const double wireU = nearest * wirePitchCm;
          const double rho = std::hypot(uE1 - wireU, vE1);
          if (rho < 1.5 * wireRadiusCm) {
            ++seedCollected;
            ++counters.collectedElectrons;
            counters.activeWires.insert(nearest);
          }
        }

        for (const auto& birth : gIonBirths) {
          if (birth.t >= observationNs) continue;

          std::array<double, 3> x =
              {birth.u, birth.v, birth.w};
          double t = birth.t;

          if (!na6p::mwpc::AdvanceIonRK4(
                  physicalSensor, x, t,
                  observationNs, ionRkDtNs)) {
            ++counters.ionTransportFailures;
            continue;
          }

          readout.AddSegment(
              response, physicalField,
              +na6p::mwpc::ElementaryChargeFc,
              birth.u, birth.v, birth.w,
              x[0], x[1], x[2]);
        }

        seedOut
            << iev << "," << clusterIndex << "," << seedIndex << ","
            << 10. * xe << "," << 10. * ye << "," << 10. * ze << ","
            << te << ","
            << nAvalancheElectrons << "," << nAvalancheIons << ","
            << gIonBirths.size() << "," << nEndpoints << ","
            << seedCollected << "\n";

        ++seedIndex;
      }

      ++clusterIndex;
      if (truncated) break;
    }

    const bool readoutValid = readout.Finalize(response);

    const double uResidualMm =
        readoutValid ? response.uCogMm - 10. * u0Cm
                     : std::numeric_limits<double>::quiet_NaN();
    const double wResidualMm =
        readoutValid ? response.wCogMm - 10. * w0Cm
                     : std::numeric_limits<double>::quiet_NaN();

    eventOut
        << iev << "," << randomSeed << ","
        << particle << "," << momentumGeV << ","
        << 10. * u0Cm << "," << 10. * w0Cm << ","
        << thetaUDeg << "," << thetaWDeg << ","
        << 10. * uStartCm << "," << 10. * vStartCm << ","
        << 10. * wStartCm << ","
        << 10. * uExitCm << "," << 10. * vExitCm << ","
        << 10. * wExitCm << ","
        << counters.clusters << ","
        << counters.primaryElectrons << ","
        << counters.seedsProcessed << ","
        << (truncated ? 1 : 0) << ","
        << counters.zeroAvalanches << ","
        << counters.energyLossEv << ","
        << counters.avalancheElectrons << ","
        << counters.avalancheIons << ","
        << counters.recordedIonBirths << ","
        << counters.ionTransportFailures << ","
        << counters.lateElectronSegments << ","
        << counters.collectedElectrons << ","
        << counters.activeWires.size() << ","
        << response.qCathodeMinusFc << ","
        << response.qCathodePlusFc << ","
        << trueXPlusMm << "," << trueXMinusMm << ","
        << response.xPlusAlphaCogMm << ","
        << response.xMinusAlphaCogMm << ","
        << response.uCogMm << "," << response.wCogMm << ","
        << uResidualMm << "," << wResidualMm << ","
        << (readoutValid ? 1 : 0) << "\n";

    const auto& stripIds = readout.StripIds();
    for (std::size_t j = 0; j < stripIds.size(); ++j) {
      const double centerMm =
          10. * stripIds[j] * readout.StripPitchCm();

      stripOut
          << iev << ",minus," << stripIds[j] << ","
          << centerMm << ","
          << response.qMinusFc[j] << ","
          << std::abs(response.qMinusFc[j]) << "\n";

      stripOut
          << iev << ",plus," << stripIds[j] << ","
          << centerMm << ","
          << response.qPlusFc[j] << ","
          << std::abs(response.qPlusFc[j]) << "\n";
    }

    std::cout
        << "event " << iev
        << " : clusters=" << counters.clusters
        << " primary_e=" << counters.primaryElectrons
        << " seeds=" << counters.seedsProcessed
        << (truncated ? " [TRUNCATED]" : "")
        << " avalanche_ions=" << counters.avalancheIons
        << " active_wires=" << counters.activeWires.size()
        << " ion_fail=" << counters.ionTransportFailures
        << " readout=" << (readoutValid ? "valid" : "INVALID");

    if (readoutValid) {
      std::cout
          << "  CoG(u,w)=("
          << response.uCogMm << ", "
          << response.wCogMm << ") mm"
          << " residual=("
          << uResidualMm << ", "
          << wResidualMm << ") mm";
    }
    std::cout << "\n";
  }

  std::cout
      << "\nWrote:\n"
      << "  " << eventFile << "\n"
      << "  " << stripFile << "\n"
      << "  " << seedFile << "\n"
      << "\nMODEL BOUNDARY:\n"
      << "  - full Heed muon + many microscopic avalanches: yes\n"
      << "  - finite-time electron/ion Shockley-Ramo sum: yes\n"
      << "  - Stage-B hybrid opposite-cathode stereo sharing: yes\n"
      << "  - electronics/noise/threshold/shaping: no\n"
      << "  - calibrated absolute detector gain/resolution: no\n";

  return 0;
}

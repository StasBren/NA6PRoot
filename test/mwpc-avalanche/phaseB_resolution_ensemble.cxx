#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"

namespace {

constexpr double ElementaryChargeFc = 1.602176634e-4;

struct IonBirth {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
};

std::vector<IonBirth> gIonBirths;

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
    if (arg.rfind(prefix, 0) == 0) {
      return arg.substr(prefix.size());
    }
  }
  return defaultValue;
}

std::vector<double> ParseDoubles(const std::string& text) {
  std::vector<double> values;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item.empty()) continue;
    values.push_back(std::stod(item));
  }
  return values;
}

std::vector<int> ParseInts(const std::string& text) {
  std::vector<int> values;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (item.empty()) continue;
    values.push_back(std::stoi(item));
  }
  return values;
}

std::string StripLabel(const int i) {
  if (i < 0) return "strip_m" + std::to_string(-i);
  if (i > 0) return "strip_p" + std::to_string(i);
  return "strip_0";
}

bool IonVelocity(Garfield::Sensor& sensor,
                 const std::array<double, 3>& x,
                 std::array<double, 3>& velocity) {
  double ex = 0., ey = 0., ez = 0., potential = 0.;
  Garfield::Medium* medium = nullptr;
  int status = 0;
  sensor.ElectricField(x[0], x[1], x[2],
                       ex, ey, ez, potential, medium, status);
  if (status != 0 || !medium) return false;

  double bx = 0., by = 0., bz = 0.;
  sensor.MagneticField(x[0], x[1], x[2], bx, by, bz, status);
  if (status != 0) {
    bx = by = bz = 0.;
  }

  return medium->IonVelocity(ex, ey, ez, bx, by, bz,
                             velocity[0], velocity[1], velocity[2]);
}

bool AdvanceIonRK4(Garfield::Sensor& sensor,
                   std::array<double, 3>& x,
                   double& t,
                   const double targetTime,
                   const double maxDtNs) {
  if (targetTime <= t) return true;

  while (t < targetTime) {
    const double h = std::min(maxDtNs, targetTime - t);

    std::array<double, 3> k1{}, k2{}, k3{}, k4{};
    if (!IonVelocity(sensor, x, k1)) return false;

    std::array<double, 3> x2 = x;
    for (int j = 0; j < 3; ++j) x2[j] += 0.5 * h * k1[j];
    if (!IonVelocity(sensor, x2, k2)) return false;

    std::array<double, 3> x3 = x;
    for (int j = 0; j < 3; ++j) x3[j] += 0.5 * h * k2[j];
    if (!IonVelocity(sensor, x3, k3)) return false;

    std::array<double, 3> x4 = x;
    for (int j = 0; j < 3; ++j) x4[j] += h * k3[j];
    if (!IonVelocity(sensor, x4, k4)) return false;

    for (int j = 0; j < 3; ++j) {
      x[j] += h * (k1[j] + 2. * k2[j] + 2. * k3[j] + k4[j]) / 6.;
    }
    t += h;

    if (!sensor.IsInside(x[0], x[1], x[2])) return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  // Stage-B ensemble study:
  //   true position w0
  //      -> independent real microscopic avalanches
  //      -> early induced charges on strips
  //      -> three-strip CoG reconstruction
  //      -> event-by-event residual distribution.
  //
  // To make O(100) events x O(10) positions practical, positive ions are
  // transported only up to the largest requested early observation window.
  // Their finite-time induced charge is evaluated from Shockley-Ramo endpoint
  // differences.  A fixed-step RK4 integrates the same Garfield ion-velocity
  // model used by DriftLineRKF.  This truncated transport must be validated
  // against the full single-event DriftLineRKF signal before production use.

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
  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 5);

  const std::string side =
      ReadStringArg(argc, argv, "--side", "minus");
  const double seedUCm =
      0.1 * ReadArg(argc, argv, "--seed-u-mm", 0.0);
  const double seedDistanceMm =
      ReadArg(argc, argv, "--seed-distance-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  auto positionsMm = ParseDoubles(
      ReadStringArg(argc, argv, "--positions-mm",
                    "-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8"));
  auto windowsNs = ParseDoubles(
      ReadStringArg(argc, argv, "--windows-ns",
                    "25,50,100,200"));
  auto recoStrips = ParseInts(
      ReadStringArg(argc, argv, "--reco-strips", "-1,0,1"));

  const int eventsPerPosition =
      ReadIntArg(argc, argv, "--events-per-position", 50);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 30000);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);
  const int minIonsForReco =
      ReadIntArg(argc, argv, "--min-ions-for-reco", 1);
  const double ionRkDtNs =
      ReadArg(argc, argv, "--ion-rk-dt-ns", 1.0);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputFile =
      ReadStringArg(argc, argv, "--output",
                    "stageB_resolution_ensemble.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      positionsMm.empty() || windowsNs.empty() || recoStrips.empty() ||
      eventsPerPosition < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || minIonsForReco < 0 ||
      ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-B resolution-ensemble parameters.\n";
    return 2;
  }

  if (side != "minus" && side != "plus") {
    std::cerr << "--side must be 'minus' or 'plus'.\n";
    return 2;
  }

  std::sort(windowsNs.begin(), windowsNs.end());
  if (windowsNs.front() <= 0.) {
    std::cerr << "Observation windows must be positive.\n";
    return 2;
  }

  const double sign = side == "plus" ? +1. : -1.;
  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const double selectedGapCm = std::abs(planeY);
  const double seedVCm = sign * 0.1 * seedDistanceMm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;

  if (std::abs(seedVCm) <= wireRadiusCm ||
      std::abs(seedVCm) >= selectedGapCm) {
    std::cerr << "Seed must be outside the wire and inside the gas gap.\n";
    return 2;
  }

  // ------------------------------------------------------------------
  // Physical gas / chamber.
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
    return 3;
  }

  Garfield::ComponentAnalyticField driftField;
  driftField.SetMedium(&gas);

  constexpr int halfNumberOfWires = 4;
  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    driftField.AddWire(i * wirePitchCm, 0.,
                       wireDiameterCm, hv, "anode");
  }
  driftField.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  driftField.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  driftField.SetMagneticField(bTesla, 0., 0.);

  // ------------------------------------------------------------------
  // Ideal strip weighting geometry, kept identical to previous Stage-B tests.
  // ------------------------------------------------------------------
  Garfield::ComponentAnalyticField weighting;
  weighting.AddPlaneY(0., 1., "weighting_back");
  weighting.AddPlaneY(planeY, 0., "weighting_front");

  std::vector<int> stripIds;
  std::vector<std::string> labels;
  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double centerW = i * stripPitchCm;
    const double wMin = centerW - 0.5 * stripWidthCm;
    const double wMax = centerW + 0.5 * stripWidthCm;
    const std::string label = StripLabel(i);

    weighting.AddStripOnPlaneY(
        'x', planeY, wMin, wMax, label, selectedGapCm);
    stripIds.push_back(i);
    labels.push_back(label);
  }

  Garfield::Sensor sensor;
  sensor.AddComponent(&driftField);
  for (const auto& label : labels) {
    sensor.AddElectrode(&weighting, label);
  }

  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 1.5) * stripPitchCm;
  sensor.SetArea(-uExtent, -gapMinusCm, -wExtent,
                 +uExtent, +gapPlusCm, +wExtent);

  // Map reconstruction-strip ids to entries in the full strip vector.
  std::vector<std::size_t> recoIndices;
  for (const int k : recoStrips) {
    auto it = std::find(stripIds.begin(), stripIds.end(), k);
    if (it == stripIds.end()) {
      std::cerr << "Reconstruction strip " << k
                << " is outside the configured strip range.\n";
      return 2;
    }
    recoIndices.push_back(
        static_cast<std::size_t>(std::distance(stripIds.begin(), it)));
  }

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out << "true_w_mm,event_index,random_seed,window_ns,"
      << "avalanche_ok,garfield_electrons,garfield_ions,"
      << "recorded_ion_births,electron_end_max_ns,"
      << "ion_birth_mean_w_mm,ion_birth_sigma_w_mm,"
      << "ion_transport_failures,reconstructable,"
      << "w_hat_mm,residual_mm,total_cluster_amplitude_fC";
  for (const int k : recoStrips) {
    out << ",Q_strip_" << k << "_fC"
        << ",A_strip_" << k << "_fC"
        << ",fraction_strip_" << k;
  }
  out << "\n";

  std::cout << std::fixed << std::setprecision(3);
  std::cout
      << "\n=== STAGE B: ENSEMBLE POSITION / RESOLUTION STUDY ===\n"
      << "positions [mm]       : ";
  for (const auto w : positionsMm) std::cout << w << " ";
  std::cout
      << "\nevents per position : " << eventsPerPosition
      << "\nwindows [ns]         : ";
  for (const auto t : windowsNs) std::cout << t << " ";
  std::cout
      << "\nreco strips          : ";
  for (const auto k : recoStrips) std::cout << k << " ";
  std::cout
      << "\nion RK4 dt           : " << ionRkDtNs << " ns"
      << "\noutput               : " << outputFile << "\n\n";

  // ------------------------------------------------------------------
  // Ensemble loop.
  //
  // Use the same EVENT seed set at every true position. This is a common-
  // random-numbers design: event j at all w0 values shares the same microscopic
  // history translated relative to the strips. Across j, events are independent.
  // ------------------------------------------------------------------
  for (const double trueWMm : positionsMm) {
    int nReco = 0;
    double sumResidual = 0.;
    double sumResidual2 = 0.;

    for (int iev = 0; iev < eventsPerPosition; ++iev) {
      const int randomSeed = baseSeed + iev;

      Garfield::RandomEngineRoot engine(
          static_cast<unsigned int>(randomSeed));
      Garfield::Random::SetEngine(engine);

      gIonBirths.clear();

      Garfield::AvalancheMicroscopic avalanche(&sensor);
      avalanche.EnableAvalancheSizeLimit(
          static_cast<unsigned int>(avalancheLimit));
      avalanche.EnableSignalCalculation(false);
      avalanche.SetUserHandleIonisation(RecordIonisation);

      const double seedWCm = 0.1 * trueWMm;
      const bool avalancheOk =
          avalanche.AvalancheElectron(
              seedUCm, seedVCm, seedWCm, 0.,
              seedEnergyEv, 0., 0., 0.);

      int nElectrons = 0;
      int nIons = 0;
      avalanche.GetAvalancheSize(nElectrons, nIons);

      // Electron induced charge from all microscopic trajectory endpoints.
      std::vector<double> qElectron(labels.size(), 0.);
      double electronEndMaxNs = 0.;

      const std::size_t nEndpoints =
          avalanche.GetNumberOfElectronEndpoints();
      for (std::size_t ie = 0; ie < nEndpoints; ++ie) {
        double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
        double u1 = 0., v1 = 0., w1 = 0., t1 = 0., e1 = 0.;
        int status = 0;
        avalanche.GetElectronEndpoint(
            ie, u0, v0, w0, t0, e0,
            u1, v1, w1, t1, e1, status);

        electronEndMaxNs = std::max(electronEndMaxNs, t1);

        for (std::size_t j = 0; j < labels.size(); ++j) {
          const double phi0 =
              weighting.WeightingPotential(u0, v0, w0, labels[j]);
          const double phi1 =
              weighting.WeightingPotential(u1, v1, w1, labels[j]);
          qElectron[j] +=
              (-ElementaryChargeFc) * (phi1 - phi0);
        }
      }

      double meanBirthW = 0.;
      double meanBirthW2 = 0.;
      for (const auto& b : gIonBirths) {
        const double wmm = 10. * b.w;
        meanBirthW += wmm;
        meanBirthW2 += wmm * wmm;
      }
      if (!gIonBirths.empty()) {
        meanBirthW /= static_cast<double>(gIonBirths.size());
        meanBirthW2 /= static_cast<double>(gIonBirths.size());
      }
      const double sigmaBirthW =
          std::sqrt(std::max(
              0., meanBirthW2 - meanBirthW * meanBirthW));

      // qIonByWindow[iWindow][iStrip]
      std::vector<std::vector<double>> qIonByWindow(
          windowsNs.size(),
          std::vector<double>(labels.size(), 0.));

      int ionTransportFailures = 0;

      for (const auto& birth : gIonBirths) {
        std::array<double, 3> x = {birth.u, birth.v, birth.w};
        double t = birth.t;

        std::vector<double> phiBirth(labels.size(), 0.);
        for (std::size_t j = 0; j < labels.size(); ++j) {
          phiBirth[j] =
              weighting.WeightingPotential(
                  birth.u, birth.v, birth.w, labels[j]);
        }

        bool transportOk = true;
        for (std::size_t iw = 0; iw < windowsNs.size(); ++iw) {
          const double T = windowsNs[iw];

          if (T > t && transportOk) {
            transportOk =
                AdvanceIonRK4(sensor, x, t, T, ionRkDtNs);
            if (!transportOk) ++ionTransportFailures;
          }

          if (!transportOk || T <= birth.t) continue;

          for (std::size_t j = 0; j < labels.size(); ++j) {
            const double phi =
                weighting.WeightingPotential(
                    x[0], x[1], x[2], labels[j]);
            qIonByWindow[iw][j] +=
                ElementaryChargeFc * (phi - phiBirth[j]);
          }
        }
      }

      for (std::size_t iw = 0; iw < windowsNs.size(); ++iw) {
        const double T = windowsNs[iw];

        const bool electronComplete =
            electronEndMaxNs <= T + 1.e-9;
        const bool gainEnough =
            nIons >= minIonsForReco;
        const bool transportGood =
            ionTransportFailures == 0;

        std::vector<double> qTotal(labels.size(), 0.);
        for (std::size_t j = 0; j < labels.size(); ++j) {
          qTotal[j] = qElectron[j] + qIonByWindow[iw][j];
        }

        double amplitudeSum = 0.;
        double weightedSum = 0.;
        for (std::size_t ir = 0; ir < recoIndices.size(); ++ir) {
          const std::size_t j = recoIndices[ir];
          const double amp = std::abs(qTotal[j]);
          const double centerMm =
              10. * stripIds[j] * stripPitchCm;
          amplitudeSum += amp;
          weightedSum += centerMm * amp;
        }

        const bool reconstructable =
            avalancheOk && gainEnough && transportGood &&
            electronComplete && amplitudeSum > 0.;

        const double wHat =
            reconstructable ? weightedSum / amplitudeSum : 0.;
        const double residual =
            reconstructable ? wHat - trueWMm : 0.;

        if (reconstructable && std::abs(T - 100.) < 1.e-9) {
          ++nReco;
          sumResidual += residual;
          sumResidual2 += residual * residual;
        }

        out
            << trueWMm << ","
            << iev << ","
            << randomSeed << ","
            << T << ","
            << (avalancheOk ? 1 : 0) << ","
            << nElectrons << ","
            << nIons << ","
            << gIonBirths.size() << ","
            << electronEndMaxNs << ","
            << meanBirthW << ","
            << sigmaBirthW << ","
            << ionTransportFailures << ","
            << (reconstructable ? 1 : 0) << ","
            << (reconstructable ? wHat : std::numeric_limits<double>::quiet_NaN()) << ","
            << (reconstructable ? residual : std::numeric_limits<double>::quiet_NaN()) << ","
            << amplitudeSum;

        for (std::size_t ir = 0; ir < recoIndices.size(); ++ir) {
          const std::size_t j = recoIndices[ir];
          const double q = qTotal[j];
          const double amp = std::abs(q);
          const double frac =
              amplitudeSum > 0. ? amp / amplitudeSum : 0.;
          out << "," << q << "," << amp << "," << frac;
        }
        out << "\n";
      }
    }

    const double mean =
        nReco > 0 ? sumResidual / nReco : 0.;
    const double variance =
        nReco > 1
            ? std::max(0., sumResidual2 / nReco - mean * mean)
            : 0.;
    const double sigma = std::sqrt(variance);

    std::cout
        << "w0=" << std::setw(6) << trueWMm << " mm"
        << "  reco@100ns=" << std::setw(3) << nReco
        << "/" << eventsPerPosition
        << "  bias=" << std::setw(8) << 1000. * mean << " um"
        << "  sigma=" << std::setw(8) << 1000. * sigma << " um\n";
  }

  std::cout
      << "\nNOTE: sigma above is the one-seed, ideal-strip, no-threshold "
         "intrinsic spread at 100 ns.\n"
      << "NOTE: zero/small-gain events are counted as unreconstructed, not "
         "resampled.\n"
      << "NOTE: this is not yet the full muon-chamber resolution; Heed "
         "primary-ionisation statistics and electronics are still absent.\n"
      << "=============================================================\n";

  return 0;
}

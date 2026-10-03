#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
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

struct ElectronEndpoint {
  double u0 = 0.;
  double v0 = 0.;
  double w0 = 0.;
  double u1 = 0.;
  double v1 = 0.;
  double w1 = 0.;
};

struct IonSegment {
  double u0 = 0.;
  double v0 = 0.;
  double w0 = 0.;
  double u1 = 0.;
  double v1 = 0.;
  double w1 = 0.;
};

struct AvalancheEvent {
  int eventIndex = 0;
  int randomSeed = 0;
  int garfieldElectrons = 0;
  int garfieldIons = 0;
  double electronEndMaxNs = 0.;
  std::vector<ElectronEndpoint> electrons;
  std::vector<IonSegment> ions;
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
    if (!item.empty()) values.push_back(std::stod(item));
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
  if (status != 0) bx = by = bz = 0.;

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
  // Stage B2: readout-geometry response at alpha = 0.
  //
  // We isolate strip pitch/width effects by generating an ensemble of real
  // microscopic avalanches ONCE in the physical MWPC field, transporting their
  // positive ions only to one early observation time, and then re-evaluating
  // Shockley-Ramo induced charge for many ideal strip geometries.
  //
  // Because the chamber field is translationally invariant along w, the same
  // stored avalanche can be shifted by w0 relative to the strips. This removes
  // avalanche-to-avalanche differences when comparing geometries.

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

  const std::string side =
      ReadStringArg(argc, argv, "--side", "minus");
  const double seedUCm =
      0.1 * ReadArg(argc, argv, "--seed-u-mm", 0.0);
  const double seedDistanceMm =
      ReadArg(argc, argv, "--seed-distance-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  auto pitchesMm = ParseDoubles(
      ReadStringArg(argc, argv, "--strip-pitches-mm",
                    "1.5,1.7,1.9,2.1,2.3"));
  auto widthsMm = ParseDoubles(
      ReadStringArg(argc, argv, "--strip-widths-mm",
                    "1.0,1.3,1.5,1.7,2.0"));
  auto xiValues = ParseDoubles(
      ReadStringArg(argc, argv, "--xi-values",
                    "-0.45,-0.30,-0.15,0,0.15,0.30,0.45"));

  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 6);
  const int events =
      ReadIntArg(argc, argv, "--events", 30);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 41000);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);
  const double observationNs =
      ReadArg(argc, argv, "--observation-ns", 100.0);
  const double ionRkDtNs =
      ReadArg(argc, argv, "--ion-rk-dt-ns", 1.0);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputFile =
      ReadStringArg(argc, argv, "--output",
                    "stageB_strip_geometry_scan.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      pitchesMm.empty() || widthsMm.empty() || xiValues.empty() ||
      halfStrips < 2 || events < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-B strip-geometry scan parameters.\n";
    return 2;
  }

  if (side != "minus" && side != "plus") {
    std::cerr << "--side must be 'minus' or 'plus'.\n";
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

  for (const auto p : pitchesMm) {
    if (p <= 0.) {
      std::cerr << "All strip pitches must be positive.\n";
      return 2;
    }
  }
  for (const auto w : widthsMm) {
    if (w <= 0.) {
      std::cerr << "All strip widths must be positive.\n";
      return 2;
    }
  }

  // ------------------------------------------------------------------
  // Physical chamber and gas.
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

  Garfield::Sensor physicalSensor;
  physicalSensor.AddComponent(&driftField);

  const double maxPitchCm =
      0.1 * *std::max_element(pitchesMm.begin(), pitchesMm.end());
  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 2.5) * maxPitchCm;

  physicalSensor.SetArea(-uExtent, -gapMinusCm, -wExtent,
                         +uExtent, +gapPlusCm, +wExtent);

  // ------------------------------------------------------------------
  // Generate one reusable ensemble of physical avalanche histories.
  // ------------------------------------------------------------------
  std::vector<AvalancheEvent> eventBank;
  eventBank.reserve(events);

  int zeroOrFailed = 0;
  int ionTransportFailedEvents = 0;

  for (int iev = 0; iev < events; ++iev) {
    const int randomSeed = baseSeed + iev;

    Garfield::RandomEngineRoot engine(
        static_cast<unsigned int>(randomSeed));
    Garfield::Random::SetEngine(engine);

    gIonBirths.clear();

    Garfield::AvalancheMicroscopic avalanche(&physicalSensor);
    avalanche.EnableAvalancheSizeLimit(
        static_cast<unsigned int>(avalancheLimit));
    avalanche.EnableSignalCalculation(false);
    avalanche.SetUserHandleIonisation(RecordIonisation);

    const bool avalancheOk =
        avalanche.AvalancheElectron(
            seedUCm, seedVCm, 0., 0.,
            seedEnergyEv, 0., 0., 0.);

    int nElectrons = 0;
    int nIons = 0;
    avalanche.GetAvalancheSize(nElectrons, nIons);

    if (!avalancheOk || nIons <= 0) {
      ++zeroOrFailed;
      continue;
    }

    AvalancheEvent evt;
    evt.eventIndex = iev;
    evt.randomSeed = randomSeed;
    evt.garfieldElectrons = nElectrons;
    evt.garfieldIons = nIons;

    const std::size_t nEndpoints =
        avalanche.GetNumberOfElectronEndpoints();
    evt.electrons.reserve(nEndpoints);

    for (std::size_t ie = 0; ie < nEndpoints; ++ie) {
      double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
      double u1 = 0., v1 = 0., w1 = 0., t1 = 0., e1 = 0.;
      int status = 0;

      avalanche.GetElectronEndpoint(
          ie, u0, v0, w0, t0, e0,
          u1, v1, w1, t1, e1, status);

      evt.electronEndMaxNs =
          std::max(evt.electronEndMaxNs, t1);
      evt.electrons.push_back(
          {u0, v0, w0, u1, v1, w1});
    }

    // At 100 ns the electron avalanche is already complete in our baseline.
    // If not, skip the event rather than approximating a partial electron path.
    if (evt.electronEndMaxNs > observationNs) {
      ++zeroOrFailed;
      continue;
    }

    bool transportGood = true;
    evt.ions.reserve(gIonBirths.size());

    for (const auto& birth : gIonBirths) {
      if (birth.t >= observationNs) continue;

      std::array<double, 3> x =
          {birth.u, birth.v, birth.w};
      double t = birth.t;

      if (!AdvanceIonRK4(
              physicalSensor, x, t, observationNs, ionRkDtNs)) {
        transportGood = false;
        break;
      }

      evt.ions.push_back(
          {birth.u, birth.v, birth.w,
           x[0], x[1], x[2]});
    }

    if (!transportGood) {
      ++ionTransportFailedEvents;
      continue;
    }

    eventBank.push_back(std::move(evt));
  }

  if (eventBank.empty()) {
    std::cerr << "No usable positive-gain avalanches were generated.\n";
    return 4;
  }

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out
      << "strip_pitch_mm,strip_width_mm,fill_factor,"
      << "xi,true_w_mm,event_index,random_seed,garfield_ions,"
      << "observation_ns,total_abs_signal_all_fC,"
      << "three_strip_abs_signal_fC,three_strip_capture_fraction,"
      << "neighbor_fraction_three_strip,central_fraction_three_strip,"
      << "left_right_asymmetry";

  for (int k = -2; k <= 2; ++k) {
    out << ",Q_strip_" << k << "_fC"
        << ",A_strip_" << k << "_fC";
  }
  out << "\n";

  std::cout << std::fixed << std::setprecision(3);
  std::cout
      << "\n=== STAGE B2: STRIP PITCH x WIDTH SCAN (alpha = 0) ===\n"
      << "side                    : " << side << "\n"
      << "observation time        : " << observationNs << " ns\n"
      << "requested events        : " << events << "\n"
      << "usable positive events  : " << eventBank.size() << "\n"
      << "zero/failed events      : " << zeroOrFailed << "\n"
      << "ion-transport failures  : " << ionTransportFailedEvents << "\n"
      << "pitch grid [mm]         : ";
  for (const auto p : pitchesMm) std::cout << p << " ";
  std::cout << "\nwidth grid [mm]         : ";
  for (const auto w : widthsMm) std::cout << w << " ";
  std::cout << "\nxi = w0/p               : ";
  for (const auto xi : xiValues) std::cout << xi << " ";
  std::cout << "\noutput                  : " << outputFile << "\n\n";

  // ------------------------------------------------------------------
  // Geometry loop. Invalid width > pitch combinations are skipped.
  // ------------------------------------------------------------------
  int nGeometries = 0;

  for (const double pitchMm : pitchesMm) {
    const double pitchCm = 0.1 * pitchMm;

    for (const double widthMm : widthsMm) {
      if (widthMm > pitchMm + 1.e-12) continue;

      const double widthCm = 0.1 * widthMm;
      const double fillFactor = widthMm / pitchMm;

      Garfield::ComponentAnalyticField weighting;
      weighting.AddPlaneY(0., 1., "weighting_back");
      weighting.AddPlaneY(planeY, 0., "weighting_front");

      std::vector<int> stripIds;
      std::vector<std::string> labels;

      for (int k = -halfStrips; k <= halfStrips; ++k) {
        const double centerW = k * pitchCm;
        const double wMin = centerW - 0.5 * widthCm;
        const double wMax = centerW + 0.5 * widthCm;
        const std::string label = StripLabel(k);

        weighting.AddStripOnPlaneY(
            'x', planeY, wMin, wMax, label, selectedGapCm);

        stripIds.push_back(k);
        labels.push_back(label);
      }

      auto indexOf = [&](const int k) -> std::size_t {
        auto it = std::find(stripIds.begin(), stripIds.end(), k);
        return static_cast<std::size_t>(
            std::distance(stripIds.begin(), it));
      };

      const std::size_t im1 = indexOf(-1);
      const std::size_t i0 = indexOf(0);
      const std::size_t ip1 = indexOf(+1);

      ++nGeometries;

      for (const double xi : xiValues) {
        const double shiftCm = xi * pitchCm;
        const double trueWMm = xi * pitchMm;

        for (const auto& evt : eventBank) {
          std::vector<double> q(labels.size(), 0.);

          // Electron contribution, complete by the fixed observation time.
          for (const auto& e : evt.electrons) {
            for (std::size_t j = 0; j < labels.size(); ++j) {
              const double phi0 =
                  weighting.WeightingPotential(
                      e.u0, e.v0, e.w0 + shiftCm, labels[j]);
              const double phi1 =
                  weighting.WeightingPotential(
                      e.u1, e.v1, e.w1 + shiftCm, labels[j]);

              q[j] +=
                  (-ElementaryChargeFc) * (phi1 - phi0);
            }
          }

          // Positive-ion contribution from birth to the fixed early time.
          for (const auto& ion : evt.ions) {
            for (std::size_t j = 0; j < labels.size(); ++j) {
              const double phi0 =
                  weighting.WeightingPotential(
                      ion.u0, ion.v0, ion.w0 + shiftCm, labels[j]);
              const double phi1 =
                  weighting.WeightingPotential(
                      ion.u1, ion.v1, ion.w1 + shiftCm, labels[j]);

              q[j] +=
                  ElementaryChargeFc * (phi1 - phi0);
            }
          }

          std::vector<double> a(labels.size(), 0.);
          double allAbs = 0.;
          for (std::size_t j = 0; j < labels.size(); ++j) {
            a[j] = std::abs(q[j]);
            allAbs += a[j];
          }

          const double a3 = a[im1] + a[i0] + a[ip1];
          const double neighborFraction =
              a3 > 0. ? (a[im1] + a[ip1]) / a3 : 0.;
          const double centralFraction =
              a3 > 0. ? a[i0] / a3 : 0.;
          const double asymmetry =
              a3 > 0. ? (a[ip1] - a[im1]) / a3 : 0.;
          const double capture =
              allAbs > 0. ? a3 / allAbs : 0.;

          out
              << pitchMm << ","
              << widthMm << ","
              << fillFactor << ","
              << xi << ","
              << trueWMm << ","
              << evt.eventIndex << ","
              << evt.randomSeed << ","
              << evt.garfieldIons << ","
              << observationNs << ","
              << allAbs << ","
              << a3 << ","
              << capture << ","
              << neighborFraction << ","
              << centralFraction << ","
              << asymmetry;

          for (int k = -2; k <= 2; ++k) {
            const std::size_t j = indexOf(k);
            out << "," << q[j] << "," << a[j];
          }
          out << "\n";
        }
      }

      std::cout
          << "done: pitch=" << std::setw(4) << pitchMm
          << " mm  width=" << std::setw(4) << widthMm
          << " mm  fill=" << std::setw(5) << fillFactor << "\n";
    }
  }

  std::cout
      << "\nvalid geometries scanned : " << nGeometries << "\n"
      << "NOTE: this is an ideal alpha=0 strip plane on the "
      << (side == "minus" ? "2-mm" : "4-mm")
      << " cathode side.\n"
      << "NOTE: only positive-gain avalanches enter normalized sharing "
         "observables; zero-gain probability is a separate issue.\n"
      << "NOTE: no threshold, electronics shaping, tilt, second strip family, "
         "or real finite cathode geometry is included yet.\n"
      << "=============================================================\n";

  return 0;
}

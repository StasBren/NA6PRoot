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
constexpr double Pi = 3.14159265358979323846;

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

std::vector<int> ParseInts(const std::string& text) {
  std::vector<int> values;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) values.push_back(std::stoi(item));
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
  // Stage B3a: one ideal tilted strip family.
  //
  // Coordinates:
  //   u = across anode wires
  //   v = chamber normal
  //   w = along anode wires
  //
  // A +alpha strip family runs at angle alpha relative to +u. Its coordinate
  // perpendicular to the strips is
  //
  //   x_alpha = w cos(alpha) - u sin(alpha).
  //
  // We generate a physical microscopic avalanche once around the central wire
  // (u_wire = 0), transport its positive ions to a fixed early observation
  // time, then reuse the exact same avalanche on neighboring wire positions by
  // translating it by integer multiples of the 4-mm wire pitch. This is the
  // cleanest test of the ideal periodic geometry: any shift in the strip
  // response should follow x_alpha.
  //
  // The readout weighting problem remains a 2-D infinite-strip problem. Tilt
  // is implemented exactly for this ideal geometry by rotating each avalanche
  // point into x_alpha before evaluating the strip weighting potential.

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
  const double seedDistanceMm =
      ReadArg(argc, argv, "--seed-distance-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  const double stripPitchMm =
      ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthMm =
      ReadArg(argc, argv, "--strip-width-mm", 1.7);

  auto tanAlphas = ParseDoubles(
      ReadStringArg(argc, argv, "--tan-alphas",
                    "0,0.05,0.10"));
  auto wireIndices = ParseInts(
      ReadStringArg(argc, argv, "--wire-indices",
                    "-2,-1,0,1,2"));
  auto xiWValues = ParseDoubles(
      ReadStringArg(argc, argv, "--xi-w-values",
                    "-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5"));

  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 7);
  const int events =
      ReadIntArg(argc, argv, "--events", 30);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 52000);
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
                    "stageB_tilt_scan.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      stripPitchMm <= 0. || stripWidthMm <= 0. ||
      stripWidthMm > stripPitchMm + 1.e-12 ||
      tanAlphas.empty() || wireIndices.empty() || xiWValues.empty() ||
      halfStrips < 2 || events < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-B tilt-scan parameters.\n";
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

  // ------------------------------------------------------------------
  // Physical chamber and gas: unchanged from Stage B2.
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

  const double stripPitchCm = 0.1 * stripPitchMm;
  const double stripWidthCm = 0.1 * stripWidthMm;

  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 3.0) * stripPitchCm;

  physicalSensor.SetArea(-uExtent, -gapMinusCm, -wExtent,
                         +uExtent, +gapPlusCm, +wExtent);

  // ------------------------------------------------------------------
  // Generate one reusable central-wire microscopic avalanche ensemble.
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
            0., seedVCm, 0., 0.,
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

    // At T=100 ns the electron avalanche is complete in our current setup.
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

  // ------------------------------------------------------------------
  // Ideal strip weighting cell in its own local frame:
  // x' = coordinate along the strip (irrelevant),
  // z' = coordinate perpendicular to the strip.
  // ------------------------------------------------------------------
  Garfield::ComponentAnalyticField weighting;
  weighting.AddPlaneY(0., 1., "weighting_back");
  weighting.AddPlaneY(planeY, 0., "weighting_front");

  std::vector<int> stripIds;
  std::vector<std::string> labels;

  for (int k = -halfStrips; k <= halfStrips; ++k) {
    const double center = k * stripPitchCm;
    const double zMin = center - 0.5 * stripWidthCm;
    const double zMax = center + 0.5 * stripWidthCm;
    const std::string label = StripLabel(k);

    weighting.AddStripOnPlaneY(
        'x', planeY, zMin, zMax, label, selectedGapCm);

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

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out
      << "tan_alpha,alpha_deg,wire_index,wire_u_mm,"
      << "xi_w,w_shift_mm,projected_x_mm,projected_eta,"
      << "event_index,random_seed,garfield_ions,observation_ns,"
      << "total_abs_signal_all_fC,three_strip_abs_signal_fC,"
      << "three_strip_capture_fraction,neighbor_fraction_three_strip,"
      << "central_fraction_three_strip,left_right_asymmetry";

  for (int k = -2; k <= 2; ++k) {
    out << ",Q_strip_" << k << "_fC"
        << ",A_strip_" << k << "_fC";
  }
  out << "\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout
      << "\n=== STAGE B3a: ONE TILTED STRIP FAMILY ===\n"
      << "side                     : " << side << "\n"
      << "strip pitch / width      : "
      << stripPitchMm << " / " << stripWidthMm << " mm\n"
      << "wire pitch               : "
      << 10. * wirePitchCm << " mm\n"
      << "observation time         : "
      << observationNs << " ns\n"
      << "requested avalanches     : " << events << "\n"
      << "usable positive events   : " << eventBank.size() << "\n"
      << "zero/failed events       : " << zeroOrFailed << "\n"
      << "ion-transport failures   : "
      << ionTransportFailedEvents << "\n"
      << "tan(alpha) values        : ";
  for (const auto x : tanAlphas) std::cout << x << " ";
  std::cout << "\nwire indices             : ";
  for (const auto x : wireIndices) std::cout << x << " ";
  std::cout << "\nxi_w = w0/p values       : ";
  for (const auto x : xiWValues) std::cout << x << " ";
  std::cout << "\noutput                   : "
            << outputFile << "\n\n";

  // ------------------------------------------------------------------
  // Re-evaluate the SAME microscopic avalanche ensemble on translated
  // wire positions and for several strip tilt angles.
  // ------------------------------------------------------------------
  for (const double tanAlpha : tanAlphas) {
    const double alpha = std::atan(tanAlpha);
    const double cosA = std::cos(alpha);
    const double sinA = std::sin(alpha);
    const double alphaDeg = alpha * 180. / Pi;

    for (const int wireIndex : wireIndices) {
      const double uShiftCm = wireIndex * wirePitchCm;
      const double wireUMm = 10. * uShiftCm;

      for (const double xiW : xiWValues) {
        const double wShiftCm = xiW * stripPitchCm;
        const double wShiftMm = 10. * wShiftCm;

        // Nominal position of the translated avalanche in the coordinate
        // perpendicular to the tilted strips.
        const double projectedXCm =
            wShiftCm * cosA - uShiftCm * sinA;
        const double projectedXMm = 10. * projectedXCm;
        const double projectedEta =
            projectedXCm / stripPitchCm;

        for (const auto& evt : eventBank) {
          std::vector<double> q(labels.size(), 0.);

          for (const auto& e : evt.electrons) {
            const double u0 = e.u0 + uShiftCm;
            const double u1 = e.u1 + uShiftCm;
            const double w0 = e.w0 + wShiftCm;
            const double w1 = e.w1 + wShiftCm;

            const double z0 = w0 * cosA - u0 * sinA;
            const double z1 = w1 * cosA - u1 * sinA;

            for (std::size_t j = 0; j < labels.size(); ++j) {
              const double phi0 =
                  weighting.WeightingPotential(
                      0., e.v0, z0, labels[j]);
              const double phi1 =
                  weighting.WeightingPotential(
                      0., e.v1, z1, labels[j]);

              q[j] +=
                  (-ElementaryChargeFc) * (phi1 - phi0);
            }
          }

          for (const auto& ion : evt.ions) {
            const double u0 = ion.u0 + uShiftCm;
            const double u1 = ion.u1 + uShiftCm;
            const double w0 = ion.w0 + wShiftCm;
            const double w1 = ion.w1 + wShiftCm;

            const double z0 = w0 * cosA - u0 * sinA;
            const double z1 = w1 * cosA - u1 * sinA;

            for (std::size_t j = 0; j < labels.size(); ++j) {
              const double phi0 =
                  weighting.WeightingPotential(
                      0., ion.v0, z0, labels[j]);
              const double phi1 =
                  weighting.WeightingPotential(
                      0., ion.v1, z1, labels[j]);

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
              << tanAlpha << ","
              << alphaDeg << ","
              << wireIndex << ","
              << wireUMm << ","
              << xiW << ","
              << wShiftMm << ","
              << projectedXMm << ","
              << projectedEta << ","
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
    }

    std::cout
        << "done tan(alpha)=" << tanAlpha
        << "  alpha=" << alphaDeg << " deg\n";
  }

  std::cout
      << "\nNOTE: wire-index translation reuses the same central-wire "
         "microscopic avalanche.\n"
      << "      This intentionally isolates the ideal periodic readout "
         "geometry from avalanche statistics.\n"
      << "NOTE: x_alpha = w cos(alpha) - u sin(alpha).\n"
      << "NOTE: one ideal strip family only; no second stereo family, "
         "threshold, electronics shaping, or finite PCB geometry yet.\n"
      << "============================================================\n";

  return 0;
}

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
#include "Garfield/DriftLineRKF.hh"
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
  double u0 = 0., v0 = 0., w0 = 0.;
  double u1 = 0., v1 = 0., w1 = 0.;
};

struct IonSegment {
  double u0 = 0., v0 = 0., w0 = 0.;
  double u1 = 0., v1 = 0., w1 = 0.;
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
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
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

std::string StripLabel(const std::string& prefix, const int k) {
  if (k < 0) return prefix + "_m" + std::to_string(-k);
  if (k > 0) return prefix + "_p" + std::to_string(k);
  return prefix + "_0";
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

double EndpointSignal(Garfield::ComponentAnalyticField& field,
                      const std::string& label,
                      const double qFc,
                      const double u0, const double v0, const double w0,
                      const double u1, const double v1, const double w1) {
  const double phi0 = field.WeightingPotential(u0, v0, w0, label);
  const double phi1 = field.WeightingPotential(u1, v1, w1, label);
  return qFc * (phi1 - phi0);
}

void ExportTrajectories(const std::string& prefix,
                        Garfield::AvalancheMicroscopic& avalanche,
                        Garfield::Sensor& sensor,
                        const std::vector<IonBirth>& ionBirths,
                        const int maxElectrons,
                        const int maxIons) {
  if (prefix.empty()) return;

  std::ofstream out(prefix + "_trajectories.csv");
  if (!out) {
    std::cerr << "WARNING: could not open trajectory output.\n";
    return;
  }
  out << "carrier,track_id,point_id,t_ns,u_mm,v_mm,w_mm\n";

  const std::size_t nElectronTracks =
      avalanche.GetNumberOfElectronEndpoints();
  const std::size_t nE =
      std::min<std::size_t>(nElectronTracks,
                            static_cast<std::size_t>(std::max(0, maxElectrons)));

  for (std::size_t ie = 0; ie < nE; ++ie) {
    const std::size_t np =
        avalanche.GetNumberOfElectronDriftLinePoints(ie);
    for (std::size_t ip = 0; ip < np; ++ip) {
      double u = 0., v = 0., w = 0., t = 0.;
      avalanche.GetElectronDriftLinePoint(
          u, v, w, t, ip, ie);
      out << "electron," << ie << "," << ip << ","
          << t << "," << 10. * u << "," << 10. * v << ","
          << 10. * w << "\n";
    }
  }

  if (maxIons <= 0 || ionBirths.empty()) return;

  const std::size_t nWanted =
      std::min<std::size_t>(ionBirths.size(),
                            static_cast<std::size_t>(maxIons));
  const double stride =
      static_cast<double>(ionBirths.size()) /
      static_cast<double>(nWanted);

  Garfield::DriftLineRKF drift(&sensor);
  drift.SetMaximumStepSize(0.005);  // 50 um in cm units.

  for (std::size_t ii = 0; ii < nWanted; ++ii) {
    const std::size_t ib =
        std::min<std::size_t>(
            ionBirths.size() - 1,
            static_cast<std::size_t>(std::floor(ii * stride)));
    const auto& birth = ionBirths[ib];

    if (!drift.DriftIon(birth.u, birth.v, birth.w, birth.t)) continue;

    const std::size_t np = drift.GetNumberOfDriftLinePoints();
    for (std::size_t ip = 0; ip < np; ++ip) {
      double u = 0., v = 0., w = 0., t = 0.;
      drift.GetDriftLinePoint(ip, u, v, w, t);
      out << "ion," << ii << "," << ip << ","
          << t << "," << 10. * u << "," << 10. * v << ","
          << 10. * w << "\n";
    }
  }

  std::cout << "trajectory output        : "
            << prefix << "_trajectories.csv\n";
}

}  // namespace

int main(int argc, char** argv) {
  // Stage B3c.2: opposite-cathode stereo readout.
  //
  // Physical drift geometry:
  //   anode wires at v = 0
  //   minus cathode at v = -2 mm (default)
  //   plus  cathode at v = +4 mm (default)
  //
  // Stereo assignment:
  //   +alpha strip family on the MINUS cathode
  //   -alpha strip family on the PLUS cathode
  //
  // Important modelling point:
  //   - total cathode coupling is computed with the full analytic MWPC cell,
  //     including the discrete anode wires;
  //   - local strip sharing is computed with Garfield's analytic planar-strip
  //     weighting solution for the corresponding cathode gap;
  //   - the local strip charges are then rescaled so that their signed sum
  //     equals the wire-aware complete-cathode signal.
  //
  // This is therefore a controlled HYBRID model. It includes the physical
  // 2/4-mm cathode asymmetry and wire-aware total coupling, while retaining
  // the ideal planar approximation for segmentation within each cathode.

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

  const std::string seedSide =
      ReadStringArg(argc, argv, "--seed-side", "minus");
  const double seedDistanceMm =
      ReadArg(argc, argv, "--seed-distance-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  const double tanAlpha =
      ReadArg(argc, argv, "--tan-alpha", 0.10);
  const double stripPitchMm =
      ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthMm =
      ReadArg(argc, argv, "--strip-width-mm", 1.7);

  auto wireIndices = ParseInts(
      ReadStringArg(argc, argv, "--wire-indices", "-2,-1,0,1,2"));
  auto xiWValues = ParseDoubles(
      ReadStringArg(argc, argv, "--xi-w-values",
                    "-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5"));

  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 7);
  const int events =
      ReadIntArg(argc, argv, "--events", 300);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 91000);
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
                    "stageB_opposite_cathode_stereo.csv");

  const std::string trajectoryPrefix =
      ReadStringArg(argc, argv, "--trajectory-prefix", "");
  const int trajectoryMaxElectrons =
      ReadIntArg(argc, argv, "--trajectory-max-electrons", 120);
  const int trajectoryMaxIons =
      ReadIntArg(argc, argv, "--trajectory-max-ions", 120);

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      tanAlpha <= 0. || stripPitchMm <= 0. || stripWidthMm <= 0. ||
      stripWidthMm > stripPitchMm + 1.e-12 ||
      wireIndices.empty() || xiWValues.empty() ||
      halfStrips < 2 || events < 1 || baseSeed < 0 ||
      avalancheLimit < 1 || observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-B3c.2 parameters.\n";
    return 2;
  }
  if (seedSide != "minus" && seedSide != "plus") {
    std::cerr << "--seed-side must be 'minus' or 'plus'.\n";
    return 2;
  }

  const double seedSign = seedSide == "plus" ? +1. : -1.;
  const double seedVCm = seedSign * 0.1 * seedDistanceMm;
  const double seedGapCm =
      seedSide == "plus" ? gapPlusCm : gapMinusCm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;
  if (std::abs(seedVCm) <= wireRadiusCm ||
      std::abs(seedVCm) >= seedGapCm) {
    std::cerr << "Seed must be outside the wire and inside the selected gap.\n";
    return 2;
  }

  const double alpha = std::atan(tanAlpha);
  const double cosA = std::cos(alpha);
  const double sinA = std::sin(alpha);
  const double alphaDeg = alpha * 180. / Pi;

  // ------------------------------------------------------------------
  // Physical MWPC field.
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

  Garfield::ComponentAnalyticField physicalField;
  physicalField.SetMedium(&gas);

  constexpr int halfNumberOfWires = 6;
  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    physicalField.AddWire(
        i * wirePitchCm, 0., wireDiameterCm, hv, "anode");
  }
  physicalField.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  physicalField.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  physicalField.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor physicalSensor;
  physicalSensor.AddComponent(&physicalField);

  const double stripPitchCm = 0.1 * stripPitchMm;
  const double stripWidthCm = 0.1 * stripWidthMm;
  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 4.0) * stripPitchCm;

  physicalSensor.SetArea(
      -uExtent, -gapMinusCm, -wExtent,
      +uExtent, +gapPlusCm, +wExtent);

  // ------------------------------------------------------------------
  // One ideal segmented weighting model per physical cathode.
  // Minus plane carries +alpha; plus plane carries -alpha.
  // ------------------------------------------------------------------
  Garfield::ComponentAnalyticField weightingMinus;
  weightingMinus.AddPlaneY(0., 1., "minus_back");
  weightingMinus.AddPlaneY(-gapMinusCm, 0., "minus_front");

  Garfield::ComponentAnalyticField weightingPlus;
  weightingPlus.AddPlaneY(0., 1., "plus_back");
  weightingPlus.AddPlaneY(+gapPlusCm, 0., "plus_front");

  std::vector<int> stripIds;
  std::vector<std::string> labelsMinus;
  std::vector<std::string> labelsPlus;

  for (int k = -halfStrips; k <= halfStrips; ++k) {
    const double center = k * stripPitchCm;
    const double sMin = center - 0.5 * stripWidthCm;
    const double sMax = center + 0.5 * stripWidthCm;

    const std::string lm = StripLabel("minus_strip", k);
    const std::string lp = StripLabel("plus_strip", k);

    weightingMinus.AddStripOnPlaneY(
        'x', -gapMinusCm, sMin, sMax, lm, gapMinusCm);
    weightingPlus.AddStripOnPlaneY(
        'x', +gapPlusCm, sMin, sMax, lp, gapPlusCm);

    stripIds.push_back(k);
    labelsMinus.push_back(lm);
    labelsPlus.push_back(lp);
  }

  // ------------------------------------------------------------------
  // Generate the microscopic avalanche bank once around the central wire.
  // ------------------------------------------------------------------
  std::vector<AvalancheEvent> eventBank;
  eventBank.reserve(events);

  int zeroOrFailed = 0;
  int ionTransportFailed = 0;
  bool trajectoryExported = false;

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

    if (!trajectoryExported && !trajectoryPrefix.empty()) {
      ExportTrajectories(
          trajectoryPrefix, avalanche, physicalSensor, gIonBirths,
          trajectoryMaxElectrons, trajectoryMaxIons);
      trajectoryExported = true;
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
      ++ionTransportFailed;
      continue;
    }

    eventBank.push_back(std::move(evt));
  }

  if (eventBank.empty()) {
    std::cerr << "No usable avalanches were generated.\n";
    return 4;
  }

  // ------------------------------------------------------------------
  // Scan wire index and along-wire avalanche position.
  // ------------------------------------------------------------------
  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out
      << "tan_alpha,alpha_deg,seed_side,wire_index,wire_u_mm,"
      << "xi_w,w_shift_mm,"
      << "x_plusalpha_true_mm,x_minusalpha_true_mm,"
      << "event_index,random_seed,garfield_ions,observation_ns,"
      << "q_cathode_minus_full_fC,q_cathode_plus_full_fC,"
      << "q_minus_local_sum_fC,q_plus_local_sum_fC,"
      << "minus_scale_factor,plus_scale_factor";

  for (const int k : stripIds) {
    out << ",Q_minus_strip_" << k << "_fC"
        << ",A_minus_strip_" << k << "_fC";
  }
  for (const int k : stripIds) {
    out << ",Q_plus_strip_" << k << "_fC"
        << ",A_plus_strip_" << k << "_fC";
  }
  out << "\n";

  int written = 0;
  int scaleFailures = 0;

  for (const int wireIndex : wireIndices) {
    const double uShiftCm = wireIndex * wirePitchCm;
    const double wireUMm = 10. * uShiftCm;

    for (const double xiW : xiWValues) {
      const double wShiftCm = xiW * stripPitchCm;
      const double wShiftMm = 10. * wShiftCm;

      const double xPlusAlphaCm =
          wShiftCm * cosA - uShiftCm * sinA;
      const double xMinusAlphaCm =
          wShiftCm * cosA + uShiftCm * sinA;

      for (const auto& evt : eventBank) {
        double qCathMinus = 0.;
        double qCathPlus = 0.;

        std::vector<double> qMinus(labelsMinus.size(), 0.);
        std::vector<double> qPlus(labelsPlus.size(), 0.);

        for (const auto& e : evt.electrons) {
          const double u0 = e.u0 + uShiftCm;
          const double u1 = e.u1 + uShiftCm;
          const double w0 = e.w0 + wShiftCm;
          const double w1 = e.w1 + wShiftCm;

          qCathMinus += EndpointSignal(
              physicalField, "cathode_minus", -ElementaryChargeFc,
              u0, e.v0, w0, u1, e.v1, w1);
          qCathPlus += EndpointSignal(
              physicalField, "cathode_plus", -ElementaryChargeFc,
              u0, e.v0, w0, u1, e.v1, w1);

          const double zMinus0 = w0 * cosA - u0 * sinA;
          const double zMinus1 = w1 * cosA - u1 * sinA;
          const double zPlus0 = w0 * cosA + u0 * sinA;
          const double zPlus1 = w1 * cosA + u1 * sinA;

          for (std::size_t j = 0; j < stripIds.size(); ++j) {
            qMinus[j] += EndpointSignal(
                weightingMinus, labelsMinus[j], -ElementaryChargeFc,
                0., e.v0, zMinus0, 0., e.v1, zMinus1);
            qPlus[j] += EndpointSignal(
                weightingPlus, labelsPlus[j], -ElementaryChargeFc,
                0., e.v0, zPlus0, 0., e.v1, zPlus1);
          }
        }

        for (const auto& ion : evt.ions) {
          const double u0 = ion.u0 + uShiftCm;
          const double u1 = ion.u1 + uShiftCm;
          const double w0 = ion.w0 + wShiftCm;
          const double w1 = ion.w1 + wShiftCm;

          qCathMinus += EndpointSignal(
              physicalField, "cathode_minus", +ElementaryChargeFc,
              u0, ion.v0, w0, u1, ion.v1, w1);
          qCathPlus += EndpointSignal(
              physicalField, "cathode_plus", +ElementaryChargeFc,
              u0, ion.v0, w0, u1, ion.v1, w1);

          const double zMinus0 = w0 * cosA - u0 * sinA;
          const double zMinus1 = w1 * cosA - u1 * sinA;
          const double zPlus0 = w0 * cosA + u0 * sinA;
          const double zPlus1 = w1 * cosA + u1 * sinA;

          for (std::size_t j = 0; j < stripIds.size(); ++j) {
            qMinus[j] += EndpointSignal(
                weightingMinus, labelsMinus[j], +ElementaryChargeFc,
                0., ion.v0, zMinus0, 0., ion.v1, zMinus1);
            qPlus[j] += EndpointSignal(
                weightingPlus, labelsPlus[j], +ElementaryChargeFc,
                0., ion.v0, zPlus0, 0., ion.v1, zPlus1);
          }
        }

        double qMinusLocalSum = 0.;
        double qPlusLocalSum = 0.;
        for (std::size_t j = 0; j < stripIds.size(); ++j) {
          qMinusLocalSum += qMinus[j];
          qPlusLocalSum += qPlus[j];
        }

        constexpr double MinDenom = 1.e-18;
        if (std::abs(qMinusLocalSum) < MinDenom ||
            std::abs(qPlusLocalSum) < MinDenom) {
          ++scaleFailures;
          continue;
        }

        const double minusScale =
            qCathMinus / qMinusLocalSum;
        const double plusScale =
            qCathPlus / qPlusLocalSum;

        for (auto& q : qMinus) q *= minusScale;
        for (auto& q : qPlus) q *= plusScale;

        out
            << tanAlpha << "," << alphaDeg << "," << seedSide << ","
            << wireIndex << "," << wireUMm << ","
            << xiW << "," << wShiftMm << ","
            << 10. * xPlusAlphaCm << ","
            << 10. * xMinusAlphaCm << ","
            << evt.eventIndex << "," << evt.randomSeed << ","
            << evt.garfieldIons << "," << observationNs << ","
            << qCathMinus << "," << qCathPlus << ","
            << qMinusLocalSum << "," << qPlusLocalSum << ","
            << minusScale << "," << plusScale;

        for (const double q : qMinus) {
          out << "," << q << "," << std::abs(q);
        }
        for (const double q : qPlus) {
          out << "," << q << "," << std::abs(q);
        }
        out << "\n";

        ++written;
      }
    }
  }

  std::cout << std::fixed << std::setprecision(6)
            << "\n=== STAGE B3c.2: OPPOSITE-CATHODE STEREO HYBRID ===\n"
            << "seed side                : " << seedSide << "\n"
            << "gaps (-/+)               : "
            << 10. * gapMinusCm << " / "
            << 10. * gapPlusCm << " mm\n"
            << "stereo families          : +alpha on minus, -alpha on plus\n"
            << "tan(alpha), alpha        : "
            << tanAlpha << ", " << alphaDeg << " deg\n"
            << "strip pitch / width      : "
            << stripPitchMm << " / " << stripWidthMm << " mm\n"
            << "wire pitch / diameter    : "
            << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "observation time         : "
            << observationNs << " ns\n"
            << "usable avalanche bank    : "
            << eventBank.size() << " / " << events << "\n"
            << "ion transport failures   : "
            << ionTransportFailed << "\n"
            << "rows written             : "
            << written << "\n"
            << "scale failures           : "
            << scaleFailures << "\n"
            << "output                    : "
            << outputFile << "\n\n"
            << "MODEL: wire-aware COMPLETE-cathode coupling + ideal planar\n"
            << "       segmented strip sharing, normalized side-by-side.\n"
            << "       This is more physical than B3b but is not yet a full\n"
            << "       numerical 3-D weighting-field solution for PCB strips.\n"
            << "==============================================================\n";

  return written > 0 ? 0 : 4;
}

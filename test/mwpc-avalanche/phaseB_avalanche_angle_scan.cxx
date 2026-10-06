#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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

double WrapDeg(double x) {
  while (x <= -180.) x += 360.;
  while (x > 180.) x -= 360.;
  return x;
}

}  // namespace

int main(int argc, char** argv) {
  // Stage B3c control requested after the PI discussion:
  // scan the direction from which a seed electron approaches one anode wire.
  //
  // Coordinates:
  //   u = across wires
  //   v = chamber normal
  //   w = along wires
  //
  // Seed convention:
  //   theta = 0 deg   -> +u
  //   theta = 90 deg  -> +v (towards the 4-mm gap)
  //   theta = 180 deg -> -u
  //   theta = 270 deg -> -v (towards the 2-mm gap)
  //
  // The electron is seeded on a circle around the central wire and Garfield
  // drifts it self-consistently into the high-field avalanche region.
  //
  // The output tests two distinct questions:
  //   (1) Does the avalanche ionisation remain localised near the incoming
  //       azimuth around the wire?
  //   (2) How does the induced complete-cathode signal split between v=-2 mm
  //       and v=+4 mm as that azimuth is changed?

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

  const double seedRadiusMm =
      ReadArg(argc, argv, "--seed-radius-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  const auto thetaValues = ParseDoubles(
      ReadStringArg(argc, argv, "--theta-deg",
                    "0,30,60,90,120,150,180,210,240,270,300,330"));

  const int eventsPerAngle =
      ReadIntArg(argc, argv, "--events-per-angle", 30);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 97000);
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
                    "stageB_avalanche_angle_scan.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      seedRadiusMm <= 0. || seedEnergyEv <= 0. ||
      thetaValues.empty() || eventsPerAngle < 1 ||
      baseSeed < 0 || avalancheLimit < 1 ||
      observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid angle-scan parameters.\n";
    return 2;
  }

  const double wireRadiusMm = 0.5 * 1.e-3 *
      ReadArg(argc, argv, "--wire-diam-um", 30.0);
  if (seedRadiusMm <= wireRadiusMm ||
      seedRadiusMm >= 0.5 * 10. * wirePitchCm ||
      seedRadiusMm >= 10. * gapMinusCm ||
      seedRadiusMm >= 10. * gapPlusCm) {
    std::cerr
        << "Seed radius must be outside the wire and inside both gaps, "
        << "and smaller than half the wire pitch.\n";
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

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);

  constexpr int halfNumberOfWires = 6;
  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    field.AddWire(i * wirePitchCm, 0.,
                  wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);

  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  sensor.SetArea(
      -uExtent, -gapMinusCm, -1.0,
      +uExtent, +gapPlusCm, +1.0);

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out
      << "theta_seed_deg,seed_u_mm,seed_v_mm,"
      << "event_index,random_seed,garfield_electrons,garfield_ions,"
      << "electron_end_max_ns,observation_ns,"
      << "birth_mean_u_um,birth_mean_v_um,birth_mean_radius_um,"
      << "birth_mean_theta_deg,birth_theta_sigma_deg,"
      << "birth_theta_minus_seed_deg,"
      << "q_minus_electron_fC,q_minus_ion_fC,q_minus_total_fC,"
      << "q_plus_electron_fC,q_plus_ion_fC,q_plus_total_fC,"
      << "abs_minus_fraction_of_two_cathodes,"
      << "abs_plus_fraction_of_two_cathodes\n";

  int requested = 0;
  int written = 0;
  int failed = 0;
  int ionTransportFailed = 0;

  for (std::size_t itheta = 0; itheta < thetaValues.size(); ++itheta) {
    const double thetaDeg = thetaValues[itheta];
    const double theta = thetaDeg * Pi / 180.;

    const double seedUCm =
        0.1 * seedRadiusMm * std::cos(theta);
    const double seedVCm =
        0.1 * seedRadiusMm * std::sin(theta);

    for (int iev = 0; iev < eventsPerAngle; ++iev) {
      ++requested;
      const int randomSeed =
          baseSeed + static_cast<int>(itheta) * 10000 + iev;

      Garfield::RandomEngineRoot engine(
          static_cast<unsigned int>(randomSeed));
      Garfield::Random::SetEngine(engine);

      gIonBirths.clear();

      Garfield::AvalancheMicroscopic avalanche(&sensor);
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

      if (!avalancheOk || nIons <= 0 || gIonBirths.empty()) {
        ++failed;
        continue;
      }

      std::vector<ElectronEndpoint> electrons;
      const std::size_t nEndpoints =
          avalanche.GetNumberOfElectronEndpoints();
      electrons.reserve(nEndpoints);

      double electronEndMaxNs = 0.;
      for (std::size_t ie = 0; ie < nEndpoints; ++ie) {
        double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
        double u1 = 0., v1 = 0., w1 = 0., t1 = 0., e1 = 0.;
        int status = 0;
        avalanche.GetElectronEndpoint(
            ie, u0, v0, w0, t0, e0,
            u1, v1, w1, t1, e1, status);
        electronEndMaxNs = std::max(electronEndMaxNs, t1);
        electrons.push_back({u0, v0, w0, u1, v1, w1});
      }

      if (electronEndMaxNs > observationNs) {
        ++failed;
        continue;
      }

      std::vector<IonSegment> ions;
      ions.reserve(gIonBirths.size());
      bool transportGood = true;

      for (const auto& birth : gIonBirths) {
        if (birth.t >= observationNs) continue;

        std::array<double, 3> x = {birth.u, birth.v, birth.w};
        double t = birth.t;

        if (!AdvanceIonRK4(sensor, x, t, observationNs, ionRkDtNs)) {
          transportGood = false;
          break;
        }

        ions.push_back(
            {birth.u, birth.v, birth.w,
             x[0], x[1], x[2]});
      }

      if (!transportGood) {
        ++ionTransportFailed;
        continue;
      }

      // Circular statistics of avalanche ionisation birth positions.
      double sumU = 0., sumV = 0., sumR = 0.;
      double sumCos = 0., sumSin = 0.;
      for (const auto& birth : gIonBirths) {
        const double r = std::hypot(birth.u, birth.v);
        const double a = std::atan2(birth.v, birth.u);
        sumU += birth.u;
        sumV += birth.v;
        sumR += r;
        sumCos += std::cos(a);
        sumSin += std::sin(a);
      }

      const double nb = static_cast<double>(gIonBirths.size());
      const double meanU = sumU / nb;
      const double meanV = sumV / nb;
      const double meanR = sumR / nb;
      const double meanTheta =
          std::atan2(sumSin, sumCos);
      const double resultant =
          std::hypot(sumCos, sumSin) / nb;
      const double sigmaTheta =
          resultant > 0.
              ? std::sqrt(std::max(0., -2. * std::log(resultant)))
              : Pi;

      const double meanThetaDeg = meanTheta * 180. / Pi;
      const double sigmaThetaDeg = sigmaTheta * 180. / Pi;
      const double deltaThetaDeg =
          WrapDeg(meanThetaDeg - thetaDeg);

      double qMinusE = 0., qMinusI = 0.;
      double qPlusE = 0., qPlusI = 0.;

      for (const auto& e : electrons) {
        qMinusE += EndpointSignal(
            field, "cathode_minus", -ElementaryChargeFc,
            e.u0, e.v0, e.w0, e.u1, e.v1, e.w1);
        qPlusE += EndpointSignal(
            field, "cathode_plus", -ElementaryChargeFc,
            e.u0, e.v0, e.w0, e.u1, e.v1, e.w1);
      }

      for (const auto& ion : ions) {
        qMinusI += EndpointSignal(
            field, "cathode_minus", +ElementaryChargeFc,
            ion.u0, ion.v0, ion.w0, ion.u1, ion.v1, ion.w1);
        qPlusI += EndpointSignal(
            field, "cathode_plus", +ElementaryChargeFc,
            ion.u0, ion.v0, ion.w0, ion.u1, ion.v1, ion.w1);
      }

      const double qMinus = qMinusE + qMinusI;
      const double qPlus = qPlusE + qPlusI;
      const double denom = std::abs(qMinus) + std::abs(qPlus);
      const double fMinus =
          denom > 0. ? std::abs(qMinus) / denom : 0.;
      const double fPlus =
          denom > 0. ? std::abs(qPlus) / denom : 0.;

      out
          << thetaDeg << ","
          << 10. * seedUCm << "," << 10. * seedVCm << ","
          << iev << "," << randomSeed << ","
          << nElectrons << "," << nIons << ","
          << electronEndMaxNs << "," << observationNs << ","
          << 1.e4 * meanU << "," << 1.e4 * meanV << ","
          << 1.e4 * meanR << ","
          << meanThetaDeg << "," << sigmaThetaDeg << ","
          << deltaThetaDeg << ","
          << qMinusE << "," << qMinusI << "," << qMinus << ","
          << qPlusE << "," << qPlusI << "," << qPlus << ","
          << fMinus << "," << fPlus << "\n";

      ++written;
    }
  }

  std::cout << std::fixed << std::setprecision(4)
            << "\n=== STAGE B3c: AVALANCHE AROUND-WIRE ANGLE SCAN ===\n"
            << "gaps (-/+)               : "
            << 10. * gapMinusCm << " / "
            << 10. * gapPlusCm << " mm\n"
            << "wire pitch / diameter    : "
            << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "seed radius              : "
            << seedRadiusMm << " mm\n"
            << "observation time         : "
            << observationNs << " ns\n"
            << "angles                    : ";
  for (const auto x : thetaValues) std::cout << x << " ";
  std::cout
      << "\nevents per angle          : "
      << eventsPerAngle << "\n"
      << "requested                  : "
      << requested << "\n"
      << "written                    : "
      << written << "\n"
      << "failed/zero gain           : "
      << failed << "\n"
      << "ion transport failures     : "
      << ionTransportFailed << "\n"
      << "output                     : "
      << outputFile << "\n"
      << "===========================================================\n";

  return written > 0 ? 0 : 4;
}

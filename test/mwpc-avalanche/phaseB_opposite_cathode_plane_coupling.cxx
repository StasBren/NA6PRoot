#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
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

}  // namespace

int main(int argc, char** argv) {
  // Stage B3c control:
  // Use the ACTUAL analytic MWPC cell (discrete anode wires plus both cathode
  // planes) to calculate the weighting potential of each COMPLETE cathode.
  //
  // This deliberately comes before segmented opposite-cathode stereo strips.
  // Garfield's analytic strip weighting formula is a parallel-plate formula
  // and does not include screening by the discrete anode wires.  The complete
  // cathode weighting potential, in contrast, belongs to the same analytic
  // cell that contains the real wires.  It is therefore the clean control for
  // asking how much one avalanche couples to each physical cathode.
  //
  // Coordinates:
  //   u = across wires, v = chamber normal, w = along wires.

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

  const int events =
      ReadIntArg(argc, argv, "--events", 300);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 83000);
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
                    "stageB_opposite_cathode_plane_coupling.csv");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      events < 1 || baseSeed < 0 || avalancheLimit < 1 ||
      observationNs <= 0. || ionRkDtNs <= 0.) {
    std::cerr << "Invalid Stage-B3c parameters.\n";
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
    std::cerr << "Seed must lie outside the anode wire and inside the gas gap.\n";
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
  const double wExtent = 2.0;  // 20 mm; complete planes are w-invariant.
  sensor.SetArea(-uExtent, -gapMinusCm, -wExtent,
                 +uExtent, +gapPlusCm, +wExtent);

  std::ofstream out(outputFile);
  if (!out) {
    std::cerr << "Could not open output file " << outputFile << "\n";
    return 4;
  }

  out
      << "event_index,random_seed,seed_side,garfield_electrons,garfield_ions,"
      << "electron_end_max_ns,observation_ns,"
      << "q_minus_electron_fC,q_minus_ion_fC,q_minus_total_fC,"
      << "q_plus_electron_fC,q_plus_ion_fC,q_plus_total_fC,"
      << "q_anode_electron_fC,q_anode_ion_fC,q_anode_total_fC,"
      << "abs_minus_fraction_of_two_cathodes,abs_plus_fraction_of_two_cathodes\n";

  int zeroOrFailed = 0;
  int ionTransportFailed = 0;
  int electronTooLate = 0;
  int written = 0;

  for (int iev = 0; iev < events; ++iev) {
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
      ++electronTooLate;
      continue;
    }

    std::vector<IonSegment> ions;
    ions.reserve(gIonBirths.size());
    bool transportGood = true;

    for (const auto& birth : gIonBirths) {
      if (birth.t >= observationNs) continue;

      std::array<double, 3> x =
          {birth.u, birth.v, birth.w};
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

    double qMinusE = 0., qMinusI = 0.;
    double qPlusE = 0., qPlusI = 0.;
    double qAnodeE = 0., qAnodeI = 0.;

    for (const auto& e : electrons) {
      qMinusE += EndpointSignal(
          field, "cathode_minus", -ElementaryChargeFc,
          e.u0, e.v0, e.w0, e.u1, e.v1, e.w1);
      qPlusE += EndpointSignal(
          field, "cathode_plus", -ElementaryChargeFc,
          e.u0, e.v0, e.w0, e.u1, e.v1, e.w1);
      qAnodeE += EndpointSignal(
          field, "anode", -ElementaryChargeFc,
          e.u0, e.v0, e.w0, e.u1, e.v1, e.w1);
    }

    for (const auto& ion : ions) {
      qMinusI += EndpointSignal(
          field, "cathode_minus", +ElementaryChargeFc,
          ion.u0, ion.v0, ion.w0, ion.u1, ion.v1, ion.w1);
      qPlusI += EndpointSignal(
          field, "cathode_plus", +ElementaryChargeFc,
          ion.u0, ion.v0, ion.w0, ion.u1, ion.v1, ion.w1);
      qAnodeI += EndpointSignal(
          field, "anode", +ElementaryChargeFc,
          ion.u0, ion.v0, ion.w0, ion.u1, ion.v1, ion.w1);
    }

    const double qMinus = qMinusE + qMinusI;
    const double qPlus = qPlusE + qPlusI;
    const double qAnode = qAnodeE + qAnodeI;

    const double cathAbs = std::abs(qMinus) + std::abs(qPlus);
    const double fMinus =
        cathAbs > 0. ? std::abs(qMinus) / cathAbs : 0.;
    const double fPlus =
        cathAbs > 0. ? std::abs(qPlus) / cathAbs : 0.;

    out
        << iev << "," << randomSeed << "," << seedSide << ","
        << nElectrons << "," << nIons << ","
        << electronEndMaxNs << "," << observationNs << ","
        << qMinusE << "," << qMinusI << "," << qMinus << ","
        << qPlusE << "," << qPlusI << "," << qPlus << ","
        << qAnodeE << "," << qAnodeI << "," << qAnode << ","
        << fMinus << "," << fPlus << "\n";

    ++written;
  }

  std::cout << std::fixed << std::setprecision(6)
            << "\n=== STAGE B3c: WIRE-AWARE OPPOSITE-CATHODE COUPLING ===\n"
            << "seed side                : " << seedSide << "\n"
            << "gaps (-/+)               : "
            << 10. * gapMinusCm << " / " << 10. * gapPlusCm << " mm\n"
            << "wire pitch / diameter    : "
            << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "anode voltage            : " << hv << " V\n"
            << "observation time         : " << observationNs << " ns\n"
            << "requested events         : " << events << "\n"
            << "written usable events    : " << written << "\n"
            << "zero/failed              : " << zeroOrFailed << "\n"
            << "electron too late        : " << electronTooLate << "\n"
            << "ion transport failures   : " << ionTransportFailed << "\n"
            << "output                    : " << outputFile << "\n\n"
            << "NOTE: cathode weighting potentials are calculated by the full\n"
            << "      analytic MWPC cell containing the discrete anode wires.\n"
            << "NOTE: this is a complete-cathode coupling control, not yet a\n"
            << "      segmented stereo-strip reconstruction.\n"
            << "==============================================================\n";

  return written > 0 ? 0 : 4;
}

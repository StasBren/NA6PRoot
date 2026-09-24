#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return std::atof(argv[i + 1]);
  }
  return defaultValue;
}

int ReadIntArg(const int argc, char** argv, const std::string& key,
               const int defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return std::atoi(argv[i + 1]);
  }
  return defaultValue;
}

double Quantile(std::vector<int> values, const double q) {
  if (values.empty()) return 0.;
  std::sort(values.begin(), values.end());
  const double pos = q * static_cast<double>(values.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(pos));
  const auto hi = static_cast<std::size_t>(std::ceil(pos));
  if (lo == hi) return static_cast<double>(values[lo]);
  const double f = pos - static_cast<double>(lo);
  return (1. - f) * values[lo] + f * values[hi];
}

}  // namespace

int main(int argc, char** argv) {
  const double pitchCm =
      0.1 * ReadArg(argc, argv, "--pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);
  const double bTesla = ReadArg(argc, argv, "--b", 0.0);
  const int events = ReadIntArg(argc, argv, "--events", 1000);
  const double x0Cm =
      0.1 * ReadArg(argc, argv, "--x0-mm", 1.0);
  const double y0Cm =
      0.1 * ReadArg(argc, argv, "--y0-mm", 0.75 * 10. * gapPlusCm);

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. || events <= 0) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (y0Cm <= -gapMinusCm || y0Cm >= gapPlusCm) {
    std::cerr << "Starting y0 is outside the gas gap.\n";
    return 2;
  }

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);

  constexpr int halfNumberOfWires = 4;
  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    field.AddWire(i * pitchCm, 0., wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");

  // Current nominal mapping:
  // detector B_y -> chamber-local u -> Garfield x.
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double xExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -5.0,
                  xExtent, +gapPlusCm, 5.0);

  const double wireRadius = 0.5 * wireDiameterCm;
  constexpr unsigned int avalancheLimit = 200000;

  std::ofstream out("gain_fluctuations.csv");
  out << "event,avalanche_electrons,avalanche_ions,collected_on_wire,"
         "attached,other_endpoints,mean_collection_time_ns\n";

  std::vector<int> collectedValues;
  collectedValues.reserve(events);

  double sumGain = 0.;
  double sumGain2 = 0.;
  int zeroCollected = 0;
  int sizeLimitLike = 0;

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "\n=== PHASE A GAIN-FLUCTUATION ENSEMBLE ===\n"
            << "events              : " << events << "\n"
            << "start (u,v)         : (" << 10. * x0Cm << ", "
            << 10. * y0Cm << ") mm\n"
            << "wire pitch          : " << 10. * pitchCm << " mm\n"
            << "wire diameter       : " << 1.e4 * wireDiameterCm << " um\n"
            << "cathode gaps        : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "anode voltage       : " << hv << " V\n"
            << "detector B_y        : " << bTesla << " T\n"
            << "-----------------------------------------------\n";

  for (int event = 0; event < events; ++event) {
    // Construct a fresh AvalancheMicroscopic object for every event so that
    // all endpoint containers and avalanche counters are guaranteed to refer
    // only to this event. The gas/field/sensor objects are shared.
    Garfield::AvalancheMicroscopic avalanche;
    avalanche.SetSensor(&sensor);
    avalanche.EnableAvalancheSizeLimit(avalancheLimit);

    const bool ok =
        avalanche.AvalancheElectron(x0Cm, y0Cm, 0., 0., 0.1);

    int ne = 0;
    int ni = 0;
    avalanche.GetAvalancheSize(ne, ni);

    int nCollected = 0;
    int nAttached = 0;
    int nOther = 0;
    double sumT = 0.;
    int nTimed = 0;

    const auto nEndpoints = avalanche.GetNumberOfElectronEndpoints();
    for (std::size_t i = 0; i < nEndpoints; ++i) {
      double xa, ya, za, ta, ea;
      double x1, y1, z1, t1, e1;
      int status = 0;
      avalanche.GetElectronEndpoint(i, xa, ya, za, ta, ea,
                                    x1, y1, z1, t1, e1, status);

      const int nearest =
          static_cast<int>(std::lround(x1 / pitchCm));
      const double wireX = nearest * pitchCm;
      const double r = std::hypot(x1 - wireX, y1);

      if (r < 1.5 * wireRadius) {
        ++nCollected;
        if (t1 >= 0.) {
          sumT += t1;
          ++nTimed;
        }
      } else if (status == -7) {
        ++nAttached;
      } else {
        ++nOther;
      }
    }

    if (!ok || nCollected == 0) ++zeroCollected;
    if (ne >= static_cast<int>(avalancheLimit) ||
        nEndpoints >= avalancheLimit) {
      ++sizeLimitLike;
    }

    const double meanT = nTimed > 0 ? sumT / nTimed : 0.;

    collectedValues.push_back(nCollected);
    sumGain += static_cast<double>(nCollected);
    sumGain2 += static_cast<double>(nCollected) *
                static_cast<double>(nCollected);

    out << event << "," << ne << "," << ni << ","
        << nCollected << "," << nAttached << "," << nOther << ","
        << meanT << "\n";

    if (event < 10) {
      std::cout << "event " << std::setw(4) << event
                << " : produced_e=" << std::setw(6) << ne
                << "  collected=" << std::setw(6) << nCollected
                << "  attached=" << std::setw(5) << nAttached
                << "  other=" << std::setw(4) << nOther
                << "  <t>=" << std::setw(8) << meanT << " ns\n";
    }
  }

  const double mean = sumGain / static_cast<double>(events);
  const double variance =
      events > 1
          ? (sumGain2 - static_cast<double>(events) * mean * mean) /
                static_cast<double>(events - 1)
          : 0.;
  const double sigma = std::sqrt(std::max(0., variance));
  const double cv = mean > 0. ? sigma / mean : 0.;

  std::cout << "-----------------------------------------------\n"
            << "effective gain = electrons collected on anode wire\n"
            << "mean gain           : " << mean << "\n"
            << "sample sigma        : " << sigma << "\n"
            << "sigma / mean        : " << cv << "\n"
            << "median              : " << Quantile(collectedValues, 0.50) << "\n"
            << "10% quantile        : " << Quantile(collectedValues, 0.10) << "\n"
            << "90% quantile        : " << Quantile(collectedValues, 0.90) << "\n"
            << "zero-collected frac : "
            << static_cast<double>(zeroCollected) / events << "\n"
            << "size-limit-like evt : " << sizeLimitLike << "\n"
            << "wrote               : gain_fluctuations.csv\n"
            << "===============================================\n";

  return 0;
}

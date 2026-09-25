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

std::string ReadStringArg(const int argc, char** argv, const std::string& key,
                          const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return argv[i + 1];
  }
  return defaultValue;
}

std::vector<double> ParseList(const std::string& csv) {
  std::vector<double> values;
  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (token.empty()) continue;
    values.push_back(std::stod(token));
  }
  return values;
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
  const int trials = ReadIntArg(argc, argv, "--trials", 2000);
  const double x0Cm =
      0.1 * ReadArg(argc, argv, "--x0-mm", 1.0);
  const double y0Cm =
      0.1 * ReadArg(argc, argv, "--y0-mm", 0.75 * 10. * gapPlusCm);
  const std::string bListText =
      ReadStringArg(argc, argv, "--b-list-t", "0,0.25,0.5,0.75,1.0");
  const std::string outName =
      ReadStringArg(argc, argv, "--output", "magnetic_drift_scan.csv");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. || trials <= 0) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (x0Cm < -0.5 * pitchCm || x0Cm > 0.5 * pitchCm) {
    std::cerr << "Starting x0 is outside the central wire cell [-p/2,+p/2].\n";
    return 2;
  }
  if (y0Cm <= -gapMinusCm || y0Cm >= gapPlusCm) {
    std::cerr << "Starting y0 is outside the gas gap.\n";
    return 2;
  }

  std::vector<double> bValues;
  try {
    bValues = ParseList(bListText);
  } catch (const std::exception& e) {
    std::cerr << "Could not parse --b-list-t: " << e.what() << "\n";
    return 2;
  }
  if (bValues.empty()) {
    std::cerr << "No magnetic-field values supplied.\n";
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

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double xExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -5.0,
                  xExtent, +gapPlusCm, 5.0);

  const double wireRadius = 0.5 * wireDiameterCm;

  std::ofstream out(outName);
  out << "B_T,trials,wire_m1,wire_0,wire_p1,attachment,other,"
         "mean_w_mm,sigma_w_mm,mean_t_ns\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== PHASE A MAGNETIC-DRIFT SCAN (NO AVALANCHE) ===\n"
            << "start (u,v,w) = (" << 10. * x0Cm << ", "
            << 10. * y0Cm << ", 0) mm\n"
            << "pitch = " << 10. * pitchCm
            << " mm, trials/B = " << trials << "\n"
            << "mapping: detector B_y -> local u -> Garfield x\n"
            << "observable: final Garfield z = local w = along-wire displacement\n\n"
            << "B[T]   wire-1  wire0  wire+1  attach  other"
            << "   <w>[mm]  sigma_w[mm]  <t>[ns]\n";

  for (const double bTesla : bValues) {
    // Nominal detector mapping:
    // global vertical MNP33 field (+Y) -> chamber-local +u -> Garfield +x.
    field.SetMagneticField(bTesla, 0., 0.);

    Garfield::AvalancheMicroscopic drift;
    drift.SetSensor(&sensor);

    int nM1 = 0, n0 = 0, nP1 = 0;
    int nAttach = 0, nOther = 0;
    double sumW = 0.;
    double sumW2 = 0.;
    double sumT = 0.;
    int nCollected = 0;
    int nTimed = 0;

    for (int it = 0; it < trials; ++it) {
      drift.DriftElectron(x0Cm, y0Cm, 0., 0., 0.1);

      const auto nEndpoints = drift.GetNumberOfElectronEndpoints();
      if (nEndpoints < 1) {
        ++nOther;
        continue;
      }

      double xa, ya, za, ta, ea;
      double x1, y1, z1, t1, e1;
      int status = 0;
      drift.GetElectronEndpoint(0, xa, ya, za, ta, ea,
                                x1, y1, z1, t1, e1, status);

      const int nearest =
          static_cast<int>(std::lround(x1 / pitchCm));
      const double wireX = nearest * pitchCm;
      const double r = std::hypot(x1 - wireX, y1);

      if (r < 1.5 * wireRadius) {
        if (nearest == -1) ++nM1;
        else if (nearest == 0) ++n0;
        else if (nearest == +1) ++nP1;
        else ++nOther;

        // Garfield z is chamber-local w, i.e. along the anode wire.
        sumW += z1;
        sumW2 += z1 * z1;
        ++nCollected;

        if (t1 >= 0.) {
          sumT += t1;
          ++nTimed;
        }
      } else if (status == -7) {
        ++nAttach;
      } else {
        ++nOther;
      }
    }

    const double meanW = nCollected > 0 ? sumW / nCollected : 0.;
    const double varW =
        nCollected > 1
            ? (sumW2 - nCollected * meanW * meanW) / (nCollected - 1)
            : 0.;
    const double sigmaW = std::sqrt(std::max(0., varW));
    const double meanT = nTimed > 0 ? sumT / nTimed : 0.;

    out << bTesla << "," << trials << ","
        << nM1 << "," << n0 << "," << nP1 << ","
        << nAttach << "," << nOther << ","
        << 10. * meanW << "," << 10. * sigmaW << ","
        << meanT << "\n";

    std::cout << std::setw(4) << bTesla << "   "
              << std::setw(6) << nM1 << "  "
              << std::setw(5) << n0 << "  "
              << std::setw(6) << nP1 << "  "
              << std::setw(6) << nAttach << "  "
              << std::setw(5) << nOther << "   "
              << std::setw(8) << 10. * meanW << "  "
              << std::setw(11) << 10. * sigmaW << "  "
              << std::setw(8) << meanT << "\n";
  }

  std::cout << "\nWrote " << outName << "\n";
  return 0;
}

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

std::vector<double> ParseMmListToCm(const std::string& csv) {
  std::vector<double> values;
  if (csv.empty()) return values;

  std::stringstream ss(csv);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (token.empty()) continue;
    values.push_back(0.1 * std::stod(token));
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
  const double bTesla = ReadArg(argc, argv, "--b", 0.0);
  const int trials = ReadIntArg(argc, argv, "--trials", 200);
  const int nSteps = ReadIntArg(argc, argv, "--steps", 17);
  const double y0Cm =
      0.1 * ReadArg(argc, argv, "--y0-mm", 0.75 * 10. * gapPlusCm);
  const std::string x0ListText =
      ReadStringArg(argc, argv, "--x0-list-mm", "");
  const std::string outName =
      ReadStringArg(argc, argv, "--output", "collection_scan.csv");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      trials <= 0 || (x0ListText.empty() && nSteps < 2)) {
    std::cerr << "Invalid scan parameters.\n";
    return 2;
  }
  if (y0Cm <= -gapMinusCm || y0Cm >= gapPlusCm) {
    std::cerr << "Starting y0 is outside the gas gap.\n";
    return 2;
  }

  std::vector<double> x0ValuesCm;
  try {
    x0ValuesCm = ParseMmListToCm(x0ListText);
  } catch (const std::exception& e) {
    std::cerr << "Could not parse --x0-list-mm: " << e.what() << "\n";
    return 2;
  }

  if (!x0ValuesCm.empty()) {
    const double halfPitch = 0.5 * pitchCm;
    for (const double x0 : x0ValuesCm) {
      if (x0 < -halfPitch || x0 > halfPitch) {
        std::cerr << "Requested x0 = " << 10. * x0
                  << " mm lies outside the central wire cell [-p/2,+p/2].\n";
        return 2;
      }
    }
  } else {
    const double xMin = -0.5 * pitchCm;
    const double xMax = +0.5 * pitchCm;
    x0ValuesCm.reserve(nSteps);
    for (int ix = 0; ix < nSteps; ++ix) {
      const double f = static_cast<double>(ix) / (nSteps - 1);
      x0ValuesCm.push_back(xMin + f * (xMax - xMin));
    }
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

  // Detector B_y -> Garfield B_x.
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double xExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -5.0,
                  xExtent, +gapPlusCm, 5.0);

  Garfield::AvalancheMicroscopic drift;
  drift.SetSensor(&sensor);

  std::ofstream out(outName);
  out << "x0_mm,trials,wire_m1,wire_0,wire_p1,attachment,left_area,"
         "other,mean_t_ns\n";

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "\n=== PHASE A COLLECTION SCAN (NO AVALANCHE) ===\n"
            << "pitch = " << 10. * pitchCm << " mm, y0 = "
            << 10. * y0Cm << " mm, trials/x = " << trials
            << ", B_y = " << bTesla << " T\n"
            << "scan mode = "
            << (x0ListText.empty() ? "uniform across one pitch"
                                   : "explicit x0 list")
            << "\n"
            << "x0_mm    wire-1    wire0    wire+1    attach    other\n";

  const double wireRadius = 0.5 * wireDiameterCm;

  for (const double x0 : x0ValuesCm) {

    int nM1 = 0, n0 = 0, nP1 = 0;
    int nAttach = 0, nLeftArea = 0, nOther = 0;
    double sumT = 0.;
    int nTimed = 0;

    for (int it = 0; it < trials; ++it) {
      drift.DriftElectron(x0, y0Cm, 0., 0., 0.1);

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

      if (t1 >= 0.) {
        sumT += t1;
        ++nTimed;
      }

      const int nearest =
          static_cast<int>(std::lround(x1 / pitchCm));
      const double wireX = nearest * pitchCm;
      const double r = std::hypot(x1 - wireX, y1);

      if (r < 1.5 * wireRadius) {
        if (nearest == -1) ++nM1;
        else if (nearest == 0) ++n0;
        else if (nearest == +1) ++nP1;
        else ++nOther;
      } else if (status == -7) {
        ++nAttach;
      } else if (status == -1) {
        ++nLeftArea;
      } else {
        ++nOther;
      }
    }

    const double meanT = nTimed > 0 ? sumT / nTimed : 0.;
    out << 10. * x0 << "," << trials << ","
        << nM1 << "," << n0 << "," << nP1 << ","
        << nAttach << "," << nLeftArea << "," << nOther << ","
        << meanT << "\n";

    std::cout << std::setw(6) << 10. * x0 << "    "
              << std::setw(6) << nM1 << "    "
              << std::setw(5) << n0 << "    "
              << std::setw(6) << nP1 << "    "
              << std::setw(6) << nAttach << "    "
              << std::setw(5) << (nLeftArea + nOther) << "\n";
  }

  std::cout << "\nWrote " << outName << "\n";
  return 0;
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {

// Garfield analytic wire cells use wires parallel to the local z axis.
// For the NA60+/DiCE chamber we use the following local mapping:
//
//   Garfield x  = detector Y  (vertical, across the horizontal wires)
//   Garfield y  = detector Z  (beam / cathode-gap direction)
//   Garfield z  = detector X  (horizontal, along the wires)
//
// Therefore the physical MNP33 field +Y maps to +x in this local cell.

struct Config {
  // All geometry values are runtime-configurable because the chamber design
  // (especially the cathode gaps and readout topology) is not frozen.
  double wirePitchCm = 0.4;          // default: 4 mm
  double wireDiameterCm = 0.003;     // default: 30 um
  double gapMinusCm = 0.2;           // cathode at y = -gapMinus
  double gapPlusCm = 0.4;            // cathode at y = +gapPlus
  double anodeVoltageV = 1800.;      // Prototype-3 starting point
  double magneticFieldT = 0.;        // Phase-A first run: B = 0
  int halfNumberOfWires = 4;         // 9 wires total

  // Initial electron position in Garfield local coordinates.
  double x0Cm = 0.10;
  double y0Cm = 0.30;
  double z0Cm = 0.0;
  double t0Ns = 0.0;
  double e0Ev = 0.10;

  unsigned int avalancheLimit = 200000;
};

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return std::atof(argv[i + 1]);
  }
  return defaultValue;
}

void PrintUsage(const char* exe) {
  std::cout
      << "Usage: " << exe << " [options]\n"
      << "  --hv V              anode voltage [V] (default 1800)\n"
      << "  --b T               detector B_y [T] (default 0)\n"
      << "  --gap-minus-mm mm    wire-to-minus-cathode gap [mm] (default 2)\n"
      << "  --gap-plus-mm mm     wire-to-plus-cathode gap [mm] (default 4)\n"
      << "  --pitch-mm mm        wire pitch [mm] (default 4)\n"
      << "  --wire-diam-um um    wire diameter [um] (default 30)\n"
      << "  --x0 cm              initial Garfield-x = detector-Y [cm]\n"
      << "  --y0 cm              initial Garfield-y = detector-Z [cm]\n"
      << "\n"
      << "The default 2+4 mm gaps reproduce the Prototype-3 starting point,\n"
      << "but they are NOT treated as the final NA60+/DiCE geometry.\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--help") {
    PrintUsage(argv[0]);
    return 0;
  }

  Config cfg;
  cfg.anodeVoltageV = ReadArg(argc, argv, "--hv", cfg.anodeVoltageV);
  cfg.magneticFieldT = ReadArg(argc, argv, "--b", cfg.magneticFieldT);
  cfg.gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 10. * cfg.gapMinusCm);
  cfg.gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 10. * cfg.gapPlusCm);
  cfg.wirePitchCm =
      0.1 * ReadArg(argc, argv, "--pitch-mm", 10. * cfg.wirePitchCm);
  cfg.wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um",
                     1.e4 * cfg.wireDiameterCm);
  cfg.x0Cm = ReadArg(argc, argv, "--x0", cfg.x0Cm);
  cfg.y0Cm = ReadArg(argc, argv, "--y0", cfg.y0Cm);

  if (cfg.gapMinusCm <= 0. || cfg.gapPlusCm <= 0. ||
      cfg.wirePitchCm <= 0. || cfg.wireDiameterCm <= 0.) {
    std::cerr << "FAIL: gaps, pitch and wire diameter must be positive.\n";
    return 2;
  }

  // Gas model. The prototype measurements used Ar/CO2 70:30.
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);

  // 2D analytic MWPC cell.
  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);

  for (int i = -cfg.halfNumberOfWires; i <= cfg.halfNumberOfWires; ++i) {
    const double x = i * cfg.wirePitchCm;
    field.AddWire(x, 0., cfg.wireDiameterCm, cfg.anodeVoltageV, "anode");
  }

  const double cathodeMinusY = -cfg.gapMinusCm;
  const double cathodePlusY = +cfg.gapPlusCm;
  field.AddPlaneY(cathodeMinusY, 0., "cathode_minus");
  field.AddPlaneY(cathodePlusY, 0., "cathode_plus");

  // Detector B is vertical (+Y). With the mapping above this is Garfield +x.
  field.SetMagneticField(cfg.magneticFieldT, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);

  const double xExtent =
      (cfg.halfNumberOfWires + 0.5) * cfg.wirePitchCm;
  sensor.SetArea(-xExtent, cathodeMinusY, -0.1,
                  xExtent, cathodePlusY, 0.1);

  Garfield::AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor);
  avalanche.EnableAvalancheSizeLimit(cfg.avalancheLimit);

  std::cout << std::fixed << std::setprecision(5);
  std::cout << "\n========== MWPC PHASE A: SINGLE ELECTRON ==========\n"
            << "gas                  : Ar/CO2 70:30\n"
            << "wire diameter        : " << 1.e4 * cfg.wireDiameterCm
            << " um\n"
            << "wire pitch           : " << 10. * cfg.wirePitchCm
            << " mm\n"
            << "cathode gaps         : "
            << 10. * cfg.gapMinusCm << " + "
            << 10. * cfg.gapPlusCm << " mm\n"
            << "anode voltage        : " << cfg.anodeVoltageV << " V\n"
            << "detector B_y         : " << cfg.magneticFieldT << " T\n"
            << "start (Garf x,y,z)   : (" << cfg.x0Cm << ", "
            << cfg.y0Cm << ", " << cfg.z0Cm << ") cm\n"
            << "mapping              : Garf(x,y,z) = Det(Y,Z,X)\n"
            << "----------------------------------------------------\n";

  avalanche.AvalancheElectron(cfg.x0Cm, cfg.y0Cm, cfg.z0Cm,
                              cfg.t0Ns, cfg.e0Ev);

  int ne = 0;
  int ni = 0;
  avalanche.GetAvalancheSize(ne, ni);

  const std::size_t nEndpoints = avalanche.GetNumberOfElectronEndpoints();
  const double wireRadius = 0.5 * cfg.wireDiameterCm;
  std::size_t endpointsOnWire = 0;

  double meanFinalX = 0.;
  double meanFinalY = 0.;

  for (std::size_t i = 0; i < nEndpoints; ++i) {
    double x0, y0, z0, t0, e0;
    double x1, y1, z1, t1, e1;
    int status = 0;
    avalanche.GetElectronEndpoint(i, x0, y0, z0, t0, e0,
                                  x1, y1, z1, t1, e1, status);

    meanFinalX += x1;
    meanFinalY += y1;

    const int nearest =
        static_cast<int>(std::lround(x1 / cfg.wirePitchCm));
    const double wireX = nearest * cfg.wirePitchCm;
    const double r = std::hypot(x1 - wireX, y1);
    if (r < 1.5 * wireRadius) ++endpointsOnWire;

    if (i < 8) {
      std::cout << "endpoint " << std::setw(4) << i
                << " : (" << x1 << ", " << y1 << ", " << z1
                << ") cm  t=" << t1 << " ns  E=" << e1
                << " eV  status=" << status << "\n";
    }
  }

  if (nEndpoints > 0) {
    meanFinalX /= static_cast<double>(nEndpoints);
    meanFinalY /= static_cast<double>(nEndpoints);
  }

  std::cout << "----------------------------------------------------\n"
            << "avalanche electrons  : " << ne << "\n"
            << "avalanche ions       : " << ni << "\n"
            << "electron endpoints   : " << nEndpoints << "\n"
            << "endpoints on wires   : " << endpointsOnWire << "\n"
            << "mean final (x,y)     : (" << meanFinalX << ", "
            << meanFinalY << ") cm\n"
            << "====================================================\n";

  if (ne == 0 || nEndpoints == 0) {
    std::cerr << "FAIL: no transported avalanche electrons.\n";
    return 2;
  }

  std::cout << "PASS: electron transport/avalanche completed.\n";
  return 0;
}

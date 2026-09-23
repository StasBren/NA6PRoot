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
  double wirePitchCm = 0.4;        // 4 mm
  double wireDiameterCm = 0.003;   // 30 um
  double cathodeMinusCm = -0.2;    // 2 mm from wire plane
  double cathodePlusCm = 0.4;      // 4 mm from wire plane
  double anodeVoltageV = 1800.;    // Prototype-3 starting point
  double magneticFieldT = 0.;      // Phase-A first run: B = 0
  int halfNumberOfWires = 4;       // 9 wires total

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
      << "  --hv V        anode voltage [V] (default 1800)\n"
      << "  --b T         detector B_y [T] (default 0)\n"
      << "  --x0 cm       initial Garfield-x = detector-Y [cm]\n"
      << "  --y0 cm       initial Garfield-y = detector-Z [cm]\n"
      << "\n"
      << "Phase-A geometry: 30 um wires, 4 mm pitch, Ar/CO2 70:30,\n"
      << "2+4 mm cathode gaps. Wires are parallel to detector X.\n";
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
  cfg.x0Cm = ReadArg(argc, argv, "--x0", cfg.x0Cm);
  cfg.y0Cm = ReadArg(argc, argv, "--y0", cfg.y0Cm);

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

  field.AddPlaneY(cfg.cathodeMinusCm, 0., "cathode_minus");
  field.AddPlaneY(cfg.cathodePlusCm, 0., "cathode_plus");

  // Detector B is vertical (+Y). With the mapping above this is Garfield +x.
  field.SetMagneticField(cfg.magneticFieldT, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);

  const double xExtent =
      (cfg.halfNumberOfWires + 0.5) * cfg.wirePitchCm;
  sensor.SetArea(-xExtent, cfg.cathodeMinusCm, -0.1,
                  xExtent, cfg.cathodePlusCm, 0.1);

  Garfield::AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor);
  avalanche.EnableAvalancheSizeLimit(cfg.avalancheLimit);

  std::cout << std::fixed << std::setprecision(5);
  std::cout << "\n========== MWPC PHASE A: SINGLE ELECTRON ==========\n"
            << "gas                  : Ar/CO2 70:30\n"
            << "wire diameter        : " << 10. * cfg.wireDiameterCm
            << " mm\n"
            << "wire pitch           : " << 10. * cfg.wirePitchCm
            << " mm\n"
            << "cathode gaps         : "
            << -10. * cfg.cathodeMinusCm << " + "
            << 10. * cfg.cathodePlusCm << " mm\n"
            << "anode voltage        : " << cfg.anodeVoltageV << " V\n"
            << "detector B_y         : " << cfg.magneticFieldT << " T\n"
            << "start (Garf x,y,z)   : (" << cfg.x0Cm << ", "
            << cfg.y0Cm << ", " << cfg.z0Cm << ") cm\n"
            << "mapping              : Garf(x,y,z) = Det(Y,Z,X)\n"
            << "----------------------------------------------------\n";

  avalanche.AvalancheElectron(cfg.x0Cm, cfg.y0Cm, cfg.z0Cm,
                              cfg.t0Ns, cfg.e0Ev);

  unsigned int ne = 0;
  unsigned int ni = 0;
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

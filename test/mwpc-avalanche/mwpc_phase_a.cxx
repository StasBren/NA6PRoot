#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/ComponentAnalyticField.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Sensor.hh>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Phase-A microscopic MWPC sandbox.
//
// IMPORTANT coordinate mapping:
//   Garfield x = detector/global Y  (across the horizontal anode wires)
//   Garfield y = detector/global Z  (through the gas gap / beam direction locally)
//   Garfield z = detector/global X  (along the horizontal wires)
//
// ComponentAnalyticField treats the wires as parallel to Garfield z, so this
// mapping makes them horizontal in the NA60+/DiCE detector convention.

namespace {

struct Config {
  double pitchCm = 0.4;          // 4 mm
  double wireDiameterCm = 0.003; // 30 um
  double gapMinusCm = 0.2;       // Prototype-3: 2 mm
  double gapPlusCm = 0.4;        // Prototype-3: 4 mm
  double anodeVoltageV = 1800.;
  double magneticFieldT = 0.;
  int halfWireCount = 3;         // total wires = 2*halfWireCount+1
  double startYcm = 0.35;        // start below the + cathode
  double initialElectronEnergyEv = 0.1;
  unsigned int avalancheSizeLimit = 200000;
};

int nearestWireIndex(const double x, const double pitch) {
  return static_cast<int>(std::lround(x / pitch));
}

} // namespace

int main(int argc, char** argv) {
  Config cfg;
  if (argc > 1) cfg.magneticFieldT = std::stod(argv[1]);
  if (argc > 2) cfg.anodeVoltageV = std::stod(argv[2]);

  std::cout << "NA60+/DiCE MWPC Phase-A microscopic test\n"
            << "  gas              = Ar/CO2 70:30\n"
            << "  wire pitch       = " << 10. * cfg.pitchCm << " mm\n"
            << "  wire diameter    = " << 10. * cfg.wireDiameterCm << " mm\n"
            << "  cathode gaps     = " << 10. * cfg.gapMinusCm << " + "
            << 10. * cfg.gapPlusCm << " mm\n"
            << "  anode voltage    = " << cfg.anodeVoltageV << " V\n"
            << "  physical B_y     = " << cfg.magneticFieldT << " T\n\n";

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);

  for (int i = -cfg.halfWireCount; i <= cfg.halfWireCount; ++i) {
    field.AddWire(i * cfg.pitchCm, 0., cfg.wireDiameterCm,
                  cfg.anodeVoltageV, "anode");
  }

  field.AddPlaneY(-cfg.gapMinusCm, 0., "cathodeMinus");
  field.AddPlaneY(+cfg.gapPlusCm, 0., "cathodePlus");

  // Physical MNP33 field is vertical (+global Y). Under the mapping above,
  // this is +Garfield x.
  field.SetMagneticField(cfg.magneticFieldT, 0., 0.);

  Garfield::Sensor sensor(&field);
  const double xExtent = (cfg.halfWireCount + 0.5) * cfg.pitchCm;
  sensor.SetArea(-xExtent, -cfg.gapMinusCm, -10.,
                 +xExtent, +cfg.gapPlusCm, +10.);

  Garfield::AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor);
  avalanche.EnableMagneticField(std::abs(cfg.magneticFieldT) > 0.);
  avalanche.EnableAvalancheSizeLimit(cfg.avalancheSizeLimit);

  std::ofstream csv("phase_a_results.csv");
  if (!csv) {
    throw std::runtime_error("Cannot open phase_a_results.csv");
  }
  csv << "start_x_cm,start_y_cm,B_T,HV_V,ne,ni,n_endpoints,"
         "dominant_wire_index,dominant_wire_fraction\n";

  // Scan one wire pitch. By translational symmetry this is the elementary
  // across-wire problem we ultimately care about.
  const std::vector<double> starts = {
      -0.19, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.19};

  std::cout << std::fixed << std::setprecision(4);
  for (const double x0 : starts) {
    avalanche.AvalancheElectron(x0, cfg.startYcm, 0.,
                                0., cfg.initialElectronEnergyEv);

    int ne = 0, ni = 0;
    avalanche.GetAvalancheSize(ne, ni);

    const std::size_t nEnd = avalanche.GetNumberOfElectronEndpoints();
    std::map<int, std::size_t> endpointWireCounts;

    for (std::size_t i = 0; i < nEnd; ++i) {
      double xe0 = 0., ye0 = 0., ze0 = 0., te0 = 0., ee0 = 0.;
      double xe1 = 0., ye1 = 0., ze1 = 0., te1 = 0., ee1 = 0.;
      int status = 0;
      avalanche.GetElectronEndpoint(i, xe0, ye0, ze0, te0, ee0,
                                    xe1, ye1, ze1, te1, ee1, status);
      endpointWireCounts[nearestWireIndex(xe1, cfg.pitchCm)]++;
    }

    int dominantWire = 999;
    std::size_t dominantCount = 0;
    for (const auto& [wire, count] : endpointWireCounts) {
      if (count > dominantCount) {
        dominantCount = count;
        dominantWire = wire;
      }
    }

    const double dominantFraction =
        nEnd > 0 ? static_cast<double>(dominantCount) / nEnd : 0.;

    std::cout << "start x=" << std::setw(7) << x0 << " cm"
              << "  avalanche e-=" << std::setw(7) << ne
              << "  ions=" << std::setw(7) << ni
              << "  endpoints=" << std::setw(7) << nEnd
              << "  dominant wire=" << std::setw(3) << dominantWire
              << "  fraction=" << dominantFraction << "\n";

    csv << x0 << ',' << cfg.startYcm << ',' << cfg.magneticFieldT << ','
        << cfg.anodeVoltageV << ',' << ne << ',' << ni << ',' << nEnd << ','
        << dominantWire << ',' << dominantFraction << '\n';
  }

  std::cout << "\nWrote phase_a_results.csv\n";
  return 0;
}

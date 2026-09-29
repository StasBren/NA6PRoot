#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Medium.hh"
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

}  // namespace

int main(int argc, char** argv) {
  // Garfield coordinates are the chamber-local coordinates used throughout
  // the Phase-A sandbox:
  //   x = u  : across the anode wires
  //   y = v  : chamber normal
  //   z = w  : along the wires
  //
  // This program visualises one *real microscopic Garfield avalanche*.
  // It stores the drift-line points of many avalanche electrons together with
  // their Garfield times. A separate Python renderer turns the event into a
  // GIF. It also exports the actual ComponentAnalyticField potential/field on
  // a u-v grid, so the GIF background is not a hand-drawn approximation.

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

  const double u0Cm = 0.1 * ReadArg(argc, argv, "--u-mm", 0.0);
  const double v0Cm = 0.1 * ReadArg(argc, argv, "--v-mm", 1.0);
  const double w0Cm = 0.1 * ReadArg(argc, argv, "--w-mm", 0.0);
  const double e0Ev = ReadArg(argc, argv, "--energy-ev", 0.1);

  const int halfNumberOfWires =
      ReadIntArg(argc, argv, "--half-wires", 4);
  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 5000);
  const int maxLines =
      ReadIntArg(argc, argv, "--max-lines", 0);
  const int maxPointsPerLine =
      ReadIntArg(argc, argv, "--max-points-per-line", 120);

  const double fieldHalfWidthCm =
      0.1 * ReadArg(argc, argv, "--field-half-width-mm", 3.0);
  const double fieldVHalfWidthCm =
      0.1 * ReadArg(argc, argv, "--field-v-half-width-mm", 0.0);
  const int fieldUSteps =
      ReadIntArg(argc, argv, "--field-u-steps", 181);
  const int fieldVSteps =
      ReadIntArg(argc, argv, "--field-v-steps", 181);

  const std::string prefix =
      ReadStringArg(argc, argv, "--output-prefix", "avalanche_vis");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      e0Ev <= 0. || halfNumberOfWires < 1 ||
      avalancheLimit <= 0 || maxLines < 0 ||
      maxPointsPerLine < 2 || fieldHalfWidthCm <= 0. ||
      fieldVHalfWidthCm < 0. ||
      fieldUSteps < 3 || fieldVSteps < 3) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (v0Cm <= -gapMinusCm || v0Cm >= gapPlusCm) {
    std::cerr << "Seed v is outside the gas gap.\n";
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

  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    field.AddWire(i * pitchCm, 0., wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double uExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-uExtent, -gapMinusCm, -2.0,
                  uExtent, +gapPlusCm, 2.0);

  // ----------------------------------------------------------------------
  // Export the *actual Garfield analytic field* on a 2-D u-v grid.
  // ElectricField(..., potential, medium, status) returns both E and V.
  // ----------------------------------------------------------------------
  const std::string fieldFile = prefix + "_field.csv";
  std::ofstream fieldOut(fieldFile);
  fieldOut << "u_mm,v_mm,ex_Vcm,ey_Vcm,ez_Vcm,E_Vcm,potential_V,status\n";

  const double fieldVMin =
      fieldVHalfWidthCm > 0.
          ? std::max(-gapMinusCm, -fieldVHalfWidthCm)
          : -gapMinusCm;
  const double fieldVMax =
      fieldVHalfWidthCm > 0.
          ? std::min(+gapPlusCm, +fieldVHalfWidthCm)
          : +gapPlusCm;

  for (int iv = 0; iv < fieldVSteps; ++iv) {
    const double fv =
        static_cast<double>(iv) / static_cast<double>(fieldVSteps - 1);
    const double v =
        fieldVMin + fv * (fieldVMax - fieldVMin);

    for (int iu = 0; iu < fieldUSteps; ++iu) {
      const double fu =
          static_cast<double>(iu) / static_cast<double>(fieldUSteps - 1);
      const double u =
          -fieldHalfWidthCm + 2. * fieldHalfWidthCm * fu;

      double ex = 0., ey = 0., ez = 0., potential = 0.;
      Garfield::Medium* medium = nullptr;
      int status = 0;
      field.ElectricField(u, v, 0., ex, ey, ez,
                          potential, medium, status);
      const double emag = std::sqrt(ex * ex + ey * ey + ez * ez);

      fieldOut << 10. * u << "," << 10. * v << ","
               << ex << "," << ey << "," << ez << ","
               << emag << "," << potential << "," << status << "\n";
    }
  }

  // ----------------------------------------------------------------------
  // One microscopic single-electron avalanche.
  // ----------------------------------------------------------------------
  Garfield::AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor);
  avalanche.EnableAvalancheSizeLimit(
      static_cast<unsigned int>(avalancheLimit));
  avalanche.EnableDriftLines();

  const bool ok =
      avalanche.AvalancheElectron(u0Cm, v0Cm, w0Cm, 0., e0Ev,
                                  0., 0., 0.);

  int ne = 0;
  int ni = 0;
  avalanche.GetAvalancheSize(ne, ni);

  const std::size_t nEndpoints =
      avalanche.GetNumberOfElectronEndpoints();

  // Select drift-line indices evenly across the complete avalanche if the
  // event contains more trajectories than we want to render.
  std::vector<std::size_t> linesToStore;
  const bool storeAllLines = (maxLines == 0);
  const std::size_t nStore =
      storeAllLines
          ? nEndpoints
          : std::min<std::size_t>(
                nEndpoints, static_cast<std::size_t>(maxLines));
  linesToStore.reserve(nStore);

  if (storeAllLines ||
      nEndpoints <= static_cast<std::size_t>(maxLines)) {
    for (std::size_t i = 0; i < nEndpoints; ++i) {
      linesToStore.push_back(i);
    }
  } else if (nStore == 1) {
    linesToStore.push_back(0);
  } else {
    // This mode is useful only for lightweight previews. For a physically
    // continuous presentation GIF prefer --max-lines 0, which stores every
    // avalanche-electron trajectory and avoids hiding parent branches.
    for (std::size_t k = 0; k < nStore; ++k) {
      const double x =
          static_cast<double>(k) *
          static_cast<double>(nEndpoints - 1) /
          static_cast<double>(nStore - 1);
      linesToStore.push_back(
          static_cast<std::size_t>(std::llround(x)));
    }
  }

  const std::string pathFile = prefix + "_paths.csv";
  std::ofstream pathOut(pathFile);
  pathOut << "electron,point,u_mm,v_mm,w_mm,t_ns,status\n";

  const std::string endpointFile = prefix + "_endpoints.csv";
  std::ofstream endpointOut(endpointFile);
  endpointOut << "electron,u0_mm,v0_mm,w0_mm,t0_ns,e0_eV,"
                 "u1_mm,v1_mm,w1_mm,t1_ns,e1_eV,status\n";

  double tMin = std::numeric_limits<double>::max();
  double tMax = -std::numeric_limits<double>::max();

  long long nAttached = 0;
  long long nOther = 0;
  long long nNearWire = 0;

  const double wireRadius = 0.5 * wireDiameterCm;

  for (std::size_t ie = 0; ie < nEndpoints; ++ie) {
    double xa = 0., ya = 0., za = 0., ta = 0., ea = 0.;
    double x1 = 0., y1 = 0., z1 = 0., t1 = 0., e1 = 0.;
    int status = 0;

    avalanche.GetElectronEndpoint(
        ie, xa, ya, za, ta, ea,
        x1, y1, z1, t1, e1, status);

    endpointOut << ie << ","
                << 10. * xa << "," << 10. * ya << "," << 10. * za << ","
                << ta << "," << ea << ","
                << 10. * x1 << "," << 10. * y1 << "," << 10. * z1 << ","
                << t1 << "," << e1 << "," << status << "\n";

    const int nearest =
        static_cast<int>(std::lround(x1 / pitchCm));
    const double wireX = nearest * pitchCm;
    const double r = std::hypot(x1 - wireX, y1);

    if (r < 1.5 * wireRadius) {
      ++nNearWire;
    } else if (status == -7) {
      ++nAttached;
    } else {
      ++nOther;
    }

    tMin = std::min(tMin, ta);
    tMax = std::max(tMax, t1);
  }

  for (const std::size_t ie : linesToStore) {
    const std::size_t nPoints =
        avalanche.GetNumberOfElectronDriftLinePoints(ie);
    if (nPoints == 0) continue;

    const std::size_t stride =
        std::max<std::size_t>(
            1, (nPoints + maxPointsPerLine - 1) /
                   static_cast<std::size_t>(maxPointsPerLine));

    int storedPoint = 0;
    int status = 0;

    double xa = 0., ya = 0., za = 0., ta = 0., ea = 0.;
    double x1 = 0., y1 = 0., z1 = 0., t1 = 0., e1 = 0.;
    avalanche.GetElectronEndpoint(
        ie, xa, ya, za, ta, ea,
        x1, y1, z1, t1, e1, status);

    for (std::size_t ip = 0; ip < nPoints; ip += stride) {
      double x = 0., y = 0., z = 0., t = 0.;
      avalanche.GetElectronDriftLinePoint(
          x, y, z, t, ip, ie);
      pathOut << ie << "," << storedPoint++ << ","
              << 10. * x << "," << 10. * y << "," << 10. * z << ","
              << t << "," << status << "\n";
    }

    // Always include the last point.
    double xl = 0., yl = 0., zl = 0., tl = 0.;
    avalanche.GetElectronDriftLinePoint(
        xl, yl, zl, tl, nPoints - 1, ie);
    pathOut << ie << "," << storedPoint << ","
            << 10. * xl << "," << 10. * yl << "," << 10. * zl << ","
            << tl << "," << status << "\n";
  }

  const bool sizeLimitReached =
      ne >= avalancheLimit ||
      nEndpoints >= static_cast<std::size_t>(avalancheLimit);

  const std::string metaFile = prefix + "_meta.txt";
  std::ofstream metaOut(metaFile);
  metaOut << std::setprecision(12)
          << "gas=Ar/CO2 70/30\n"
          << "wire_pitch_mm=" << 10. * pitchCm << "\n"
          << "wire_diameter_um=" << 1.e4 * wireDiameterCm << "\n"
          << "gap_minus_mm=" << 10. * gapMinusCm << "\n"
          << "gap_plus_mm=" << 10. * gapPlusCm << "\n"
          << "anode_voltage_V=" << hv << "\n"
          << "B_T=" << bTesla << "\n"
          << "seed_u_mm=" << 10. * u0Cm << "\n"
          << "seed_v_mm=" << 10. * v0Cm << "\n"
          << "seed_w_mm=" << 10. * w0Cm << "\n"
          << "seed_energy_eV=" << e0Ev << "\n"
          << "field_half_width_mm=" << 10. * fieldHalfWidthCm << "\n"
          << "field_v_min_mm=" << 10. * fieldVMin << "\n"
          << "field_v_max_mm=" << 10. * fieldVMax << "\n"
          << "avalanche_ok=" << (ok ? 1 : 0) << "\n"
          << "avalanche_electrons=" << ne << "\n"
          << "avalanche_ions=" << ni << "\n"
          << "avalanche_size_limit=" << avalancheLimit << "\n"
          << "size_limit_reached=" << (sizeLimitReached ? 1 : 0) << "\n"
          << "electron_endpoints=" << nEndpoints << "\n"
          << "stored_drift_lines=" << linesToStore.size() << "\n"
          << "near_wire_endpoints=" << nNearWire << "\n"
          << "attached_endpoints=" << nAttached << "\n"
          << "other_endpoints=" << nOther << "\n"
          << "t_min_ns=" << tMin << "\n"
          << "t_max_ns=" << tMax << "\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== MICROSCOPIC AVALANCHE VISUALISATION ===\n"
            << "seed (u,v,w)       : (" << 10. * u0Cm << ", "
            << 10. * v0Cm << ", " << 10. * w0Cm << ") mm\n"
            << "seed energy        : " << e0Ev << " eV\n"
            << "wire pitch        : " << 10. * pitchCm << " mm\n"
            << "wire diameter     : " << 1.e4 * wireDiameterCm << " um\n"
            << "gaps              : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "anode voltage     : " << hv << " V\n"
            << "B                 : " << bTesla << " T\n"
            << "avalanche e-      : " << ne << "\n"
            << "avalanche ions    : " << ni << "\n"
            << "size limit        : " << avalancheLimit << "\n"
            << "limit reached     : "
            << (sizeLimitReached ? "YES" : "no") << "\n"
            << "electron paths    : " << nEndpoints << "\n"
            << "stored paths      : " << linesToStore.size() << "\n"
            << "near-wire ends    : " << nNearWire << "\n"
            << "attached ends     : " << nAttached << "\n"
            << "other ends        : " << nOther << "\n"
            << "time span         : " << (tMax - tMin) << " ns\n\n"
            << (sizeLimitReached
                    ? "WARNING: avalanche hit the configured size limit; "
                      "rerun before using this event as a complete avalanche.\n\n"
                    : "")
            << "Wrote:\n"
            << "  " << metaFile << "\n"
            << "  " << fieldFile << "\n"
            << "  " << endpointFile << "\n"
            << "  " << pathFile << "\n";

  return 0;
}

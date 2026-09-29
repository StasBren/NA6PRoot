#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

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
  // Garfield coordinates are identified with chamber-local coordinates:
  //   x = u  : across the anode wires
  //   y = v  : chamber normal / drift direction
  //   z = w  : along the anode wires
  //
  // This executable is meant for VISUALISATION, not detector production.
  // It generates ONE Heed muon event and stores:
  //   - Heed ionisation-cluster positions,
  //   - individual low-energy conduction-electron seed positions,
  //   - the actual microscopic Garfield drift-line points for a configurable
  //     subset of those electrons.
  //
  // Avalanche multiplication is deliberately not included here. The goal is
  // to make the Heed -> microscopic transport handoff visible first.

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
  const double u0Cm = 0.1 * ReadArg(argc, argv, "--u0-mm", 1.0);
  const double momentumEv =
      1.e9 * ReadArg(argc, argv, "--momentum-gev", 10.0);

  const int maxVisualElectrons =
      ReadIntArg(argc, argv, "--max-electrons", 40);
  const int maxPathPoints =
      ReadIntArg(argc, argv, "--max-path-points", 250);
  const std::string prefix =
      ReadStringArg(argc, argv, "--output-prefix", "heed_vis");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      momentumEv <= 0. || maxVisualElectrons <= 0 ||
      maxPathPoints < 2) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (u0Cm < -0.5 * pitchCm || u0Cm > 0.5 * pitchCm) {
    std::cerr << "Starting u is outside the central wire cell [-p/2,+p/2].\n";
    return 2;
  }

  // ------------------------------------------------------------------------
  // Gas and ideal MWPC field: same Phase-A reference setup used elsewhere.
  // ------------------------------------------------------------------------
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

  // Nominal NA60 mapping used in the Phase-A sandbox:
  // global vertical B -> local +u -> Garfield +x.
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double uExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-uExtent, -gapMinusCm, -5.0,
                  uExtent, +gapPlusCm, 5.0);

  // ------------------------------------------------------------------------
  // Heed: one muon crossing from +v to -v.
  // ------------------------------------------------------------------------
  Garfield::TrackHeed heed;
  heed.SetSensor(&sensor);
  heed.SetParticle("mu-");
  heed.SetMomentum(momentumEv);
  heed.EnableDeltaElectronTransport();

  const double marginCm = std::min(0.01, 0.05 * gapPlusCm);
  const double vStartCm = gapPlusCm - marginCm;
  heed.NewTrack(u0Cm, vStartCm, 0., 0.,
                0., -1., 0.);

  const std::string metaName = prefix + "_meta.txt";
  const std::string clusterName = prefix + "_clusters.csv";
  const std::string electronName = prefix + "_electrons.csv";
  const std::string pathName = prefix + "_paths.csv";

  std::ofstream metaOut(metaName);
  std::ofstream clusterOut(clusterName);
  std::ofstream electronOut(electronName);
  std::ofstream pathOut(pathName);

  if (!metaOut || !clusterOut || !electronOut || !pathOut) {
    std::cerr << "Could not open one or more output files.\n";
    return 3;
  }

  metaOut << std::setprecision(12)
          << "gas=Ar/CO2 70/30\n"
          << "particle=mu-\n"
          << "momentum_GeV=" << momentumEv * 1.e-9 << "\n"
          << "wire_pitch_mm=" << 10. * pitchCm << "\n"
          << "wire_diameter_um=" << 1.e4 * wireDiameterCm << "\n"
          << "gap_minus_mm=" << 10. * gapMinusCm << "\n"
          << "gap_plus_mm=" << 10. * gapPlusCm << "\n"
          << "anode_voltage_V=" << hv << "\n"
          << "B_T=" << bTesla << "\n"
          << "muon_u_mm=" << 10. * u0Cm << "\n"
          << "muon_v_start_mm=" << 10. * vStartCm << "\n"
          << "muon_v_end_mm=" << -10. * gapMinusCm << "\n"
          << "muon_w_mm=0\n";

  clusterOut << "cluster,u_mm,v_mm,w_mm,t_ns,electrons,energy_transfer_eV\n";
  electronOut << "electron,cluster,u_mm,v_mm,w_mm,t_ns,heed_energy_eV,"
                 "heed_dx,heed_dy,heed_dz,drifted\n";
  pathOut << "electron,cluster,point,u_mm,v_mm,w_mm,t_ns,status\n";

  int nClusters = 0;
  int nConduction = 0;
  int nVisual = 0;
  double totalEnergyEv = 0.;

  double xc = 0., yc = 0., zc = 0., tc = 0.;
  double ec = 0., extra = 0.;
  int nc = 0;

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== HEED / GARFIELD MUON VISUALISATION EVENT ===\n"
            << "particle            : mu-\n"
            << "momentum            : " << momentumEv * 1.e-9 << " GeV/c\n"
            << "gas                 : Ar/CO2 70/30\n"
            << "wire pitch          : " << 10. * pitchCm << " mm\n"
            << "gas gap             : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "anode voltage       : " << hv << " V\n"
            << "B                   : " << bTesla << " T\n"
            << "muon start u        : " << 10. * u0Cm << " mm\n"
            << "max drifted e-      : " << maxVisualElectrons << "\n"
            << "max stored path pts : " << maxPathPoints << " / electron\n\n";

  while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    const int clusterId = nClusters++;
    totalEnergyEv += ec;

    clusterOut << clusterId << ","
               << 10. * xc << "," << 10. * yc << "," << 10. * zc << ","
               << tc << "," << nc << "," << ec << "\n";

    for (int ie = 0; ie < nc; ++ie) {
      const int electronId = nConduction++;

      double xe = 0., ye = 0., ze = 0., te = 0.;
      double ee = 0., dxe = 0., dye = 0., dze = 0.;
      if (!heed.GetElectron(ie, xe, ye, ze, te,
                            ee, dxe, dye, dze)) {
        continue;
      }

      const bool doDrift = nVisual < maxVisualElectrons;
      electronOut << electronId << "," << clusterId << ","
                  << 10. * xe << "," << 10. * ye << "," << 10. * ze << ","
                  << te << "," << ee << ","
                  << dxe << "," << dye << "," << dze << ","
                  << (doDrift ? 1 : 0) << "\n";

      if (!doDrift) continue;
      ++nVisual;

      // For low-energy conduction electrons returned after Heed delta-electron
      // transport, position/time are the meaningful handoff. We start the
      // microscopic Garfield drift at 0.1 eV with zero direction; Garfield
      // then randomises the initial direction and follows individual gas
      // collisions.
      Garfield::AvalancheMicroscopic drift;
      drift.SetSensor(&sensor);
      drift.EnableDriftLines();
      const bool ok = drift.DriftElectron(xe, ye, ze, te, 0.1, 0., 0., 0.);

      if (!ok || drift.GetNumberOfElectronEndpoints() < 1) continue;

      double xa = 0., ya = 0., za = 0., ta = 0., ea = 0.;
      double x1 = 0., y1 = 0., z1 = 0., t1 = 0., e1 = 0.;
      int status = 0;
      drift.GetElectronEndpoint(0, xa, ya, za, ta, ea,
                                x1, y1, z1, t1, e1, status);

      const std::size_t nPoints =
          drift.GetNumberOfElectronDriftLinePoints(0);
      if (nPoints == 0) continue;

      const std::size_t stride =
          std::max<std::size_t>(1, (nPoints + maxPathPoints - 1) /
                                      maxPathPoints);

      int storedPoint = 0;
      for (std::size_t ip = 0; ip < nPoints; ip += stride) {
        double x = 0., y = 0., z = 0., t = 0.;
        drift.GetElectronDriftLinePoint(x, y, z, t, ip, 0);
        pathOut << electronId << "," << clusterId << "," << storedPoint++ << ","
                << 10. * x << "," << 10. * y << "," << 10. * z << ","
                << t << "," << status << "\n";
      }

      // Make sure the physical endpoint is represented even if the stride did
      // not land exactly on the final stored point.
      double xl = 0., yl = 0., zl = 0., tl = 0.;
      drift.GetElectronDriftLinePoint(
          xl, yl, zl, tl, nPoints - 1, 0);
      pathOut << electronId << "," << clusterId << "," << storedPoint << ","
              << 10. * xl << "," << 10. * yl << "," << 10. * zl << ","
              << tl << "," << status << "\n";
    }
  }

  metaOut << "n_clusters=" << nClusters << "\n"
          << "n_conduction_electrons=" << nConduction << "\n"
          << "n_visualised_electrons=" << nVisual << "\n"
          << "total_energy_transfer_eV=" << totalEnergyEv << "\n";

  std::cout << "clusters             : " << nClusters << "\n"
            << "conduction electrons : " << nConduction << "\n"
            << "drifted for GIF      : " << nVisual << "\n"
            << "total dE             : " << totalEnergyEv << " eV\n\n"
            << "Wrote:\n"
            << "  " << metaName << "\n"
            << "  " << clusterName << "\n"
            << "  " << electronName << "\n"
            << "  " << pathName << "\n";

  return 0;
}

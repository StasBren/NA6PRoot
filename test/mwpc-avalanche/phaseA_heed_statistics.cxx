#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

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

std::string ReadStringArg(const int argc, char** argv,
                          const std::string& key,
                          const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return argv[i + 1];
  }
  return defaultValue;
}

}  // namespace

int main(int argc, char** argv) {
  // Fast Heed-only ensemble for primary-ionisation statistics.
  //
  // This intentionally stops BEFORE Garfield electron drift/avalanche.
  // For every muon track we store:
  //   - number of Heed ionisation clusters,
  //   - total number of conduction electrons,
  //   - total energy transferred in the gas.
  //
  // For every cluster we store:
  //   - number of conduction electrons in that cluster,
  //   - energy transferred to that cluster,
  //   - cluster position/time.
  //
  // This lets us study separately:
  //   1) total ionisation electrons per muon;
  //   2) cluster-size fluctuations;
  //   3) dE_track vs N_e(track);
  //   4) dE_cluster vs N_e(cluster).

  const double pitchCm =
      0.1 * ReadArg(argc, argv, "--pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);
  const double u0Cm = 0.1 * ReadArg(argc, argv, "--u0-mm", 1.0);
  const double momentumEv =
      1.e9 * ReadArg(argc, argv, "--momentum-gev", 10.0);
  const double thetaDeg = ReadArg(argc, argv, "--theta-deg", 0.0);
  const double phiDeg = ReadArg(argc, argv, "--phi-deg", 0.0);
  const int tracks = ReadIntArg(argc, argv, "--tracks", 5000);
  const std::string prefix =
      ReadStringArg(argc, argv, "--output-prefix", "heed_stats");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      momentumEv <= 0. || tracks <= 0 ||
      std::abs(thetaDeg) >= 80.) {
    std::cerr << "Invalid input parameters.\n";
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
  for (int iw = -halfNumberOfWires; iw <= halfNumberOfWires; ++iw) {
    field.AddWire(iw * pitchCm, 0., wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double uExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-uExtent, -gapMinusCm, -5.0,
                  uExtent, +gapPlusCm, 5.0);

  Garfield::TrackHeed heed;
  heed.SetSensor(&sensor);
  heed.SetParticle("mu-");
  heed.SetMomentum(momentumEv);
  heed.EnableDeltaElectronTransport();

  const double pi = std::acos(-1.0);
  const double theta = thetaDeg * pi / 180.;
  const double phi = phiDeg * pi / 180.;
  const double dirU = std::sin(theta) * std::cos(phi);
  const double dirV = -std::cos(theta);
  const double dirW = std::sin(theta) * std::sin(phi);

  const double marginCm = std::min(0.01, 0.05 * gapPlusCm);
  const double vStartCm = gapPlusCm - marginCm;

  const std::string trackName = prefix + "_tracks.csv";
  const std::string clusterName = prefix + "_clusters.csv";

  std::ofstream trackOut(trackName);
  std::ofstream clusterOut(clusterName);
  if (!trackOut || !clusterOut) {
    std::cerr << "Could not open output files.\n";
    return 3;
  }

  trackOut << "track,clusters,electrons,total_dE_eV,"
              "mean_cluster_electrons,max_cluster_electrons,"
              "max_cluster_dE_eV\n";
  clusterOut << "track,cluster,electrons,energy_transfer_eV,"
                "u_mm,v_mm,w_mm,t_ns\n";

  long long totalClusters = 0;
  long long totalElectrons = 0;
  double totalEnergyEv = 0.;

  for (int itrack = 0; itrack < tracks; ++itrack) {
    heed.NewTrack(
        u0Cm, vStartCm, 0., 0.,
        dirU, dirV, dirW);

    int nClusters = 0;
    int nElectrons = 0;
    int maxClusterElectrons = 0;
    double dETrack = 0.;
    double maxClusterDE = 0.;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    double ec = 0., extra = 0.;
    int nc = 0;

    while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      clusterOut
          << itrack << "," << nClusters << ","
          << nc << "," << ec << ","
          << 10. * xc << "," << 10. * yc << "," << 10. * zc << ","
          << tc << "\n";

      ++nClusters;
      nElectrons += nc;
      dETrack += ec;
      maxClusterElectrons = std::max(maxClusterElectrons, nc);
      maxClusterDE = std::max(maxClusterDE, ec);
    }

    const double meanClusterElectrons =
        nClusters > 0
            ? static_cast<double>(nElectrons) /
                  static_cast<double>(nClusters)
            : 0.;

    trackOut
        << itrack << ","
        << nClusters << ","
        << nElectrons << ","
        << dETrack << ","
        << meanClusterElectrons << ","
        << maxClusterElectrons << ","
        << maxClusterDE << "\n";

    totalClusters += nClusters;
    totalElectrons += nElectrons;
    totalEnergyEv += dETrack;

    if ((itrack + 1) % 500 == 0 || itrack + 1 == tracks) {
      std::cout << "completed " << (itrack + 1)
                << " / " << tracks << " tracks\n";
    }
  }

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "\n=== HEED PRIMARY-IONISATION STATISTICS ===\n"
            << "tracks                    : " << tracks << "\n"
            << "particle                  : mu-\n"
            << "momentum                  : "
            << momentumEv * 1.e-9 << " GeV/c\n"
            << "gas                       : Ar/CO2 70/30\n"
            << "gas gap                   : "
            << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "theta / phi               : "
            << thetaDeg << " / " << phiDeg << " deg\n"
            << "mean clusters / track     : "
            << static_cast<double>(totalClusters) / tracks << "\n"
            << "mean electrons / track    : "
            << static_cast<double>(totalElectrons) / tracks << "\n"
            << "mean dE / track           : "
            << totalEnergyEv / tracks << " eV\n"
            << "global dE / electron      : "
            << (totalElectrons > 0
                    ? totalEnergyEv / static_cast<double>(totalElectrons)
                    : 0.)
            << " eV/e-\n"
            << "\nWrote:\n"
            << "  " << trackName << "\n"
            << "  " << clusterName << "\n";

  return 0;
}

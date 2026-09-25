#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
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
  const double x0Cm = 0.1 * ReadArg(argc, argv, "--x0-mm", 1.0);
  const double momentumEv =
      1.e9 * ReadArg(argc, argv, "--momentum-gev", 10.0);
  const int tracks = ReadIntArg(argc, argv, "--tracks", 3);
  const int maxPrintedClusters =
      ReadIntArg(argc, argv, "--print-clusters", 12);

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      tracks <= 0 || momentumEv <= 0.) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (x0Cm < -0.5 * pitchCm || x0Cm > 0.5 * pitchCm) {
    std::cerr << "Starting u is outside the central wire cell [-p/2,+p/2].\n";
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

  // Nominal mapping:
  // detector/global +Y magnetic field -> chamber-local +u -> Garfield +x.
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double xExtent = (halfNumberOfWires + 0.5) * pitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -5.0,
                  xExtent, +gapPlusCm, 5.0);

  Garfield::TrackHeed heed;
  heed.SetSensor(&sensor);
  heed.SetParticle("mu-");
  heed.SetMomentum(momentumEv);
  // Keep the default Heed delta-electron transport explicit:
  // GetElectron then returns low-energy conduction electrons from the
  // ionisation cascade, which are the appropriate seeds for drift.
  heed.EnableDeltaElectronTransport();

  const double wireRadius = 0.5 * wireDiameterCm;

  std::ofstream trackOut("heed_track_summary.csv");
  trackOut << "track,clusters,conduction_electrons,energy_loss_eV,"
              "wire_m1,wire_0,wire_p1,attachment,other,active_wires,"
              "mean_drift_time_ns\n";

  std::ofstream clusterOut("heed_clusters.csv");
  clusterOut << "track,cluster,u_mm,v_mm,w_mm,t_ns,"
                "electrons,energy_transfer_eV\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== PHASE A HEED MUON TRACK INSPECTION ===\n"
            << "tracks             : " << tracks << "\n"
            << "particle           : mu-\n"
            << "momentum           : " << momentumEv * 1.e-9 << " GeV/c\n"
            << "start u            : " << 10. * x0Cm << " mm\n"
            << "track direction    : local -v (through the gas gap)\n"
            << "gas gap            : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "wire pitch         : " << 10. * pitchCm << " mm\n"
            << "B_y                : " << bTesla << " T\n"
            << "Heed output        : ionisation clusters + conduction electrons\n"
            << "Garfield step      : drift each conduction electron to a wire"
               " (no avalanche yet)\n\n";

  for (int itrack = 0; itrack < tracks; ++itrack) {
    // Start just inside the +v cathode and propagate toward -v.
    // Garfield coordinates are (x,y,z) = (u,v,w).
    const double marginCm = std::min(0.01, 0.05 * gapPlusCm);
    const double yStartCm = gapPlusCm - marginCm;

    heed.NewTrack(x0Cm, yStartCm, 0., 0.,
                  0., -1., 0.);

    int nClusters = 0;
    int nElectrons = 0;
    double energyLossEv = 0.;

    int nM1 = 0, n0 = 0, nP1 = 0;
    int nAttach = 0, nOther = 0;
    double sumDriftTime = 0.;
    int nTimed = 0;
    std::set<int> activeWires;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    double ec = 0., extra = 0.;
    int nc = 0;

    while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      clusterOut << itrack << "," << nClusters << ","
                 << 10. * xc << "," << 10. * yc << "," << 10. * zc << ","
                 << tc << "," << nc << "," << ec << "\n";

      if (itrack < 3 && nClusters < maxPrintedClusters) {
        std::cout << "track " << itrack
                  << " cluster " << std::setw(3) << nClusters
                  << " : (u,v,w)=(" << std::setw(7) << 10. * xc << ", "
                  << std::setw(7) << 10. * yc << ", "
                  << std::setw(7) << 10. * zc << ") mm"
                  << "  ne=" << std::setw(3) << nc
                  << "  dE=" << std::setw(9) << ec << " eV\n";
      }

      ++nClusters;
      nElectrons += nc;
      energyLossEv += ec;

      for (int ie = 0; ie < nc; ++ie) {
        double xe = 0., ye = 0., ze = 0., te = 0.;
        double ee = 0., dxe = 0., dye = 0., dze = 0.;
        if (!heed.GetElectron(ie, xe, ye, ze, te,
                              ee, dxe, dye, dze)) {
          ++nOther;
          continue;
        }

        // With Heed delta-electron transport enabled, GetElectron gives
        // low-energy conduction-electron positions. The returned energy/
        // direction are not physically meaningful for these electrons, so
        // start microscopic drift at a small thermal-scale energy and let
        // Garfield randomise the initial direction.
        Garfield::AvalancheMicroscopic drift;
        drift.SetSensor(&sensor);
        drift.DriftElectron(xe, ye, ze, te, 0.1);

        if (drift.GetNumberOfElectronEndpoints() < 1) {
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
          activeWires.insert(nearest);
          if (nearest == -1) ++nM1;
          else if (nearest == 0) ++n0;
          else if (nearest == +1) ++nP1;
          else ++nOther;

          if (t1 >= te) {
            sumDriftTime += t1 - te;
            ++nTimed;
          }
        } else if (status == -7) {
          ++nAttach;
        } else {
          ++nOther;
        }
      }
    }

    const double meanDriftTime =
        nTimed > 0 ? sumDriftTime / nTimed : 0.;

    trackOut << itrack << "," << nClusters << ","
             << nElectrons << "," << energyLossEv << ","
             << nM1 << "," << n0 << "," << nP1 << ","
             << nAttach << "," << nOther << ","
             << activeWires.size() << "," << meanDriftTime << "\n";

    std::cout << "\ntrack " << itrack << " summary:"
              << " clusters=" << nClusters
              << ", conduction e-=" << nElectrons
              << ", dE=" << energyLossEv << " eV"
              << ", wire[-1,0,+1]=[" << nM1 << "," << n0 << "," << nP1 << "]"
              << ", attached=" << nAttach
              << ", other=" << nOther
              << ", active wires=" << activeWires.size()
              << ", <drift t>=" << meanDriftTime << " ns\n\n";
  }

  std::cout << "Wrote heed_track_summary.csv and heed_clusters.csv\n";
  return 0;
}

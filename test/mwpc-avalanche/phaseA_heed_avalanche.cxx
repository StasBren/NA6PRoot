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

std::string ReadStringArg(const int argc, char** argv, const std::string& key,
                          const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return argv[i + 1];
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
  const int maxPrintedSeeds =
      ReadIntArg(argc, argv, "--print-seeds", 8);
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix", "heed_avalanche");

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
  heed.EnableDeltaElectronTransport();

  const double wireRadius = 0.5 * wireDiameterCm;
  constexpr unsigned int avalancheLimit = 200000;

  const std::string trackFile = outputPrefix + "_track_summary.csv";
  const std::string seedFile = outputPrefix + "_seed_summary.csv";

  std::ofstream trackOut(trackFile);
  trackOut << "track,clusters,primary_electrons,energy_loss_eV,"
              "seed_to_m1,seed_to_0,seed_to_p1,seed_multi_wire,seed_zero,"
              "avalanche_electrons,wire_m1_charge_e,wire_0_charge_e,"
              "wire_p1_charge_e,attached_endpoints,other_endpoints,"
              "active_wires,total_collected_e,total_collected_fC,"
              "mean_collected_time_ns\n";

  std::ofstream seedOut(seedFile);
  seedOut << "track,seed,seed_u_mm,seed_v_mm,seed_w_mm,"
             "avalanche_electrons,avalanche_ions,collected_e,"
             "wire_m1_charge_e,wire_0_charge_e,wire_p1_charge_e,"
             "attached_endpoints,other_endpoints,mean_collected_time_ns\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== PHASE A HEED + FULL AVALANCHE TRACK TEST ===\n"
            << "tracks             : " << tracks << "\n"
            << "particle           : mu-\n"
            << "momentum           : " << momentumEv * 1.e-9 << " GeV/c\n"
            << "start u            : " << 10. * x0Cm << " mm\n"
            << "track direction    : local -v\n"
            << "gas gap            : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "wire pitch         : " << 10. * pitchCm << " mm\n"
            << "B_y                : " << bTesla << " T\n"
            << "chain              : Heed ionisation -> conduction e-"
               " -> microscopic avalanche -> wire charge\n\n";

  for (int itrack = 0; itrack < tracks; ++itrack) {
    const double marginCm = std::min(0.01, 0.05 * gapPlusCm);
    const double yStartCm = gapPlusCm - marginCm;

    heed.NewTrack(x0Cm, yStartCm, 0., 0.,
                  0., -1., 0.);

    int nClusters = 0;
    int nPrimaryElectrons = 0;
    double energyLossEv = 0.;

    long long totalAvalancheElectrons = 0;
    long long qM1 = 0, q0 = 0, qP1 = 0;
    int seedToM1 = 0, seedTo0 = 0, seedToP1 = 0;
    int seedMultiWire = 0, seedZero = 0;
    long long nAttachedEndpoints = 0;
    long long nOtherEndpoints = 0;
    long long totalCollected = 0;
    double sumCollectedTime = 0.;
    long long nTimedCollected = 0;
    std::set<int> activeWires;

    int seedIndex = 0;

    double xc = 0., yc = 0., zc = 0., tc = 0.;
    double ec = 0., extra = 0.;
    int nc = 0;

    while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      ++nClusters;
      nPrimaryElectrons += nc;
      energyLossEv += ec;

      for (int ie = 0; ie < nc; ++ie) {
        double xe = 0., ye = 0., ze = 0., te = 0.;
        double ee = 0., dxe = 0., dye = 0., dze = 0.;
        if (!heed.GetElectron(ie, xe, ye, ze, te,
                              ee, dxe, dye, dze)) {
          ++nOtherEndpoints;
          ++seedIndex;
          continue;
        }

        Garfield::AvalancheMicroscopic avalanche;
        avalanche.SetSensor(&sensor);
        avalanche.EnableAvalancheSizeLimit(avalancheLimit);
        avalanche.AvalancheElectron(xe, ye, ze, te, 0.1);

        int ne = 0;
        int ni = 0;
        avalanche.GetAvalancheSize(ne, ni);
        totalAvalancheElectrons += ne;

        int seedCollected = 0;
        int seedAttached = 0;
        int seedOther = 0;
        long long seedQM1 = 0, seedQ0 = 0, seedQP1 = 0;
        double seedTimeSum = 0.;
        int seedTimed = 0;

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
            ++seedCollected;
            ++totalCollected;
            activeWires.insert(nearest);

            if (nearest == -1) {
              ++qM1;
              ++seedQM1;
            } else if (nearest == 0) {
              ++q0;
              ++seedQ0;
            } else if (nearest == +1) {
              ++qP1;
              ++seedQP1;
            } else {
              ++nOtherEndpoints;
            }

            if (t1 >= te) {
              const double dt = t1 - te;
              seedTimeSum += dt;
              ++seedTimed;
              sumCollectedTime += dt;
              ++nTimedCollected;
            }
          } else if (status == -7) {
            ++seedAttached;
            ++nAttachedEndpoints;
          } else {
            ++seedOther;
            ++nOtherEndpoints;
          }
        }

        const int seedWireCount =
            (seedQM1 > 0 ? 1 : 0) +
            (seedQ0 > 0 ? 1 : 0) +
            (seedQP1 > 0 ? 1 : 0);
        if (seedWireCount == 0) {
          ++seedZero;
        } else if (seedWireCount > 1) {
          ++seedMultiWire;
        } else if (seedQM1 > 0) {
          ++seedToM1;
        } else if (seedQ0 > 0) {
          ++seedTo0;
        } else if (seedQP1 > 0) {
          ++seedToP1;
        }

        const double seedMeanTime =
            seedTimed > 0 ? seedTimeSum / seedTimed : 0.;

        seedOut << itrack << "," << seedIndex << ","
                << 10. * xe << "," << 10. * ye << "," << 10. * ze << ","
                << ne << "," << ni << ","
                << seedCollected << ","
                << seedQM1 << "," << seedQ0 << "," << seedQP1 << ","
                << seedAttached << "," << seedOther << ","
                << seedMeanTime << "\n";

        if (itrack < 3 && seedIndex < maxPrintedSeeds) {
          std::cout << "track " << itrack
                    << " seed " << std::setw(3) << seedIndex
                    << " : (u,v,w)=(" << std::setw(7) << 10. * xe << ", "
                    << std::setw(7) << 10. * ye << ", "
                    << std::setw(7) << 10. * ze << ") mm"
                    << "  avalanche_e=" << std::setw(6) << ne
                    << "  collected=" << std::setw(6) << seedCollected
                    << "  attached=" << std::setw(4) << seedAttached
                    << "\n";
        }

        ++seedIndex;
      }
    }

    const double meanCollectedTime =
        nTimedCollected > 0
            ? sumCollectedTime / static_cast<double>(nTimedCollected)
            : 0.;
    const double totalCollectedFc =
        static_cast<double>(totalCollected) * 1.602176634e-4;

    trackOut << itrack << "," << nClusters << ","
             << nPrimaryElectrons << "," << energyLossEv << ","
             << seedToM1 << "," << seedTo0 << "," << seedToP1 << ","
             << seedMultiWire << "," << seedZero << ","
             << totalAvalancheElectrons << ","
             << qM1 << "," << q0 << "," << qP1 << ","
             << nAttachedEndpoints << "," << nOtherEndpoints << ","
             << activeWires.size() << "," << totalCollected << ","
             << totalCollectedFc << "," << meanCollectedTime << "\n";

    std::cout << "\ntrack " << itrack << " summary:"
              << " clusters=" << nClusters
              << ", primary e-=" << nPrimaryElectrons
              << ", dE=" << energyLossEv << " eV"
              << ", seed wire[-1,0,+1]=[" << seedToM1 << ","
              << seedTo0 << "," << seedToP1 << "]"
              << ", seed multi=" << seedMultiWire
              << ", seed zero=" << seedZero
              << ", avalanche e-=" << totalAvalancheElectrons
              << ", wire charge[-1,0,+1]=[" << qM1 << ","
              << q0 << "," << qP1 << "] e-"
              << ", total collected=" << totalCollected
              << " e- (" << totalCollectedFc << " fC)"
              << ", attached endpoints=" << nAttachedEndpoints
              << ", active wires=" << activeWires.size()
              << ", <arrival dt>=" << meanCollectedTime << " ns\n\n";
  }

  std::cout << "Wrote " << trackFile << " and " << seedFile << "\n";
  return 0;
}

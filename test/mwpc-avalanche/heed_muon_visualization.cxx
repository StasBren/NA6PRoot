#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

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

struct ElectronSeed {
  int electronId = -1;
  int clusterId = -1;
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;
  double heedEnergy = 0.;
  double dx = 0.;
  double dy = 0.;
  double dz = 0.;
};

std::vector<int> SelectFirst(const int n, const int maxShow) {
  std::vector<int> out;
  const int nTake = std::min(n, maxShow);
  out.reserve(nTake);
  for (int i = 0; i < nTake; ++i) out.push_back(i);
  return out;
}

std::vector<int> SelectRandom(const int n, const int maxShow,
                              const unsigned int seed) {
  std::vector<int> out(n);
  for (int i = 0; i < n; ++i) out[i] = i;
  std::mt19937 rng(seed);
  std::shuffle(out.begin(), out.end(), rng);
  if (static_cast<int>(out.size()) > maxShow) out.resize(maxShow);
  return out;
}

std::vector<int> SelectBalanced(
    const std::vector<ElectronSeed>& seeds,
    const int maxShow,
    const int maxPerCluster,
    const unsigned int rngSeed) {
  std::unordered_map<int, std::vector<int>> byCluster;
  std::vector<int> clusterOrder;

  for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
    const int cid = seeds[i].clusterId;
    if (byCluster.find(cid) == byCluster.end()) {
      clusterOrder.push_back(cid);
    }
    byCluster[cid].push_back(i);
  }

  std::mt19937 rng(rngSeed);
  for (auto& kv : byCluster) {
    std::shuffle(kv.second.begin(), kv.second.end(), rng);
  }

  std::vector<int> selected;
  selected.reserve(std::min<int>(maxShow, seeds.size()));

  // Round-robin over clusters. This avoids one large ionisation cluster
  // monopolising the presentation sample.
  for (int round = 0; round < maxPerCluster; ++round) {
    bool added = false;
    for (const int cid : clusterOrder) {
      auto& indices = byCluster[cid];
      if (round >= static_cast<int>(indices.size())) continue;
      selected.push_back(indices[round]);
      added = true;
      if (static_cast<int>(selected.size()) >= maxShow) return selected;
    }
    if (!added) break;
  }

  // If the per-cluster cap leaves unused display slots, fill them from the
  // remaining electrons without biasing toward early clusters.
  if (static_cast<int>(selected.size()) < maxShow) {
    std::vector<char> already(seeds.size(), 0);
    for (const int idx : selected) already[idx] = 1;

    std::vector<int> leftovers;
    leftovers.reserve(seeds.size() - selected.size());
    for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
      if (!already[i]) leftovers.push_back(i);
    }
    std::shuffle(leftovers.begin(), leftovers.end(), rng);

    for (const int idx : leftovers) {
      selected.push_back(idx);
      if (static_cast<int>(selected.size()) >= maxShow) break;
    }
  }

  return selected;
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
  //   - all individual low-energy conduction-electron seed positions,
  //   - actual microscopic Garfield drift-line points for a balanced subset.
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
  // Incident-angle convention:
  //   theta = 0 deg : normal incidence along local -v.
  //   theta > 0     : tilt away from -v.
  //   phi = 0 deg   : tilt in the +u direction (visible in the u-v GIF).
  //   phi = 90 deg  : tilt in the +w direction (along the anode wires).
  const double thetaDeg = ReadArg(argc, argv, "--theta-deg", 0.0);
  const double phiDeg = ReadArg(argc, argv, "--phi-deg", 0.0);
  const double momentumEv =
      1.e9 * ReadArg(argc, argv, "--momentum-gev", 10.0);

  const int maxVisualElectrons =
      ReadIntArg(argc, argv, "--max-electrons", 40);
  const int maxPerCluster =
      ReadIntArg(argc, argv, "--max-per-cluster", 3);
  const int rngSeed =
      ReadIntArg(argc, argv, "--rng-seed", 12345);
  const int maxPathPoints =
      ReadIntArg(argc, argv, "--max-path-points", 250);
  const std::string sampleMode =
      ReadStringArg(argc, argv, "--sample-mode", "balanced");
  const std::string prefix =
      ReadStringArg(argc, argv, "--output-prefix", "heed_vis");

  if (pitchCm <= 0. || wireDiameterCm <= 0. ||
      gapMinusCm <= 0. || gapPlusCm <= 0. ||
      momentumEv <= 0. || std::abs(thetaDeg) >= 80. ||
      maxVisualElectrons <= 0 ||
      maxPerCluster <= 0 || maxPathPoints < 2) {
    std::cerr << "Invalid input parameters.\n";
    return 2;
  }
  if (u0Cm < -0.5 * pitchCm || u0Cm > 0.5 * pitchCm) {
    std::cerr << "Starting u is outside the central wire cell [-p/2,+p/2].\n";
    return 2;
  }
  if (sampleMode != "balanced" &&
      sampleMode != "random" &&
      sampleMode != "first") {
    std::cerr << "--sample-mode must be balanced, random, or first.\n";
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

  const double pi = std::acos(-1.0);
  const double theta = thetaDeg * pi / 180.;
  const double phi = phiDeg * pi / 180.;
  const double dirU = std::sin(theta) * std::cos(phi);
  const double dirV = -std::cos(theta);
  const double dirW = std::sin(theta) * std::sin(phi);

  // Expected intercept at the lower cathode, useful for rendering the
  // projected muon line. Heed itself is still given the direction cosines.
  const double vEndCm = -gapMinusCm;
  const double flightToLowerCm = (vStartCm - vEndCm) / std::cos(theta);
  const double uEndCm = u0Cm + flightToLowerCm * dirU;
  const double wEndCm = flightToLowerCm * dirW;

  heed.NewTrack(u0Cm, vStartCm, 0., 0.,
                dirU, dirV, dirW);

  const std::string metaName = prefix + "_meta.txt";
  const std::string clusterName = prefix + "_clusters.csv";
  const std::string electronName = prefix + "_electrons.csv";
  const std::string pathName = prefix + "_paths.csv";

  std::ofstream metaOut(metaName);
  std::ofstream clusterOut(clusterName);

  if (!metaOut || !clusterOut) {
    std::cerr << "Could not open metadata/cluster output files.\n";
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
          << "muon_u_start_mm=" << 10. * u0Cm << "\n"
          << "muon_u_end_mm=" << 10. * uEndCm << "\n"
          << "muon_v_start_mm=" << 10. * vStartCm << "\n"
          << "muon_v_end_mm=" << 10. * vEndCm << "\n"
          << "muon_w_start_mm=0\n"
          << "muon_w_end_mm=" << 10. * wEndCm << "\n"
          << "theta_deg=" << thetaDeg << "\n"
          << "phi_deg=" << phiDeg << "\n"
          << "dir_u=" << dirU << "\n"
          << "dir_v=" << dirV << "\n"
          << "dir_w=" << dirW << "\n"
          << "sample_mode=" << sampleMode << "\n"
          << "max_per_cluster=" << maxPerCluster << "\n"
          << "rng_seed=" << rngSeed << "\n";

  clusterOut << "cluster,u_mm,v_mm,w_mm,t_ns,electrons,energy_transfer_eV\n";

  int nClusters = 0;
  int nConduction = 0;
  double totalEnergyEv = 0.;
  std::vector<ElectronSeed> seeds;

  double xc = 0., yc = 0., zc = 0., tc = 0.;
  double ec = 0., extra = 0.;
  int nc = 0;

  // First pass: let Heed generate the complete event and store all
  // conduction-electron seed positions. No Garfield drift yet.
  while (heed.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    const int clusterId = nClusters++;
    totalEnergyEv += ec;

    clusterOut << clusterId << ","
               << 10. * xc << "," << 10. * yc << "," << 10. * zc << ","
               << tc << "," << nc << "," << ec << "\n";

    for (int ie = 0; ie < nc; ++ie) {
      double xe = 0., ye = 0., ze = 0., te = 0.;
      double ee = 0., dxe = 0., dye = 0., dze = 0.;
      if (!heed.GetElectron(ie, xe, ye, ze, te,
                            ee, dxe, dye, dze)) {
        continue;
      }

      ElectronSeed seed;
      seed.electronId = nConduction++;
      seed.clusterId = clusterId;
      seed.x = xe;
      seed.y = ye;
      seed.z = ze;
      seed.t = te;
      seed.heedEnergy = ee;
      seed.dx = dxe;
      seed.dy = dye;
      seed.dz = dze;
      seeds.push_back(seed);
    }
  }

  // Choose the subset for presentation AFTER seeing the whole Heed event.
  // This is the key change relative to the first visualiser: a large early
  // cluster can no longer monopolise the displayed Garfield trajectories.
  std::vector<int> selected;
  if (sampleMode == "balanced") {
    selected = SelectBalanced(
        seeds, maxVisualElectrons, maxPerCluster,
        static_cast<unsigned int>(rngSeed));
  } else if (sampleMode == "random") {
    selected = SelectRandom(
        static_cast<int>(seeds.size()), maxVisualElectrons,
        static_cast<unsigned int>(rngSeed));
  } else {
    selected = SelectFirst(
        static_cast<int>(seeds.size()), maxVisualElectrons);
  }

  std::vector<char> isSelected(seeds.size(), 0);
  for (const int idx : selected) {
    if (idx >= 0 && idx < static_cast<int>(isSelected.size())) {
      isSelected[idx] = 1;
    }
  }

  std::ofstream electronOut(electronName);
  std::ofstream pathOut(pathName);
  if (!electronOut || !pathOut) {
    std::cerr << "Could not open electron/path output files.\n";
    return 3;
  }

  electronOut << "electron,cluster,u_mm,v_mm,w_mm,t_ns,heed_energy_eV,"
                 "heed_dx,heed_dy,heed_dz,drifted\n";
  pathOut << "electron,cluster,point,u_mm,v_mm,w_mm,t_ns,status\n";

  for (int idx = 0; idx < static_cast<int>(seeds.size()); ++idx) {
    const auto& seed = seeds[idx];
    electronOut << seed.electronId << "," << seed.clusterId << ","
                << 10. * seed.x << "," << 10. * seed.y << "," << 10. * seed.z
                << "," << seed.t << "," << seed.heedEnergy << ","
                << seed.dx << "," << seed.dy << "," << seed.dz << ","
                << (isSelected[idx] ? 1 : 0) << "\n";
  }

  // Second pass: Garfield microscopic drift only for the selected,
  // presentation-friendly subset.
  int nVisualWithPath = 0;
  for (const int idx : selected) {
    if (idx < 0 || idx >= static_cast<int>(seeds.size())) continue;
    const auto& seed = seeds[idx];

    Garfield::AvalancheMicroscopic drift;
    drift.SetSensor(&sensor);
    drift.EnableDriftLines();

    const bool ok = drift.DriftElectron(
        seed.x, seed.y, seed.z, seed.t,
        0.1, 0., 0., 0.);

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
        std::max<std::size_t>(
            1, (nPoints + maxPathPoints - 1) / maxPathPoints);

    int storedPoint = 0;
    for (std::size_t ip = 0; ip < nPoints; ip += stride) {
      double x = 0., y = 0., z = 0., t = 0.;
      drift.GetElectronDriftLinePoint(x, y, z, t, ip, 0);
      pathOut << seed.electronId << "," << seed.clusterId << ","
              << storedPoint++ << ","
              << 10. * x << "," << 10. * y << "," << 10. * z << ","
              << t << "," << status << "\n";
    }

    double xl = 0., yl = 0., zl = 0., tl = 0.;
    drift.GetElectronDriftLinePoint(
        xl, yl, zl, tl, nPoints - 1, 0);
    pathOut << seed.electronId << "," << seed.clusterId << ","
            << storedPoint << ","
            << 10. * xl << "," << 10. * yl << "," << 10. * zl << ","
            << tl << "," << status << "\n";

    ++nVisualWithPath;
  }

  metaOut << "n_clusters=" << nClusters << "\n"
          << "n_conduction_electrons=" << nConduction << "\n"
          << "n_selected_electrons=" << selected.size() << "\n"
          << "n_visualised_electrons=" << nVisualWithPath << "\n"
          << "total_energy_transfer_eV=" << totalEnergyEv << "\n";

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "\n=== HEED / GARFIELD MUON VISUALISATION EVENT ===\n"
            << "particle             : mu-\n"
            << "momentum             : " << momentumEv * 1.e-9 << " GeV/c\n"
            << "gas                  : Ar/CO2 70/30\n"
            << "wire pitch           : " << 10. * pitchCm << " mm\n"
            << "gas gap              : " << 10. * gapMinusCm << " + "
            << 10. * gapPlusCm << " mm\n"
            << "anode voltage        : " << hv << " V\n"
            << "B                    : " << bTesla << " T\n"
            << "muon start u         : " << 10. * u0Cm << " mm\n"
            << "incident theta       : " << thetaDeg
            << " deg from local -v\n"
            << "incident phi         : " << phiDeg
            << " deg (0 -> +u, 90 -> +w)\n"
            << "projected lower hit  : u=" << 10. * uEndCm
            << " mm, w=" << 10. * wEndCm << " mm\n"
            << "clusters             : " << nClusters << "\n"
            << "conduction electrons : " << nConduction << "\n"
            << "sampling mode        : " << sampleMode << "\n"
            << "selected e-          : " << selected.size() << "\n"
            << "Garfield paths       : " << nVisualWithPath << "\n"
            << "total dE             : " << totalEnergyEv << " eV\n\n"
            << "Wrote:\n"
            << "  " << metaName << "\n"
            << "  " << clusterName << "\n"
            << "  " << electronName << "\n"
            << "  " << pathName << "\n";

  return 0;
}

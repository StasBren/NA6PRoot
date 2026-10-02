#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
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

std::string ReadStringArg(const int argc, char** argv,
                          const std::string& key,
                          const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == key) return argv[i + 1];
  }
  return defaultValue;
}

std::string StripLabel(const int i) {
  if (i < 0) return "strip_m" + std::to_string(-i);
  if (i > 0) return "strip_p" + std::to_string(i);
  return "strip_0";
}

struct PairStart {
  double u = 0.;
  double v = 0.;
  double w = 0.;
};

struct ElectronEnd {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
  int status = -999;
  bool success = false;
};

struct IonEnd {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
  int status = -999;
  bool success = false;
};

void ConfigureElectron(Garfield::AvalancheMicroscopic& electron) {
  // Controlled one-electron transport: DriftElectron follows the seed
  // electron without transporting avalanche secondaries. We use microscopic
  // tracking because the starts are only a few microns from the anode wire.
  electron.EnableSignalCalculation(true);
  electron.UseWeightingPotential(true);
}

void ConfigureIon(Garfield::DriftLineRKF& ion,
                  const int averagingOrder,
                  const double maxStepCm) {
  ion.EnableSignalCalculation(true);
  ion.UseWeightingPotential(true);
  ion.SetSignalAveragingOrder(
      static_cast<std::size_t>(averagingOrder));
  ion.SetMaximumStepSize(maxStepCm);
}

}  // namespace

int main(int argc, char** argv) {
  // Local chamber coordinates:
  //   Garfield x = u : across anode wires
  //   Garfield y = v : chamber normal
  //   Garfield z = w : along anode wires
  //
  // This is a deliberately synthetic Stage-B model. Each sampled point
  // represents one electron-positive-ion pair born close to the anode wire.
  // The spatial cloud is identical in spirit to phaseB_ion_cloud_signal.cxx.
  // It is NOT yet a calibrated microscopic avalanche distribution.

  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double wirePitchCm =
      0.1 * ReadArg(argc, argv, "--wire-pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);

  const double stripPitchCm =
      0.1 * ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthCm =
      0.1 * ReadArg(argc, argv, "--strip-width-mm", 1.7);
  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 5);

  const std::string side =
      ReadStringArg(argc, argv, "--side", "minus");

  const double wireRadiusMm = 0.5 * 10. * wireDiameterCm;
  const double defaultOffsetMm = wireRadiusMm + 0.005;
  const double cloudU0Mm =
      ReadArg(argc, argv, "--cloud-u0-mm", 0.0);
  const double cloudW0Mm =
      ReadArg(argc, argv, "--cloud-w0-mm", 0.0);
  const double cloudOffsetMm =
      ReadArg(argc, argv, "--cloud-offset-mm", defaultOffsetMm);

  const double sigmaUUm =
      ReadArg(argc, argv, "--sigma-u-um", 5.0);
  const double sigmaVUm =
      ReadArg(argc, argv, "--sigma-v-um", 5.0);
  const double sigmaWUm =
      ReadArg(argc, argv, "--sigma-w-um", 200.0);

  const int nPairs =
      ReadIntArg(argc, argv, "--pairs", 500);
  const double pairWeight =
      ReadArg(argc, argv, "--pair-weight", 1.0);
  const int randomSeed =
      ReadIntArg(argc, argv, "--seed", 12345);

  const double electronEnergyEv =
      ReadArg(argc, argv, "--electron-energy-ev", 0.1);
  const double electronDtNs =
      ReadArg(argc, argv, "--electron-dt-ns", 0.001);
  const double electronTmaxNs =
      ReadArg(argc, argv, "--electron-tmax-ns", 0.10);

  const double ionDtNs =
      ReadArg(argc, argv, "--ion-dt-ns", 5.0);
  const double tMaxUs =
      ReadArg(argc, argv, "--tmax-us", 30.0);
  const double maxStepMm =
      ReadArg(argc, argv, "--max-step-mm", 0.05);
  const int signalAveragingOrder =
      ReadIntArg(argc, argv, "--signal-averaging-order", 2);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "stageB_ei_cloud_2mm");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      nPairs < 1 || pairWeight <= 0. ||
      sigmaUUm < 0. || sigmaVUm < 0. || sigmaWUm < 0. ||
      electronEnergyEv <= 0. || electronDtNs <= 0. ||
      electronTmaxNs <= 0. || ionDtNs <= 0. ||
      tMaxUs <= 0. || maxStepMm <= 0. ||
      signalAveragingOrder < 1) {
    std::cerr << "Invalid Stage-B electron-ion cloud parameters.\n";
    return 2;
  }

  if (side != "plus" && side != "minus") {
    std::cerr << "--side must be 'plus' or 'minus'.\n";
    return 2;
  }

  const double sign = side == "plus" ? +1. : -1.;
  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const double gapCm = std::abs(planeY);
  const double cloudVCentreMm = sign * cloudOffsetMm;

  if (cloudOffsetMm <= wireRadiusMm || 0.1 * cloudOffsetMm >= gapCm) {
    std::cerr << "Cloud centre must lie outside the wire and inside "
                 "the selected gas gap.\n";
    return 2;
  }

  // ------------------------------------------------------------------
  // 1) Physical MWPC drift field.
  // ------------------------------------------------------------------
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  if (!gas.LoadIonMobility(mobilityFile)) {
    std::cerr << "Could not load ion mobility file: "
              << mobilityFile << "\n";
    return 3;
  }

  Garfield::ComponentAnalyticField driftField;
  driftField.SetMedium(&gas);

  constexpr int halfNumberOfWires = 4;
  for (int i = -halfNumberOfWires; i <= halfNumberOfWires; ++i) {
    driftField.AddWire(i * wirePitchCm, 0.,
                       wireDiameterCm, hv, "anode");
  }
  driftField.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  driftField.AddPlaneY(+gapPlusCm, 0., "cathode_plus");

  // ------------------------------------------------------------------
  // 2) Same ideal strip-weighting model used by the ions-only baseline.
  // ------------------------------------------------------------------
  Garfield::ComponentAnalyticField weighting;
  weighting.AddPlaneY(0., 1., "weighting_back");
  weighting.AddPlaneY(planeY, 0., "weighting_front");

  std::vector<std::string> labels;
  labels.reserve(2 * halfStrips + 1);
  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double centerW = i * stripPitchCm;
    const double wMin = centerW - 0.5 * stripWidthCm;
    const double wMax = centerW + 0.5 * stripWidthCm;
    const std::string label = StripLabel(i);
    weighting.AddStripOnPlaneY('x', planeY, wMin, wMax, label, gapCm);
    labels.push_back(label);
  }

  Garfield::Sensor sensor;
  sensor.AddComponent(&driftField);
  for (const auto& label : labels) {
    sensor.AddElectrode(&weighting, label);
  }

  const double xExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double zExtent =
      (halfStrips + 1.5) * stripPitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -zExtent,
                 xExtent, +gapPlusCm, +zExtent);

  constexpr double elementaryChargeFc = 1.602176634e-4;
  const double wireRadiusCm = 0.5 * wireDiameterCm;

  // ------------------------------------------------------------------
  // 3) Sample the pair cloud ONCE.
  //
  // The random-number sequence and acceptance logic intentionally match the
  // ions-only generator so the same seed/parameters give the same starts.
  // ------------------------------------------------------------------
  std::mt19937 rng(static_cast<std::uint32_t>(randomSeed));
  std::normal_distribution<double> gaussU(0., 1.e-4 * sigmaUUm);
  std::normal_distribution<double> gaussV(0., 1.e-4 * sigmaVUm);
  std::normal_distribution<double> gaussW(0., 1.e-4 * sigmaWUm);

  const double cloudU0Cm = 0.1 * cloudU0Mm;
  const double cloudV0Cm = 0.1 * cloudVCentreMm;
  const double cloudW0Cm = 0.1 * cloudW0Mm;

  std::vector<PairStart> starts;
  starts.reserve(nPairs);

  constexpr int maxAttemptsPerPair = 10000;
  for (int iPair = 0; iPair < nPairs; ++iPair) {
    PairStart p;
    bool accepted = false;

    for (int attempt = 0; attempt < maxAttemptsPerPair; ++attempt) {
      p.u = cloudU0Cm + gaussU(rng);
      p.v = cloudV0Cm + sign * gaussV(rng);
      p.w = cloudW0Cm + gaussW(rng);

      const double rho = std::hypot(p.u, p.v);
      const bool outsideWire = rho > wireRadiusCm;
      const bool onSelectedSide = sign * p.v > 0.;
      const bool insideGap = std::abs(p.v) < gapCm;

      if (outsideWire && onSelectedSide && insideGap) {
        accepted = true;
        break;
      }
    }

    if (!accepted) {
      std::cerr << "Failed to sample a valid pair after "
                << maxAttemptsPerPair << " attempts.\n";
      return 4;
    }
    starts.push_back(p);
  }

  std::vector<std::vector<double>> phiStart(
      starts.size(), std::vector<double>(labels.size(), 0.));
  for (std::size_t ip = 0; ip < starts.size(); ++ip) {
    for (std::size_t j = 0; j < labels.size(); ++j) {
      phiStart[ip][j] = weighting.WeightingPotential(
          starts[ip].u, starts[ip].v, starts[ip].w, labels[j]);
    }
  }

  // ------------------------------------------------------------------
  // 4) Fine electron run.
  //
  // This gives the real raw prompt-electron current profile. We do NOT force
  // it onto the 5-ns ion grid because that would destroy its peak current.
  // For finite charge windows >= 50 ns, the full electron charge can simply
  // be added to the ion integral after this pulse has ended.
  // ------------------------------------------------------------------
  const int nElectronBins =
      static_cast<int>(std::ceil(electronTmaxNs / electronDtNs));
  sensor.SetTimeWindow(0., electronDtNs, nElectronBins);

  Garfield::AvalancheMicroscopic electron(&sensor);
  ConfigureElectron(electron);

  std::vector<ElectronEnd> electronEnds(starts.size());
  std::vector<std::vector<double>> phiElectronEnd(
      starts.size(), std::vector<double>(labels.size(), 0.));

  int nElectronSuccessful = 0;
  double maxElectronEndTimeNs = 0.;

  for (std::size_t ip = 0; ip < starts.size(); ++ip) {
    const auto& p = starts[ip];

    // Initial microscopic direction points radially toward the anode-wire
    // centre in the u-v plane.
    const double rho = std::hypot(p.u, p.v);
    const double du = rho > 0. ? -p.u / rho : 0.;
    const double dv = rho > 0. ? -p.v / rho : -sign;

    const bool ok = electron.DriftElectron(
        p.u, p.v, p.w, 0.,
        electronEnergyEv, du, dv, 0.);

    ElectronEnd end;
    end.success = ok;

    if (electron.GetNumberOfElectronEndpoints() > 0) {
      double us = 0., vs = 0., ws = 0., ts = 0.;
      double es = electronEnergyEv, ee = electronEnergyEv;
      electron.GetElectronEndpoint(
          0, us, vs, ws, ts, es,
          end.u, end.v, end.w, end.t, ee, end.status);
    } else {
      end.u = p.u;
      end.v = p.v;
      end.w = p.w;
      end.t = 0.;
      end.status = -999;
    }

    if (ok) ++nElectronSuccessful;
    maxElectronEndTimeNs = std::max(maxElectronEndTimeNs, end.t);
    electronEnds[ip] = end;

    for (std::size_t j = 0; j < labels.size(); ++j) {
      phiElectronEnd[ip][j] = weighting.WeightingPotential(
          end.u, end.v, end.w, labels[j]);
    }
  }

  const std::string electronWaveFile =
      outputPrefix + "_electron_waveforms.csv";
  std::ofstream eWaveOut(electronWaveFile);
  eWaveOut << "time_ns";
  for (const auto& label : labels) {
    eWaveOut << "," << label << "_electron_fC_per_ns";
  }
  eWaveOut << "\n";

  std::vector<double> qElectronFine(labels.size(), 0.);
  std::vector<double> peakElectron(labels.size(), 0.);

  for (int ibin = 0; ibin < nElectronBins; ++ibin) {
    eWaveOut << (ibin + 0.5) * electronDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          pairWeight * sensor.GetElectronSignal(labels[j], ibin);
      eWaveOut << "," << current;
      qElectronFine[j] += current * electronDtNs;
      peakElectron[j] =
          std::max(peakElectron[j], std::abs(current));
    }
    eWaveOut << "\n";
  }

  // ------------------------------------------------------------------
  // 5) Full ion run using the SAME sampled starts.
  // ------------------------------------------------------------------
  sensor.ClearSignal();

  const int nIonBins =
      static_cast<int>(std::ceil(tMaxUs * 1000. / ionDtNs));
  sensor.SetTimeWindow(0., ionDtNs, nIonBins);

  Garfield::DriftLineRKF ion(&sensor);
  ConfigureIon(ion, signalAveragingOrder, 0.1 * maxStepMm);

  std::vector<IonEnd> ionEnds(starts.size());
  std::vector<std::vector<double>> phiIonEnd(
      starts.size(), std::vector<double>(labels.size(), 0.));

  int nIonSuccessful = 0;
  double maxIonEndTimeNs = 0.;

  for (std::size_t ip = 0; ip < starts.size(); ++ip) {
    const auto& p = starts[ip];

    const bool ok = ion.DriftIon(p.u, p.v, p.w, 0.);

    IonEnd end;
    end.success = ok;
    ion.GetEndPoint(end.u, end.v, end.w, end.t, end.status);

    if (ok) ++nIonSuccessful;
    maxIonEndTimeNs = std::max(maxIonEndTimeNs, end.t);
    ionEnds[ip] = end;

    for (std::size_t j = 0; j < labels.size(); ++j) {
      phiIonEnd[ip][j] = weighting.WeightingPotential(
          end.u, end.v, end.w, labels[j]);
    }
  }

  const std::string ionWaveFile =
      outputPrefix + "_ion_waveforms.csv";
  std::ofstream iWaveOut(ionWaveFile);
  iWaveOut << "time_ns";
  for (const auto& label : labels) {
    iWaveOut << "," << label << "_ion_fC_per_ns";
  }
  iWaveOut << "\n";

  std::vector<double> qIon(labels.size(), 0.);
  std::vector<double> peakIon(labels.size(), 0.);

  for (int ibin = 0; ibin < nIonBins; ++ibin) {
    iWaveOut << (ibin + 0.5) * ionDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          pairWeight * sensor.GetIonSignal(labels[j], ibin);
      iWaveOut << "," << current;
      qIon[j] += current * ionDtNs;
      peakIon[j] = std::max(peakIon[j], std::abs(current));
    }
    iWaveOut << "\n";
  }

  // ------------------------------------------------------------------
  // 6) Endpoint Ramo checks and per-pair diagnostic file.
  // ------------------------------------------------------------------
  std::vector<double> qElectronRamo(labels.size(), 0.);
  std::vector<double> qIonRamo(labels.size(), 0.);

  for (std::size_t ip = 0; ip < starts.size(); ++ip) {
    for (std::size_t j = 0; j < labels.size(); ++j) {
      qElectronRamo[j] += pairWeight * (-elementaryChargeFc) *
          (phiElectronEnd[ip][j] - phiStart[ip][j]);
      qIonRamo[j] += pairWeight * (+elementaryChargeFc) *
          (phiIonEnd[ip][j] - phiStart[ip][j]);
    }
  }

  const std::string pairsFile =
      outputPrefix + "_pairs.csv";
  std::ofstream pairsOut(pairsFile);
  pairsOut
      << "pair,u0_mm,v0_mm,w0_mm,"
      << "electron_u1_mm,electron_v1_mm,electron_w1_mm,"
      << "electron_t_ns,electron_status,electron_success,"
      << "ion_u1_mm,ion_v1_mm,ion_w1_mm,"
      << "ion_t_ns,ion_status,ion_success\n";

  for (std::size_t ip = 0; ip < starts.size(); ++ip) {
    const auto& p = starts[ip];
    const auto& e = electronEnds[ip];
    const auto& q = ionEnds[ip];

    pairsOut
        << ip << ","
        << 10. * p.u << "," << 10. * p.v << "," << 10. * p.w << ","
        << 10. * e.u << "," << 10. * e.v << "," << 10. * e.w << ","
        << e.t << "," << e.status << "," << (e.success ? 1 : 0) << ","
        << 10. * q.u << "," << 10. * q.v << "," << 10. * q.w << ","
        << q.t << "," << q.status << "," << (q.success ? 1 : 0)
        << "\n";
  }

  const std::string summaryFile =
      outputPrefix + "_summary.csv";
  std::ofstream sumOut(summaryFile);
  sumOut
      << "strip,center_w_mm,"
      << "Qe_fine_fC,Qi_full_fC,Qtotal_full_fC,"
      << "Ipeak_e_fC_per_ns,Ipeak_i_fC_per_ns,"
      << "Qe_ramo_fC,Qi_ramo_fC,Qtotal_ramo_fC\n";

  std::cout << std::scientific << std::setprecision(6);
  std::cout
      << "\n=== STAGE B: SYNTHETIC ELECTRON + ION CLOUD ===\n"
      << "side                    : " << side << "\n"
      << "gas                     : Ar/CO2 70:30\n"
      << "wire pitch / diameter   : " << 10. * wirePitchCm
      << " mm / " << 1.e4 * wireDiameterCm << " um\n"
      << "anode voltage           : " << hv << " V\n"
      << "strip pitch / width     : " << 10. * stripPitchCm
      << " / " << 10. * stripWidthCm << " mm\n"
      << "cloud centre (u,v,w) mm : (" << cloudU0Mm << ", "
      << cloudVCentreMm << ", " << cloudW0Mm << ")\n"
      << "cloud sigma (u,v,w) um  : (" << sigmaUUm << ", "
      << sigmaVUm << ", " << sigmaWUm << ")\n"
      << "sampled pairs           : " << nPairs << "\n"
      << "pair weight             : " << pairWeight << "\n"
      << "random seed             : " << randomSeed << "\n"
      << "electron energy         : " << electronEnergyEv << " eV\n"
      << "electron success        : " << nElectronSuccessful
      << " / " << nPairs << "\n"
      << "max electron end time   : " << maxElectronEndTimeNs
      << " ns\n"
      << "electron window         : " << electronTmaxNs
      << " ns, dt=" << electronDtNs << " ns\n"
      << "ion success             : " << nIonSuccessful
      << " / " << nPairs << "\n"
      << "max ion end time        : " << maxIonEndTimeNs
      << " ns\n"
      << "ion window              : " << tMaxUs
      << " us, dt=" << ionDtNs << " ns\n\n"
      << "NOTE: the pair cloud is a controlled synthetic Stage-B model, "
         "not yet a calibrated Garfield avalanche distribution.\n"
      << "NOTE: electron and ion waveforms intentionally use separate "
         "time grids. Finite windows >= the electron pulse duration are "
         "combined in post-processing by adding the full electron charge.\n\n"
      << "strip   Qe[fC]          Qi[fC]          Qtot[fC]"
         "        Ipeak_e         Ipeak_i\n";

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const std::size_t j =
        static_cast<std::size_t>(i + halfStrips);

    const double qTotal = qElectronFine[j] + qIon[j];
    const double qTotalRamo = qElectronRamo[j] + qIonRamo[j];

    sumOut
        << i << "," << 10. * i * stripPitchCm << ","
        << qElectronFine[j] << "," << qIon[j] << ","
        << qTotal << ","
        << peakElectron[j] << "," << peakIon[j] << ","
        << qElectronRamo[j] << "," << qIonRamo[j] << ","
        << qTotalRamo << "\n";

    std::cout
        << std::setw(5) << i << "   "
        << std::setw(14) << qElectronFine[j] << "   "
        << std::setw(14) << qIon[j] << "   "
        << std::setw(14) << qTotal << "   "
        << std::setw(14) << peakElectron[j] << "   "
        << std::setw(14) << peakIon[j] << "\n";
  }

  if (maxElectronEndTimeNs > electronTmaxNs) {
    std::cout
        << "\nWARNING: at least one electron extends beyond the fine "
           "signal window. Increase --electron-tmax-ns before interpreting "
           "electron charge.\n";
  }
  if (maxIonEndTimeNs > tMaxUs * 1000.) {
    std::cout
        << "\nWARNING: at least one ion extends beyond the ion signal "
           "window. Increase --tmax-us before interpreting full charge.\n";
  }

  std::cout
      << "\nWrote:\n"
      << "  " << electronWaveFile << "\n"
      << "  " << ionWaveFile << "\n"
      << "  " << summaryFile << "\n"
      << "  " << pairsFile << "\n"
      << "======================================================\n";

  return (nElectronSuccessful == nPairs && nIonSuccessful == nPairs) ? 0 : 5;
}

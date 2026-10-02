#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {

struct IonBirth {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
};

std::vector<IonBirth> gIonBirths;

void RecordIonisation(const double x, const double y, const double z,
                      const double t, const int, const int,
                      Garfield::Medium*) {
  gIonBirths.push_back({x, y, z, t});
}

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

struct IonEnd {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
  int status = -999;
  bool success = false;
};

}  // namespace

int main(int argc, char** argv) {
  // Local coordinates:
  //   Garfield x = u : across anode wires
  //   Garfield y = v : chamber normal
  //   Garfield z = w : along anode wires
  //
  // Stage-B goal:
  //   one REAL microscopic Garfield avalanche from one seed electron,
  //   followed by explicit drift of every positive avalanche ion whose
  //   creation point/time is recorded from Garfield's ionisation callback.
  //
  // This is the replacement for the synthetic electron+ion cloud.
  // It is still a one-seed avalanche, NOT yet a full Heed muon event.

  const double gapMinusCm =
      0.1 * ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusCm =
      0.1 * ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double wirePitchCm =
      0.1 * ReadArg(argc, argv, "--wire-pitch-mm", 4.0);
  const double wireDiameterCm =
      1.e-4 * ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double hv = ReadArg(argc, argv, "--hv", 1800.0);
  const double bTesla = ReadArg(argc, argv, "--b-tesla", 0.0);

  const double stripPitchCm =
      0.1 * ReadArg(argc, argv, "--strip-pitch-mm", 1.7);
  const double stripWidthCm =
      0.1 * ReadArg(argc, argv, "--strip-width-mm", 1.7);
  const int halfStrips =
      ReadIntArg(argc, argv, "--half-strips", 5);

  const std::string side =
      ReadStringArg(argc, argv, "--side", "minus");
  const double seedUCm =
      0.1 * ReadArg(argc, argv, "--seed-u-mm", 0.0);
  const double seedWCm =
      0.1 * ReadArg(argc, argv, "--seed-w-mm", 0.0);
  const double seedDistanceMm =
      ReadArg(argc, argv, "--seed-distance-mm", 1.0);
  const double seedEnergyEv =
      ReadArg(argc, argv, "--seed-energy-ev", 0.1);

  const int avalancheLimit =
      ReadIntArg(argc, argv, "--avalanche-limit", 50000);

  const double electronDtNs =
      ReadArg(argc, argv, "--electron-dt-ns", 0.02);
  const double electronTmaxNs =
      ReadArg(argc, argv, "--electron-tmax-ns", 200.0);

  const double ionDtNs =
      ReadArg(argc, argv, "--ion-dt-ns", 5.0);
  const double ionTmaxUs =
      ReadArg(argc, argv, "--ion-tmax-us", 30.0);
  const double maxIonStepMm =
      ReadArg(argc, argv, "--max-ion-step-mm", 0.05);
  const int signalAveragingOrder =
      ReadIntArg(argc, argv, "--signal-averaging-order", 2);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "stageB_real_avalanche");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      seedDistanceMm <= 0. || seedEnergyEv <= 0. ||
      avalancheLimit < 1 ||
      electronDtNs <= 0. || electronTmaxNs <= 0. ||
      ionDtNs <= 0. || ionTmaxUs <= 0. ||
      maxIonStepMm <= 0. || signalAveragingOrder < 1) {
    std::cerr << "Invalid Stage-B real-avalanche parameters.\n";
    return 2;
  }

  if (side != "minus" && side != "plus") {
    std::cerr << "--side must be 'minus' or 'plus'.\n";
    return 2;
  }

  const double sign = side == "plus" ? +1. : -1.;
  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const double selectedGapCm = std::abs(planeY);
  const double seedVCm = sign * 0.1 * seedDistanceMm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;

  if (std::abs(seedVCm) <= wireRadiusCm ||
      std::abs(seedVCm) >= selectedGapCm) {
    std::cerr << "Seed must be outside the anode wire and inside the "
                 "selected gas gap.\n";
    return 2;
  }

  // ------------------------------------------------------------------
  // 1) Physical MWPC field and gas.
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
  driftField.SetMagneticField(bTesla, 0., 0.);

  // ------------------------------------------------------------------
  // 2) Same ideal strip-weighting model used in the earlier Stage-B tests.
  //
  // Limitation retained deliberately for continuity:
  // the strip weighting problem is the analytic plane-strip approximation,
  // not yet a full 3-D wire + segmented-cathode field map.
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
    weighting.AddStripOnPlaneY(
        'x', planeY, wMin, wMax, label, selectedGapCm);
    labels.push_back(label);
  }

  Garfield::Sensor sensor;
  sensor.AddComponent(&driftField);
  for (const auto& label : labels) {
    sensor.AddElectrode(&weighting, label);
  }

  const double uExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double wExtent =
      (halfStrips + 1.5) * stripPitchCm;
  sensor.SetArea(-uExtent, -gapMinusCm, -wExtent,
                 +uExtent, +gapPlusCm, +wExtent);

  constexpr double elementaryChargeFc = 1.602176634e-4;

  // ------------------------------------------------------------------
  // 3) REAL microscopic Garfield electron avalanche.
  // ------------------------------------------------------------------
  const int nElectronBins =
      static_cast<int>(std::ceil(electronTmaxNs / electronDtNs));
  sensor.SetTimeWindow(0., electronDtNs, nElectronBins);

  gIonBirths.clear();

  Garfield::AvalancheMicroscopic avalanche(&sensor);
  avalanche.EnableAvalancheSizeLimit(
      static_cast<unsigned int>(avalancheLimit));
  avalanche.EnableSignalCalculation(true);
  avalanche.UseWeightingPotential(true);
  avalanche.SetUserHandleIonisation(RecordIonisation);

  // Zero direction is intentional: Garfield samples a random initial
  // direction for a gas electron when the supplied vector has zero norm.
  const bool avalancheOk =
      avalanche.AvalancheElectron(
          seedUCm, seedVCm, seedWCm, 0.,
          seedEnergyEv, 0., 0., 0.);

  int nElectrons = 0;
  int nIonsGarfield = 0;
  avalanche.GetAvalancheSize(nElectrons, nIonsGarfield);

  const std::size_t nElectronEndpoints =
      avalanche.GetNumberOfElectronEndpoints();

  std::vector<double> qElectron(labels.size(), 0.);
  std::vector<double> qElectronRamo(labels.size(), 0.);
  std::vector<double> peakElectron(labels.size(), 0.);

  double maxElectronEndTimeNs = 0.;
  double minElectronStartTimeNs = 1.e99;
  double maxElectronStartTimeNs = 0.;

  // Endpoint Shockley-Ramo sum over every microscopic electron trajectory.
  for (std::size_t ie = 0; ie < nElectronEndpoints; ++ie) {
    double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
    double u1 = 0., v1 = 0., w1 = 0., t1 = 0., e1 = 0.;
    int status = 0;
    avalanche.GetElectronEndpoint(
        ie, u0, v0, w0, t0, e0,
        u1, v1, w1, t1, e1, status);

    minElectronStartTimeNs = std::min(minElectronStartTimeNs, t0);
    maxElectronStartTimeNs = std::max(maxElectronStartTimeNs, t0);
    maxElectronEndTimeNs = std::max(maxElectronEndTimeNs, t1);

    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double phi0 =
          weighting.WeightingPotential(u0, v0, w0, labels[j]);
      const double phi1 =
          weighting.WeightingPotential(u1, v1, w1, labels[j]);
      qElectronRamo[j] +=
          (-elementaryChargeFc) * (phi1 - phi0);
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

  for (int ibin = 0; ibin < nElectronBins; ++ibin) {
    eWaveOut << (ibin + 0.5) * electronDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          sensor.GetElectronSignal(labels[j], ibin);
      eWaveOut << "," << current;
      qElectron[j] += current * electronDtNs;
      peakElectron[j] =
          std::max(peakElectron[j], std::abs(current));
    }
    eWaveOut << "\n";
  }

  // Store microscopic electron endpoints.  The starting points of secondary
  // trajectories show where/when the avalanche develops.
  const std::string electronEndpointsFile =
      outputPrefix + "_electron_endpoints.csv";
  std::ofstream endpointOut(electronEndpointsFile);
  endpointOut
      << "electron,u0_mm,v0_mm,w0_mm,t0_ns,e0_eV,"
      << "u1_mm,v1_mm,w1_mm,t1_ns,e1_eV,status\n";

  for (std::size_t ie = 0; ie < nElectronEndpoints; ++ie) {
    double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
    double u1 = 0., v1 = 0., w1 = 0., t1 = 0., e1 = 0.;
    int status = 0;
    avalanche.GetElectronEndpoint(
        ie, u0, v0, w0, t0, e0,
        u1, v1, w1, t1, e1, status);
    endpointOut
        << ie << ","
        << 10. * u0 << "," << 10. * v0 << "," << 10. * w0 << ","
        << t0 << "," << e0 << ","
        << 10. * u1 << "," << 10. * v1 << "," << 10. * w1 << ","
        << t1 << "," << e1 << "," << status << "\n";
  }

  // ------------------------------------------------------------------
  // 4) Explicit drift of every positive avalanche ion.
  //
  // Garfield's ionisation callback records the actual microscopic creation
  // point/time.  This preserves the real avalanche chronology and spatial
  // distribution instead of imposing a Gaussian cloud by hand.
  // ------------------------------------------------------------------
  sensor.ClearSignal();

  const int nIonBins =
      static_cast<int>(std::ceil(ionTmaxUs * 1000. / ionDtNs));
  sensor.SetTimeWindow(0., ionDtNs, nIonBins);

  Garfield::DriftLineRKF ion(&sensor);
  ion.EnableSignalCalculation(true);
  ion.UseWeightingPotential(true);
  ion.SetSignalAveragingOrder(
      static_cast<std::size_t>(signalAveragingOrder));
  ion.SetMaximumStepSize(0.1 * maxIonStepMm);

  std::vector<IonEnd> ionEnds(gIonBirths.size());
  std::vector<double> qIonRamo(labels.size(), 0.);

  int nIonSuccessful = 0;
  double maxIonEndTimeNs = 0.;
  double meanIonBirthW = 0.;
  double meanIonBirthW2 = 0.;

  for (std::size_t ii = 0; ii < gIonBirths.size(); ++ii) {
    const auto& birth = gIonBirths[ii];

    meanIonBirthW += 10. * birth.w;
    meanIonBirthW2 += 100. * birth.w * birth.w;

    const bool ok =
        ion.DriftIon(birth.u, birth.v, birth.w, birth.t);

    IonEnd end;
    end.success = ok;
    ion.GetEndPoint(end.u, end.v, end.w, end.t, end.status);
    ionEnds[ii] = end;

    if (ok) ++nIonSuccessful;
    maxIonEndTimeNs = std::max(maxIonEndTimeNs, end.t);

    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double phi0 =
          weighting.WeightingPotential(
              birth.u, birth.v, birth.w, labels[j]);
      const double phi1 =
          weighting.WeightingPotential(
              end.u, end.v, end.w, labels[j]);
      qIonRamo[j] +=
          elementaryChargeFc * (phi1 - phi0);
    }
  }

  if (!gIonBirths.empty()) {
    meanIonBirthW /= static_cast<double>(gIonBirths.size());
    meanIonBirthW2 /= static_cast<double>(gIonBirths.size());
  }
  const double sigmaIonBirthWMm =
      std::sqrt(std::max(
          0., meanIonBirthW2 - meanIonBirthW * meanIonBirthW));

  const std::string ionBirthFile =
      outputPrefix + "_ion_births.csv";
  std::ofstream ionBirthOut(ionBirthFile);
  ionBirthOut
      << "ion,u0_mm,v0_mm,w0_mm,t0_ns,"
      << "u1_mm,v1_mm,w1_mm,t1_ns,status,success\n";

  for (std::size_t ii = 0; ii < gIonBirths.size(); ++ii) {
    const auto& b = gIonBirths[ii];
    const auto& e = ionEnds[ii];
    ionBirthOut
        << ii << ","
        << 10. * b.u << "," << 10. * b.v << "," << 10. * b.w << ","
        << b.t << ","
        << 10. * e.u << "," << 10. * e.v << "," << 10. * e.w << ","
        << e.t << "," << e.status << "," << (e.success ? 1 : 0)
        << "\n";
  }

  std::vector<double> qIon(labels.size(), 0.);
  std::vector<double> peakIon(labels.size(), 0.);

  const std::string ionWaveFile =
      outputPrefix + "_ion_waveforms.csv";
  std::ofstream iWaveOut(ionWaveFile);
  iWaveOut << "time_ns";
  for (const auto& label : labels) {
    iWaveOut << "," << label << "_ion_fC_per_ns";
  }
  iWaveOut << "\n";

  for (int ibin = 0; ibin < nIonBins; ++ibin) {
    iWaveOut << (ibin + 0.5) * ionDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          sensor.GetIonSignal(labels[j], ibin);
      iWaveOut << "," << current;
      qIon[j] += current * ionDtNs;
      peakIon[j] =
          std::max(peakIon[j], std::abs(current));
    }
    iWaveOut << "\n";
  }

  // ------------------------------------------------------------------
  // 5) Compact strip and event summaries.
  // ------------------------------------------------------------------
  const std::string stripSummaryFile =
      outputPrefix + "_strip_summary.csv";
  std::ofstream stripOut(stripSummaryFile);
  stripOut
      << "strip,center_w_mm,"
      << "Qe_fC,Qi_fC,Qtotal_fC,"
      << "Qe_ramo_fC,Qi_ramo_fC,Qtotal_ramo_fC,"
      << "Ipeak_e_fC_per_ns,Ipeak_i_fC_per_ns\n";

  std::cout << std::scientific << std::setprecision(6);
  std::cout
      << "\n=== STAGE B: REAL MICROSCOPIC GARFIELD AVALANCHE ===\n"
      << "side                       : " << side << "\n"
      << "gas                        : Ar/CO2 70:30\n"
      << "seed (u,v,w) mm            : ("
      << 10. * seedUCm << ", " << 10. * seedVCm << ", "
      << 10. * seedWCm << ")\n"
      << "seed energy                : " << seedEnergyEv << " eV\n"
      << "B(local u)                 : " << bTesla << " T\n"
      << "wire pitch / diameter      : "
      << 10. * wirePitchCm << " mm / "
      << 1.e4 * wireDiameterCm << " um\n"
      << "anode voltage              : " << hv << " V\n"
      << "strip pitch / width        : "
      << 10. * stripPitchCm << " / "
      << 10. * stripWidthCm << " mm\n"
      << "avalanche success          : "
      << (avalancheOk ? "yes" : "no") << "\n"
      << "Garfield avalanche e / ions: "
      << nElectrons << " / " << nIonsGarfield << "\n"
      << "electron trajectories      : "
      << nElectronEndpoints << "\n"
      << "recorded ionisation births : "
      << gIonBirths.size() << "\n"
      << "ion drift success          : "
      << nIonSuccessful << " / " << gIonBirths.size() << "\n"
      << "electron start-time range  : "
      << minElectronStartTimeNs << " .. "
      << maxElectronStartTimeNs << " ns\n"
      << "max electron end time      : "
      << maxElectronEndTimeNs << " ns\n"
      << "ion-birth sigma_w          : "
      << sigmaIonBirthWMm << " mm\n"
      << "max ion end time           : "
      << maxIonEndTimeNs << " ns\n"
      << "electron signal grid       : dt="
      << electronDtNs << " ns, tmax="
      << electronTmaxNs << " ns\n"
      << "ion signal grid            : dt="
      << ionDtNs << " ns, tmax="
      << ionTmaxUs << " us\n\n"
      << "strip    Qe[fC]          Qi[fC]          Qtot[fC]"
         "        Ipeak_e         Ipeak_i\n";

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const std::size_t j =
        static_cast<std::size_t>(i + halfStrips);
    const double qTotal = qElectron[j] + qIon[j];
    const double qTotalRamo =
        qElectronRamo[j] + qIonRamo[j];

    stripOut
        << i << "," << 10. * i * stripPitchCm << ","
        << qElectron[j] << "," << qIon[j] << "," << qTotal << ","
        << qElectronRamo[j] << "," << qIonRamo[j] << ","
        << qTotalRamo << ","
        << peakElectron[j] << "," << peakIon[j] << "\n";

    std::cout
        << std::setw(5) << i << "   "
        << std::setw(14) << qElectron[j] << "   "
        << std::setw(14) << qIon[j] << "   "
        << std::setw(14) << qTotal << "   "
        << std::setw(14) << peakElectron[j] << "   "
        << std::setw(14) << peakIon[j] << "\n";
  }

  const std::string eventSummaryFile =
      outputPrefix + "_event_summary.csv";
  std::ofstream eventOut(eventSummaryFile);
  eventOut
      << "avalanche_success,garfield_electrons,garfield_ions,"
      << "electron_trajectories,recorded_ion_births,"
      << "successful_ion_drifts,"
      << "electron_start_min_ns,electron_start_max_ns,"
      << "electron_end_max_ns,ion_birth_sigma_w_mm,"
      << "ion_end_max_ns,avalanche_limit\n";
  eventOut
      << (avalancheOk ? 1 : 0) << ","
      << nElectrons << "," << nIonsGarfield << ","
      << nElectronEndpoints << "," << gIonBirths.size() << ","
      << nIonSuccessful << ","
      << minElectronStartTimeNs << ","
      << maxElectronStartTimeNs << ","
      << maxElectronEndTimeNs << ","
      << sigmaIonBirthWMm << ","
      << maxIonEndTimeNs << ","
      << avalancheLimit << "\n";

  if (nElectrons >= avalancheLimit) {
    std::cout
        << "\nWARNING: avalanche reached the configured size limit. "
           "Increase --avalanche-limit before treating the gain/shape "
           "as an unconstrained avalanche.\n";
  }
  if (maxElectronEndTimeNs > electronTmaxNs) {
    std::cout
        << "\nWARNING: electron avalanche extends beyond the electron "
           "signal window. Increase --electron-tmax-ns.\n";
  }
  if (maxIonEndTimeNs > ionTmaxUs * 1000.) {
    std::cout
        << "\nWARNING: ion drift extends beyond the ion signal window. "
           "Increase --ion-tmax-us.\n";
  }
  if (static_cast<int>(gIonBirths.size()) != nIonsGarfield) {
    std::cout
        << "\nNOTE: recorded ionisation callbacks (" << gIonBirths.size()
        << ") differ from Garfield ion count (" << nIonsGarfield
        << "). Inspect this before interpreting absolute ion charge.\n";
  }

  std::cout
      << "\nWrote:\n"
      << "  " << electronWaveFile << "\n"
      << "  " << ionWaveFile << "\n"
      << "  " << electronEndpointsFile << "\n"
      << "  " << ionBirthFile << "\n"
      << "  " << stripSummaryFile << "\n"
      << "  " << eventSummaryFile << "\n"
      << "=======================================================\n";

  return avalancheOk ? 0 : 5;
}

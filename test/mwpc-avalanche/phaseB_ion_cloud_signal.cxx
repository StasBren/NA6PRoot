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

}  // namespace

int main(int argc, char** argv) {
  // Local chamber coordinates and Garfield coordinates:
  //   Garfield x = local u : across anode wires
  //   Garfield y = local v : chamber normal
  //   Garfield z = local w : along anode wires

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
      ReadStringArg(argc, argv, "--side", "plus");

  // Cloud centre. The default radial offset is wire radius + 5 um,
  // as in the converged B1a single-ion test.
  const double wireRadiusMm = 0.5 * 1.e1 * wireDiameterCm;
  const double defaultOffsetMm = wireRadiusMm + 0.005;
  const double cloudU0Mm = ReadArg(argc, argv, "--cloud-u0-mm", 0.0);
  const double cloudW0Mm = ReadArg(argc, argv, "--cloud-w0-mm", 0.0);
  const double cloudOffsetMm =
      ReadArg(argc, argv, "--cloud-offset-mm", defaultOffsetMm);

  // Compact synthetic cloud widths. These are reference values only, not
  // yet calibrated to the real avalanche distribution.
  const double sigmaUUm = ReadArg(argc, argv, "--sigma-u-um", 5.0);
  const double sigmaVUm = ReadArg(argc, argv, "--sigma-v-um", 5.0);
  const double sigmaWUm = ReadArg(argc, argv, "--sigma-w-um", 50.0);

  const int nIons = ReadIntArg(argc, argv, "--ions", 200);
  const double ionWeight = ReadArg(argc, argv, "--ion-weight", 1.0);
  const int randomSeed = ReadIntArg(argc, argv, "--seed", 12345);

  const double dtNs = ReadArg(argc, argv, "--dt-ns", 5.0);
  const double tMaxUs = ReadArg(argc, argv, "--tmax-us", 300.0);
  const double maxStepMm = ReadArg(argc, argv, "--max-step-mm", 0.05);
  const int signalAveragingOrder =
      ReadIntArg(argc, argv, "--signal-averaging-order", 2);

  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "phaseB_ion_cloud");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      nIons < 1 || ionWeight <= 0. ||
      sigmaUUm < 0. || sigmaVUm < 0. || sigmaWUm < 0. ||
      dtNs <= 0. || tMaxUs <= 0. || maxStepMm <= 0. ||
      signalAveragingOrder < 1) {
    std::cerr << "Invalid Phase-B1b ion-cloud parameters.\n";
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

  // As in B1a: Ar+ mobility in Ar is used as a timing approximation.
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
  // 2) Strip weighting-field component for the selected cathode side.
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

  // ------------------------------------------------------------------
  // 3) Sensor and converged B1a numerical settings.
  // ------------------------------------------------------------------
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

  const int nBins =
      static_cast<int>(std::ceil(tMaxUs * 1000. / dtNs));
  sensor.SetTimeWindow(0., dtNs, nBins);

  Garfield::DriftLineRKF ion(&sensor);
  ion.EnableSignalCalculation(true);
  ion.UseWeightingPotential(true);
  ion.SetSignalAveragingOrder(
      static_cast<std::size_t>(signalAveragingOrder));
  ion.SetMaximumStepSize(0.1 * maxStepMm);  // mm -> cm

  // ------------------------------------------------------------------
  // 4) Generate a reproducible compact ion cloud near the anode.
  // ------------------------------------------------------------------
  std::mt19937 rng(static_cast<std::uint32_t>(randomSeed));
  std::normal_distribution<double> gaussU(0., 1.e-4 * sigmaUUm);
  std::normal_distribution<double> gaussV(0., 1.e-4 * sigmaVUm);
  std::normal_distribution<double> gaussW(0., 1.e-4 * sigmaWUm);

  const double cloudU0Cm = 0.1 * cloudU0Mm;
  const double cloudV0Cm = 0.1 * cloudVCentreMm;
  const double cloudW0Cm = 0.1 * cloudW0Mm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;

  const std::string ionsFile = outputPrefix + "_ions.csv";
  std::ofstream ionsOut(ionsFile);
  ionsOut << "ion,u0_mm,v0_mm,w0_mm,u1_mm,v1_mm,w1_mm,"
             "drift_time_ns,status,success\n";

  std::vector<double> ramoSum(labels.size(), 0.);
  int nAccepted = 0;
  int nSuccessful = 0;
  double maxDriftTimeNs = 0.;

  constexpr int maxAttemptsPerIon = 10000;
  constexpr double elementaryChargeFc = 1.602176634e-4;

  for (int iIon = 0; iIon < nIons; ++iIon) {
    double u0 = 0., v0 = 0., w0 = 0.;
    bool accepted = false;

    for (int attempt = 0; attempt < maxAttemptsPerIon; ++attempt) {
      u0 = cloudU0Cm + gaussU(rng);
      // Sigma-v is sampled along the selected outward normal.
      v0 = cloudV0Cm + sign * gaussV(rng);
      w0 = cloudW0Cm + gaussW(rng);

      const double rho = std::hypot(u0, v0);
      const bool outsideWire = rho > wireRadiusCm;
      const bool onSelectedSide = sign * v0 > 0.;
      const bool insideGap = std::abs(v0) < gapCm;

      if (outsideWire && onSelectedSide && insideGap) {
        accepted = true;
        break;
      }
    }

    if (!accepted) {
      std::cerr << "Failed to sample a valid cloud ion after "
                << maxAttemptsPerIon << " attempts.\n";
      return 4;
    }
    ++nAccepted;

    // Record the weighting-potential start values before drifting.
    std::vector<double> phiStart(labels.size(), 0.);
    for (std::size_t j = 0; j < labels.size(); ++j) {
      phiStart[j] =
          weighting.WeightingPotential(u0, v0, w0, labels[j]);
    }

    const bool ok = ion.DriftIon(u0, v0, w0, 0.);

    double u1 = 0., v1 = 0., w1 = 0., t1 = 0.;
    int status = 0;
    ion.GetEndPoint(u1, v1, w1, t1, status);

    if (ok) ++nSuccessful;
    maxDriftTimeNs = std::max(maxDriftTimeNs, t1);

    ionsOut << iIon << ","
            << 10. * u0 << "," << 10. * v0 << "," << 10. * w0 << ","
            << 10. * u1 << "," << 10. * v1 << "," << 10. * w1 << ","
            << t1 << "," << status << "," << (ok ? 1 : 0) << "\n";

    // Endpoint Shockley-Ramo prediction for the total cloud signal.
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double phiEnd =
          weighting.WeightingPotential(u1, v1, w1, labels[j]);
      ramoSum[j] += ionWeight * elementaryChargeFc *
                    (phiEnd - phiStart[j]);
    }
  }

  // ------------------------------------------------------------------
  // 5) Save accumulated cloud waveforms and strip summary.
  // ------------------------------------------------------------------
  const std::string waveformFile =
      outputPrefix + "_waveforms.csv";
  const std::string summaryFile =
      outputPrefix + "_summary.csv";

  std::ofstream waveOut(waveformFile);
  waveOut << "time_ns";
  for (const auto& label : labels) {
    waveOut << "," << label << "_fC_per_ns";
  }
  waveOut << "\n";

  std::vector<double> qInt(labels.size(), 0.);
  std::vector<double> peakAbs(labels.size(), 0.);

  for (int ibin = 0; ibin < nBins; ++ibin) {
    const double t = (ibin + 0.5) * dtNs;
    waveOut << t;

    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          ionWeight * sensor.GetIonSignal(labels[j], ibin);
      waveOut << "," << current;
      qInt[j] += current * dtNs;
      peakAbs[j] = std::max(peakAbs[j], std::abs(current));
    }
    waveOut << "\n";
  }

  double totalAbsCharge = 0.;
  for (const double q : qInt) totalAbsCharge += std::abs(q);

  std::ofstream sumOut(summaryFile);
  sumOut << "strip,center_w_mm,integrated_ion_signal_fC,"
            "peak_abs_current_fC_per_ns,ramo_endpoint_fC,"
            "fraction_of_abs_integrated_charge\n";

  std::cout << std::scientific << std::setprecision(6);
  std::cout << "\n=== PHASE B1b: COMPACT AVALANCHE-ION CLOUD ===\n"
            << "side                    : " << side << "\n"
            << "gas                     : Ar/CO2 70:30\n"
            << "ion mobility model      : " << mobilityFile
            << "  [timing approximation]\n"
            << "wire pitch / diameter   : " << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "anode voltage           : " << hv << " V\n"
            << "strip pitch / width     : " << 10. * stripPitchCm << " / "
            << 10. * stripWidthCm << " mm\n"
            << "cloud centre (u,v,w) mm : (" << cloudU0Mm << ", "
            << cloudVCentreMm << ", " << cloudW0Mm << ")\n"
            << "cloud sigma (u,v,w) um  : (" << sigmaUUm << ", "
            << sigmaVUm << ", " << sigmaWUm << ")\n"
            << "sampled ions            : " << nIons << "\n"
            << "ion weight              : " << ionWeight << "\n"
            << "effective ion count     : " << nIons * ionWeight << "\n"
            << "random seed             : " << randomSeed << "\n"
            << "accepted / successful   : " << nAccepted << " / "
            << nSuccessful << "\n"
            << "max ion drift time      : " << maxDriftTimeNs
            << " ns\n"
            << "signal window           : " << tMaxUs
            << " us, dt=" << dtNs << " ns\n"
            << "RKF max step            : " << maxStepMm << " mm\n"
            << "signal avg. order       : " << signalAveragingOrder
            << "\n\n"
            << "NOTE: cloud widths are a synthetic reference model, not yet"
               " a calibrated Garfield avalanche distribution.\n\n"
            << "strip   center_w[mm]       Q_signal[fC]       |I|_peak[fC/ns]"
               "       Q_Ramo[fC]       |Q| fraction\n";

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const std::size_t j = static_cast<std::size_t>(i + halfStrips);
    const double fraction =
        totalAbsCharge > 0. ? std::abs(qInt[j]) / totalAbsCharge : 0.;

    sumOut << i << "," << 10. * i * stripPitchCm << ","
           << qInt[j] << "," << peakAbs[j] << ","
           << ramoSum[j] << "," << fraction << "\n";

    std::cout << std::setw(5) << i << "   "
              << std::setw(12) << 10. * i * stripPitchCm << "   "
              << std::setw(16) << qInt[j] << "   "
              << std::setw(18) << peakAbs[j] << "   "
              << std::setw(16) << ramoSum[j] << "   "
              << std::setw(12) << fraction << "\n";
  }

  if (maxDriftTimeNs > tMaxUs * 1000.) {
    std::cout << "\nWARNING: at least one ion drift exceeds the signal "
                 "window. Increase --tmax-us before interpreting charge.\n";
  }

  std::cout << "\nWrote " << waveformFile
            << ", " << summaryFile
            << " and " << ionsFile << "\n"
            << "=====================================================\n";

  return nSuccessful == nIons ? 0 : 5;
}

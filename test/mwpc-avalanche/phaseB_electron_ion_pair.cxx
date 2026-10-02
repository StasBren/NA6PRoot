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

struct DriftPoint {
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
};

std::vector<DriftPoint> ExtractDriftLine(const Garfield::DriftLineRKF& drift) {
  const std::size_t n = drift.GetNumberOfDriftLinePoints();
  std::vector<DriftPoint> points;
  points.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    DriftPoint p;
    drift.GetDriftLinePoint(i, p.u, p.v, p.w, p.t);
    points.push_back(p);
  }
  return points;
}

std::vector<DriftPoint> ExtractElectronDriftLine(
    const Garfield::AvalancheMicroscopic& drift) {
  std::vector<DriftPoint> points;
  if (drift.GetNumberOfElectronEndpoints() == 0) return points;

  const std::size_t n = drift.GetNumberOfElectronDriftLinePoints(0);
  points.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    DriftPoint p;
    drift.GetElectronDriftLinePoint(p.u, p.v, p.w, p.t, i, 0);
    points.push_back(p);
  }
  return points;
}

void WriteDriftLine(const std::string& filename,
                    const std::vector<DriftPoint>& points) {
  std::ofstream out(filename);
  out << "index,u_mm,v_mm,w_mm,time_ns,step_mm,dt_ns\n";
  for (std::size_t i = 0; i < points.size(); ++i) {
    double dsMm = 0.;
    double dtNs = 0.;
    if (i > 0) {
      const double du = points[i].u - points[i - 1].u;
      const double dv = points[i].v - points[i - 1].v;
      const double dw = points[i].w - points[i - 1].w;
      dsMm = 10. * std::sqrt(du * du + dv * dv + dw * dw);
      dtNs = points[i].t - points[i - 1].t;
    }
    out << i << "," << 10. * points[i].u << ","
        << 10. * points[i].v << ","
        << 10. * points[i].w << ","
        << points[i].t << "," << dsMm << "," << dtNs << "\n";
  }
}

void ConfigureElectronDrift(Garfield::AvalancheMicroscopic& drift) {
  // Microscopic tracking is used for the electron because the pair is created
  // only a few microns from the anode-wire surface, where the macroscopic RKF
  // drift-line integrator can fail. DriftElectron follows only the seed
  // electron (no secondary-electron transport), which is exactly what we want
  // for this controlled one-electron + one-ion building block.
  drift.EnableSignalCalculation(true);
  drift.UseWeightingPotential(true);
  drift.EnableDriftLines(true);
}

void ConfigureIonDrift(Garfield::DriftLineRKF& drift,
                       const int averagingOrder) {
  drift.EnableSignalCalculation(true);
  drift.UseWeightingPotential(true);
  drift.SetSignalAveragingOrder(
      static_cast<std::size_t>(averagingOrder));
  drift.SetMaximumStepSize();
}

}  // namespace

int main(int argc, char** argv) {
  // Local chamber coordinates:
  //   Garfield x = u : across anode wires
  //   Garfield y = v : chamber normal
  //   Garfield z = w : along anode wires

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
  const double u0Cm =
      0.1 * ReadArg(argc, argv, "--u0-mm", 0.0);
  const double w0Cm =
      0.1 * ReadArg(argc, argv, "--w0-mm", 0.0);

  // The electron-ion pair is created at one common point just outside the
  // anode wire. Default: wire radius + 5 um.
  const double defaultOffsetMm = 0.5 * 10. * wireDiameterCm + 0.005;
  const double pairOffsetMm =
      ReadArg(argc, argv, "--pair-offset-mm", defaultOffsetMm);
  const double electronEnergyEv =
      ReadArg(argc, argv, "--electron-energy-ev", 0.1);

  // Two time resolutions:
  //  - fine electron window resolves the prompt electron pulse;
  //  - coarse window follows the slow ion for the full drift.
  const double electronDtNs =
      ReadArg(argc, argv, "--electron-dt-ns", 0.001);
  const double electronTmaxNs =
      ReadArg(argc, argv, "--electron-tmax-ns", 0.10);
  const double ionDtNs =
      ReadArg(argc, argv, "--ion-dt-ns", 50.0);
  const double tMaxUs =
      ReadArg(argc, argv, "--tmax-us", 150.0);

  const int signalAveragingOrder =
      ReadIntArg(argc, argv, "--signal-averaging-order", 2);
  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "phaseB_electron_ion_pair");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      electronEnergyEv <= 0. ||
      electronDtNs <= 0. || electronTmaxNs <= 0. ||
      ionDtNs <= 0. || tMaxUs <= 0. ||
      signalAveragingOrder < 1) {
    std::cerr << "Invalid Phase-B electron-ion-pair parameters.\n";
    return 2;
  }
  if (side != "plus" && side != "minus") {
    std::cerr << "--side must be 'plus' or 'minus'.\n";
    return 2;
  }

  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const std::string planeReadoutLabel =
      side == "plus" ? "cathode_plus" : "cathode_minus";
  const double gapCm = std::abs(planeY);
  const double sign = side == "plus" ? +1. : -1.;
  const double v0Cm = sign * 0.1 * pairOffsetMm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;

  if (std::abs(v0Cm) <= wireRadiusCm) {
    std::cerr << "Pair starting point lies inside/on the anode wire.\n"
              << "wire radius = " << 10. * wireRadiusCm
              << " mm, |v0| = " << 10. * std::abs(v0Cm) << " mm\n";
    return 2;
  }
  if (std::abs(v0Cm) >= gapCm) {
    std::cerr << "Pair starting point lies outside the selected gas gap.\n";
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

  // Positive-ion timing remains approximate, as in the single-ion test.
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
  // 2) Weighting field for the readout strips on the selected side.
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
    weighting.AddStripOnPlaneY('x', planeY, wMin, wMax,
                               label, gapCm);
    labels.push_back(label);
  }

  Garfield::Sensor sensor;
  sensor.AddComponent(&driftField);
  for (const auto& label : labels) {
    sensor.AddElectrode(&weighting, label);
  }

  // Wire-aware control observable.
  //
  // ComponentAnalyticField's analytic strip weighting formula is a
  // plane-condenser solution and does not impose the perpendicular anode-wire
  // surface as a boundary. The weighting field of the COMPLETE cathode plane,
  // however, is computed by the same analytic cell that contains the actual
  // anode wires. We therefore record it as a control for the near-wire
  // electron signal. On the anode-wire surface its weighting potential must
  // approach zero.
  sensor.AddElectrode(&driftField, planeReadoutLabel);

  const double xExtent =
      (halfNumberOfWires + 0.5) * wirePitchCm;
  const double zExtent =
      (halfStrips + 1.5) * stripPitchCm;
  sensor.SetArea(-xExtent, -gapMinusCm, -zExtent,
                  xExtent, +gapPlusCm, +zExtent);

  const int nElectronBins =
      static_cast<int>(std::ceil(electronTmaxNs / electronDtNs));
  const int nIonBins =
      static_cast<int>(std::ceil(tMaxUs * 1000. / ionDtNs));

  constexpr double elementaryChargeFc = 1.602176634e-4;

  // ------------------------------------------------------------------
  // 3) Fine-time electron-only run.
  // ------------------------------------------------------------------
  sensor.SetTimeWindow(0., electronDtNs, nElectronBins);

  Garfield::AvalancheMicroscopic electron(&sensor);
  ConfigureElectronDrift(electron);

  // Point the initial electron direction toward the anode wire. Subsequent
  // motion is microscopic and includes collisions in the gas.
  const double electronDirectionV = -sign;
  const bool electronOk =
      electron.DriftElectron(u0Cm, v0Cm, w0Cm, 0.,
                             electronEnergyEv,
                             0., electronDirectionV, 0.);

  double ue = u0Cm, ve = v0Cm, we = w0Cm, te = 0.;
  double eStart = electronEnergyEv, eEnd = electronEnergyEv;
  int electronStatus = -3;
  if (electron.GetNumberOfElectronEndpoints() > 0) {
    double us = 0., vs = 0., ws = 0., ts = 0.;
    electron.GetElectronEndpoint(0, us, vs, ws, ts, eStart,
                                 ue, ve, we, te, eEnd, electronStatus);
  }
  const auto electronPoints = ExtractElectronDriftLine(electron);

  std::vector<std::vector<double>> electronFine(
      labels.size(), std::vector<double>(nElectronBins, 0.));
  std::vector<double> qElectronFine(labels.size(), 0.);
  std::vector<double> planeElectronFine(nElectronBins, 0.);
  double qPlaneElectronFine = 0.;

  for (int ibin = 0; ibin < nElectronBins; ++ibin) {
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          sensor.GetElectronSignal(labels[j], ibin);
      electronFine[j][ibin] = current;
      qElectronFine[j] += current * electronDtNs;
    }
    planeElectronFine[ibin] =
        sensor.GetElectronSignal(planeReadoutLabel, ibin);
    qPlaneElectronFine += planeElectronFine[ibin] * electronDtNs;
  }

  const std::string electronWaveFile =
      outputPrefix + "_electron_waveforms.csv";
  std::ofstream eWaveOut(electronWaveFile);
  eWaveOut << "time_ns";
  for (const auto& label : labels) {
    eWaveOut << "," << label << "_electron_fC_per_ns";
  }
  eWaveOut << ",readout_plane_electron_fC_per_ns\n";
  for (int ibin = 0; ibin < nElectronBins; ++ibin) {
    eWaveOut << (ibin + 0.5) * electronDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      eWaveOut << "," << electronFine[j][ibin];
    }
    eWaveOut << "," << planeElectronFine[ibin] << "\n";
  }

  WriteDriftLine(outputPrefix + "_electron_driftline.csv",
                 electronPoints);

  // ------------------------------------------------------------------
  // 4) Coarse-time ion-only run.
  // ------------------------------------------------------------------
  sensor.ClearSignal();
  sensor.SetTimeWindow(0., ionDtNs, nIonBins);

  Garfield::DriftLineRKF ion(&sensor);
  ConfigureIonDrift(ion, signalAveragingOrder);

  const bool ionOk =
      ion.DriftIon(u0Cm, v0Cm, w0Cm, 0.);

  double ui = 0., vi = 0., wi = 0., ti = 0.;
  int ionStatus = 0;
  ion.GetEndPoint(ui, vi, wi, ti, ionStatus);
  const auto ionPoints = ExtractDriftLine(ion);

  std::vector<std::vector<double>> ionOnly(
      labels.size(), std::vector<double>(nIonBins, 0.));
  std::vector<double> qIonOnly(labels.size(), 0.);
  std::vector<double> planeIonOnly(nIonBins, 0.);

  for (int ibin = 0; ibin < nIonBins; ++ibin) {
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current =
          sensor.GetIonSignal(labels[j], ibin);
      ionOnly[j][ibin] = current;
      qIonOnly[j] += current * ionDtNs;
    }
    planeIonOnly[ibin] =
        sensor.GetIonSignal(planeReadoutLabel, ibin);
  }

  const std::string ionWaveFile =
      outputPrefix + "_ion_waveforms.csv";
  std::ofstream iWaveOut(ionWaveFile);
  iWaveOut << "time_ns";
  for (const auto& label : labels) {
    iWaveOut << "," << label << "_ion_fC_per_ns";
  }
  iWaveOut << ",readout_plane_ion_fC_per_ns\n";
  for (int ibin = 0; ibin < nIonBins; ++ibin) {
    iWaveOut << (ibin + 0.5) * ionDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      iWaveOut << "," << ionOnly[j][ibin];
    }
    iWaveOut << "," << planeIonOnly[ibin] << "\n";
  }

  WriteDriftLine(outputPrefix + "_ion_driftline.csv",
                 ionPoints);

  // ------------------------------------------------------------------
  // 5) Electron + ion together on ONE coarse time grid.
  //
  // Sensor accumulates drift-line signals. Query electron and ion components
  // separately and also GetSignal = their sum. This is an explicit linearity
  // check. The prompt electron is intentionally under-resolved here; its
  // detailed shape is stored by the fine-time run above.
  // ------------------------------------------------------------------
  sensor.ClearSignal();
  sensor.SetTimeWindow(0., ionDtNs, nIonBins);

  Garfield::AvalancheMicroscopic electronCombined(&sensor);
  ConfigureElectronDrift(electronCombined);
  const bool electronCombinedOk =
      electronCombined.DriftElectron(u0Cm, v0Cm, w0Cm, 0.,
                                     electronEnergyEv,
                                     0., electronDirectionV, 0.);

  Garfield::DriftLineRKF ionCombined(&sensor);
  ConfigureIonDrift(ionCombined, signalAveragingOrder);
  const bool ionCombinedOk =
      ionCombined.DriftIon(u0Cm, v0Cm, w0Cm, 0.);

  const std::string combinedWaveFile =
      outputPrefix + "_combined_waveforms.csv";
  std::ofstream cWaveOut(combinedWaveFile);
  cWaveOut << "time_ns";
  for (const auto& label : labels) {
    cWaveOut << "," << label << "_electron_fC_per_ns"
             << "," << label << "_ion_fC_per_ns"
             << "," << label << "_total_fC_per_ns";
  }
  cWaveOut << ",readout_plane_electron_fC_per_ns"
           << ",readout_plane_ion_fC_per_ns"
           << ",readout_plane_total_fC_per_ns\n";

  std::vector<double> qElectronCoarse(labels.size(), 0.);
  std::vector<double> qIonCoarse(labels.size(), 0.);
  std::vector<double> qTotalCoarse(labels.size(), 0.);
  double maxLinearityResidual = 0.;
  double qPlaneElectronCoarse = 0.;
  double qPlaneIonCoarse = 0.;
  double qPlaneTotalCoarse = 0.;

  for (int ibin = 0; ibin < nIonBins; ++ibin) {
    cWaveOut << (ibin + 0.5) * ionDtNs;
    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double ie =
          sensor.GetElectronSignal(labels[j], ibin);
      const double ii =
          sensor.GetIonSignal(labels[j], ibin);
      const double it =
          sensor.GetSignal(labels[j], ibin);

      cWaveOut << "," << ie << "," << ii << "," << it;

      qElectronCoarse[j] += ie * ionDtNs;
      qIonCoarse[j] += ii * ionDtNs;
      qTotalCoarse[j] += it * ionDtNs;
      maxLinearityResidual =
          std::max(maxLinearityResidual, std::abs(it - ie - ii));
    }
    const double planeIe =
        sensor.GetElectronSignal(planeReadoutLabel, ibin);
    const double planeIi =
        sensor.GetIonSignal(planeReadoutLabel, ibin);
    const double planeIt =
        sensor.GetSignal(planeReadoutLabel, ibin);
    cWaveOut
        << "," << planeIe
        << "," << planeIi
        << "," << planeIt
        << "\n";
    qPlaneElectronCoarse += planeIe * ionDtNs;
    qPlaneIonCoarse += planeIi * ionDtNs;
    qPlaneTotalCoarse += planeIt * ionDtNs;
  }

  // ------------------------------------------------------------------
  // 6) Endpoint Shockley-Ramo check for each strip.
  // ------------------------------------------------------------------
  const std::string summaryFile =
      outputPrefix + "_summary.csv";
  std::ofstream summaryOut(summaryFile);
  summaryOut
      << "strip,center_w_mm,phi_start,phi_electron_end,phi_ion_end,"
      << "q_electron_expected_fC,q_ion_expected_fC,q_total_expected_fC,"
      << "q_electron_fine_fC,q_electron_coarse_fC,"
      << "q_ion_coarse_fC,q_total_coarse_fC\n";

  std::cout << std::fixed << std::setprecision(9);
  std::cout
      << "\n=== PHASE B: ONE ELECTRON + ONE POSITIVE ION ===\n"
      << "side                    : " << side << "\n"
      << "gas                     : Ar/CO2 70:30\n"
      << "pair start (u,v,w)      : (" << 10. * u0Cm << ", "
      << 10. * v0Cm << ", " << 10. * w0Cm << ") mm\n"
      << "wire pitch / diameter   : " << 10. * wirePitchCm
      << " mm / " << 1.e4 * wireDiameterCm << " um\n"
      << "anode voltage           : " << hv << " V\n"
      << "strip pitch / width     : " << 10. * stripPitchCm
      << " / " << 10. * stripWidthCm << " mm\n"
      << "electron initial energy : " << electronEnergyEv << " eV\n"
      << "electron fine window    : 0 .. " << electronTmaxNs
      << " ns, dt=" << electronDtNs << " ns\n"
      << "ion/combined window     : 0 .. " << tMaxUs
      << " us, dt=" << ionDtNs << " ns\n"
      << "electron end (u,v,w,t)  : (" << 10. * ue << ", "
      << 10. * ve << ", " << 10. * we << ") mm, "
      << te << " ns, status=" << electronStatus
      << ", success=" << (electronOk ? "yes" : "no") << "\n"
      << "ion end (u,v,w,t)       : (" << 10. * ui << ", "
      << 10. * vi << ", " << 10. * wi << ") mm, "
      << ti << " ns, status=" << ionStatus
      << ", success=" << (ionOk ? "yes" : "no") << "\n"
      << "combined re-run success : electron="
      << (electronCombinedOk ? "yes" : "no")
      << ", ion=" << (ionCombinedOk ? "yes" : "no") << "\n"
      << "max |Itot-Ie-Ii|         : "
      << maxLinearityResidual << " fC/ns\n\n"
      << "strip   phi0        Qe_exp[fC]    Qi_exp[fC]"
         "    Qtot_exp[fC]  Qtot_Garfield[fC]\n";

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const std::size_t j =
        static_cast<std::size_t>(i + halfStrips);
    const std::string label = StripLabel(i);

    const double phi0 =
        weighting.WeightingPotential(u0Cm, v0Cm, w0Cm, label);
    const double phiElectronEnd =
        weighting.WeightingPotential(ue, ve, we, label);
    const double phiIonEnd =
        weighting.WeightingPotential(ui, vi, wi, label);

    // Signal convention used by Garfield:
    // Q_signal = q [phi_w(final) - phi_w(initial)].
    const double qElectronExpected =
        -elementaryChargeFc * (phiElectronEnd - phi0);
    const double qIonExpected =
        +elementaryChargeFc * (phiIonEnd - phi0);
    const double qTotalExpected =
        qElectronExpected + qIonExpected;

    summaryOut
        << i << "," << 10. * i * stripPitchCm << ","
        << phi0 << "," << phiElectronEnd << "," << phiIonEnd << ","
        << qElectronExpected << "," << qIonExpected << ","
        << qTotalExpected << ","
        << qElectronFine[j] << "," << qElectronCoarse[j] << ","
        << qIonCoarse[j] << "," << qTotalCoarse[j] << "\n";

    std::cout << std::setw(5) << i << "   "
              << std::setw(10) << phi0 << "   "
              << std::setw(12) << qElectronExpected << "   "
              << std::setw(12) << qIonExpected << "   "
              << std::setw(13) << qTotalExpected << "   "
              << std::setw(16) << qTotalCoarse[j] << "\n";
  }

  // Wire-aware complete-cathode control. Unlike the strip approximation,
  // this weighting potential contains the actual anode-wire boundary.
  const double phiPlaneStart =
      driftField.WeightingPotential(u0Cm, v0Cm, w0Cm, planeReadoutLabel);
  const double phiPlaneElectronEnd =
      driftField.WeightingPotential(ue, ve, we, planeReadoutLabel);
  const double phiPlaneIonEnd =
      driftField.WeightingPotential(ui, vi, wi, planeReadoutLabel);

  const double qPlaneElectronExpected =
      -elementaryChargeFc * (phiPlaneElectronEnd - phiPlaneStart);
  const double qPlaneIonExpected =
      +elementaryChargeFc * (phiPlaneIonEnd - phiPlaneStart);
  const double qPlaneTotalExpected =
      qPlaneElectronExpected + qPlaneIonExpected;

  const std::string planeSummaryFile =
      outputPrefix + "_plane_control_summary.csv";
  std::ofstream planeSummaryOut(planeSummaryFile);
  planeSummaryOut
      << "label,phi_start,phi_electron_end,phi_ion_end,"
      << "q_electron_expected_fC,q_ion_expected_fC,q_total_expected_fC,"
      << "q_electron_fine_fC,q_electron_coarse_fC,"
      << "q_ion_coarse_fC,q_total_coarse_fC\n";
  planeSummaryOut
      << planeReadoutLabel << ","
      << phiPlaneStart << "," << phiPlaneElectronEnd << ","
      << phiPlaneIonEnd << ","
      << qPlaneElectronExpected << "," << qPlaneIonExpected << ","
      << qPlaneTotalExpected << ","
      << qPlaneElectronFine << "," << qPlaneElectronCoarse << ","
      << qPlaneIonCoarse << "," << qPlaneTotalCoarse << "\n";

  std::cout
      << "\nWIRE-AWARE COMPLETE-CATHODE CONTROL (" << planeReadoutLabel
      << ")\n"
      << "  phi(start)          = " << phiPlaneStart << "\n"
      << "  phi(electron end)   = " << phiPlaneElectronEnd
      << "   [should approach 0 on anode wire]\n"
      << "  phi(ion end)        = " << phiPlaneIonEnd
      << "   [should approach 1 on selected cathode]\n"
      << "  Qe expected         = " << qPlaneElectronExpected << " fC\n"
      << "  Qi expected         = " << qPlaneIonExpected << " fC\n"
      << "  Qtotal expected     = " << qPlaneTotalExpected << " fC\n"
      << "  Qtotal Garfield     = " << qPlaneTotalCoarse << " fC\n";

  if (te > electronTmaxNs) {
    std::cout
        << "\nWARNING: electron drift time exceeds the fine signal window. "
        << "Increase --electron-tmax-ns.\n";
  }
  if (ti > tMaxUs * 1000.) {
    std::cout
        << "\nWARNING: ion drift time exceeds the coarse signal window. "
        << "Increase --tmax-us.\n";
  }

  std::cout
      << "\nWrote:\n"
      << "  " << electronWaveFile << "\n"
      << "  " << ionWaveFile << "\n"
      << "  " << combinedWaveFile << "\n"
      << "  " << outputPrefix << "_electron_driftline.csv\n"
      << "  " << outputPrefix << "_ion_driftline.csv\n"
      << "  " << summaryFile << "\n"
      << "  " << planeSummaryFile << "\n"
      << "\nInterpretation target:\n"
      << "  electron = prompt component; ion = slow component;\n"
      << "  total signal must equal their linear sum.\n"
      << "===================================================\n";

  return (electronOk && ionOk &&
          electronCombinedOk && ionCombinedOk) ? 0 : 4;
}

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
  const double u0Cm =
      0.1 * ReadArg(argc, argv, "--u0-mm", 0.0);
  const double w0Cm =
      0.1 * ReadArg(argc, argv, "--w0-mm", 0.0);

  // Start the ion just outside the anode-wire surface on the selected side.
  // Default is wire radius + 5 um.
  const double defaultOffsetMm = 0.5 * 1.e1 * wireDiameterCm + 0.005;
  const double ionOffsetMm =
      ReadArg(argc, argv, "--ion-offset-mm", defaultOffsetMm);

  const double dtNs =
      ReadArg(argc, argv, "--dt-ns", 50.0);
  const double tMaxUs =
      ReadArg(argc, argv, "--tmax-us", 100.0);
  // Negative value keeps Garfield's geometry-based automatic step limit.
  // Positive values are interpreted in mm and converted to cm below.
  const double maxStepMm =
      ReadArg(argc, argv, "--max-step-mm", -1.0);
  const int signalAveragingOrder =
      ReadIntArg(argc, argv, "--signal-averaging-order", 2);
  const std::string mobilityFile =
      ReadStringArg(argc, argv, "--ion-mobility",
                    "IonMobility_Ar+_Ar.txt");
  const std::string outputPrefix =
      ReadStringArg(argc, argv, "--output-prefix",
                    "phaseB_single_ion");
  const std::string verifyStrip =
      ReadStringArg(argc, argv, "--verify-strip", "strip_0");

  if (gapMinusCm <= 0. || gapPlusCm <= 0. ||
      wirePitchCm <= 0. || wireDiameterCm <= 0. ||
      stripPitchCm <= 0. || stripWidthCm <= 0. ||
      stripWidthCm > stripPitchCm || halfStrips < 1 ||
      dtNs <= 0. || tMaxUs <= 0. || signalAveragingOrder < 1) {
    std::cerr << "Invalid Phase-B1 single-ion parameters.\n";
    return 2;
  }
  if (side != "plus" && side != "minus") {
    std::cerr << "--side must be 'plus' or 'minus'.\n";
    return 2;
  }

  const double planeY = side == "plus" ? gapPlusCm : -gapMinusCm;
  const double gapCm = std::abs(planeY);
  const double sign = side == "plus" ? +1. : -1.;
  const double v0Cm = sign * 0.1 * ionOffsetMm;

  const double wireRadiusCm = 0.5 * wireDiameterCm;
  if (std::abs(v0Cm) <= wireRadiusCm) {
    std::cerr << "Ion starting point lies inside/on the anode wire.\n"
              << "wire radius = " << 10. * wireRadiusCm
              << " mm, |v0| = " << 10. * std::abs(v0Cm) << " mm\n";
    return 2;
  }
  if (std::abs(v0Cm) >= gapCm) {
    std::cerr << "Ion starting point lies outside the selected gas gap.\n";
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

  // Garfield++ does not compute positive-ion mobility from Magboltz.
  // For this first signal test we use the distributed Ar+ in Ar table
  // as an explicit approximation. Timing is therefore provisional.
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
  // 2) Dedicated strip weighting-field component for the selected side.
  // ------------------------------------------------------------------
  Garfield::ComponentAnalyticField weighting;

  // Auxiliary voltages only make a valid analytic cell. They are not the
  // physical detector bias and do not set the weighting-electrode voltage.
  weighting.AddPlaneY(0., 1., "weighting_back");
  weighting.AddPlaneY(planeY, 0., "weighting_front");

  std::vector<std::string> labels;
  labels.reserve(2 * halfStrips + 1);
  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const double centerW = i * stripPitchCm;
    const double wMin = centerW - 0.5 * stripWidthCm;
    const double wMax = centerW + 0.5 * stripWidthCm;
    const std::string label = StripLabel(i);

    // Strip runs along local u (Garfield x), and is segmented in local w
    // (Garfield z).
    weighting.AddStripOnPlaneY('x', planeY, wMin, wMax,
                               label, gapCm);
    labels.push_back(label);
  }

  if (std::find(labels.begin(), labels.end(), verifyStrip) == labels.end()) {
    std::cerr << "--verify-strip label '" << verifyStrip
              << "' is not among the configured readout strips.\n";
    return 2;
  }

  // ------------------------------------------------------------------
  // 3) Sensor: use the physical field for transport and weighting
  //    component for Shockley-Ramo signal calculation.
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

  // ------------------------------------------------------------------
  // 4) Drift one positive ion and calculate the induced strip signals.
  // ------------------------------------------------------------------
  Garfield::DriftLineRKF ion(&sensor);
  ion.EnableSignalCalculation(true);
  ion.UseWeightingPotential(true);
  ion.SetSignalAveragingOrder(
      static_cast<std::size_t>(signalAveragingOrder));
  if (maxStepMm > 0.) {
    ion.SetMaximumStepSize(0.1 * maxStepMm);  // mm -> cm
  } else {
    ion.SetMaximumStepSize();
  }

  const bool ok = ion.DriftIon(u0Cm, v0Cm, w0Cm, 0.);

  double u1 = 0., v1 = 0., w1 = 0., t1 = 0.;
  int status = 0;
  ion.GetEndPoint(u1, v1, w1, t1, status);

  // ------------------------------------------------------------------
  // 5) Save current waveforms and integrated induced charge.
  //    Garfield current units are fC/ns, so current * dt gives fC.
  // ------------------------------------------------------------------
  const std::string waveformFile =
      outputPrefix + "_waveforms.csv";
  const std::string summaryFile =
      outputPrefix + "_summary.csv";
  const std::string driftLineFile =
      outputPrefix + "_driftline.csv";
  const std::string ramoCheckFile =
      outputPrefix + "_ramo_check.csv";

  std::ofstream waveOut(waveformFile);
  waveOut << "time_ns";
  for (const auto& label : labels) waveOut << "," << label << "_fC_per_ns";
  waveOut << "\n";

  std::vector<double> qInt(labels.size(), 0.);
  std::vector<double> peakAbs(labels.size(), 0.);

  for (int ibin = 0; ibin < nBins; ++ibin) {
    const double t = (ibin + 0.5) * dtNs;
    waveOut << t;

    for (std::size_t j = 0; j < labels.size(); ++j) {
      const double current = sensor.GetIonSignal(labels[j], ibin);
      waveOut << "," << current;
      qInt[j] += current * dtNs;
      peakAbs[j] = std::max(peakAbs[j], std::abs(current));
    }
    waveOut << "\n";
  }

  // Cache the drift line so it can be used both for the ordinary drift-line
  // CSV and for an independent local Shockley-Ramo diagnostic.
  struct DriftPoint {
    double u = 0.;
    double v = 0.;
    double w = 0.;
    double t = 0.;
  };

  const std::size_t nDriftPoints = ion.GetNumberOfDriftLinePoints();
  std::vector<DriftPoint> driftPoints;
  driftPoints.reserve(nDriftPoints);
  for (std::size_t ip = 0; ip < nDriftPoints; ++ip) {
    DriftPoint p;
    ion.GetDriftLinePoint(ip, p.u, p.v, p.w, p.t);
    driftPoints.push_back(p);
  }

  std::ofstream driftOut(driftLineFile);
  driftOut << "index,u_mm,v_mm,w_mm,time_ns,step_mm,dt_ns\n";
  for (std::size_t ip = 0; ip < driftPoints.size(); ++ip) {
    double dsMm = 0.;
    double dtns = 0.;
    if (ip > 0) {
      const double dx = driftPoints[ip].u - driftPoints[ip - 1].u;
      const double dy = driftPoints[ip].v - driftPoints[ip - 1].v;
      const double dz = driftPoints[ip].w - driftPoints[ip - 1].w;
      dsMm = 10. * std::sqrt(dx * dx + dy * dy + dz * dz);
      dtns = driftPoints[ip].t - driftPoints[ip - 1].t;
    }

    driftOut << ip << "," << 10. * driftPoints[ip].u << ","
             << 10. * driftPoints[ip].v << ","
             << 10. * driftPoints[ip].w << ","
             << driftPoints[ip].t << "," << dsMm << "," << dtns << "\n";
  }

  // ------------------------------------------------------------------
  // 6) Direct Shockley-Ramo verification on each RKF trajectory segment.
  //
  // This is independent of Sensor::GetIonSignal:
  //   - query the real drift field directly from driftField,
  //   - query the weighting field directly from weighting,
  //   - estimate v from consecutive trajectory points,
  //   - form i = -q v . Ew at the segment midpoint.
  //
  // The minus sign matches the signal-current convention used here:
  // Q_ind = q [phi_w(final) - phi_w(initial)].
  // Since Ew = -grad(phi_w), dQ_ind/dt = -q v.Ew.
  // ------------------------------------------------------------------
  constexpr double elementaryChargeFc = 1.602176634e-4;

  std::ofstream ramoOut(ramoCheckFile);
  ramoOut
      << "segment,t_mid_ns,u_mid_mm,v_mid_mm,w_mid_mm,"
      << "vx_cm_per_ns,vy_cm_per_ns,vz_cm_per_ns,speed_mm_per_us,"
      << "Ereal_x_V_per_cm,Ereal_y_V_per_cm,Ereal_z_V_per_cm,"
      << "Ereal_mag_V_per_cm,field_status,"
      << "Ew_x_per_cm,Ew_y_per_cm,Ew_z_per_cm,Ew_mag_per_cm,"
      << "phi_start,phi_end,"
      << "i_ramo_field_fC_per_ns,i_dphi_fC_per_ns\n";

  for (std::size_t ip = 1; ip < driftPoints.size(); ++ip) {
    const auto& a = driftPoints[ip - 1];
    const auto& b = driftPoints[ip];
    const double dt = b.t - a.t;
    if (dt <= 0.) continue;

    const double um = 0.5 * (a.u + b.u);
    const double vm = 0.5 * (a.v + b.v);
    const double wm = 0.5 * (a.w + b.w);
    const double tm = 0.5 * (a.t + b.t);

    const double vx = (b.u - a.u) / dt;
    const double vy = (b.v - a.v) / dt;
    const double vz = (b.w - a.w) / dt;
    const double speedMmPerUs =
        10000. * std::sqrt(vx * vx + vy * vy + vz * vz);

    double ex = 0., ey = 0., ez = 0.;
    Garfield::Medium* medium = nullptr;
    int fieldStatus = 0;
    driftField.ElectricField(um, vm, wm, ex, ey, ez, medium, fieldStatus);
    const double eMag = std::sqrt(ex * ex + ey * ey + ez * ez);

    double wx = 0., wy = 0., wz = 0.;
    weighting.WeightingField(um, vm, wm, wx, wy, wz, verifyStrip);
    const double ewMag = std::sqrt(wx * wx + wy * wy + wz * wz);

    const double phiA =
        weighting.WeightingPotential(a.u, a.v, a.w, verifyStrip);
    const double phiB =
        weighting.WeightingPotential(b.u, b.v, b.w, verifyStrip);

    // q [fC] * v [cm/ns] * Ew [1/cm] = fC/ns.
    const double iRamoField =
        -elementaryChargeFc * (vx * wx + vy * wy + vz * wz);
    const double iFromDeltaPhi =
        elementaryChargeFc * (phiB - phiA) / dt;

    ramoOut << (ip - 1) << "," << tm << ","
            << 10. * um << "," << 10. * vm << "," << 10. * wm << ","
            << vx << "," << vy << "," << vz << ","
            << speedMmPerUs << ","
            << ex << "," << ey << "," << ez << "," << eMag << ","
            << fieldStatus << ","
            << wx << "," << wy << "," << wz << "," << ewMag << ","
            << phiA << "," << phiB << ","
            << iRamoField << "," << iFromDeltaPhi << "\n";
  }

  std::ofstream sumOut(summaryFile);
  sumOut << "strip,center_w_mm,integrated_ion_signal_fC,"
            "peak_abs_current_fC_per_ns,phi_start,phi_end,"
            "ramo_delta_phi_fC\n";

  std::cout << std::fixed << std::setprecision(7);
  std::cout << "\n=== PHASE B1a: SINGLE-ION SHOCKLEY-RAMO SIGNAL ===\n"
            << "side                 : " << side << "\n"
            << "gas                  : Ar/CO2 70:30\n"
            << "ion mobility model   : " << mobilityFile
            << "  [timing approximation]\n"
            << "wire pitch / diameter: " << 10. * wirePitchCm << " mm / "
            << 1.e4 * wireDiameterCm << " um\n"
            << "anode voltage        : " << hv << " V\n"
            << "strip pitch / width  : " << 10. * stripPitchCm << " / "
            << 10. * stripWidthCm << " mm\n"
            << "ion start (u,v,w)    : (" << 10. * u0Cm << ", "
            << 10. * v0Cm << ", " << 10. * w0Cm << ") mm\n"
            << "ion end   (u,v,w)    : (" << 10. * u1 << ", "
            << 10. * v1 << ", " << 10. * w1 << ") mm\n"
            << "drift status         : " << status
            << "   success=" << (ok ? "yes" : "no") << "\n"
            << "ion drift time       : " << t1 << " ns\n"
            << "signal window        : 0 .. " << tMaxUs
            << " us, dt=" << dtNs << " ns\n"
            << "RKF max step         : ";
  if (maxStepMm > 0.) {
    std::cout << maxStepMm << " mm (explicit)\n";
  } else {
    std::cout << "Garfield automatic geometry limit\n";
  }
  std::cout << "signal avg. order    : " << signalAveragingOrder << "\n"
            << "drift-line points    : " << nDriftPoints << "\n"
            << "direct Ramo check    : " << verifyStrip << "\n\n"
            << "strip   center_w[mm]   Q_signal[fC]   |I|_peak[fC/ns]"
               "   q*dphi[fC]\n";

  for (int i = -halfStrips; i <= halfStrips; ++i) {
    const std::size_t j = static_cast<std::size_t>(i + halfStrips);
    const std::string label = StripLabel(i);

    const double phiStart =
        weighting.WeightingPotential(u0Cm, v0Cm, w0Cm, label);
    const double phiEnd =
        weighting.WeightingPotential(u1, v1, w1, label);

    // For one positive ion, q = +e.
    const double qRamo =
        elementaryChargeFc * (phiEnd - phiStart);

    sumOut << i << "," << 10. * i * stripPitchCm << ","
           << qInt[j] << "," << peakAbs[j] << ","
           << phiStart << "," << phiEnd << ","
           << qRamo << "\n";

    std::cout << std::setw(5) << i << "   "
              << std::setw(12) << 10. * i * stripPitchCm << "   "
              << std::setw(12) << qInt[j] << "   "
              << std::setw(16) << peakAbs[j] << "   "
              << std::setw(12) << qRamo << "\n";
  }

  if (t1 > tMaxUs * 1000.) {
    std::cout << "\nWARNING: ion drift time exceeds the signal window. "
                 "Increase --tmax-us before interpreting integrated charge.\n";
  }

  std::cout << "\nWrote " << waveformFile
            << ", " << summaryFile
            << ", " << driftLineFile
            << " and " << ramoCheckFile << "\n"
            << "Garfield signal unit: fC/ns; integrated values above are fC.\n"
            << "Direct Ramo check uses Garfield WeightingField independently "
               "of GetIonSignal.\n"
            << "===========================================================\n";

  return ok ? 0 : 4;
}

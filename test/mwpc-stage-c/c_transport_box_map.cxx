#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"

#include "../mwpc-avalanche/MWPCCoordinateFrame.h"

namespace {

constexpr double Pi = 3.14159265358979323846;

double ReadArg(const int argc, char** argv, const std::string& key,
               const double defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return std::atof(argv[i + 1]);
    if (arg.rfind(prefix, 0) == 0) {
      return std::atof(arg.substr(prefix.size()).c_str());
    }
  }
  return defaultValue;
}

int ReadIntArg(const int argc, char** argv, const std::string& key,
               const int defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return std::atoi(argv[i + 1]);
    if (arg.rfind(prefix, 0) == 0) {
      return std::atoi(arg.substr(prefix.size()).c_str());
    }
  }
  return defaultValue;
}

std::string ReadStringArg(const int argc, char** argv,
                          const std::string& key,
                          const std::string& defaultValue) {
  const std::string prefix = key + "=";
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == key && i + 1 < argc) return argv[i + 1];
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
  }
  return defaultValue;
}

double Clamp(const double x, const double lo, const double hi) {
  return std::max(lo, std::min(hi, x));
}

struct LastStep {
  bool valid = false;
  double u = 0.;
  double v = 0.;
  double w = 0.;
  double t = 0.;
  double energy = 0.;
  double du = 0.;
  double dv = 0.;
  double dw = 0.;
};

LastStep gLastStep;

void RecordStep(const double x, const double y, const double z,
                const double t, const double e,
                const double dx, const double dy, const double dz,
                const bool hole, const int /*band*/) {
  if (hole) return;
  gLastStep.valid = true;
  gLastStep.u = x;
  gLastStep.v = y;
  gLastStep.w = z;
  gLastStep.t = t;
  gLastStep.energy = e;
  gLastStep.du = dx;
  gLastStep.dv = dy;
  gLastStep.dw = dz;
}

}  // namespace

int main(int argc, char** argv) {
  // ------------------------------------------------------------------
  // Stage-C transport response map.
  //
  // Divide one wire cell in (u,v) into boxes. From every box launch many
  // single electrons and record:
  //   * destination anode wire,
  //   * drift time,
  //   * hit position around the wire (phi_hit),
  //   * microscopic direction immediately before collection,
  //   * incidence angle with respect to the inward surface normal.
  //
  // The map deliberately stops before avalanche gain and cathode readout.
  // It is intended as a reusable transport kernel for later GEANT cluster
  // distributions / digitisation.
  // ------------------------------------------------------------------

  const double gapMinusMm =
      ReadArg(argc, argv, "--gap-minus-mm", 2.0);
  const double gapPlusMm =
      ReadArg(argc, argv, "--gap-plus-mm", 4.0);
  const double wirePitchMm =
      ReadArg(argc, argv, "--wire-pitch-mm", 4.0);
  const double wireDiameterUm =
      ReadArg(argc, argv, "--wire-diam-um", 30.0);
  const double hv =
      ReadArg(argc, argv, "--hv", 1800.0);
  const double bTesla =
      ReadArg(argc, argv, "--b-tesla", 0.0);
  const int halfWires =
      ReadIntArg(argc, argv, "--half-wires", 4);

  const double uMinMm =
      ReadArg(argc, argv, "--u-min-mm", -2.0);
  const double uMaxMm =
      ReadArg(argc, argv, "--u-max-mm", +2.0);
  const double vMinMm =
      ReadArg(argc, argv, "--v-min-mm", -gapMinusMm);
  const double vMaxMm =
      ReadArg(argc, argv, "--v-max-mm", +gapPlusMm);
  const double duMm =
      ReadArg(argc, argv, "--du-mm", 0.2);
  const double dvMm =
      ReadArg(argc, argv, "--dv-mm", 0.2);

  const int electronsPerBox =
      ReadIntArg(argc, argv, "--electrons-per-box", 50);
  const int baseSeed =
      ReadIntArg(argc, argv, "--base-seed", 170001);
  const double startEnergyEv =
      ReadArg(argc, argv, "--start-energy-ev", 0.1);
  const double w0Mm =
      ReadArg(argc, argv, "--w0-mm", 0.0);
  const double wHalfSpanMm =
      ReadArg(argc, argv, "--w-half-span-mm", 20.0);

  // Optional u-bin slicing is useful for resumable shell runs.
  // Indices are inclusive; -1 means "all".
  const int iuMinRequested =
      ReadIntArg(argc, argv, "--iu-min", -1);
  const int iuMaxRequested =
      ReadIntArg(argc, argv, "--iu-max", -1);

  const std::string output =
      ReadStringArg(argc, argv, "--output", "c_transport_box_electrons.csv");

  if (gapMinusMm <= 0. || gapPlusMm <= 0. ||
      wirePitchMm <= 0. || wireDiameterUm <= 0. ||
      halfWires < 2 || uMaxMm <= uMinMm || vMaxMm <= vMinMm ||
      duMm <= 0. || dvMm <= 0. || electronsPerBox < 1 ||
      baseSeed < 0 || startEnergyEv <= 0. || wHalfSpanMm <= 0.) {
    std::cerr << "Invalid Stage-C transport-map parameters.\n";
    return 2;
  }

  const double nUExact = (uMaxMm - uMinMm) / duMm;
  const double nVExact = (vMaxMm - vMinMm) / dvMm;
  const int nU = static_cast<int>(std::llround(nUExact));
  const int nV = static_cast<int>(std::llround(nVExact));
  if (std::abs(nUExact - nU) > 1.e-9 ||
      std::abs(nVExact - nV) > 1.e-9 ||
      nU < 1 || nV < 1) {
    std::cerr
        << "Grid ranges must be integer multiples of du/dv.\n";
    return 2;
  }

  int iuFirst = 0;
  int iuLast = nU - 1;
  if (iuMinRequested >= 0) iuFirst = iuMinRequested;
  if (iuMaxRequested >= 0) iuLast = iuMaxRequested;
  if (iuFirst < 0 || iuLast >= nU || iuLast < iuFirst) {
    std::cerr << "Invalid iu range. nU=" << nU << "\n";
    return 2;
  }

  const auto frame = na6p::mwpc::MakeNominalNA60Frame();
  const auto localFromGlobalX = frame.GlobalToLocalVector({1., 0., 0.});
  const auto localFromGlobalY = frame.GlobalToLocalVector({0., 1., 0.});
  const auto localFromGlobalZ = frame.GlobalToLocalVector({0., 0., 1.});
  constexpr double frameTol = 1.e-12;
  if (std::abs(localFromGlobalX.w - 1.) > frameTol ||
      std::abs(localFromGlobalY.u - 1.) > frameTol ||
      std::abs(localFromGlobalZ.v - 1.) > frameTol) {
    std::cerr << "Nominal ChamberFrame no longer matches Stage-C convention.\n";
    return 3;
  }

  const double gapMinusCm = 0.1 * gapMinusMm;
  const double gapPlusCm = 0.1 * gapPlusMm;
  const double wirePitchCm = 0.1 * wirePitchMm;
  const double wireDiameterCm = 1.e-4 * wireDiameterUm;
  const double wireRadiusCm = 0.5 * wireDiameterCm;
  const double w0Cm = 0.1 * w0Mm;
  const double wHalfSpanCm = 0.1 * wHalfSpanMm;

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);
  for (int i = -halfWires; i <= halfWires; ++i) {
    field.AddWire(i * wirePitchCm, 0., wireDiameterCm, hv, "anode");
  }
  field.AddPlaneY(-gapMinusCm, 0., "cathode_minus");
  field.AddPlaneY(+gapPlusCm, 0., "cathode_plus");
  field.SetMagneticField(bTesla, 0., 0.);

  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  const double uExtent =
      (halfWires + 0.5) * wirePitchCm;
  sensor.SetArea(
      -uExtent, -gapMinusCm, -wHalfSpanCm,
      +uExtent, +gapPlusCm, +wHalfSpanCm);

  // Optional NEW pure microscopic DRIFT-LINE collection mode for the
  // PI's P(u(t)),P(v(t)),P(w(t)) plots conditional on true native muon angle.
  // Uniform response-map output and its historical behaviour are untouched.
  const std::string traceStarts =
      ReadStringArg(argc, argv, "--trace-starts", "");
  if (!traceStarts.empty()) {
    const std::string traceOutput = ReadStringArg(
        argc, argv, "--trace-output", "c_centered_native_drift_lines.csv");
    const int traceCollisionSteps = ReadIntArg(
        argc, argv, "--trace-collision-steps", 20);
    if (std::abs(gapMinusMm - 3.) > 1.e-12 ||
        std::abs(gapPlusMm - 3.) > 1.e-12 ||
        std::abs(wirePitchMm - 4.) > 1.e-12 ||
        std::abs(bTesla) > 1.e-12 ||
        traceCollisionSteps < 1) {
      std::cerr << "Native track drift snapshots demand centred 3+3mm, "
                   "4mm pitch, B=0, positive collision sample stride.\n";
      return 21;
    }
    std::ifstream source(traceStarts);
    if (!source) {
      std::cerr << "Cannot read native Heed starting electron CSV: "
                << traceStarts << "\n";
      return 22;
    }
    std::ofstream trace(traceOutput);
    if (!trace) {
      std::cerr << "Cannot open drift-line output CSV: "
                << traceOutput << "\n";
      return 23;
    }
    std::string header;
    std::getline(source, header);
    if (header != "electron_id,u_cell_mm,v_birth_mm,w_start_mm,seed") {
      std::cerr << "Incorrect native electron start CSV header. "
                   "Expected numeric ID,u_cell,v,w,seed.\n";
      return 24;
    }
    trace << std::setprecision(12)
          << "electron_id,point_index,t_rel_birth_ns,"
             "u_cell_mm,v_mm,w_mm,end_status,number_points\n";
    std::string line;
    long long nElectronStarts = 0;
    long long nMicroPoints = 0;
    while (std::getline(source, line)) {
      if (line.empty()) continue;
      std::istringstream row(line);
      std::string part;
      std::vector<std::string> v;
      while (std::getline(row, part, ',')) v.push_back(part);
      if (v.size() != 5) {
        std::cerr << "Bad native Heed seed row: " << line << "\n";
        return 25;
      }
      try {
        const long long id = std::stoll(v[0]);
        const double u0mm = std::stod(v[1]);
        const double v0mm = std::stod(v[2]);
        const double w0mm = std::stod(v[3]);
        const long long seed = std::stoll(v[4]);
        if (id < 0 || seed < 0 || seed > 2147483647LL ||
            !std::isfinite(u0mm) || !std::isfinite(v0mm) ||
            !std::isfinite(w0mm) ||
            std::abs(u0mm) > 2.001 ||
            std::abs(v0mm) > 3.001 ||
            std::abs(w0mm) > 1.e-5) {
          std::cerr << "Native drift starts must be within centered "
                       "wire cell and gas half gap, with w0=0.\n";
          return 26;
        }
        Garfield::RandomEngineRoot engine(
            static_cast<unsigned int>(seed));
        Garfield::Random::SetEngine(engine);
        Garfield::AvalancheMicroscopic drift(&sensor);
        drift.EnableSignalCalculation(false);
        drift.EnableDriftLines();
        drift.SetCollisionSteps(traceCollisionSteps);
        drift.DriftElectron(0.1*u0mm,0.1*v0mm,0.,0.,
                            startEnergyEv,0.,0.,0.);
        const auto endCount=drift.GetNumberOfElectronEndpoints();
        int status=-999;
        if (endCount) {
          double uStart=0.,vStart=0.,wStart=0.,tStart=0.,eStart=0.;
          double uEnd=0.,vEnd=0.,wEnd=0.,tEnd=0.,eEnd=0.;
          drift.GetElectronEndpoint(0,uStart,vStart,wStart,tStart,eStart,
                  uEnd,vEnd,wEnd,tEnd,eEnd,status);
        }
        if (!endCount) {
          trace << id << ",0,0," << u0mm << "," << v0mm
                << ",0,-999,1\n";
          ++nMicroPoints;
        } else {
          const auto count=drift.GetNumberOfElectronDriftLinePoints(0);
          for (std::size_t k=0;k<count;++k) {
            double uu=0.,vv=0.,ww=0.,tt=0.;
            drift.GetElectronDriftLinePoint(uu,vv,ww,tt,k,0);
            trace << id << "," << k << "," << tt << ","
                  << 10.*uu << "," << 10.*vv << "," << 10.*ww
                  << "," << status << "," << count << "\n";
            ++nMicroPoints;
          }
        }
        ++nElectronStarts;
      } catch (const std::exception& e) {
        std::cerr << "Failed native Garfield electron seed/trace parse: "
                  << e.what() << "\n";
        return 27;
      }
    }
    if (!nElectronStarts) {
      std::cerr << "No native Heed electrons to drift.\n";
      return 28;
    }
    std::cout << "PASS: centered native-angle microscopic Garfield "
              << nElectronStarts << " electrons, "
              << nMicroPoints << " REAL drift-line points -> "
              << traceOutput << "\n";
    return 0;
  }

  std::ofstream out(output);
  if (!out) {
    std::cerr << "Could not open output file: " << output << "\n";
    return 4;
  }

  out << std::setprecision(12);
  out
      << "iu,iv,box_id,u_center_mm,v_center_mm,"
      << "u_start_mm,v_start_mm,w_start_mm,"
      << "sample,random_seed,drift_ok,collected,status,"
      << "wire_index,wire_u_mm,"
      << "u_end_mm,v_end_mm,w_end_mm,drift_time_ns,"
      << "phi_hit_deg,cos_phi_hit,sin_phi_hit,hit_side_u,"
      << "incidence_normal_deg,"
      << "last_step_u_mm,last_step_v_mm,last_step_w_mm,"
      << "last_step_energy_eV,last_dir_u,last_dir_v,last_dir_w\n";

  std::cout << std::fixed << std::setprecision(4)
            << "\n=== STAGE C: ELECTRON TRANSPORT BOX MAP ===\n"
            << "cell u range              : "
            << uMinMm << " .. " << uMaxMm << " mm\n"
            << "gas v range               : "
            << vMinMm << " .. " << vMaxMm << " mm\n"
            << "grid                       : "
            << nU << " x " << nV
            << " boxes (du=" << duMm
            << " mm, dv=" << dvMm << " mm)\n"
            << "u bins in this invocation : "
            << iuFirst << " .. " << iuLast << "\n"
            << "electrons / box           : "
            << electronsPerBox << "\n"
            << "wire pitch / diameter     : "
            << wirePitchMm << " mm / "
            << wireDiameterUm << " um\n"
            << "gap (-,+)                 : "
            << gapMinusMm << ", " << gapPlusMm << " mm\n"
            << "B(local u)                : "
            << bTesla << " T\n"
            << "start energy              : "
            << startEnergyEv << " eV\n\n";

  long long nLaunched = 0;
  long long nCollected = 0;

  for (int iu = iuFirst; iu <= iuLast; ++iu) {
    const double uLoMm = uMinMm + iu * duMm;
    const double uHiMm = uLoMm + duMm;
    const double uCenterMm = 0.5 * (uLoMm + uHiMm);

    std::cout << "u-bin " << iu << "/" << (nU - 1)
              << "  center=" << uCenterMm << " mm\n";

    for (int iv = 0; iv < nV; ++iv) {
      const double vLoMm = vMinMm + iv * dvMm;
      const double vHiMm = vLoMm + dvMm;
      const double vCenterMm = 0.5 * (vLoMm + vHiMm);
      const int boxId = iu * nV + iv;

      for (int isample = 0; isample < electronsPerBox; ++isample) {
        const int randomSeed =
            baseSeed + boxId * electronsPerBox + isample;

        Garfield::RandomEngineRoot randomEngine(
            static_cast<unsigned int>(randomSeed));
        Garfield::Random::SetEngine(randomEngine);

        const double ru = Garfield::RndmUniform();
        const double rv = Garfield::RndmUniform();
        double uStartMm = uLoMm + ru * duMm;
        double vStartMm = vLoMm + rv * dvMm;

        auto insideAnyWire = [&](const double uMm, const double vMm) {
          const int iw =
              static_cast<int>(std::lround(uMm / wirePitchMm));
          if (std::abs(iw) > halfWires) return false;
          const double duWireMm = uMm - iw * wirePitchMm;
          const double rMm = std::hypot(duWireMm, vMm);
          return rMm <= 0.5e-3 * wireDiameterUm;
        };

        if (insideAnyWire(uStartMm, vStartMm)) {
          uStartMm = uCenterMm;
          vStartMm = vCenterMm;
        }
        if (insideAnyWire(uStartMm, vStartMm)) continue;

        const double uStartCm = 0.1 * uStartMm;
        const double vStartCm = 0.1 * vStartMm;

        gLastStep = LastStep{};

        Garfield::AvalancheMicroscopic drift(&sensor);
        drift.EnableSignalCalculation(false);
        drift.SetUserHandleStep(RecordStep);

        const bool driftOk =
            drift.DriftElectron(
                uStartCm, vStartCm, w0Cm,
                0., startEnergyEv, 0., 0., 0.);

        double u0 = 0., v0 = 0., w0 = 0., t0 = 0., e0 = 0.;
        double u1 = std::numeric_limits<double>::quiet_NaN();
        double v1 = std::numeric_limits<double>::quiet_NaN();
        double w1 = std::numeric_limits<double>::quiet_NaN();
        double t1 = std::numeric_limits<double>::quiet_NaN();
        double e1 = std::numeric_limits<double>::quiet_NaN();
        int status = 0;

        const std::size_t nEndpoints =
            drift.GetNumberOfElectronEndpoints();
        if (nEndpoints > 0) {
          drift.GetElectronEndpoint(
              0,
              u0, v0, w0, t0, e0,
              u1, v1, w1, t1, e1,
              status);
        }

        bool collected = false;
        int wireIndex = 999999;
        double wireUMm =
            std::numeric_limits<double>::quiet_NaN();
        double phiHitDeg =
            std::numeric_limits<double>::quiet_NaN();
        double cosPhi =
            std::numeric_limits<double>::quiet_NaN();
        double sinPhi =
            std::numeric_limits<double>::quiet_NaN();
        int hitSideU = 0;
        double incidenceNormalDeg =
            std::numeric_limits<double>::quiet_NaN();

        if (nEndpoints > 0 &&
            std::isfinite(u1) && std::isfinite(v1)) {
          wireIndex =
              static_cast<int>(std::lround(u1 / wirePitchCm));
          if (std::abs(wireIndex) <= halfWires) {
            const double wireUCm = wireIndex * wirePitchCm;
            const double radialU = u1 - wireUCm;
            const double radialV = v1;
            const double rho = std::hypot(radialU, radialV);

            if (rho < 1.5 * wireRadiusCm && status != -7) {
              collected = true;
              wireUMm = 10. * wireUCm;
              phiHitDeg =
                  std::atan2(radialV, radialU) * 180. / Pi;
              if (rho > 0.) {
                cosPhi = radialU / rho;
                sinPhi = radialV / rho;
              }
              if (radialU > 0.) {
                hitSideU = +1;
              } else if (radialU < 0.) {
                hitSideU = -1;
              }

              if (gLastStep.valid && rho > 0.) {
                const double dNorm =
                    std::sqrt(gLastStep.du * gLastStep.du +
                              gLastStep.dv * gLastStep.dv +
                              gLastStep.dw * gLastStep.dw);
                if (dNorm > 0.) {
                  const double inwardU = -radialU / rho;
                  const double inwardV = -radialV / rho;
                  const double c =
                      Clamp(
                          (gLastStep.du * inwardU +
                           gLastStep.dv * inwardV) / dNorm,
                          -1., 1.);
                  incidenceNormalDeg =
                      std::acos(c) * 180. / Pi;
                }
              }
            }
          }
        }

        ++nLaunched;
        if (collected) ++nCollected;

        out
            << iu << "," << iv << "," << boxId << ","
            << uCenterMm << "," << vCenterMm << ","
            << uStartMm << "," << vStartMm << "," << w0Mm << ","
            << isample << "," << randomSeed << ","
            << (driftOk ? 1 : 0) << ","
            << (collected ? 1 : 0) << ","
            << status << ","
            << wireIndex << "," << wireUMm << ","
            << (10. * u1) << "," << (10. * v1) << ","
            << (10. * w1) << "," << (t1 - t0) << ","
            << phiHitDeg << "," << cosPhi << "," << sinPhi << ","
            << hitSideU << ","
            << incidenceNormalDeg << ","
            << (gLastStep.valid ? 10. * gLastStep.u :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? 10. * gLastStep.v :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? 10. * gLastStep.w :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? gLastStep.energy :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? gLastStep.du :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? gLastStep.dv :
                std::numeric_limits<double>::quiet_NaN())
            << ","
            << (gLastStep.valid ? gLastStep.dw :
                std::numeric_limits<double>::quiet_NaN())
            << "\n";
      }
    }
  }

  std::cout
      << "\nWrote: " << output << "\n"
      << "launched  : " << nLaunched << "\n"
      << "collected : " << nCollected << "\n"
      << "fraction  : "
      << (nLaunched > 0
              ? static_cast<double>(nCollected) / nLaunched
              : 0.)
      << "\n"
      << "\nMODEL BOUNDARY:\n"
      << "  - one-wire-cell (u,v) transport kernel: yes\n"
      << "  - single-electron microscopic drift: yes\n"
      << "  - destination wire / drift time: yes\n"
      << "  - wire hit position and incidence angle: yes\n"
      << "  - avalanche multiplication: intentionally no\n"
      << "  - GEANT cluster weights: not applied here\n"
      << "  - cathode-strip weighting/electronics: no\n";

  return 0;
}

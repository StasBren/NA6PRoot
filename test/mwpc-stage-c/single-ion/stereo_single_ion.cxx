// NA60+ P3 noFW single-positive-ion cathode-strip benchmark.
// Transport: Garfield++ ComponentAnalyticField + DriftLineRKF.
// Weighting: existing 3D Gmsh/Elmer solution imported with ComponentElmer.
// No Geant4/Heed simulation, avalanche gain, electronics or Python analysis.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {
namespace fs = std::filesystem;
constexpr double kElementaryChargeFc = 1.602176634e-4;
struct Opt {
  fs::path mapDir, outDir, mobility;
  double hv = 1800., gapMinus = 3., gapPlus = 3.;
  double pitch = 3., wireDiameterUm = 30.;
  int halfWires = 2;
  double u0 = 0., v0 = 0.04, w0 = 0.2, t0 = 0.;
  double tmaxUs = 250.;
  std::vector<double> dts{10., 5., 2.};
  std::vector<double> steps{0.04, 0.02};
};
std::string Arg(int argc, char** argv, const std::string& k,
                const std::string& fallback = "") {
  for (int j = 1; j + 1 < argc; ++j) {
    if (k == argv[j]) return argv[j+1];
  }
  return fallback;
}
double Number(int argc, char** argv, const std::string& k, double v) {
  const auto s = Arg(argc, argv, k);
  return s.empty() ? v : std::stod(s);
}
std::vector<double> Numbers(const std::string& s) {
  std::vector<double> result;
  std::size_t first = 0;
  while (first < s.size()) {
    const auto last = s.find(',', first);
    const auto part = s.substr(first, last == std::string::npos ?
                                      last : last-first);
    if (part.empty()) throw std::runtime_error("Empty numeric-list entry");
    result.push_back(std::stod(part));
    if (last == std::string::npos) break;
    first = last+1;
  }
  if (result.empty()) throw std::runtime_error("Empty numeric list");
  for (const double x : result)
    if (!std::isfinite(x) || x <= 0.)
      throw std::runtime_error("Numeric grid must be positive and finite");
  return result;
}
void RequireFile(const fs::path& p) {
  if (!fs::exists(p) || !fs::is_regular_file(p))
    throw std::runtime_error("Missing required file: "+p.string());
}
double Phi(Garfield::ComponentElmer& map, double u, double v,
           double w, const std::string& label) {
  const double p = map.WeightingPotential(u, v, w, label);
  if (!std::isfinite(p) || p < -0.1 || p > 1.1)
    throw std::runtime_error("Invalid FEM weighting potential for "+label);
  return p;
}
struct Result {
  double dt, step, ue, ve, we, te;
  int status;
  bool ok, complete, arrived;
  double qPlus, qMinus, ramoPlus, ramoMinus;
  double errPlus, errMinus;
};
double Tolerance(double qRamo) {
  // Absolute floor of 0.005 elementary charges; 3% otherwise.
  return std::max(0.03 * std::abs(qRamo), 0.005*kElementaryChargeFc);
}
bool CloseCharge(double measured, double predicted) {
  return std::abs(measured-predicted) <= Tolerance(predicted);
}
}  // namespace

int main(int argc, char** argv) {
  try {
    Opt o;
    o.mapDir = Arg(argc, argv, "--map-dir");
    o.outDir = Arg(argc, argv, "--out-dir");
    o.mobility = Arg(argc, argv, "--mobility");
    o.hv = Number(argc, argv, "--hv", o.hv);
    o.gapMinus = Number(argc, argv, "--gap-minus-mm", o.gapMinus);
    o.gapPlus = Number(argc, argv, "--gap-plus-mm", o.gapPlus);
    o.pitch = Number(argc, argv, "--wire-pitch-mm", o.pitch);
    o.wireDiameterUm = Number(argc, argv, "--wire-diameter-um", o.wireDiameterUm);
    o.halfWires = static_cast<int>(Number(argc, argv, "--half-wires", o.halfWires));
    o.u0 = Number(argc, argv, "--u-mm", o.u0);
    o.v0 = Number(argc, argv, "--v-mm", o.v0);
    o.w0 = Number(argc, argv, "--w-mm", o.w0);
    o.t0 = Number(argc, argv, "--t0-ns", o.t0);
    o.tmaxUs = Number(argc, argv, "--tmax-us", o.tmaxUs);
    const auto dtArg = Arg(argc, argv, "--dt-ns");
    const auto stepArg = Arg(argc, argv, "--steps-mm");
    if (!dtArg.empty()) o.dts = Numbers(dtArg);
    if (!stepArg.empty()) o.steps = Numbers(stepArg);
    if (o.mapDir.empty() || o.outDir.empty() || o.mobility.empty())
      throw std::runtime_error("Require --map-dir, --out-dir and --mobility");
    if (!(o.hv > 0 && o.gapMinus > 0 && o.gapPlus > 0 &&
          o.pitch > 0 && o.wireDiameterUm > 0 && o.halfWires >= 1 &&
          o.tmaxUs > 0 && o.tmaxUs <= 5000 && o.t0 >= 0))
      throw std::runtime_error("Invalid geometry, voltage or time range");
    for (auto p : {o.mobility, o.mapDir/"chamber"/"mesh.header",
                   o.mapDir/"chamber"/"mesh.nodes",
                   o.mapDir/"chamber"/"mesh.elements",
                   o.mapDir/"dielectrics.dat",
                   o.mapDir/"chamber"/"drift.result",
                   o.mapDir/"chamber"/"weight_minus_0.result",
                   o.mapDir/"chamber"/"weight_plus_0.result"}) RequireFile(p);
    // The runner provides a fresh isolated directory. Refuse to overwrite it.
    if (!fs::is_directory(o.outDir))
      throw std::runtime_error("Output directory must exist: "+o.outDir.string());
    if (fs::exists(o.outDir/"single_ion_summary.csv"))
      throw std::runtime_error("Refusing to overwrite single-ion results");

    const double wireR = 0.5e-4 * o.wireDiameterUm;  // um -> cm
    const double pitchCm = 0.1*o.pitch;
    const double uCm = 0.1*o.u0, vCm = 0.1*o.v0, wCm = 0.1*o.w0;
    const double nearestU = std::round(uCm/pitchCm)*pitchCm;
    if (std::hypot(uCm-nearestU, vCm) <= wireR ||
        vCm <= -0.1*o.gapMinus || vCm >= 0.1*o.gapPlus)
      throw std::runtime_error("Ion birth lies inside wire/cathode metal");
    if (std::abs(o.u0) >= (o.halfWires + 0.5)*o.pitch)
      throw std::runtime_error("Start outside analytic wire patch");

    Garfield::MediumMagboltz gas;
    gas.SetComposition("ar", 70., "co2", 30.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.Initialise(false);
    if (!gas.LoadIonMobility(o.mobility.string()))
      throw std::runtime_error("Could not load requested ion mobility");

    Garfield::ComponentAnalyticField realField;
    realField.SetMedium(&gas);
    for (int k = -o.halfWires; k <= o.halfWires; ++k) {
      realField.AddWire(k*pitchCm, 0., 2.*wireR, o.hv, "anode");
    }
    realField.AddPlaneY(-0.1*o.gapMinus, 0., "cathode_minus");
    realField.AddPlaneY(+0.1*o.gapPlus, 0., "cathode_plus");
    realField.SetMagneticField(0., 0., 0.);

    const auto chamber = o.mapDir/"chamber";
    Garfield::ComponentElmer weighting(
        (chamber/"mesh.header").string(), (chamber/"mesh.elements").string(),
        (chamber/"mesh.nodes").string(), (o.mapDir/"dielectrics.dat").string(),
        (chamber/"drift.result").string(), "mm");
    const std::string minus = "minus_strip_0", plus = "plus_strip_0";
    if (!weighting.SetWeightingField((chamber/"weight_minus_0.result").string(), minus) ||
        !weighting.SetWeightingField((chamber/"weight_plus_0.result").string(), plus))
      throw std::runtime_error("Cannot load stereo cathode weighting maps");

    Garfield::Sensor sensor;
    sensor.AddComponent(&realField); // Real drift from Garfield analytic MWPC.
    sensor.AddElectrode(&weighting, minus); // 3D FEM only for signals.
    sensor.AddElectrode(&weighting, plus);
    sensor.SetArea(-(o.halfWires+0.5)*pitchCm, -0.1*o.gapMinus,
                    -1., (o.halfWires+0.5)*pitchCm, 0.1*o.gapPlus, 1.);

    const double pm0 = Phi(weighting,uCm,vCm,wCm,minus);
    const double pp0 = Phi(weighting,uCm,vCm,wCm,plus);
    // Electrode checks are near (not exactly on) a boundary: the map
    // interpolator is defined in the gas, not in the copper.
    const double nearPlus = Phi(weighting,0.,0.1*(o.gapPlus-0.02),0.,plus);
    const double nearMinus = Phi(weighting,0.,-0.1*(o.gapMinus-0.02),0.,minus);
    const double farPlus = Phi(weighting,0.,-0.1*(o.gapMinus-0.02),0.,plus);
    const double farMinus = Phi(weighting,0.,0.1*(o.gapPlus-0.02),0.,minus);
    const bool fieldGate = nearPlus > 0.6 && nearMinus > 0.6 &&
                           std::abs(farPlus) < 0.1 && std::abs(farMinus) < 0.1;
    if (!fieldGate) throw std::runtime_error(
        "3D FEM map failed elementary strip/opposite-cathode boundary gate");

    std::ofstream summary(o.outDir/"single_ion_summary.csv");
    if (!summary) throw std::runtime_error("Cannot write summary");
    summary << std::setprecision(15)
       << "step_index,dt_index,dt_ns,max_step_mm,"
          "u_start_mm,v_start_mm,w_start_mm,t_start_ns,"
          "u_end_mm,v_end_mm,w_end_mm,t_end_ns,status,drift_ok,"
          "complete_window,cathode_arrival,"
          "minus_phi_start,minus_phi_end,plus_phi_start,plus_phi_end,"
          "minus_Q_sensor_fC,minus_Q_ramo_fC,minus_abs_error_fC,"
          "plus_Q_sensor_fC,plus_Q_ramo_fC,plus_abs_error_fC\n";

    bool gates = true;
    std::vector<Result> all;
    for (std::size_t js = 0; js < o.steps.size(); ++js) {
      for (std::size_t jd = 0; jd < o.dts.size(); ++jd) {
        const double dt = o.dts[jd], step = o.steps[js];
        const double nExact = std::ceil((1000.*o.tmaxUs-o.t0)/dt);
        if (nExact < 2 || nExact > 250000)
          throw std::runtime_error("Time window requires 2..250000 bins");
        const int nBins = static_cast<int>(nExact);
        sensor.SetTimeWindow(0., dt, nBins);
        sensor.ClearSignal();

        Garfield::DriftLineRKF ion(&sensor);
        ion.EnableSignalCalculation(true);
        ion.UseWeightingPotential(true);
        ion.EnableDiffusion(false); // Deterministic numerical-convergence test.
        ion.SetSignalAveragingOrder(2);
        ion.SetMaximumStepSize(0.1*step); // step mm -> cm

        const bool ok = ion.DriftIon(uCm, vCm, wCm, o.t0);
        double ue=0.,ve=0.,we=0.,te=0.;
        int status=-999;
        ion.GetEndPoint(ue,ve,we,te,status);
        const bool complete = te <= dt*nBins && std::isfinite(te) && te>o.t0;
        const bool arrived = std::abs(std::abs(ve)-0.1*(ve>0?o.gapPlus:o.gapMinus))
                                  < 0.01;  // within 0.1 mm of cathode
        const double pm1 = Phi(weighting,ue,ve,we,minus);
        const double pp1 = Phi(weighting,ue,ve,we,plus);
        const double qmRamo = kElementaryChargeFc*(pm1-pm0);
        const double qpRamo = kElementaryChargeFc*(pp1-pp0);

        std::ofstream wave(o.outDir/("wave_"+std::to_string(js)+"_"+
                                      std::to_string(jd)+".csv"));
        if (!wave) throw std::runtime_error("Cannot write waveform");
        wave << std::setprecision(15)
             << "time_ns,minus_ion_fC_per_ns,plus_ion_fC_per_ns,"
                "minus_cumulative_fC,plus_cumulative_fC\n";
        double qm=0.,qp=0.;
        for (int k=0;k<nBins;++k) {
          const double im=sensor.GetIonSignal(minus,k);
          const double ip=sensor.GetIonSignal(plus,k);
          if (!std::isfinite(im) || !std::isfinite(ip))
            throw std::runtime_error("Non-finite Garfield strip current");
          qm+=im*dt;
          qp+=ip*dt;
          wave << (k+0.5)*dt << "," << im << "," << ip
               << "," << qm << "," << qp << "\n";
        }
        Result result{dt,step,10.*ue,10.*ve,10.*we,te,
                      status,ok,complete,arrived,qm,qp,qmRamo,qpRamo,
                      std::abs(qm-qmRamo),std::abs(qp-qpRamo)};
        all.push_back(result);
        const bool pass = ok && complete && arrived &&
                   CloseCharge(qm,qmRamo) && CloseCharge(qp,qpRamo);
        gates &= pass;
        summary << js << "," << jd << "," << dt << "," << step << ","
                << o.u0 << "," << o.v0 << "," << o.w0 << "," << o.t0 << ","
                << 10.*ue << "," << 10.*ve << "," << 10.*we << "," << te << ","
                << status << "," << int(ok) << "," << int(complete) << ","
                << int(arrived) << "," << pm0 << "," << pm1 << ","
                << pp0 << "," << pp1 << ","
                << qm << "," << qmRamo << "," << result.errMinus << ","
                << qp << "," << qpRamo << "," << result.errPlus << "\n";
        std::cout << std::setprecision(9) << "dt=" << dt << " ns, step="
                  << step << " mm: ion end v=" << 10.*ve << " mm, t="
                  << te << " ns, Q_minus=" << qm << " (Ramo " << qmRamo
                  << "), Q_plus=" << qp << " (Ramo " << qpRamo
                  << "), gate=" << (pass?"PASS":"FAIL") << "\n";
      }
    }
    // Drift-step convergence: compare the same finest time bin across steps.
    const auto finest = std::min_element(o.dts.begin(),o.dts.end())-o.dts.begin();
    for(std::size_t js=1;js<o.steps.size();++js) {
      const auto& a=all[(js-1)*o.dts.size()+finest];
      const auto& b=all[js*o.dts.size()+finest];
      const double du=std::hypot(a.ue-b.ue,a.ve-b.ve);
      const double dw=std::abs(a.we-b.we);
      const bool agree = du<0.05 && dw<0.05 &&
            std::abs(a.ramoPlus-b.ramoPlus)<Tolerance(b.ramoPlus) &&
            std::abs(a.ramoMinus-b.ramoMinus)<Tolerance(b.ramoMinus);
      std::cout<<"drift-step comparison "<<js-1<<"/"<<js
               <<" endpoint displacement="<<du<<" mm, dw="<<dw
               <<" mm; "<<(agree?"PASS":"FAIL")<<"\n";
      gates &= agree;
    }
    std::ofstream metadata(o.outDir/"single_ion_metadata.txt");
    metadata << "NA60+ Simulation / numerical one-ion cathode-strip benchmark\n"
             << "Ion: one +e charge; no avalanche, electrons or electronics\n"
             << "Real E field: Garfield ComponentAnalyticField, P3 noFW\n"
             << "Weighting potential: 3D Gmsh/Elmer, explicit grounded anode wires\n"
             << "Geometry: pitch_mm="<<o.pitch<<", wire_diam_um="<<o.wireDiameterUm
             <<", gap_minus_mm="<<o.gapMinus<<", gap_plus_mm="<<o.gapPlus<<"\n"
             << "HV_V="<<o.hv<<" (1800 V is a test-only historical baseline)\n"
             << "Ion mobility="<<o.mobility<<" (species/mixture not validated)\n"
             << "Diffusion disabled for deterministic convergence test\n"
             << "Finite FEM patch; ideal grounded inter-strip gaps; no PCB FR4\n"
             << "Error tolerance: max(3% of Ramo charge, 0.005 elementary charge)\n"
             << "Recorded window must cover full trajectory; cathode arrival required\n"
             << "GATES="<<(gates?"PASS":"FAIL")<<"\n";
    std::cout << (gates ? "PASS" : "FAIL")
              << " single-ion current/charge and convergence gates\n";
    return gates ? 0 : 5;
  } catch (const std::exception& e) {
    std::cerr << "SINGLE ION BENCHMARK FAILURE: " << e.what() << "\n";
    return 2;
  }
}

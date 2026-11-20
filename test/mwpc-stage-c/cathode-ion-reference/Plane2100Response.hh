// NA60+ 2100-V P3 noFW: real nine-wire analytic field, two WHOLE planes.
// Positive-ion motion uses a PROVISIONAL Ar+ in Ar mobility table.
// Never interpret these as physical U/V stereo strips or calibrated ion times.
#pragma once
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace NA60Plane {
constexpr double eFc=1.602176634e-4;
constexpr double pitch=0.3, wireRadius=0.0015, gap=0.3;
inline void Require(bool ok,const std::string& msg) {
  if(!ok)throw std::runtime_error("NA60+ plane response: "+msg);
}
struct Direct {
  std::array<double,2> endpointQ{},integralQ{};
  double endU=0,endV=0,endW=0,endT=0;
  int status=0;
  bool ok=false;
};
class Transport {
  Garfield::MediumMagboltz gas_;
  Garfield::ComponentAnalyticField field_;
  Garfield::Sensor sensor_;
  Garfield::DriftLineRKF drift_;
  double dt_;
  int bins_;
public:
  Transport(const std::string& mobility,double dt,int bins,double stepMm)
    : drift_(&sensor_),dt_(dt),bins_(bins) {
    Require(dt>0 && bins>=2 && bins<=10000,"invalid time grid");
    Require(stepMm>0 && stepMm<=0.5,"invalid RKF max step");
    gas_.SetComposition("ar",70.,"co2",30.);
    gas_.SetTemperature(293.15);gas_.SetPressure(760.);
    gas_.Initialise(false);
    Require(gas_.LoadIonMobility(mobility),"could not load ion mobility "+mobility);
    field_.SetMedium(&gas_);
    for(int k=-4;k<=4;++k)field_.AddWire(k*pitch,0.,0.003,2100.,"anode");
    field_.AddPlaneY(-gap,0.,"cathode_minus");
    field_.AddPlaneY(+gap,0.,"cathode_plus");
    field_.SetMagneticField(0.,0.,0.);
    field_.AddReadout("cathode_plus");
    field_.AddReadout("cathode_minus");
    sensor_.AddComponent(&field_);
    sensor_.AddElectrode(&field_,"cathode_plus");
    sensor_.AddElectrode(&field_,"cathode_minus");
    sensor_.SetArea(-1.35,-gap,-5.,1.35,gap,5.);
    sensor_.SetTimeWindow(0.,dt_,bins_);
    drift_.EnableSignalCalculation(true);
    drift_.UseWeightingPotential(true);
    // DriftLineRKF integrates mean ion drift deterministically; unlike
    // AvalancheMC it does not add stochastic diffusion at each step.
    // No EnableDiffusion API exists on DriftLineRKF in installed Garfield++.
    drift_.SetSignalAveragingOrder(2);
    drift_.SetMaximumStepSize(0.1*stepMm);
  }
  void CheckStart(double u,double v,double w,double t) const {
    Require(std::isfinite(u)&&std::isfinite(v)&&std::isfinite(w)&&std::isfinite(t),
            "nonfinite birth");
    const double nearest=std::round(u/pitch)*pitch;
    Require(std::hypot(u-nearest,v)>wireRadius &&
            std::abs(u)<1.35 && std::abs(v)<gap &&
            std::abs(w)<5 && t>=0,"birth outside gas or within metal");
    // This symmetry-axis shortcut is undefined when a birth is EXACTLY at v=0.
    Require(v!=0.,"birth exactly on v=0: endpoint-side shortcut undefined");
  }
  std::array<double,2> PredictedEndpointCharge(double u,double v,double w) {
    const double pp=field_.WeightingPotential(u,v,w,"cathode_plus");
    const double pm=field_.WeightingPotential(u,v,w,"cathode_minus");
    Require(std::isfinite(pp)&&std::isfinite(pm) &&
            pp>-0.01&&pm>-0.01&&pp<1.01&&pm<1.01,
            "nonphysical whole-plane weighting potential");
    // In perfectly symmetric P3 noFW field ion drifts to cathode sign(v).
    // IMPORTANT: this is an endpoint-side hypothesis, NOT a true ion drift.
    return {eFc*((v>0?1.:0.)-pp),eFc*((v<0?1.:0.)-pm)};
  }
  Direct Drift(double u,double v,double w,double t) {
    sensor_.ClearSignal();
    Direct out;
    out.ok=drift_.DriftIon(u,v,w,t);
    drift_.GetEndPoint(out.endU,out.endV,out.endW,out.endT,out.status);
    if(!(out.ok&&std::isfinite(out.endT)&&out.endT>t &&
         out.endT < dt_*bins_)) {
      std::ostringstream diagnostic;
      diagnostic<<std::setprecision(17)
       <<"ion drift did not satisfy transport/time-window integrity: "
       <<"transport_ok="<<out.ok
       <<" endpoint_status="<<out.status
       <<" start_cm=("<<u<<","<<v<<","<<w<<")"
       <<" start_time_ns="<<t
       <<" end_cm=("<<out.endU<<","<<out.endV<<","<<out.endW<<")"
       <<" end_time_ns="<<out.endT
       <<" duration_ns="<<(out.endT-t)
       <<" window_end_ns="<<(dt_*bins_)
       <<" valid_end_time="<<std::isfinite(out.endT)
       <<" inside_time_window="<<(std::isfinite(out.endT)&&out.endT<dt_*bins_);
      Require(false,diagnostic.str());
    }
    Require(std::abs(std::abs(out.endV)-gap)<0.001,
            "ion did not end on a physical cathode");
    const bool plus=out.endV>0;
    Require(plus==(v>0),"initial-v endpoint-side hypothesis failed");
    const double pp0=field_.WeightingPotential(u,v,w,"cathode_plus");
    const double pm0=field_.WeightingPotential(u,v,w,"cathode_minus");
    const double pp1=field_.WeightingPotential(out.endU,out.endV,out.endW,"cathode_plus");
    const double pm1=field_.WeightingPotential(out.endU,out.endV,out.endW,"cathode_minus");
    out.endpointQ={eFc*(pp1-pp0),eFc*(pm1-pm0)};
    for(int bin=0;bin<bins_;++bin) {
      out.integralQ[0]+=dt_*sensor_.GetIonSignal("cathode_plus",bin);
      out.integralQ[1]+=dt_*sensor_.GetIonSignal("cathode_minus",bin);
    }
    for(int k=0;k<2;++k) {
      const double tol=std::max(2.e-5*eFc,1.e-3*std::abs(out.endpointQ[k]));
      Require(std::abs(out.integralQ[k]-out.endpointQ[k])<tol,
              "Ramo endpoint versus integrated current failed");
    }
    return out;
  }
  double Current(int side,int bin) {
    return sensor_.GetIonSignal(side==0?"cathode_plus":"cathode_minus",bin);
  }
};
} // namespace NA60Plane

// Read-only Phase C control: real p3_noFW field, positive-ion drift,
// wire-aware weighting fields of BOTH COMPLETE cathode planes.
// Test starts (or weighted input ion births). No individual U/V strips.
// No Garfield primary/avalanche simulation in this program.
// Ar+ in pure Ar mobility is an explicit first-pass APPROXIMATION to
// transport in Ar/CO2 70:30. Quantitative arrival times are provisional.
//
// Build as independent target; never change archived microscopic physics.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {
struct Ion {double u,v,w,t,weight;};
constexpr double eFc = 1.602176634e-4;
constexpr double pitchCm=.30;
constexpr double radiusCm=.0015;
constexpr double cathodeCm=.30;

std::string Arg(int argc,char**argv,const std::string& key,
                const std::string& fallback="") {
  for(int i=1;i+1<argc;++i)
    if(argv[i]==key)return argv[i+1];
  return fallback;
}
double ArgDouble(int argc,char**argv,const std::string& key,double v) {
  const auto s=Arg(argc,argv,key);
  return s.empty()?v:std::stod(s);
}
std::vector<Ion> Starts(const std::string& filename) {
  if(filename.empty()) {
    // Eight physically allowed unit positive ions to verify upper/lower
    // mirror response. Radii are >15um, the 30um anode wire radius.
    return {{0.,+.003,0.,0.,1.},{0.,-.003,0.,0.,1.},
            {+.003,+.002,0.,0.,1.},{+.003,-.002,0.,0.,1.},
            {-.003,+.002,0.,0.,1.},{-.003,-.002,0.,0.,1.},
            {0.,+.020,0.,0.,1.},{0.,-.020,0.,0.,1.}};
  }
  std::ifstream in(filename);
  if(!in)throw std::runtime_error("Cannot read frozen cloud sample: "+filename);
  std::string h;
  std::getline(in,h);
  if(!h.empty()&&h.back()=='\r')h.pop_back();
  if(h!="u_mm,v_mm,w_mm,t_ns,weight")
    throw std::runtime_error("Ion input must have u_mm,v_mm,w_mm,t_ns,weight header");
  std::vector<Ion> out;
  std::string line;
  while(std::getline(in,line)) {
    if(line.empty())continue;
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream ss(line);
    Ion i{};
    if(!(ss>>i.u>>i.v>>i.w>>i.t>>i.weight))
      throw std::runtime_error("Malformed input ion birth CSV");
    i.u*=.1;i.v*=.1;i.w*=.1;
    if(!std::isfinite(i.u)||!std::isfinite(i.v)||!std::isfinite(i.w)||
       !std::isfinite(i.t)||!std::isfinite(i.weight)||i.weight<=0.)
      throw std::runtime_error("Invalid source ion coordinate/time/weight");
    out.push_back(i);
  }
  if(out.empty()||out.size()>5000)
    throw std::runtime_error("Provide 1..5000 selected/weighted real ion births");
  return out;
}
void CheckStart(const Ion& i) {
  const double du=i.u-std::round(i.u/pitchCm)*pitchCm;
  const double rho=std::hypot(du,i.v);
  if(rho<=radiusCm||std::abs(i.v)>=cathodeCm||
     std::abs(i.u)>1.3||std::abs(i.w)>4.8||i.t<0.)
    throw std::runtime_error("Input positive ion is inside a metal or outside p3 noFW gas patch");
}
}

int main(int argc,char**argv) {
 try {
  if(argc>1&&std::string(argv[1])=="--selftest") {
    const auto s=Starts("");
    for(const auto& i:s)CheckStart(i);
    if(s.size()!=8)throw std::runtime_error("bad default ion sample");
    std::cout<<"PASS two-cathode p3 noFW input/mirror/metal geometry checks\n";
    return 0;
  }
  const auto mobility=Arg(argc,argv,"--mobility");
  const auto source=Arg(argc,argv,"--input");
  const auto output=Arg(argc,argv,"--out-prefix","phaseC_p3_noFW_two_cathodes");
  const double dt=ArgDouble(argc,argv,"--dt-ns",100.);
  const double tmaxUs=ArgDouble(argc,argv,"--tmax-us",250.);
  if(mobility.empty()||dt<=0||tmaxUs<=0||tmaxUs>5000.)
    throw std::runtime_error("Provide --mobility FILE and positive --dt-ns/--tmax-us");
  const int nBins=static_cast<int>(std::ceil(1000.*tmaxUs/dt));
  if(nBins<2||nBins>100000)
    throw std::runtime_error("Time grid must contain 2..100000 bins");
  const auto ions=Starts(source);
  for(const auto& i:ions)CheckStart(i);
  for(const std::string name:{output+"_waveforms.csv",output+"_ions.csv",
                              output+"_metadata.txt"}) {
    if(std::ifstream(name).good())
      throw std::runtime_error("Refusing to overwrite existing result: "+name);
  }

  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar",70.,"co2",30.);
  gas.SetTemperature(293.15);gas.SetPressure(760.);
  gas.Initialise(false);
  if(!gas.LoadIonMobility(mobility))
    throw std::runtime_error("Could not load explicit ion mobility "+mobility);

  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);
  for(int k=-4;k<=4;++k)field.AddWire(k*pitchCm,0.,.003,1800.,"anode");
  field.AddPlaneY(-cathodeCm,0.,"cathode_minus");
  field.AddPlaneY(+cathodeCm,0.,"cathode_plus");
  field.SetMagneticField(0.,0.,0.);
  // Exactly these physical planes at unit weighting potential,
  // with all other conductors (including 9 anodes) at ground.
  field.AddReadout("cathode_plus");
  field.AddReadout("cathode_minus");
  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  sensor.AddElectrode(&field,"cathode_plus");
  sensor.AddElectrode(&field,"cathode_minus");
  sensor.SetArea(-1.35,-.3,-5.,1.35,.3,5.);
  sensor.SetTimeWindow(0.,dt,nBins);
  Garfield::DriftLineRKF ion(&sensor);
  ion.EnableSignalCalculation(true);
  ion.UseWeightingPotential(true);
  ion.SetSignalAveragingOrder(2);
  ion.SetMaximumStepSize();

  std::vector<double> totalPlus(nBins,0.),totalMinus(nBins,0.);
  std::ofstream points(output+"_ions.csv");
  points<<std::setprecision(13)
        <<"index,u0_mm,v0_mm,w0_mm,t0_ns,weight,"
          "u_end_mm,v_end_mm,w_end_mm,t_end_ns,status,drift_ok,"
          "plus_phi_start,plus_phi_end,minus_phi_start,minus_phi_end,"
          "plus_ramo_delta_Q_fC,minus_ramo_delta_Q_fC,"
          "plus_window_Q_fC,minus_window_Q_fC\n";
  int outOfWindow=0;
  for(std::size_t n=0;n<ions.size();++n) {
    const auto& p=ions[n];
    sensor.ClearSignal();
    const double pp0=field.WeightingPotential(p.u,p.v,p.w,"cathode_plus");
    const double pm0=field.WeightingPotential(p.u,p.v,p.w,"cathode_minus");
    const bool ok=ion.DriftIon(p.u,p.v,p.w,p.t);
    double ue=0,ve=0,we=0,te=0;
    int status=0;
    ion.GetEndPoint(ue,ve,we,te,status);
    const double pp1=field.WeightingPotential(ue,ve,we,"cathode_plus");
    const double pm1=field.WeightingPotential(ue,ve,we,"cathode_minus");
    double qp=0,qm=0;
    for(int j=0;j<nBins;++j){
      const double ip=sensor.GetIonSignal("cathode_plus",j)*p.weight;
      const double im=sensor.GetIonSignal("cathode_minus",j)*p.weight;
      totalPlus[j]+=ip;totalMinus[j]+=im;
      qp+=dt*ip;qm+=dt*im;
    }
    if(te>=dt*nBins)++outOfWindow;
    points<<n<<","<<10*p.u<<","<<10*p.v<<","<<10*p.w<<","<<p.t<<","
          <<p.weight<<","<<10*ue<<","<<10*ve<<","<<10*we<<","<<te
          <<","<<status<<","<<(ok?1:0)<<","<<pp0<<","<<pp1<<","
          <<pm0<<","<<pm1<<","
          <<eFc*p.weight*(pp1-pp0)<<","<<eFc*p.weight*(pm1-pm0)
          <<","<<qp<<","<<qm<<"\n";
    if(n%50==0)std::cout<<"Ion drift "<<n+1<<"/"<<ions.size()<<std::endl;
  }
  std::ofstream waves(output+"_waveforms.csv");
  waves<<std::setprecision(13)<<"time_ns,plus_ion_fC_per_ns,minus_ion_fC_per_ns\n";
  for(int j=0;j<nBins;++j)
    waves<<(j+.5)*dt<<","<<totalPlus[j]<<","<<totalMinus[j]<<"\n";
  std::ofstream meta(output+"_metadata.txt");
  meta<<"p3_noFW only: anode pitch 3mm, diameter 30um, gap +-3mm\n";
  meta<<"Gas Ar/CO2 70:30, 760 Torr, 293.15 K, anodes 1800V, B=0\n";
  meta<<"Positive ions only; no prompt electron current or electronics\n";
  meta<<"Ion-mobility approximation Ar+ in Ar, not calibrated Ar/CO2 mixture\n";
  meta<<"Mobility file="<<mobility<<"\n";
  meta<<"Weighting potentials EXACT same wire-aware ComponentAnalyticField\n";
  meta<<"Ion sample file="<<(source.empty()?"8 arbitrary symmetric smoke points":source)<<"\n";
  meta<<"N ions or macroions="<<ions.size()<<" time_window_us="<<tmaxUs
      <<" dt_ns="<<dt<<"\n";
  meta<<"Ions beyond chosen window="<<outOfWindow<<"\n";
  meta<<"Ramo q delta phi compared per-ion with integrated sampled Sensor signal\n";
  std::cout<<"PASS p3 noFW two WHOLE cathode ion-only signal; ions "
           <<ions.size()<<", beyond_window="<<outOfWindow<<"\n"
           <<"Saved "<<output<<"_{ions,waveforms,metadata}.*\n";
  return 0;
 } catch(const std::exception& e) {
   std::cerr<<"PHASE C ION REFERENCE FAILURE: "<<e.what()<<"\n";
   return 2;
 }
}

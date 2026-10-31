// SEPARATE experimental RKF+Polya MWPC baseline. No Heed or microscopic avalanches.
// Not a substitute for the microscopic ion-birth census / cathode signal.
// Uses the exact ORIGINAL D500 start CSV schema and the SAME frozen B=0
// analytic 9-anode/optional-75um-field-wire geometry from native B1.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {
struct Config {std::string name; double pitch, fwDiameter;};
const std::array<Config,4> configs{{{"p3_noFW",0.3,0.},{"p3_FW75",0.3,0.0075},
                                    {"p4_noFW",0.4,0.},{"p4_FW75",0.4,0.0075}}};
struct Start {std::string geom, originalSHA; int trials,box,iu,iv,sample; double u,v,w;};
std::vector<std::string> split(const std::string& s) {
  std::stringstream in(s);std::vector<std::string> z;std::string x;
  while(std::getline(in,x,','))z.push_back(x);
  return z;
}
double numeric(const std::string& s) {
  std::size_t n=0;const double v=std::stod(s,&n);
  if(n!=s.size()||!std::isfinite(v))throw std::runtime_error("Nonfinite/trailing numeric input");
  return v;
}
int integer(const std::string& s) {
  std::size_t n=0;const int v=std::stoi(s,&n);
  if(n!=s.size())throw std::runtime_error("Trailing integer input");
  return v;
}
const Config& cfg(const std::string& s){
  const auto it=std::find_if(configs.begin(),configs.end(),[&](const Config& g){return g.name==s;});
  if(it==configs.end())throw std::runtime_error("Unknown geometry: "+s);
  return *it;
}
std::vector<Start> readStarts(const std::string& path) {
  std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read: "+path);
  const std::string head=
    "lut_trials,geometry,original_lut_sha256,box_id,iu,iv,sample,random_seed,"
    "u_start_mm,v_start_mm,w_start_mm,old_drift_ok,old_collected,"
    "old_field_wire_collected,old_wire_index,old_phi_hit_deg,"
    "old_drift_time_ns,old_last_step_energy_eV";
  std::string line;if(!std::getline(in,line))throw std::runtime_error("Empty CSV");
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  if(line!=head)throw std::runtime_error("Wrong original exact D500 CSV header");
  std::vector<Start> out;
  while(std::getline(in,line)){
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line.empty())continue;
    auto f=split(line);if(f.size()!=18)throw std::runtime_error("Wrong CSV column count");
    Start z;
    z.trials=integer(f[0]);z.geom=f[1];z.originalSHA=f[2];
    z.box=integer(f[3]);z.iu=integer(f[4]);z.iv=integer(f[5]);
    z.sample=integer(f[6]);const auto originalSeed=integer(f[7]);
    z.u=numeric(f[8]);z.v=numeric(f[9]);z.w=numeric(f[10]);
    const auto& g=cfg(z.geom);const int nu=static_cast<int>(std::lround(10*g.pitch/.1));
    if(z.trials!=500 || z.sample<0 || z.sample>=500 ||
       z.iu<0||z.iu>=nu||z.iv<0||z.iv>=30||z.box!=z.iu*30+z.iv||
       originalSeed<0 || z.originalSHA.size()!=64 ||
       !std::all_of(z.originalSHA.begin(),z.originalSHA.end(),[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');}) ||
       z.u < -5*g.pitch+z.iu*.1-1e-6 || z.u > -5*g.pitch+(z.iu+1)*.1+1e-6 ||
       z.v < -3.+z.iv*.2-1e-6 || z.v> -3.+(z.iv+1)*.2+1e-6 ||
       std::abs(z.w)>1e-8) throw std::runtime_error("Invalid original D500 exact start");
    out.push_back(z);
    if(out.size()>400)throw std::runtime_error("Pilot limit 400 rows; split input into batches");
  }
  if(out.empty())throw std::runtime_error("No selected original D500 starts");
  return out;
}
bool inMetal(const Start& s) {
  const auto& g=cfg(s.geom);
  const double u=.1*s.u, v=.1*s.v;
  const double a=u-std::round(u/g.pitch)*g.pitch;
  if(std::hypot(a,v)<=.0015)return true;  // anode diameter 30um
  if(g.fwDiameter==0)return false;
  // field-wire centres lie at half-pitch from the nearest anode
  const double b=u-(std::round(u/g.pitch-.5)+.5)*g.pitch;
  return std::hypot(b,v)<=g.fwDiameter*.5;
}
void makeGas(const std::string& target,int ncoll) {
  if(ncoll<1||ncoll>10)throw std::runtime_error("ncoll must be 1..10 for a new pilot");
  if(std::filesystem::exists(target))throw std::runtime_error("Refuse to overwrite existing gas table: "+target);
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar",70.,"co2",30.);
  gas.SetTemperature(293.15);gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
  // Log-spaced E from drift to anode field strengths. B=0 only.
  gas.SetFieldGrid(100.,300000.,40,true);
  std::cout<<"Generating Ar70 CO2 30 B=0 293.15K 760Torr Magboltz gas table; ";
  std::cout<<"40 log E nodes from 100 to 300000 V/cm, ncoll="<<ncoll<<"\n"<<std::flush;
  gas.GenerateGasTable(ncoll);
  gas.WriteGasFile(target);
  if(!std::filesystem::is_regular_file(target)||std::filesystem::file_size(target)==0)
    throw std::runtime_error("Gas table not created correctly");
  std::cout<<"PASS RKF GAS TABLE WRITTEN: "<<target<<"\n";
}
void pilot(const std::vector<Start>& starts,const std::string& gasFile,
           const std::string& output,double theta,double p0) {
  if(std::filesystem::exists(output))throw std::runtime_error("Refuse to overwrite pilot output: "+output);
  if(!std::filesystem::is_regular_file(gasFile))throw std::runtime_error("Missing gas table");
  if(theta<0||theta>10||p0<0||p0>=1)throw std::runtime_error("Unphysical Pólya parameters");
  Garfield::MediumMagboltz gas;
  gas.LoadGasFile(gasFile);
  // Table itself contains the reference mixture, T and P; do not reset them.
  std::ofstream out(output);if(!out)throw std::runtime_error("Output not writable: "+output);
  out<<"geometry,box_id,sample,source_u_mm,source_v_mm,source_w_mm,physical_metal,"
        "rkf_transport_ok,rkf_mean_gain,rkf_loss,rkf_final_u_mm,rkf_final_v_mm,"
        "rkf_final_w_mm,rkf_end_status,polya_proxy_ion_births,polya_theta,p0_proxy,"
        "seconds_drift\n";
  out<<std::setprecision(16);
  for(const auto& geometry:configs){
    Garfield::ComponentAnalyticField field;
    field.SetMedium(&gas);
    for(int j=-4;j<=4;++j)field.AddWire(j*geometry.pitch,0.,.003,1800.,"anode");
    if(geometry.fwDiameter>0)
      for(int j=-4;j<4;++j)field.AddWire((j+.5)*geometry.pitch,0.,geometry.fwDiameter,0.,"field_wire");
    field.AddPlaneY(-.3,0.,"cathode_minus");field.AddPlaneY(.3,0.,"cathode_plus");
    field.SetMagneticField(0.,0.,0.);
    Garfield::Sensor sensor;
    sensor.AddComponent(&field);
    sensor.SetArea(-4.5*geometry.pitch,-.3,-5.,4.5*geometry.pitch,.3,5.);
    Garfield::DriftLineRKF drift;
    drift.SetSensor(&sensor);
    drift.EnableSignalCalculation(false);
    drift.EnableIonTail(false);
    drift.EnableAvalanche(true);
    // Do NOT call SetGainFluctuationsPolya here: Garfield applies it to
    // induced SIGNAL scaling, not necessarily the ion-birth count.
    // Instead explicitly sample a separate, named PROXY after obtaining
    // the macroscopic RKF mean gain. No ion cloud or Ramo signal claimed.
    std::mt19937_64 rng(20261031ULL);
    std::uniform_real_distribution<double> U(0.,1.);
    for(const auto& z:starts){
      if(z.geom!=geometry.name)continue;
      const bool metal=inMetal(z);
      bool good=false;
      double gain=0.,loss=0.,xx=std::numeric_limits<double>::quiet_NaN();
      double yy=xx,zz=xx,tt=xx;
      int endStatus=-999;
      double polyaProxy=0.;
      const auto start=std::chrono::steady_clock::now();
      if(!metal){
        good=drift.DriftElectron(.1*z.u,.1*z.v,.1*z.w,0.);
        if(good){
          drift.GetEndPoint(xx,yy,zz,tt,endStatus);
          gain=drift.GetGain();
          loss=drift.GetLoss();
          if(!std::isfinite(gain)||!std::isfinite(loss)||gain<0||loss<0||gain>1.e8)
            throw std::runtime_error("RKF returned invalid mean gain/loss; abort pilot");
          if(gain>0. && U(rng)>=p0){
            const double shape=theta+1.;
            std::gamma_distribution<double> polya(shape,gain/shape);
            polyaProxy=polya(rng);
          }
        }
      }
      const auto secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      // Failed nonmetal transport is NOT silently interpreted as zero-ion physics.
      // Physical metal starts get a zero proxy by explicit geometry only.
      out<<geometry.name<<","<<z.box<<","<<z.sample<<","<<z.u<<","<<z.v<<","<<z.w<<","
         <<(metal?1:0)<<","<<(good?1:0)<<",";
      if(good||metal)out<<gain<<","<<loss<<",";
      else out<<",,";
      if(good)out<<10.*xx<<","<<10.*yy<<","<<10.*zz<<","<<endStatus<<",";
      else out<<",,,,";
      if(good||metal)out<<polyaProxy;
      out<<","<<theta<<","<<p0<<","<<secs<<"\n";
    }
  }
  out.flush();if(!out)throw std::runtime_error("RKF pilot CSV write failure");
  std::cout<<"PASS RKF-POLYA PILOT OUTPUT (not yet validated versus physical ion births): "<<output<<"\n";
}
} // namespace
int main(int argc,char**argv){
  try{
    if(argc==3 && std::string(argv[1])=="--check"){
      const auto starts=readStarts(argv[2]);
      std::cout<<"PASS C++ ORIGINAL D500 exact-start parser ("<<starts.size()
               <<" rows); NO gas, DriftLineRKF, Garfield microscopic physics\n";
      return 0;
    }
    if((argc==3||argc==4) && std::string(argv[1])=="--make-gas"){
      makeGas(argv[2],argc==4?integer(argv[3]):5);
      return 0;
    }
    if(argc==7 && std::string(argv[1])=="--run"){
      const auto starts=readStarts(argv[2]);
      pilot(starts,argv[3],argv[4],numeric(argv[5]),numeric(argv[6]));
      return 0;
    }
    std::cerr<<"Usage:\n  --check OLD_D500_CSV\n  --make-gas NEW_GASFILE [ncoll]\n"
                "  --run OLD_D500_CSV AR_CO2_70_30_B0.gas OUT.csv POLYA_THETA P0_PROXY\n";
    return 2;
  }catch(const std::exception& e){std::cerr<<"RKF PILOT FAILURE: "<<e.what()<<"\n";return 2;}
}
// PI: compare microscopic electron drift, diffusion MC, deterministic RKF.
// Entirely separate from c_transport_box_map.cxx / frozen 20/500 maps.
// All methods see IDENTICAL start coordinates in the same 3+3 mm field.
// No avalanche gain: one launched electron -> at most one tracked end point.
// Outputs small CSV with per-electron timings and the EXACT frozen
// microscopic model's geometrical collection classification.
// A macroscopic Magboltz gas-table file is required for MC/RKF.
// Prepare it explicitly once (--mode prepare); startup not timed as drift.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"

namespace {
constexpr double Pi=3.14159265358979323846;
constexpr double PitchCm=0.3, WireDiamCm=.003, WireRadiusCm=.0015;
constexpr double FieldDiamCm=.0075, FieldRadiusCm=.00375;
constexpr double HalfGapCm=.3;
constexpr int HalfWires=4;

std::string Arg(int argc,char**argv,const std::string&key,
                const std::string&fallback) {
  for(int i=1;i<argc;++i) {
    if(std::string(argv[i])==key && i+1<argc)return argv[i+1];
    const std::string pref=key+"=";
    const std::string a=argv[i];
    if(a.rfind(pref,0)==0)return a.substr(pref.size());
  }
  return fallback;
}
int Integer(int argc,char**argv,const std::string&key,int fallback){
  return std::stoi(Arg(argc,argv,key,std::to_string(fallback)));
}
double Float(int argc,char**argv,const std::string&key,double fallback){
  return std::stod(Arg(argc,argv,key,std::to_string(fallback)));
}
struct Start{int iu,iv,trial;double u,v;unsigned seed;};
struct End{bool ok=false;int status=-999;double u=std::numeric_limits<double>::quiet_NaN();
 double v=std::numeric_limits<double>::quiet_NaN();
 double t=std::numeric_limits<double>::quiet_NaN();};
struct Classified{int anode=0,fieldwire=0,other=0,wire=999;double phi=std::numeric_limits<double>::quiet_NaN();};
Classified Classify(const End&e,bool useField) {
  Classified c;
  if(!std::isfinite(e.u)||!std::isfinite(e.v)){c.other=1;return c;}
  // Identical to original Stage-C microscopic endpoint proximity criterion:
  // nearest finite conductor, radius < 1.5 physical radius, status != -7.
  const int iw=static_cast<int>(std::lround(e.u/PitchCm));
  if(std::abs(iw)<=HalfWires &&
     std::hypot(e.u-iw*PitchCm,e.v)<1.5*WireRadiusCm &&
     e.status!=-7){
    c.anode=1;c.wire=iw;
    c.phi=180./Pi*std::atan2(e.v,e.u-iw*PitchCm);
    return c;
  }
  if(useField){
    const int jf=static_cast<int>(std::floor(e.u/PitchCm));
    if(jf>=-HalfWires&&jf<HalfWires &&
       std::hypot(e.u-(jf+0.5)*PitchCm,e.v)<1.5*FieldRadiusCm &&
       e.status!=-7){c.fieldwire=1;return c;}
  }
  c.other=1;return c;
}
bool InsideWire(double u,double v,bool useField) {
  int iw=static_cast<int>(std::lround(u/PitchCm));
  if(std::abs(iw)<=HalfWires && std::hypot(u-iw*PitchCm,v)<=WireRadiusCm)
    return true;
  if(useField){
    int jf=static_cast<int>(std::floor(u/PitchCm));
    if(jf>=-HalfWires&&jf<HalfWires &&
       std::hypot(u-(jf+.5)*PitchCm,v)<=FieldRadiusCm)return true;
  }
  return false;
}
std::vector<Start> Starts(int n,bool field) {
  // Sparse, stratified, intentionally NON-symmetrised coverage.
  // Include both near-anode and near-field-wire regions, both gas halves.
  const std::array<int,5> iu={1,7,14,22,28};
  const std::array<int,5> iv={1,7,14,22,28};
  std::mt19937 rng(73511);
  std::uniform_real_distribution<double> unif(0.,1.);
  std::vector<Start> points;
  for(int i:iu)for(int j:iv)for(int k=0;k<n;++k){
    double u=(-1.5+0.1*i+0.1*unif(rng))*.1;
    double v=(-3.+0.2*j+0.2*unif(rng))*.1;
    if(InsideWire(u,v,field)){
      u=(-1.5+.1*(i+.5))*.1;
      v=(-3.+.2*(j+.5))*.1;
    }
    if(InsideWire(u,v,field))throw std::runtime_error("Start in conductor");
    points.push_back({i,j,k,u,v,static_cast<unsigned>(410000+i*1000+j*100+k)});
  }
  return points;
}
void Geometry(Garfield::ComponentAnalyticField& f,Garfield::MediumMagboltz&gas,bool field){
  f.SetMedium(&gas);
  for(int i=-HalfWires;i<=HalfWires;++i)
    f.AddWire(i*PitchCm,0.,WireDiamCm,1800.,"anode");
  if(field)for(int i=-HalfWires;i<HalfWires;++i)
    f.AddWire((i+.5)*PitchCm,0.,FieldDiamCm,0.,"field_wire");
  f.AddPlaneY(-HalfGapCm,0.,"cathode_minus");
  f.AddPlaneY(+HalfGapCm,0.,"cathode_plus");
  f.SetMagneticField(0.,0.,0.);
}
void GasSetup(Garfield::MediumMagboltz&gas){
  gas.SetComposition("ar",70.,"co2",30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);
}
}

int main(int argc,char**argv){
  try{
    const std::string mode=Arg(argc,argv,"--mode","run");
    const std::string gasfile=Arg(argc,argv,"--gas-file","");
    const std::string output=Arg(argc,argv,"--output","");
    if(gasfile.empty())throw std::runtime_error("--gas-file must be provided");
    Garfield::MediumMagboltz gas;
    GasSetup(gas);
    if(mode=="prepare"){
      const int ncoll=Integer(argc,argv,"--ncoll",2);
      if(ncoll<1||ncoll>50)throw std::runtime_error("--ncoll out of bounds");
      const auto start=std::chrono::steady_clock::now();
      gas.SetFieldGrid(100.,500000.,16,true); // V/cm, logarithmic
      std::cout<<"PREPARE 70/30 Ar/CO2, 293.15 K, 760 Torr, 16 log E points "
               <<"100..500000 V/cm, B=0, Magboltz ncoll="<<ncoll
               <<" (in units of 10^7 collisions per point).\n"<<std::flush;
      gas.GenerateGasTable(ncoll,false);
      if(!gas.WriteGasFile(gasfile))throw std::runtime_error("Cannot save gas table");
      const auto dt=std::chrono::duration<double>(
         std::chrono::steady_clock::now()-start).count();
      std::cout<<"PASS GAS_TABLE wall_seconds="<<dt
               <<" destination="<<gasfile<<"\n";
      return 0;
    }
    if(mode!="run"||output.empty())throw std::runtime_error(
        "Mode run requires --output CSV; --gas-file required for MC/RKF");
    // LOAD FIRST: the microscopic cross sections still need Initialise().
    if(!gas.LoadGasFile(gasfile))throw std::runtime_error(
       "Cannot load transport table; first run --mode prepare");
    gas.SetMaxElectronEnergy(200.);
    gas.Initialise(false);
    const int trials=Integer(argc,argv,"--trials-per-box",4);
    const double stepUm=Float(argc,argv,"--mc-step-um",5.);
    if(trials<1||trials>100||stepUm<0.2||stepUm>20.)
       throw std::runtime_error("Invalid trial count/MC step (0.2..20um)");
    std::ofstream csv(output);
    if(!csv)throw std::runtime_error("Cannot write benchmark output CSV");
    csv<<std::setprecision(13)
       <<"geometry,method,mc_step_um,iu,iv,trial,random_seed,"
       <<"u_start_mm,v_start_mm,drift_ok,status,"
       <<"u_end_mm,v_end_mm,drift_time_ns,"
       <<"anode_collected,field_wire_collected,other_end,wire_index,"
       <<"phi_hit_deg,drift_wall_us\n";
    for(const bool fieldFlag: {false,true}){
      Garfield::ComponentAnalyticField field;
      Geometry(field,gas,fieldFlag);
      Garfield::Sensor sensor;
      sensor.AddComponent(&field);
      sensor.SetArea(-1.35,-.3,-2.,1.35,.3,2.);
      Garfield::AvalancheMicroscopic microscopic(&sensor);
      microscopic.EnableSignalCalculation(false);
      Garfield::AvalancheMC mc(&sensor);
      mc.SetDistanceSteps(stepUm*1.e-4);
      Garfield::DriftLineRKF rkf(&sensor);
      rkf.EnableSignalCalculation(false);
      rkf.EnableAvalanche(false);
      rkf.SetMaximumStepSize(5.e-4); // 5um fixed cap for thin wires.
      const auto points=Starts(trials,fieldFlag);
      const std::string geometry=fieldFlag?"p3_field75um":"p3_no_fw";
      for(const std::string method: {"microscopic","mc","rkf"}){
        long long totalUs=0,anodes=0;
        for(const auto& p:points){
          Garfield::RandomEngineRoot engine(p.seed);
          Garfield::Random::SetEngine(engine);
          End e;
          const auto t0=std::chrono::steady_clock::now();
          if(method=="microscopic"){
            e.ok=microscopic.DriftElectron(p.u,p.v,0.,0.,.1,0.,0.,0.);
            if(microscopic.GetNumberOfElectronEndpoints()>0){
              double u0,v0,w0,s0,e0,u1,v1,w1,s1,e1;
              microscopic.GetElectronEndpoint(0,u0,v0,w0,s0,e0,
                      u1,v1,w1,s1,e1,e.status);
              e.u=u1;e.v=v1;e.t=s1-s0;
            }
          }else if(method=="mc"){
            e.ok=mc.DriftElectron(p.u,p.v,0.,0.);
            if(mc.GetNumberOfElectronEndpoints()>0){
              double u0,v0,w0,s0,u1,v1,w1,s1;
              mc.GetElectronEndpoint(0,u0,v0,w0,s0,u1,v1,w1,s1,e.status);
              e.u=u1;e.v=v1;e.t=s1-s0;
            }
          }else{
            e.ok=rkf.DriftElectron(p.u,p.v,0.,0.);
            if(rkf.GetNumberOfDriftLinePoints()>0){
              double u1,v1,w1,s1;
              rkf.GetEndPoint(u1,v1,w1,s1,e.status);
              e.u=u1;e.v=v1;e.t=s1;
            }
          }
          const auto microseconds=std::chrono::duration_cast<std::chrono::microseconds>(
             std::chrono::steady_clock::now()-t0).count();
          totalUs+=microseconds;
          auto c=Classify(e,fieldFlag);
          anodes+=c.anode;
          csv<<geometry<<","<<method<<","<<stepUm<<","
             <<p.iu<<","<<p.iv<<","<<p.trial<<","<<p.seed<<","
             <<10*p.u<<","<<10*p.v<<","<<(e.ok?1:0)<<","<<e.status<<","
             <<10*e.u<<","<<10*e.v<<","<<e.t<<","
             <<c.anode<<","<<c.fieldwire<<","<<c.other<<","<<c.wire<<","
             <<c.phi<<","<<microseconds<<"\n";
        }
        csv.flush();
        std::cout<<geometry<<" "<<method<<": "<<points.size()
                 <<" matched starts, anode="<<anodes
                 <<" drift_us_total="<<totalUs
                 <<" mean_us="<<double(totalUs)/points.size()<<"\n"<<std::flush;
      }
    }
    std::cout<<"PASS: matched-start benchmark CSV: "<<output<<"\n";
    return 0;
  }catch(const std::exception&e){
    std::cerr<<"Benchmark error: "<<e.what()<<"\n";
    return 2;
  }
}

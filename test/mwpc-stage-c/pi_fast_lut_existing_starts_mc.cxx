// Garfield++ read-only MICRO-LUT-matched-start AvalancheMC benchmark.
// Does NOT calculate microscopic tracks or touch frozen LUTs.
// Input: independently selected and SHA-verified starts from existing 500
// microscopic trajectories/cell. Only macroscopic MC drift is run.
// Note: attachment off is a PHYSICS DIAGNOSTIC, not an approved model.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"

namespace {
constexpr double Pi=3.14159265358979323846;
constexpr double Pitch=.3, AnodeRadius=.0015, FieldRadius=.00375;
constexpr double Nan=std::numeric_limits<double>::quiet_NaN();
std::string Arg(int argc,char** argv,const std::string& key,
                const std::string& fallback="") {
  for(int i=1;i<argc;++i){
    const std::string a=argv[i];
    if(a==key && i+1<argc)return argv[i+1];
    if(a.rfind(key+"=",0)==0)return a.substr(key.size()+1);
  }
  return fallback;
}
struct Start{
  std::string geometry,split;
  int iu=0,iv=0,sample=0;
  unsigned seed=0;
  double u=0.,v=0.;
};
std::vector<std::string> Parts(const std::string& line){
  std::vector<std::string> out;
  std::istringstream ss(line);
  std::string part;
  while(std::getline(ss,part,','))out.push_back(part);
  return out;
}
bool IsFw(const std::string& key){
  if(key=="p3p0_anode_only")return false;
  if(key=="p3p0_field75p0um")return true;
  throw std::runtime_error("Unknown geometry: "+key);
}
std::vector<Start> ReadStarts(const std::string& file){
  std::ifstream in(file);
  if(!in)throw std::runtime_error("Missing selected frozen starts: "+file);
  std::string header;
  if(!std::getline(in,header))throw std::runtime_error("Empty starts CSV");
  const std::string expected=
      "geometry,split,iu,iv,sample,random_seed,u_start_mm,v_start_mm,"
      "reference_collected,reference_field_wire_collected,reference_status,"
      "reference_wire_index,reference_phi_hit_deg,reference_drift_time_ns";
  if(!header.empty()&&header.back()=='\r')header.pop_back();
  if(header!=expected)throw std::runtime_error("Unexpected selected reference CSV schema");
  std::vector<Start> starts;
  std::unordered_set<std::string> seen;
  int lineNo=1;
  std::string line;
  while(std::getline(in,line)){
    ++lineNo;
    if(line.empty())continue;
    auto parts=Parts(line);
    if(parts.size()!=14)throw std::runtime_error("Bad selected start CSV row "+std::to_string(lineNo));
    Start s;
    s.geometry=parts[0];
    IsFw(s.geometry);
    s.split=parts[1];
    if(s.split!="fit"&&s.split!="validation")
      throw std::runtime_error("Bad split label at row "+std::to_string(lineNo));
    s.iu=std::stoi(parts[2]);
    s.iv=std::stoi(parts[3]);
    s.sample=std::stoi(parts[4]);
    const long long value=std::stoll(parts[5]);
    if(value<1||value>2147483647LL)throw std::runtime_error("Bad MC seed");
    s.seed=static_cast<unsigned>(value);
    s.u=.1*std::stod(parts[6]);
    s.v=.1*std::stod(parts[7]);
    if(s.iu<0||s.iu>=30||s.iv<0||s.iv>=30||
       s.sample<0||s.sample>=500||!std::isfinite(s.u)||!std::isfinite(s.v))
      throw std::runtime_error("Invalid selected start indices/positions");
    const double uLo=(-1.5+.1*s.iu)*.1;
    const double vLo=(-3.+.2*s.iv)*.1;
    if(s.u<uLo-1.e-8||s.u>uLo+.01+1.e-8||
       s.v<vLo-1.e-8||s.v>vLo+.02+1.e-8)
      throw std::runtime_error("Selected start not inside claimed LUT cell");
    const std::string id=s.geometry+":"+std::to_string(s.iu)+":"+
      std::to_string(s.iv)+":"+std::to_string(s.sample);
    if(!seen.insert(id).second)throw std::runtime_error("Duplicate microscopic trial ID");
    starts.push_back(s);
  }
  if(starts.empty()||starts.size()>900000)
    throw std::runtime_error("Selected reference size out of bounds");
  return starts;
}
void Geometry(Garfield::ComponentAnalyticField& field,
              Garfield::MediumMagboltz& gas,bool fw){
  field.SetMedium(&gas);
  for(int i=-4;i<=4;++i)field.AddWire(i*Pitch,0.,.003,1800.,"anode");
  if(fw)for(int i=-4;i<4;++i)
    field.AddWire((i+.5)*Pitch,0.,.0075,0.,"field_wire");
  field.AddPlaneY(-.3,0.,"cathode_minus");
  field.AddPlaneY(+.3,0.,"cathode_plus");
  field.SetMagneticField(0.,0.,0.);
}
struct Classification{
  int anode=0,fw=0,other=0,wire=999999;
  double phi=Nan;
};
Classification Classify(double u,double v,int status,bool field){
  Classification result;
  if(!std::isfinite(u)||!std::isfinite(v)){
    result.other=1;return result;
  }
  const int iw=static_cast<int>(std::lround(u/Pitch));
  if(std::abs(iw)<=4&&std::hypot(u-iw*Pitch,v)<1.5*AnodeRadius&&
     status!=-7){
    result.anode=1;result.wire=iw;
    result.phi=180./Pi*std::atan2(v,u-iw*Pitch);
    return result;
  }
  if(field){
    const int jf=static_cast<int>(std::floor(u/Pitch));
    if(jf>=-4&&jf<4&&
       std::hypot(u-(jf+.5)*Pitch,v)<1.5*FieldRadius&&status!=-7){
      result.fw=1;return result;
    }
  }
  result.other=1;return result;
}
}
int main(int argc,char** argv){
 try{
  const std::string input=Arg(argc,argv,"--starts");
  const std::string gasfile=Arg(argc,argv,"--gas-file");
  const std::string output=Arg(argc,argv,"--output");
  const std::string attach=Arg(argc,argv,"--mc-attachment","off");
  const double step=std::stod(Arg(argc,argv,"--mc-step-um","5"));
  if(input.empty()||gasfile.empty()||output.empty())
    throw std::runtime_error("Required --starts FILE --gas-file FILE --output FILE");
  if(attach!="on"&&attach!="off")
    throw std::runtime_error("--mc-attachment must be on or off");
  if(step<.2||step>20.||!std::isfinite(step))
    throw std::runtime_error("Step must be 0.2..20 um");
  std::ifstream old(output);
  if(old.good())throw std::runtime_error("Output exists, refusing overwrite: "+output);
  const auto starts=ReadStarts(input);
  Garfield::MediumMagboltz gas;
  gas.SetComposition("ar",70.,"co2",30.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  if(!gas.LoadGasFile(gasfile))
    throw std::runtime_error("Cannot load original Magboltz table");
  gas.SetMaxElectronEnergy(200.);
  gas.Initialise(false);
  std::ofstream out(output);
  if(!out)throw std::runtime_error("Cannot create MC CSV");
  out<<std::setprecision(13)
     <<"geometry,split,iu,iv,sample,random_seed,u_start_mm,v_start_mm,"
       "mc_step_um,mc_attachment,drift_ok,status,"
       "u_end_mm,v_end_mm,w_end_mm,drift_time_ns,"
       "anode_collected,field_wire_collected,other_end,wire_index,"
       "phi_hit_deg,drift_wall_us\n";
  for(const std::string geom:{"p3p0_anode_only","p3p0_field75p0um"}){
    const bool fw=IsFw(geom);
    Garfield::ComponentAnalyticField field;
    Geometry(field,gas,fw);
    Garfield::Sensor sensor;
    sensor.AddComponent(&field);
    sensor.SetArea(-1.35,-.3,-2.,1.35,.3,2.);
    Garfield::AvalancheMC mc(&sensor);
    mc.SetDistanceSteps(step*1.e-4);
    mc.EnableAttachment(attach=="on");
    long long n=0,collected=0,statusAttach=0;
    double totalUs=0.;
    for(const auto& s:starts){
      if(s.geometry!=geom)continue;
      Garfield::RandomEngineRoot engine(s.seed);
      Garfield::Random::SetEngine(engine);
      const auto t0=std::chrono::steady_clock::now();
      const bool ok=mc.DriftElectron(s.u,s.v,0.,0.);
      double u=Nan,v=Nan,w=Nan,t=Nan;
      int status=-999;
      if(mc.GetNumberOfElectronEndpoints()>0){
        double u0=0.,v0=0.,w0=0.,startTime=0.,endTime=0.;
        mc.GetElectronEndpoint(0,u0,v0,w0,startTime,u,v,w,endTime,status);
        t=endTime-startTime;
      }
      const double us=1.e-3*std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now()-t0).count();
      const auto c=Classify(u,v,status,fw);
      ++n;collected+=c.anode;statusAttach+=(status==-7);totalUs+=us;
      out<<geom<<","<<s.split<<","<<s.iu<<","<<s.iv<<","
         <<s.sample<<","<<s.seed<<","<<10*s.u<<","<<10*s.v<<","
         <<step<<","<<attach<<","<<int(ok)<<","<<status<<","
         <<10*u<<","<<10*v<<","<<10*w<<","<<t<<","
         <<c.anode<<","<<c.fw<<","<<c.other<<","<<c.wire<<","
         <<c.phi<<","<<us<<"\n";
    }
    out.flush();
    std::cout<<geom<<": starts="<<n<<", anode="<<collected
             <<", attachment_status_minus7="<<statusAttach
             <<", measured_MC_drift_seconds="<<totalUs*1.e-6
             <<", mean_us_per_electron="<<(n>0?totalUs/n:0.)
             <<"\n"<<std::flush;
    if(n==0)throw std::runtime_error("No selected reference rows for geometry "+geom);
  }
  std::cout<<"PASS matched frozen LUT500 starts -> "<<output
           <<", MC step um="<<step<<", attachment="<<attach<<"\n";
  return 0;
 }catch(const std::exception& e){
  std::cerr<<"Fast-LUT-MC benchmark error: "<<e.what()<<"\n";
  return 2;
 }
}

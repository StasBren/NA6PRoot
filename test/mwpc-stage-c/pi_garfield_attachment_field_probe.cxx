// Diagnose transport-table attachment against the actual 3 mm MWPC E field.
// Reads the SAME Magboltz .gas file as the three-method benchmark.
// NO microscopic trajectories or new Magboltz table generation.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

namespace {
std::string Arg(const int argc,char** argv,const std::string& key) {
  for(int i=1;i<argc;++i){
    const std::string a=argv[i];
    if(a==key && i+1<argc)return argv[i+1];
    if(a.rfind(key+"=",0)==0)return a.substr(key.size()+1);
  }
  return "";
}
void Geometry(Garfield::ComponentAnalyticField& field,
              Garfield::MediumMagboltz& gas,const bool fw){
  field.SetMedium(&gas);
  for(int i=-4;i<=4;++i)field.AddWire(i*.3,0.,.003,1800.,"anode");
  if(fw)for(int i=-4;i<4;++i)
    field.AddWire((i+.5)*.3,0.,.0075,0.,"field_wire");
  field.AddPlaneY(-.3,0.,"cathode_minus");
  field.AddPlaneY(+.3,0.,"cathode_plus");
  field.SetMagneticField(0.,0.,0.);
}
struct Transport{
  bool ok=false;
  double eta=std::numeric_limits<double>::quiet_NaN();
  double alpha=std::numeric_limits<double>::quiet_NaN();
  double velocity=std::numeric_limits<double>::quiet_NaN();
  double dl=std::numeric_limits<double>::quiet_NaN();
  double dt=std::numeric_limits<double>::quiet_NaN();
};
Transport Coeff(Garfield::MediumMagboltz& gas,double ex,double ey,double ez) {
  Transport t;
  double vx=0.,vy=0.,vz=0.;
  const bool a=gas.ElectronAttachment(ex,ey,ez,0.,0.,0.,t.eta);
  const bool b=gas.ElectronTownsend(ex,ey,ez,0.,0.,0.,t.alpha);
  const bool c=gas.ElectronVelocity(ex,ey,ez,0.,0.,0.,vx,vy,vz);
  const bool d=gas.ElectronDiffusion(ex,ey,ez,0.,0.,0.,t.dl,t.dt);
  t.velocity=std::sqrt(vx*vx+vy*vy+vz*vz);
  t.ok=a&&b&&c&&d&&std::isfinite(t.eta)&&std::isfinite(t.alpha)
    &&std::isfinite(t.velocity)&&std::isfinite(t.dl)&&std::isfinite(t.dt);
  return t;
}
void Row(std::ofstream& o,const std::string& kind,const std::string& geom,
         double u,double v,double ex,double ey,double ez,
         Garfield::MediumMagboltz& gas){
  const double e=std::sqrt(ex*ex+ey*ey+ez*ez);
  const auto c=Coeff(gas,ex,ey,ez);
  o<<kind<<","<<geom<<","<<u<<","<<v<<","<<e<<","
   <<int(c.ok)<<","<<c.eta<<","<<c.alpha<<","<<c.velocity<<","
   <<c.dl<<","<<c.dt<<"\n";
}
}
int main(int argc,char** argv){
  try{
    const auto gasfile=Arg(argc,argv,"--gas-file");
    const auto output=Arg(argc,argv,"--output");
    if(gasfile.empty()||output.empty())
      throw std::runtime_error("Use --gas-file path --output path");
    std::ifstream gtest(gasfile);
    if(!gtest)throw std::runtime_error("Gas table absent: "+gasfile);
    std::ifstream old(output);
    if(old.good())throw std::runtime_error("Output exists: "+output);
    Garfield::MediumMagboltz gas;
    gas.SetComposition("ar",70.,"co2",30.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    if(!gas.LoadGasFile(gasfile))
      throw std::runtime_error("Cannot load gas table");
    gas.SetMaxElectronEnergy(200.);
    gas.Initialise(false);
    std::ofstream out(output);
    if(!out)throw std::runtime_error("Cannot create CSV: "+output);
    out<<std::setprecision(12);
    out<<"kind,geometry,u_mm,v_mm,E_V_cm,transport_ok,"
       <<"eta_per_cm,alpha_per_cm,drift_cm_per_ns,"
       <<"diff_long_sqrt_cm,diff_trans_sqrt_cm\n";

    for(int j=0;j<=80;++j){
      const double E=100.*std::pow(5000.,j/80.);
      Row(out,"gas_sweep","n_a_n",0.,0.,E,0.,0.,gas);
    }
    for(const bool fw: {false,true}){
      Garfield::ComponentAnalyticField field;
      Geometry(field,gas,fw);
      Garfield::Sensor sensor;
      sensor.AddComponent(&field);
      sensor.SetArea(-1.35,-.3,-2.,1.35,.3,2.);
      const std::string geom=fw?"p3_field75um":"p3_no_fw";
      // Center and near-midpoint positions sample low- and high-field
      // regions; all points are gas, not inside a wire.
      const std::vector<double> uu={0.,.025,.075,.125,.15};
      const std::vector<double> vv={-.25,-.15,-.05,-.015,
                                    .015,.05,.15,.25};
      int selected=0;
      for(const double u:uu)for(const double v:vv){
        double ex=0.,ey=0.,ez=0.;
        Garfield::Medium* medium=nullptr;
        int status=-999;
        sensor.ElectricField(u,v,0.,ex,ey,ez,medium,status);
        if(status!=0||medium==nullptr)continue;
        Row(out,"chamber",geom,10.*u,10.*v,ex,ey,ez,gas);
        ++selected;
      }
      std::cout<<geom<<" sampled valid chamber points: "
               <<selected<<"\n";
    }
    std::cout<<"PASS attachment coefficients at gas-table grid and "
             "MWPC field points (units: eta/alpha 1/cm): "
             <<output<<"\n";
    return 0;
  }catch(const std::exception& e){
    std::cerr<<"Attachment-field probe error: "<<e.what()<<"\n";
    return 2;
  }
}

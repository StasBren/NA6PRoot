// Read-only independent Garfield analytic weighting-potential reference.
// Same nine wires, cathode gaps and voltage as Plane2100Response.hh.
// No particles, gas transport, archived data, or FEM.
#include "Garfield/ComponentAnalyticField.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
int main(int argc,char** argv) {
  if(argc!=2 || std::filesystem::exists(argv[1])) {
    std::cerr<<"Usage: PROGRAM new_reference.csv\n"; return 2;
  }
  Garfield::ComponentAnalyticField field;
  for(int j=-4;j<=4;++j) field.AddWire(.3*j,0.,.003,2100.,"anode");
  field.AddPlaneY(-.3,0.,"cathode_minus");
  field.AddPlaneY(.3,0.,"cathode_plus");
  field.AddReadout("cathode_minus");
  field.AddReadout("cathode_plus");
  std::ofstream out(argv[1]);
  if(!out) return 2;
  out<<std::setprecision(17)
     <<"region,wire_index,r_over_R,u_mm,v_mm,phi_plus,phi_minus\n";
  auto row=[&](const char* region,int j,double ratio,double u,double v) {
    const auto pp=field.WeightingPotential(.1*u,.1*v,0.,"cathode_plus");
    const auto pm=field.WeightingPotential(.1*u,.1*v,0.,"cathode_minus");
    if(!std::isfinite(pp)||!std::isfinite(pm)||pp<-.05||pm<-.05||pp>1.05||pm>1.05)
      return false;
    out<<region<<","<<j<<","<<ratio<<","<<u<<","<<v<<","<<pp<<","<<pm<<"\n";
    return bool(out);
  };
  constexpr double pi=3.14159265358979323846;
  for(int j=-4;j<=4;++j)
    for(double ratio:{1.01,1.1,1.5,2.,5.,20.})
      for(int a=0;a<16;++a) {
        double theta=(a+.37)*2*pi/16;
        if(!row("shell",j,ratio,3.*j+.015*ratio*std::cos(theta),
                   .015*ratio*std::sin(theta))) return 3;
      }
  for(int j=-13;j<=13;++j)
    for(double v:{-2.95,-2.,-1.,-.5,.5,1.,2.,2.95})
      if(!row("bulk",0,0.,j,v)) return 3;
  out.close();
  if(!out) return 3;
  std::cout<<"Wrote 1080 analytic whole-plane reference points\n";
  return 0;
}

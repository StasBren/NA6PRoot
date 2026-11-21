// NA60+ P3: read-only sample of real 3D Gmsh/Elmer U/V strip
// weighting potentials and all three local gradients near anode cylinders.
// A successful import or a Ramo identity is NOT proof of FEM convergence.
// No ion drift, avalanche, current, electronics or edits to source maps.
#include "Garfield/ComponentElmer.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
namespace fs=std::filesystem;
constexpr double kPi=3.1415926535897932384626433832795;
std::string Arg(int argc,char**argv,const std::string&key){
 for(int i=1;i+1<argc;++i)if(key==argv[i])return argv[i+1];
 return "";
}
void Need(bool b,const std::string&s) {
 if(!b)throw std::runtime_error(s);
}
void NeedFile(const fs::path&p){
 Need(fs::is_regular_file(p)&&fs::file_size(p)>0,
      "missing/empty FEM map component: "+p.string());
}
struct Eval {double phi=0,eu=0,ev=0,ew=0;};
class Probe {
 Garfield::ComponentElmer map_;
 double gap_=3.,radius_=0.015,pitch_=3.;
 double Phi(double u,double v,double w,const std::string&label) {
   const double p=map_.WeightingPotential(0.1*u,0.1*v,0.1*w,label);
   Need(std::isfinite(p),"nonfinite Elmer weighting potential for "+label);
   // Large excursions are invalid, but permit finite-element ringing.
   Need(p>-.15&&p<1.15,"unphysical Elmer weighting potential for "+label);
   return p;
 }
public:
 Probe(const fs::path&dir,double gap,double radius,double pitch):
   map_((dir/"chamber/mesh.header").string(),
        (dir/"chamber/mesh.elements").string(),
        (dir/"chamber/mesh.nodes").string(),
        (dir/"dielectrics.dat").string(),
        (dir/"chamber/drift.result").string(),"mm"),
   gap_(gap),radius_(radius),pitch_(pitch) {
    Need(map_.SetWeightingField(
      (dir/"chamber/weight_plus_0.result").string(),"plus_strip_0"),
      "could not import plus weighting map");
    Need(map_.SetWeightingField(
      (dir/"chamber/weight_minus_0.result").string(),"minus_strip_0"),
      "could not import minus weighting map");
 }
 Eval At(double u,double v,double w,double h,
         const std::string&label) {
   Need(std::abs(v)+h<gap_,"gradient point too close to cathode");
   const double uwire=std::round(u/pitch_)*pitch_;
   const double rho=std::hypot(u-uwire,v);
   // Both centered +/-h probes remain strictly outside all wire metal.
   Need(rho-h>radius_,"gradient point crosses anode conductor");
   Eval q;
   q.phi=Phi(u,v,w,label);
   q.eu=(Phi(u-h,v,w,label)-Phi(u+h,v,w,label))/(2*h);
   q.ev=(Phi(u,v-h,w,label)-Phi(u,v+h,w,label))/(2*h);
   q.ew=(Phi(u,v,w-h,label)-Phi(u,v,w+h,label))/(2*h);
   Need(std::isfinite(q.eu)&&std::isfinite(q.ev)&&std::isfinite(q.ew),
        "nonfinite finite-difference weighting field");
   return q;
 }
 double Raw(double u,double v,double w,const std::string&label){
   return Phi(u,v,w,label);
 }
};
} // namespace

int main(int argc,char**argv) {
 try {
  const fs::path dir=Arg(argc,argv,"--map-dir");
  const fs::path output=Arg(argc,argv,"--out-csv");
  const auto tag=Arg(argc,argv,"--tag");
  Need(!dir.empty()&&!output.empty()&&!tag.empty(),
       "usage: --map-dir MAP --out-csv NEW.csv --tag baseline|refined");
  Need(!fs::exists(output),"refusing to overwrite existing output CSV");
  Need(fs::is_directory(dir),"FEM map directory missing");
  for(const auto&part:{"chamber/mesh.header","chamber/mesh.elements",
       "chamber/mesh.nodes","dielectrics.dat","chamber/drift.result",
       "chamber/weight_plus_0.result","chamber/weight_minus_0.result"})
    NeedFile(dir/part);
  const double gap=3.,wirePitch=3.,wireR=0.015;
  const std::array<double,5> shells{{1.1,1.5,2.,5.,20.}};
  const std::array<double,8> angles{{
    kPi/8,3*kPi/8,5*kPi/8,7*kPi/8,
    9*kPi/8,11*kPi/8,13*kPi/8,15*kPi/8}};
  // A reusable read-only survey of the same predefined points, for the
  // two different FEM mesh resolutions. No physical ion transport here.
  Probe probe(dir,gap,wireR,wirePitch);
  const double nearPlus=probe.Raw(0.,+gap-.05,0.,"plus_strip_0");
  const double farPlus =probe.Raw(0.,-gap+.05,0.,"plus_strip_0");
  const double nearMinus=probe.Raw(0.,-gap+.05,0.,"minus_strip_0");
  const double farMinus =probe.Raw(0.,+gap-.05,0.,"minus_strip_0");
  std::cout<<std::setprecision(12)
   <<"boundary_near_plus="<<nearPlus
   <<" boundary_far_plus="<<farPlus
   <<" boundary_near_minus="<<nearMinus
   <<" boundary_far_minus="<<farMinus<<"\n";
  Need(nearPlus>0.5&&nearMinus>0.5&&
       std::abs(farPlus)<0.1&&std::abs(farMinus)<0.1,
       "FEM cathode weighting sanity check failed; refusing invalid sample");
  std::ofstream out(output);
  Need(bool(out),"unable to create sampling CSV");
  out<<std::setprecision(16)
    <<"mesh_tag,side,wire_index,r_over_R,theta_index,theta_rad,w_mm,"
      "u_mm,v_mm,h_mm,phi,Eu_per_mm,Ev_per_mm,Ew_per_mm,"
      "Eu_halfstep_per_mm,Ev_halfstep_per_mm,Ew_halfstep_per_mm,"
      "gradient_step_difference_per_mm\n";
  double maxHalfStep=0.,worstMirror=0.,maxP=0.;
  int points=0;
  for(int wire:{-1,0,1})for(double ratio:shells)for(int a=0;a<8;++a)
    for(double w:{0.,1.1}) {
      const double r=ratio*wireR, theta=angles[a];
      const double u=wire*wirePitch+r*std::cos(theta);
      const double v=r*std::sin(theta);
      const double clearance=r-wireR;
      const double h=std::min(0.0015,clearance*.15);
      Need(h>0 && h<clearance*.5,"invalid gradient step near wire");
      for(const auto&side:{"plus","minus"}){
        const std::string label=std::string(side)+"_strip_0";
        const auto q=probe.At(u,v,w,h,label);
        const auto fine=probe.At(u,v,w,h*.5,label);
        const double discrepancy=std::sqrt(
         std::pow(q.eu-fine.eu,2)+std::pow(q.ev-fine.ev,2)+
         std::pow(q.ew-fine.ew,2));
        maxHalfStep=std::max(maxHalfStep,discrepancy);
        maxP=std::max(maxP,std::abs(q.phi));
        const auto mirrored=probe.Raw(-u,-v,w,
                          std::string(side=="plus"?"minus":"plus")+"_strip_0");
        worstMirror=std::max(worstMirror,std::abs(q.phi-mirrored));
        Need(std::isfinite(discrepancy),
             "nonfinite gradient finite-difference step check");
        out<<tag<<","<<side<<","<<wire<<","<<ratio<<","<<a<<","
           <<theta<<","<<w<<","<<u<<","<<v<<","<<h<<","
           <<q.phi<<","<<q.eu<<","<<q.ev<<","<<q.ew<<","
           <<fine.eu<<","<<fine.ev<<","<<fine.ew<<","
           <<discrepancy<<"\n";
        ++points;
      }
    }
  out.close();
  Need(bool(out),"error writing FEM probe CSV");
  std::cout<<"COMPLETED near-wire 3D FEM individual-strip survey: "
           <<"points="<<points
           <<" max_gradient_h_vs_h2_difference_per_mm="<<maxHalfStep
           <<" max_reflection_potential_difference="<<worstMirror
           <<" max_abs_strip_potential="<<maxP<<"\n"
           <<"This is READ-ONLY diagnostic sampling, NOT a mesh convergence pass.\n"
           <<"No ions, no avalanche, no electronics. CSV: "<<output<<"\n";
  return 0;
 }catch(const std::exception&e){
  std::cerr<<"FEM STRIP PROBE FAILED: "<<e.what()<<"\n";
  return 2;
 }
}

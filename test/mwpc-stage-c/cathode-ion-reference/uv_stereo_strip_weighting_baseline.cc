// NA60+/DiCE U/V stereo-strip geometry: preliminary electrostatic baseline.
// Source: Milov "Readout for the MS modules ..." (2026-08-30),
// Sec 4.2: 1.5-2mm proposed MuTr strip width and tan(alpha)=0.1.
// The width (default 1.7mm) is a proposed option; the strip PITCH IS NOT
// SPECIFIED by that report and MUST be provided explicitly.
//
// Exact harmonic weighting potential for ONE infinitely long strip on one
// of TWO grounded infinite parallel planar boundaries, in coordinates
// rotated to the stereo strip direction. Analytic solution of Laplace's
// equation in the bare parallel-plate slab. No anode-wire conductors,
// PCB vias/bridges, dielectric layers, finite module edges or multiplexing.
// This is an orientation/boundary-convention control, NOT the physical
// wire-inclusive MWPC U/V weighting map and NOT an induced-current solver.
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
constexpr double pi=3.1415926535897932384626433832795;
void Need(bool pass,const std::string&what) {
 if(!pass)throw std::runtime_error("UV STRIP BASELINE: "+what);
}
std::string Arg(int n,char**v,const std::string&key) {
 for(int i=1;i+1<n;++i)if(key==v[i])return v[i+1];
 return "";
}
double Num(int n,char**v,const std::string&key,double fallback){
 const auto s=Arg(n,v,key);
 return s.empty()?fallback:std::stod(s);
}
struct Geometry {
 double upperGapMm=3.,lowerGapMm=3.;
 double stripWidthMm=1.7,stripPitchMm=0.;
 double stereoTan=0.1;
 int halfIndex=8;
 double Height()const{return upperGapMm+lowerGapMm;}
 double S(bool upper,double uMm,double wMm)const {
   // u: across wire; w: along wire. Long axes of strips deviate by
   // +/-atan(stereoTan) from u axis. s: normal-to-strip coordinate.
   const double sense=upper?1.:-1.;
   return (wMm-sense*stereoTan*uMm)/std::hypot(1.,stereoTan);
 }
 double ToUpperDistance(bool upper,double vMm) const {
   // Convert to slab coordinate y: 0 at the opposite grounded plane,
   // H at this readout plane.
   return upper?vMm+lowerGapMm:upperGapMm-vMm;
 }
 double Phi(bool upper,int strip,double uMm,double vMm,double wMm) const{
   const double H=Height(),y=ToUpperDistance(upper,vMm);
   Need(y>=-1.e-10&&y<=H+1.e-10,"query outside gas slab");
   const double s=S(upper,uMm,wMm);
   const double center=strip*stripPitchMm;
   const double a=center-.5*stripWidthMm,b=center+.5*stripWidthMm;
   if(y<=0.)return 0.;
   if(y>=H)return (s>a && s<b)?1.:0.;
   const double t=std::tan(pi*y/(2.*H));
   const double phi=(std::atan(std::tanh(pi*(b-s)/(2.*H))*t)-
                     std::atan(std::tanh(pi*(a-s)/(2.*H))*t))/pi;
   Need(std::isfinite(phi)&&phi>=-1.e-10&&phi<=1.+1.e-10,
        "Laplace strip potential outside [0,1]");
   return std::clamp(phi,0.,1.);
 }
 double Sum(bool upper,double u,double v,double w)const {
   double q=0.;
   for(int j=-halfIndex;j<=halfIndex;++j)q+=Phi(upper,j,u,v,w);
   return q;
 }
};
void Validate(const Geometry&g) {
 Need(g.upperGapMm>0&&g.lowerGapMm>0&&g.upperGapMm<20&&g.lowerGapMm<20,
      "positive finite cathode gaps required");
 Need(g.stripPitchMm>0&&g.stripWidthMm>0&&
      g.stripWidthMm<=g.stripPitchMm&&g.stripPitchMm<=15,
      "strip pitch required, must exceed physical strip width");
 Need(std::abs(g.stereoTan)<=1.&&g.halfIndex>=4&&g.halfIndex<=64,
      "invalid stereo tangent or preview number of strips");
 for(bool upper:{false,true}) {
   const double thisSide=upper?g.upperGapMm:-g.lowerGapMm;
   const double opposite=upper?-g.lowerGapMm:g.upperGapMm;
   const double eps=1.e-6;
   Need(g.Phi(upper,0,0.,thisSide,0.)==1.,"unit strip surface boundary");
   Need(g.Phi(upper,1,0.,thisSide,0.)==0.,"other strip surface boundary");
   Need(g.Phi(upper,0,0.,opposite,0.)==0.,"opposite grounded plane boundary");
   Need(g.Phi(upper,0,0.,thisSide+(upper?-eps:+eps),0.)>.99,
        "strip face interior near-boundary check");
   for(double u:{-2.,0.,2.})for(double w:{-2.,0.,2.})
    for(double v:{-.5*g.lowerGapMm,0.,.5*g.upperGapMm}) {
      const double x=g.Phi(upper,0,u,v,w);
      Need(x>=0.&&x<=1.,"electrostatic maximum principle");
      Need(g.Sum(upper,u,v,w)<=1.0000001,
           "nonoverlapping strip weighting potentials exceed whole plane");
    }
 }
 if(std::abs(g.upperGapMm-g.lowerGapMm)<1.e-12) {
   for(double u:{-.5,0.,.5})for(double v:{-1.,0.,1.})
    for(double w:{-1.,0.,1.}){
      const double a=g.Phi(true,0,u,v,w);
      const double b=g.Phi(false,0,-u,-v,w);
      Need(std::abs(a-b)<1.e-12,
           "U/V stereo reflection failed in symmetric candidate geometry");
    }
 }
}
} // namespace
int main(int argc,char**argv){
 try {
  Geometry g;
  const auto out=Arg(argc,argv,"--out-prefix");
  Need(!out.empty(),"usage: --out-prefix PREFIX --strip-pitch-mm PITCH "
      "[--strip-width-mm 1.7 --gap-plus-mm 3 --gap-minus-mm 3 "
      "--stereo-tan 0.1 --half-index 8]");
  const auto pitch=Arg(argc,argv,"--strip-pitch-mm");
  Need(!pitch.empty(),
       "strip pitch NOT provided by design report: explicitly supply "
       "--strip-pitch-mm for illustrative geometry");
  g.stripPitchMm=std::stod(pitch);
  g.stripWidthMm=Num(argc,argv,"--strip-width-mm",1.7);
  g.upperGapMm=Num(argc,argv,"--gap-plus-mm",3.);
  g.lowerGapMm=Num(argc,argv,"--gap-minus-mm",3.);
  g.stereoTan=Num(argc,argv,"--stereo-tan",.1);
  g.halfIndex=(int)Num(argc,argv,"--half-index",8.);
  Validate(g);
  const auto csv=out+"_uv_bare_slab_map.csv";
  const auto info=out+"_scope.txt";
  Need(!std::filesystem::exists(csv)&&!std::filesystem::exists(info),
       "refusing to overwrite existing strip maps");
  std::ofstream result(csv),meta(info);
  Need(bool(result)&&bool(meta),"cannot write geometry output directory");
  meta<<std::setprecision(15)
    <<"NA60+/DiCE candidate two-cathode U/V stereo layout\n"
    <<"SOURCE: Milov 2026-08-30 readout design, Sec 4.2 & 6\n"
    <<"strip_width_mm="<<g.stripWidthMm<<"\n"
    <<"strip_pitch_mm="<<g.stripPitchMm
    <<" (EXPLICIT ILLUSTRATIVE INPUT; report does not fix pitch)\n"
    <<"stereo_tan="<<g.stereoTan<<" (U and V have opposite slopes)\n"
    <<"upper_gap_mm="<<g.upperGapMm<<"\n"
    <<"lower_gap_mm="<<g.lowerGapMm<<"\n"
    <<"Wires run along w; U/V axes rotate relative to transverse u.\n"
    <<"BARE SLAB DIRICHLET LAPLACE EXACT; NO ANODE WIRES in weighting field.\n"
    <<"NO BOARD PIXELS, PCB DIELECTRICS, REAL STRIP ENDS OR READOUT.\n"
    <<"NOT the physical 3D MWPC weighting solution; cannot predict real signals.\n";
  result<<std::setprecision(15)
    <<"electrode,strip_index,u_mm,v_mm,w_mm,normal_strip_coordinate_mm,"
      "weighting_potential_unitless,all_included_strip_sum,"
      "bare_whole_plane_potential\n";
  const double H=g.Height();
  for(bool upper:{true,false})for(int strip=-g.halfIndex;strip<=g.halfIndex;++strip)
   for(int iv=1;iv<=5;++iv)for(int iu=-3;iu<=3;++iu)
    for(int iw=-6;iw<=6;++iw){
      const double u=iu;
      const double v=-g.lowerGapMm+H*iv/6.;
      const double w=iw;
      const double phi=g.Phi(upper,strip,u,v,w);
      const double all=g.Sum(upper,u,v,w);
      const double whole=g.ToUpperDistance(upper,v)/H;
      result<<(upper?"U_plus":"V_minus")<<","<<strip<<","
            <<u<<","<<v<<","<<w<<","<<g.S(upper,u,w)<<","
            <<phi<<","<<all<<","<<whole<<"\n";
    }
  std::cout<<"PASS candidate U/V stereo orientation and bare-slab boundary checks\n"
           <<"map="<<csv<<" scope="<<info<<"\n"
           <<"IMPORTANT: geometry/weighting BASELINE, NOT wire-inclusive "
           <<"physical-strip weighting or signal simulation.\n";
  return 0;
 }catch(const std::exception&e){
   std::cerr<<e.what()<<"\n";return 2;
 }
}

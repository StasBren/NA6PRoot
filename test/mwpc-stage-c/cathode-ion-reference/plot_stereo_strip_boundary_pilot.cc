// Numerical/analytic validation of candidate simple two-sided stereo-strip
// coordinates and a WIRE-FREE, two-parallel-plate boundary weighting field.
// Not a physically complete cathode weighting field for an MWPC (no anodes).
// No Garfield avalanche, RKF current, or any electronics.
#include "StereoStripBoundaryApprox.hh"
#include <TCanvas.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
namespace {
namespace ss=NA60Stereo;
std::string Arg(int n,char**v,const std::string&name){
 for(int i=1;i+1<n;++i)if(name==v[i])return v[i+1];
 return "";
}
double Num(int n,char**v,const std::string&name,double fallback){
 const auto value=Arg(n,v,name);
 return value.empty()?fallback:std::stod(value);
}
void Close(double value,double expected,double tol,const std::string&msg){
 ss::Check(std::isfinite(value)&&std::abs(value-expected)<=tol,
           msg+" observed="+std::to_string(value)+
           " expected="+std::to_string(expected));
}
void UnitTests(const ss::Design&d) {
 using ss::Side;
 ss::Validate(d);
 const double g=d.gapMm,m=d.tanStereo;
 // Upper electrode boundary, opposite cathode, neighboring strips.
 Close(ss::Potential(d,Side::Plus,0,0.,+g,0.),1.,1.e-14,
       "upper strip center at electrode");
 Close(ss::Potential(d,Side::Minus,0,0.,+g,0.),0.,1.e-14,
       "lower strip cannot be excited at upper cathode");
 Close(ss::Potential(d,Side::Plus,0,0.,+g,2.*d.stripPitchMm),
       0.,1.e-14,"non-target strip region must be zero at boundary");
 // Plus and minus orientations are opposite. This identity reflects both
 // u and v; reflecting only v would incorrectly align the strip geometry.
 for(const auto& pos:
     {std::array<double,3>{0.,0.,0.},
      {2.3,0.75,-1.1},{-3.7,-1.2,0.6}}){
  const auto [u,v,w]=pos;
  Close(ss::Potential(d,Side::Plus,0,u,v,w),
        ss::Potential(d,Side::Minus,0,-u,-v,w),1.e-13,
        "two-sided mirrored potential");
 }
 Close(ss::StripCoordinateMm(Side::Plus,10.,10.*m,d),
       0.,1.e-13,"plus stereo strip orientation");
 Close(ss::StripCoordinateMm(Side::Minus,10.,-10.*m,d),
       0.,1.e-13,"minus stereo strip orientation");
 if(m>0.){
   const double u=12.3,w=-2.4;
   const double sPlus=w-m*u,sMinus=w+m*u;
   Close(ss::DecodeUmm(sPlus,sMinus,d),u,1.e-13,
         "stereo measured-coordinate decoding u");
   Close(ss::DecodeWmm(sPlus,sMinus),w,1.e-13,
         "stereo measured-coordinate decoding w");
 }
 // In the artificial zero-gap, wire-free, two-parallel-plane model, the
 // sum of ALL contiguous-strip weighting potentials is the whole-plate
 // linear potential. For a wire-aware MWPC these values differ.
 if(std::abs(d.stripPitchMm-d.stripWidthMm)<1.e-12){
  for(double v:{-1.2,0.,1.8})
   for(auto side:{Side::Plus,Side::Minus}){
    double total=0.;
    for(int i=-45;i<=45;++i)
      total+=ss::Potential(d,side,i,1.3,v,0.7);
    const double whole=(side==Side::Plus?(v+g)/(2.*g):(g-v)/(2.*g));
    Close(total,whole,1.e-8,
          "sum of ideal contiguous strips vs wire-free full-plane potential");
   }
 }
 // Harmonicity is the actual 3D Laplacian in u,v,w, including stereo slope.
 const double u=2.1,v=0.9,w=0.4,h=0.01;
 for(auto side:{Side::Plus,Side::Minus}){
  auto p=[&](double x,double y,double z){
    return ss::Potential(d,side,0,x,y,z);
  };
  const double o=p(u,v,w);
  const double lap=(p(u+h,v,w)+p(u-h,v,w)+
    p(u,v+h,w)+p(u,v-h,w)+p(u,v,w+h)+p(u,v,w-h)-6.*o)/(h*h);
  Close(lap,0.,1.e-4,
        "3D Laplace equation away from conductor boundaries");
  auto field=ss::WeightingField(d,side,0,u,v,w);
  ss::Check(std::isfinite(field.eu)&&std::isfinite(field.ev)&&
            std::isfinite(field.ew),"weighting E gradient invalid");
 }
 std::cout<<"PASS stereo projected-coordinate, reflecting, electrode boundary, "
     "3D Laplace and wire-free strip partition numerical controls\n";
}
}
int main(int argc,char**argv){
 try {
 ss::Design design;
 design.gapMm=Num(argc,argv,"--gap-mm",3.);
 design.stripWidthMm=Num(argc,argv,"--strip-width-mm",1.7);
 design.stripPitchMm=Num(argc,argv,"--strip-pitch-mm",1.7);
 design.tanStereo=Num(argc,argv,"--tan-stereo",0.1);
 const auto prefix=Arg(argc,argv,"--out-prefix");
 ss::Validate(design);
 UnitTests(design);
 if(prefix.empty()) {
   std::cout<<"PASS geometry self-test only; no output maps requested\n";
   return 0;
 }
 for(auto ext:{"_map.csv","_geometry.root","_potential.png"}){
   ss::Check(!std::filesystem::exists(prefix+ext),
             "refusing overwrite output "+prefix+ext);
 }
 std::ofstream csv(prefix+"_map.csv");
 ss::Check(bool(csv),"cannot write stereo map CSV");
 csv<<std::setprecision(15)<<"u_mm,w_mm,"
    "U_plus_strip0_potential_v_plus1p5mm,"
    "V_minus_strip0_potential_v_minus1p5mm,"
    "U_plus_strip0_potential_at_anode_plane,"
    "V_minus_strip0_potential_at_anode_plane\n";
 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
 TFile root((prefix+"_geometry.root").c_str(),"RECREATE",
       "Ideal wire-free two-sided stereo-strip geometry and Laplace controls",4);
 ss::Check(!root.IsZombie(),"cannot create stereo pilot ROOT");
 TNamed limitations("limitations",
  "NA60+/DiCE preferred candidate under discussion: 1.7mm strip width "
  "measured ALONG anode wire, tan stereo angle +/-0.1, 3+3mm gaps. "
  "This test assumes contiguous strips (pitch=width); production "
  "interstrip copper gaps, PCB dielectric, traces, finite module edges, "
  "and ALL anode wires are ABSENT from this illustrative electrostatic "
  "parallel-plate solution. Do NOT derive physical MWPC strip currents "
  "from this test yet. It is not an electronics simulation.");
 limitations.Write();
 const double range=8.;
 TH2D uPlus("U_plus_potential",
       "Upper (+v) candidate stereo strip 0, v=+1.5 mm;u across wires (mm);w along wires (mm)",
       101,-range,range,101,-range,range);
 TH2D vMinus("V_minus_potential",
       "Lower (-v) candidate stereo strip 0, v=-1.5 mm;u across wires (mm);w along wires (mm)",
       101,-range,range,101,-range,range);
 TH2D uAtAnode("U_plus_anodeplane",
       "Upper candidate strip 0, wire-plane cross-section v=0;u across wires (mm);w along wires (mm)",
       101,-range,range,101,-range,range);
 TH2D vAtAnode("V_minus_anodeplane",
       "Lower candidate strip 0, wire-plane cross-section v=0;u across wires (mm);w along wires (mm)",
       101,-range,range,101,-range,range);
 for(int iu=1;iu<=101;++iu)
  for(int iw=1;iw<=101;++iw){
   const double u=uPlus.GetXaxis()->GetBinCenter(iu);
   const double w=uPlus.GetYaxis()->GetBinCenter(iw);
   const double up=ss::Potential(design,ss::Side::Plus,0,u,+1.5,w);
   const double down=ss::Potential(design,ss::Side::Minus,0,u,-1.5,w);
   const double up0=ss::Potential(design,ss::Side::Plus,0,u,0.,w);
   const double down0=ss::Potential(design,ss::Side::Minus,0,u,0.,w);
   uPlus.SetBinContent(iu,iw,up);
   vMinus.SetBinContent(iu,iw,down);
   uAtAnode.SetBinContent(iu,iw,up0);
   vAtAnode.SetBinContent(iu,iw,down0);
   csv<<u<<","<<w<<","<<up<<","<<down<<","<<up0<<","<<down0<<"\n";
  }
 root.cd();
 for(auto h:{&uPlus,&vMinus,&uAtAnode,&vAtAnode})h->Write();
 TCanvas canvas("stereo_potential",
     "Wire-free candidate stereo strip potential, no anode screening",1500,1100);
 canvas.Divide(2,2);
 canvas.cd(1);uPlus.Draw("COLZ");
 canvas.cd(2);vMinus.Draw("COLZ");
 canvas.cd(3);uAtAnode.Draw("COLZ");
 canvas.cd(4);vAtAnode.Draw("COLZ");
 canvas.cd(0);canvas.Modified();canvas.Update();
 root.cd();canvas.Write();canvas.SaveAs((prefix+"_potential.png").c_str());
 root.Close();csv.close();
 std::cout<<"PASS candidate strip coordinate and wire-free weighting-boundary "
   "pilot. width along wire="<<design.stripWidthMm
   <<" mm, pitch along wire="<<design.stripPitchMm
   <<" mm, tanStereo="<<design.tanStereo
   <<", gap=+"<<design.gapMm<<"/-"<<design.gapMm<<" mm\n"
   <<"Outputs "<<prefix<<"_map.csv, _geometry.root, _potential.png\n"
   <<"CRITICAL NOT anode-wire-aware field, NOT a physical strip signal.\n";
 return 0;
 }catch(const std::exception&e){
  std::cerr<<"STEREO GEOMETRY PILOT FAILED: "<<e.what()<<"\n";
  return 2;
 }
}

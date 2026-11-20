// NA60+/DiCE: actual grounded thin-anode-wire shielding in a candidate
// 3D U/V stereo-strip weighting potential (no induced-current simulation).
// Analytic strip base + boundary line-element MFS/collocation correction.
// Validation: zero-potential anode surfaces, cathode Dirichlet values,
// plus/minus symmetry, strip partition, Garfield++ complete-plane
// reference, optional independent wire-segment size/span convergence.
#include "WireAwareStereoWeighting.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
namespace {
namespace ss=NA60Stereo;
namespace ws=NA60WireStrip;
std::string Arg(int n,char**v,const std::string&key) {
 for(int i=1;i+1<n;++i)if(key==v[i])return v[i+1];
 return "";
}
bool Flag(int n,char**v,const std::string&key) {
 for(int i=1;i<n;++i)if(key==v[i])return true;
 return false;
}
double Numeric(int n,char**v,const std::string&key,double def) {
 const auto value=Arg(n,v,key);
 return value.empty()?def:std::stod(value);
}
void CheckClose(double a,double b,double tol,const std::string&name) {
 ss::Check(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<tol,
           name+" mismatch: "+std::to_string(a)+" versus "+std::to_string(b));
}
struct Point {double u=0,v=0,w=0;};
const std::array<Point,6> sample{{
 {0.70,0.35,0.20}, {-1.10,-0.80,0.65}, {2.30,1.20,-0.50},
 {-2.10,-1.25,0.15}, {4.60,0.55,-0.45}, {0.65,-0.55,1.10}
}};
struct Results {
 double wirePlus=0,wireMinus=0,wireWhole=0;
 double maxReflection=0,maxPartition=0,maxWholeGarfield=0;
 double maxSymmetryU=0;
};
Garfield::ComponentAnalyticField GarfieldPlaneReference(
       const ws::Settings&g) {
 Garfield::ComponentAnalyticField field;
 for(int i=-(g.nWires-1)/2;i<=(g.nWires-1)/2;++i) {
   // Geometry uses cm in Garfield++, mm in our 3D line-element solver.
   field.AddWire(0.1*i*g.wirePitchMm,0.,
                 0.2*g.wireRadiusMm,2100.,"anode");
 }
 field.AddPlaneY(-0.1*g.strip.gapMm,0.,"cathode_minus");
 field.AddPlaneY(+0.1*g.strip.gapMm,0.,"cathode_plus");
 field.SetMagneticField(0.,0.,0.);
 field.AddReadout("cathode_plus");
 field.AddReadout("cathode_minus");
 return field;
}
Results Validate(const ws::Model&model,
      const ws::Solution&plus,const ws::Solution&minus,
      const ws::Solution&wholePlus,const ws::Solution&wholeMinus,
      std::ostream&metric) {
 const auto&g=model.Geometry();
 Results res;
 res.wirePlus=model.WireResidual(plus,3.);
 res.wireMinus=model.WireResidual(minus,3.);
 res.wireWhole=model.WireResidual(wholePlus,3.);
 // Wire *surface* is only approximated by line monopoles and angular
 // collocation. This is a measurable, deliberately non-zero residual.
 ss::Check(res.wirePlus<0.025&&res.wireMinus<0.025&&
           res.wireWhole<0.025,
           "grounded anode cylinder residual exceeds pilot tolerance");
 const double gap=g.strip.gapMm;
 for(auto side:{ss::Side::Plus,ss::Side::Minus}){
   const auto&sol=(side==ss::Side::Plus?plus:minus);
   const double vHere=(side==ss::Side::Plus?gap:-gap);
   const double vOther=-vHere;
   CheckClose(model.Potential(sol,0.,vHere,0.),1.,1.e-12,
              "selected electrode strip center Dirichlet");
   CheckClose(model.Potential(sol,0.,vOther,0.),0.,1.e-12,
              "other cathode Dirichlet");
   CheckClose(model.Potential(sol,0.,vHere,2*g.strip.stripPitchMm),
              0.,1.e-12,"neighbor electrode strip boundary");
 }
 auto gf=GarfieldPlaneReference(g);
 const auto sumTest=std::array<Point,4>{{
   {0.4,0.75,0.2},{1.5,-0.65,0.8},
   {-1.2,1.1,-0.7},{2.7,0.4,0.3}
 }};
 for(const auto&p:sample) {
   const double pUp=model.Potential(plus,p.u,p.v,p.w);
   const double pDn=model.Potential(minus,-p.u,-p.v,p.w);
   res.maxReflection=std::max(res.maxReflection,std::abs(pUp-pDn));
   const double qUp=model.Potential(wholePlus,p.u,p.v,p.w);
   const double qDn=model.Potential(wholeMinus,p.u,p.v,p.w);
   // Two full-plane weighting potentials must sum to something <=1
   // in the gas (the grounded anodes share some weighting-potential sum).
   ss::Check(pUp>-0.03&&pUp<1.03&&qUp>-0.03&&qUp<1.03 &&
             qDn>-0.03&&qDn<1.03,
             "nonphysical electrostatic weighting-potential range");
   const double garP=gf.WeightingPotential(.1*p.u,.1*p.v,.1*p.w,
                                           "cathode_plus");
   const double garM=gf.WeightingPotential(.1*p.u,.1*p.v,.1*p.w,
                                           "cathode_minus");
   ss::Check(std::isfinite(garP)&&std::isfinite(garM),
             "Garfield analytic whole-plane weighting returned nonfinite");
   res.maxWholeGarfield=std::max(res.maxWholeGarfield,
       std::max(std::abs(qUp-garP),std::abs(qDn-garM)));
   metric<<"garfield_reference,upper,u="<<p.u<<";v="<<p.v<<";w="<<p.w
         <<","<<qUp<<","<<garP<<","<<qUp-garP<<"\n";
   metric<<"garfield_reference,lower,u="<<p.u<<";v="<<p.v<<";w="<<p.w
         <<","<<qDn<<","<<garM<<","<<qDn-garM<<"\n";
 }
 // Partition by ideal contiguous upper strips: the sum of strip 3D
 // potentials MUST approximate the fully wire-screened upper plate
 // solution in the local region. No individual-strip response physics
 // is claimed outside the finite wire-source domain.
 std::vector<ws::Solution> neighbors;
 for(int strip=-11;strip<=11;++strip)
   neighbors.push_back(model.Solve(ss::Side::Plus,
        ws::Electrode::SingleStrip,strip));
 for(const auto&p:sumTest) {
   double sum=0.;
   for(const auto&strip:neighbors)
     sum+=model.Potential(strip,p.u,p.v,p.w);
   const double whole=model.Potential(wholePlus,p.u,p.v,p.w);
   res.maxPartition=std::max(res.maxPartition,std::abs(sum-whole));
   metric<<"strip_partition,upper,u="<<p.u<<";v="<<p.v<<";w="<<p.w
         <<","<<sum<<","<<whole<<","<<sum-whole<<"\n";
 }
 ss::Check(res.maxReflection<0.0001,
           "opposite cathode stereo reflection symmetry broken");
 ss::Check(res.maxPartition<0.01,
           "sum of contiguous strips does not reproduce screened whole plate");
 metric<<"quality,wire_surface_residual_plus,central,0,"
       <<res.wirePlus<<","<<res.wirePlus<<"\n";
 metric<<"quality,wire_surface_residual_minus,central,0,"
       <<res.wireMinus<<","<<res.wireMinus<<"\n";
 metric<<"quality,wire_surface_residual_whole,central,0,"
       <<res.wireWhole<<","<<res.wireWhole<<"\n";
 metric<<"quality,reflection_error,central,0,"
       <<res.maxReflection<<","<<res.maxReflection<<"\n";
 metric<<"quality,strip_sum_error,central,0,"
       <<res.maxPartition<<","<<res.maxPartition<<"\n";
 metric<<"quality,max_garfield_whole_plane_deviation,central,0,"
       <<res.maxWholeGarfield<<","<<res.maxWholeGarfield<<"\n";
 return res;
}
struct Convergence {
 double maxStepDifference=0,maxExtentDifference=0;
 int refinedUnknowns=0,wideUnknowns=0;
};
Convergence Refine(const ws::Settings&g,
                   const ws::Model&base,
                   const ws::Solution&baseU,
                   const ws::Solution&baseWhole,
                   std::ostream&metric){
 ws::Settings small=g;
 small.segmentStepMm*=0.5;
 ws::Settings wider=g;
 wider.halfSpanWmm+=4.;
 ws::Validate(small);
 ws::Validate(wider);
 ws::Model refined(small),expanded(wider);
 const auto refinedU=refined.Solve(ss::Side::Plus,
                        ws::Electrode::SingleStrip,0);
 const auto refinedWhole=refined.Solve(ss::Side::Plus,
                        ws::Electrode::WholeCathode);
 const auto expandedU=expanded.Solve(ss::Side::Plus,
                        ws::Electrode::SingleStrip,0);
 const auto expandedWhole=expanded.Solve(ss::Side::Plus,
                        ws::Electrode::WholeCathode);
 Convergence out;
 out.refinedUnknowns=refined.Unknowns();
 out.wideUnknowns=expanded.Unknowns();
 for(const auto&p:sample){
  const double oldU=base.Potential(baseU,p.u,p.v,p.w);
  const double newU=refined.Potential(refinedU,p.u,p.v,p.w);
  const double wideU=expanded.Potential(expandedU,p.u,p.v,p.w);
  const double oldWhole=base.Potential(baseWhole,p.u,p.v,p.w);
  const double newWhole=refined.Potential(refinedWhole,p.u,p.v,p.w);
  const double wideWhole=expanded.Potential(expandedWhole,p.u,p.v,p.w);
  out.maxStepDifference=std::max(out.maxStepDifference,
       std::max(std::abs(oldU-newU),std::abs(oldWhole-newWhole)));
  out.maxExtentDifference=std::max(out.maxExtentDifference,
       std::max(std::abs(oldU-wideU),std::abs(oldWhole-wideWhole)));
  metric<<"refine_step,plus_strip0,u="<<p.u<<";v="<<p.v<<";w="<<p.w
        <<","<<oldU<<","<<newU<<","<<oldU-newU<<"\n";
  metric<<"extend_w,plus_strip0,u="<<p.u<<";v="<<p.v<<";w="<<p.w
        <<","<<oldU<<","<<wideU<<","<<oldU-wideU<<"\n";
  metric<<"refine_step,whole_plus,u="<<p.u<<";v="<<p.v<<";w="<<p.w
        <<","<<oldWhole<<","<<newWhole<<","<<oldWhole-newWhole<<"\n";
  metric<<"extend_w,whole_plus,u="<<p.u<<";v="<<p.v<<";w="<<p.w
        <<","<<oldWhole<<","<<wideWhole<<","<<oldWhole-wideWhole<<"\n";
 }
 ss::Check(out.maxStepDifference<0.03&&out.maxExtentDifference<0.03,
           "large 3D line-element segmentation/w-extent dependence");
 return out;
}
}
int main(int argc,char**argv){
 try {
  ws::Settings g;
  g.strip.gapMm=Numeric(argc,argv,"--gap-mm",3.);
  g.strip.stripWidthMm=Numeric(argc,argv,"--strip-width-mm",1.7);
  g.strip.stripPitchMm=Numeric(argc,argv,"--strip-pitch-mm",1.7);
  g.strip.tanStereo=Numeric(argc,argv,"--tan-stereo",0.1);
  g.wirePitchMm=Numeric(argc,argv,"--wire-pitch-mm",3.);
  g.wireRadiusMm=Numeric(argc,argv,"--wire-radius-mm",0.015);
  g.nWires=(int)Numeric(argc,argv,"--wires",9.);
  g.halfSpanWmm=Numeric(argc,argv,"--halfspan-w-mm",10.);
  g.segmentStepMm=Numeric(argc,argv,"--segment-step-mm",1.);
  g.imagePairs=(int)Numeric(argc,argv,"--image-pairs",12.);
  const auto prefix=Arg(argc,argv,"--out-prefix");
  const bool refine=Flag(argc,argv,"--convergence");
  ws::Validate(g);
  if(!prefix.empty())for(const auto&ext:{"_checks.csv","_potential.csv",
           "_potential.root","_potential.png"}) {
    ss::Check(!std::filesystem::exists(prefix+ext),
       "refusing overwrite result "+prefix+ext);
  }
  ws::Model model(g);
  std::cout<<"3D line-element boundary matrix: unknowns="
     <<model.Unknowns()<<" image_pairs="<<g.imagePairs
     <<" per side, wire_span=[-"<<g.halfSpanWmm
     <<",+"<<g.halfSpanWmm<<"]mm\n";
  const auto up=model.Solve(ss::Side::Plus,ws::Electrode::SingleStrip,0);
  const auto up1=model.Solve(ss::Side::Plus,ws::Electrode::SingleStrip,1);
  const auto minus=model.Solve(ss::Side::Minus,ws::Electrode::SingleStrip,0);
  const auto wholeUp=model.Solve(ss::Side::Plus,ws::Electrode::WholeCathode);
  const auto wholeDown=model.Solve(ss::Side::Minus,ws::Electrode::WholeCathode);
  std::ostringstream measurements;
  measurements<<std::setprecision(15)
              <<"test,side,location,calculated,reference,difference\n";
  const auto metrics=Validate(model,up,minus,wholeUp,wholeDown,measurements);
  Convergence conv;
  if(refine)conv=Refine(g,model,up,wholeUp,measurements);
  std::cout<<std::setprecision(7)
     <<"WIRE-GROUNDED SURFACE residual U="<<metrics.wirePlus
     <<" V="<<metrics.wireMinus<<" whole="<<metrics.wireWhole<<"\n"
     <<"opposing_stereo_reflection_error="<<metrics.maxReflection
     <<" strip_sum_vs_wire_aware_plane_error="<<metrics.maxPartition<<"\n"
     <<"max_abs_whole_plane_difference_vs_Garfield_analytic="
     <<metrics.maxWholeGarfield<<"\n";
  if(refine)
    std::cout<<"CONVERGENCE step_max_deviation="<<conv.maxStepDifference
     <<" wire_span_max_deviation="<<conv.maxExtentDifference
     <<" refined_unknowns="<<conv.refinedUnknowns
     <<" extended_unknowns="<<conv.wideUnknowns<<"\n";
  if(!prefix.empty()) {
    std::ofstream checks(prefix+"_checks.csv");
    ss::Check(bool(checks),"cannot write checks CSV");
    checks<<measurements.str();checks.close();
    std::ofstream data(prefix+"_potential.csv");
    ss::Check(bool(data),"cannot write field map CSV");
    data<<std::setprecision(15)
       <<"u_mm,w_mm,v_mm,upper_strip0_potential,"
         "upper_strip1_potential,lower_strip0_potential,"
         "upper_whole_wire_aware_potential,"
         "upper_strip0_wire_free_potential\n";
    gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
    TFile root((prefix+"_potential.root").c_str(),"RECREATE",
      "NA60+ candidate stereo cathode strip weighting with thin grounded wires",4);
    ss::Check(!root.IsZombie(),"cannot make ROOT diagnostic output");
    TNamed caveat("scope",
      "3D wire-aware thin-wire line ELEMENT collocation pilot, anode cylinders "
      "grounded only approximately (surface residual quantified). "
      "Dirichlet cathodes via image Green, strip width/pitch along wires. "
      "Finite w wire correction source span and 4-angle cylinder average; "
      "needs convergence/Garfield whole-plane reference, realistic PCB/copper "
      "gaps, finite boundaries before physical strip currents. "
      "No electron/ion drift, no induced strip signals, no electronics.");
    caveat.Write();
    constexpr int n=45;
    const double rangeU=12.0,rangeW=7.0,v=0.30;
    TH2D U0("U0_3D_wire_aware",
        "Plus cathode strip 0, v=+0.30mm (wire-aware);u across wires (mm);w along wires (mm)",
        n,-rangeU,rangeU,n,-rangeW,rangeW);
    TH2D U1("U1_3D_wire_aware",
        "Plus cathode strip 1, v=+0.30mm (wire-aware);u across wires (mm);w along wires (mm)",
        n,-rangeU,rangeU,n,-rangeW,rangeW);
    TH2D V0("V0_3D_wire_aware",
        "Minus cathode strip 0, v=+0.30mm (wire-aware);u across wires (mm);w along wires (mm)",
        n,-rangeU,rangeU,n,-rangeW,rangeW);
    TH2D delta("U0_3D_wire_screening_correction",
        "Grounded anode wire correction to plus strip 0;u across wires (mm);w along wires (mm)",
        n,-rangeU,rangeU,n,-rangeW,rangeW);
    for(int i=1;i<=n;++i)for(int j=1;j<=n;++j){
      const double u=U0.GetXaxis()->GetBinCenter(i);
      const double w=U0.GetYaxis()->GetBinCenter(j);
      const double a=model.Potential(up,u,v,w);
      const double b=model.Potential(up1,u,v,w);
      const double c=model.Potential(minus,u,v,w);
      const double d=model.Potential(wholeUp,u,v,w);
      const double bare=ss::Potential(g.strip,ss::Side::Plus,0,u,v,w);
      U0.SetBinContent(i,j,a);
      U1.SetBinContent(i,j,b);
      V0.SetBinContent(i,j,c);
      delta.SetBinContent(i,j,a-bare);
      data<<u<<","<<w<<","<<v<<","<<a<<","<<b<<","<<c
          <<","<<d<<","<<bare<<"\n";
    }
    root.cd();U0.Write();U1.Write();V0.Write();delta.Write();
    TCanvas canvas("wire_aware_stereo_pilot",
      "Candidate stereo-strip weighting potentials with grounded thin-wire screening",1480,1050);
    canvas.Divide(2,2);
    canvas.cd(1);U0.Draw("COLZ");
    canvas.cd(2);U1.Draw("COLZ");
    canvas.cd(3);V0.Draw("COLZ");
    canvas.cd(4);delta.Draw("COLZ");
    canvas.cd(0);canvas.Modified();canvas.Update();
    root.cd();canvas.Write();canvas.SaveAs((prefix+"_potential.png").c_str());
    root.Close();data.close();
  }
  std::cout<<"PASS 3D thin-wire screening FIELD PILOT with boundary, "
     <<"symmetry and strip-partition checks. "
     <<"Whole-plane Garfield difference is a reported diagnostic, "
     <<"NOT yet accepted as production field accuracy.\n"
     <<"NO current, NO electron/ion drift, NO PCB dielectric, NO electronics.\n";
  return 0;
 }catch(const std::exception&e) {
  std::cerr<<"WIRE-AWARE STRIP PILOT FAILED: "<<e.what()<<"\n";
  return 2;
 }
}

// Stage C: near-wire *individual stereo-strip* gradient convergence.
// No Garfield, no avalanche drift, no electronics. Checks how actual
// U/V strip weighting fields change with thin-wire numerical parameters.
// This is an internal convergence study, NOT external physical validation.
#include "WireAwareStereoWeighting.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
namespace ws=NA60WireStrip;
namespace ss=NA60Stereo;
std::string Arg(int argc,char**argv,const std::string&k) {
 for(int i=1;i+1<argc;++i)if(k==argv[i])return argv[i+1];
 return "";
}
double Number(int argc,char**argv,const std::string&k,double def) {
 const auto a=Arg(argc,argv,k);
 return a.empty()?def:std::stod(a);
}
struct Quantity {
 double potential=0;
 ss::Field e{};
};
Quantity Eval(const ws::Model&m,const ws::Solution&s,
              double u,double v,double w,double h) {
 ss::Check(h>0&&std::abs(v)+h<m.Geometry().strip.gapMm,
           "invalid near-anode gradient step / cathode proximity");
 const auto pot=[&](double a,double b,double c){
   return m.Potential(s,a,b,c);
 };
 Quantity q;
 q.potential=pot(u,v,w);
 q.e.eu=(pot(u-h,v,w)-pot(u+h,v,w))/(2*h);
 q.e.ev=(pot(u,v-h,w)-pot(u,v+h,w))/(2*h);
 q.e.ew=(pot(u,v,w-h)-pot(u,v,w+h))/(2*h);
 ss::Check(std::isfinite(q.potential)&&std::isfinite(q.e.eu)&&
           std::isfinite(q.e.ev)&&std::isfinite(q.e.ew),
           "nonfinite individual strip weighting-gradient");
 return q;
}
double Norm2(const ss::Field&f){
 return f.eu*f.eu+f.ev*f.ev+f.ew*f.ew;
}
double Difference2(const ss::Field&a,const ss::Field&b){
 return (a.eu-b.eu)*(a.eu-b.eu)+
        (a.ev-b.ev)*(a.ev-b.ev)+
        (a.ew-b.ew)*(a.ew-b.ew);
}
struct Errors {
 int n=0;
 double maxPotentialDelta=0;
 double potentialSquared=0;
 double maxGradientDelta=0;
 double maxRelativeGradient=0;
 double gradientDeltaSquared=0;
 double baselineGradientSquared=0;
 double minBaselinePotential=std::numeric_limits<double>::infinity();
 double maxBaselineGradient=0;
};
struct Case {
 std::string name;
 ws::Settings geometry{};
};
struct Response {
 ws::Model model;
 ws::Solution upper,lower;
 explicit Response(const ws::Settings&g)
   :model(g),
     upper(model.Solve(ss::Side::Plus,ws::Electrode::SingleStrip,0)),
     lower(model.Solve(ss::Side::Minus,ws::Electrode::SingleStrip,0)){}
 const ws::Solution& Strip(int side)const{return side?lower:upper;}
};
std::vector<Case> Variations(const ws::Settings&base) {
 std::vector<Case> out;
 auto angle=base;
 angle.angularCollocationSamples=8;
 out.push_back({"angle_4to8",angle});
 auto segment=base;
 segment.segmentStepMm/=2.;
 out.push_back({"segment_half",segment});
 auto length=base;
 length.halfSpanWmm+=4.;
 out.push_back({"span_plus4mm",length});
 auto image=base;
 image.imagePairs+=8;
 out.push_back({"image_plus8",image});
 if(base.localHalfSpanWmm>0. && base.localStepMm>0.){
   auto local=base;
   local.localStepMm/=2.;
   out.push_back({"local_segment_half",local});
   auto outer=base;
   outer.segmentStepMm/=2.;
   out.push_back({"outer_segment_half",outer});
   auto span=base;
   span.localHalfSpanWmm+=1.;
   out.push_back({"local_span_plus1mm",span});
 }
 for(const auto&v:out)ws::Validate(v.geometry);
 return out;
}
}
int main(int argc,char**argv) {
 try {
  ws::Settings g;
  g.strip.gapMm=Number(argc,argv,"--gap-mm",3.);
  g.strip.stripWidthMm=Number(argc,argv,"--strip-width-mm",1.7);
  g.strip.stripPitchMm=Number(argc,argv,"--strip-pitch-mm",1.7);
  g.strip.tanStereo=Number(argc,argv,"--tan-stereo",0.1);
  g.wirePitchMm=Number(argc,argv,"--wire-pitch-mm",3.);
  g.wireRadiusMm=Number(argc,argv,"--wire-radius-mm",0.015);
  g.nWires=(int)Number(argc,argv,"--wires",9.);
  g.halfSpanWmm=Number(argc,argv,"--halfspan-w-mm",10.);
  g.segmentStepMm=Number(argc,argv,"--segment-step-mm",1.);
  g.localHalfSpanWmm=Number(argc,argv,"--local-halfspan-w-mm",0.);
  g.localStepMm=Number(argc,argv,"--local-step-mm",0.);
  g.imagePairs=(int)Number(argc,argv,"--image-pairs",12.);
  g.angularCollocationSamples=4;  // original baseline, fixed
  ws::Validate(g);
  const auto prefix=Arg(argc,argv,"--out-prefix");
  const auto onlyVariation=Arg(argc,argv,"--only-variation");
  ss::Check(!prefix.empty(),"usage: --out-prefix PREFIX [--only-variation segment_half]");
  if(!onlyVariation.empty()) {
    ss::Check(onlyVariation=="segment_half"||
              onlyVariation=="angle_4to8"||
              onlyVariation=="span_plus4mm"||
              onlyVariation=="image_plus8"||
              onlyVariation=="local_segment_half"||
              onlyVariation=="outer_segment_half"||
              onlyVariation=="local_span_plus1mm",
              "--only-variation must name an existing numerical perturbation");
  }
  for(const auto&suffix:{"_points.csv","_summary.csv"})
    ss::Check(!std::filesystem::exists(prefix+suffix),
       "refusing to overwrite convergence output "+prefix+suffix);
  auto variations=Variations(g);
  if(!onlyVariation.empty()) {
    variations.erase(std::remove_if(variations.begin(),variations.end(),
       [&](const Case&candidate){return candidate.name!=onlyVariation;}),
       variations.end());
    ss::Check(variations.size()==1,"unknown filtered convergence case");
  }
  std::cout<<"Building baseline 3D model, 4 wire-surface collocation angles, "
           <<"outer segments="<<g.segmentStepMm<<" mm"
           <<" central_halfspan="<<g.localHalfSpanWmm<<" mm"
           <<" central_step="<<g.localStepMm<<" mm ...\n";
  Response baseline(g);
  std::array<double,2> baseWireResidual{};
  for(int side=0;side<2;++side)
    baseWireResidual[side]=baseline.model.WireResidual(
          baseline.Strip(side),1.1);
  std::cout<<"baseline_unknowns="<<baseline.model.Unknowns()
           <<" baseline_U_wire_surface_residual="<<baseWireResidual[0]
           <<" baseline_V_wire_surface_residual="<<baseWireResidual[1]
           <<"\n";
  std::ofstream csv(prefix+"_points.csv"),summary(prefix+"_summary.csv");
  ss::Check(bool(csv)&&bool(summary),"cannot open strip-convergence CSVs");
  csv<<std::setprecision(16)
    <<"variation,side,wire_index,r_over_R,theta_rad,w_mm,"
      "baseline_potential,variant_potential,absolute_potential_difference,"
      "baseline_Eu_per_mm,baseline_Ev_per_mm,baseline_Ew_per_mm,"
      "variant_Eu_per_mm,variant_Ev_per_mm,variant_Ew_per_mm,"
      "gradient_absolute_difference_per_mm,"
      "baseline_gradient_norm_per_mm,local_relative_gradient_difference\n";
  summary<<std::setprecision(16)
    <<"variation,side,n_points,mean_potential_rms_difference,"
      "max_abs_potential_difference,max_abs_gradient_difference_per_mm,"
      "max_local_relative_gradient_difference,"
      "vector_gradient_relative_L2,minimum_baseline_abs_potential,"
      "max_baseline_gradient_per_mm,"
      "baseline_wire_surface_residual,"
      "variant_wire_surface_residual\n";

  // Independent U and V strip gradients; no whole-cathode reference here.
  // Each sample is above a wire cylinder at a finite radial clearance.
  // Include one w=1.1 mm location to probe stereo longitudinal derivatives.
  constexpr std::array<double,4> angles{
       ws::pi/8.,5.*ws::pi/8.,9.*ws::pi/8.,13.*ws::pi/8.};
  constexpr std::array<double,3> shells{1.10,2.,5.};
  double largestGradientL2=0.;
  for(const auto&variation:variations){
    std::cout<<"BEGIN "<<variation.name<<" ";
    std::cout.flush();
    Response test(variation.geometry);
    std::array<double,2> testWireResidual{};
    for(int side=0;side<2;++side)
      testWireResidual[side]=test.model.WireResidual(
        test.Strip(side),1.1);
    std::cout<<"variant_unknowns="<<test.model.Unknowns()
             <<" variant_U_wire_surface_residual="<<testWireResidual[0]
             <<" variant_V_wire_surface_residual="<<testWireResidual[1]
             <<" ";
    std::array<Errors,2> error;
    for(int index:{-1,0,1})
     for(double ratio:shells)
      for(double theta:angles)
       for(double w:{0.,1.1}) {
        const double r=ratio*g.wireRadiusMm;
        const double u=index*g.wirePitchMm+r*std::cos(theta);
        const double v=r*std::sin(theta);
        const double clearance=r-g.wireRadiusMm;
        const double h=std::min(0.0015,0.15*clearance);
        ss::Check(h>0&&h<clearance/2.,
            "near-anode finite-difference step crosses wire metal");
        for(int side=0;side<2;++side){
          const auto base=Eval(baseline.model,
                   baseline.Strip(side),u,v,w,h);
          const auto variant=Eval(test.model,test.Strip(side),u,v,w,h);
          const double pAbs=std::abs(base.potential-variant.potential);
          const double gradDelta=std::sqrt(Difference2(base.e,variant.e));
          const double gradNorm=std::sqrt(Norm2(base.e));
          ss::Check(std::isfinite(gradDelta)&&std::isfinite(gradNorm),
              "nonfinite strip gradient-convergence metric");
          auto& stat=error[side];
          stat.n++;
          stat.maxPotentialDelta=std::max(stat.maxPotentialDelta,pAbs);
          stat.potentialSquared+=pAbs*pAbs;
          stat.maxGradientDelta=std::max(stat.maxGradientDelta,gradDelta);
          stat.maxRelativeGradient=std::max(stat.maxRelativeGradient,
                gradDelta/std::max(gradNorm,1.e-30));
          stat.gradientDeltaSquared+=gradDelta*gradDelta;
          stat.baselineGradientSquared+=gradNorm*gradNorm;
          stat.minBaselinePotential=std::min(stat.minBaselinePotential,
                    std::abs(base.potential));
          stat.maxBaselineGradient=std::max(stat.maxBaselineGradient,
                    gradNorm);
          csv<<variation.name<<","
             <<(side?"V_minus":"U_plus")<<","
             <<index<<","<<ratio<<","<<theta<<","<<w<<","
             <<base.potential<<","<<variant.potential<<","<<pAbs<<","
             <<base.e.eu<<","<<base.e.ev<<","<<base.e.ew<<","
             <<variant.e.eu<<","<<variant.e.ev<<","<<variant.e.ew<<","
             <<gradDelta<<","<<gradNorm<<","
             <<gradDelta/std::max(gradNorm,1.e-30)<<"\n";
        }
       }
    for(int side=0;side<2;++side){
      const auto &x=error[side];
      ss::Check(x.n>0 && x.baselineGradientSquared>0,
                "missing valid near-wire strip gradients");
      const double L2=std::sqrt(x.gradientDeltaSquared/
                                x.baselineGradientSquared);
      largestGradientL2=std::max(largestGradientL2,L2);
      summary<<variation.name<<","
        <<(side?"V_minus":"U_plus")<<","<<x.n<<","
        <<std::sqrt(x.potentialSquared/x.n)<<","
        <<x.maxPotentialDelta<<","<<x.maxGradientDelta<<","
        <<x.maxRelativeGradient<<","<<L2<<","
        <<x.minBaselinePotential<<","<<x.maxBaselineGradient<<","
        <<baseWireResidual[side]<<","<<testWireResidual[side]<<"\n";
      std::cout<<(side?" V":" U")
        <<"_strip_gradient_relative_L2="<<L2
        <<(side?"\n":" ");
    }
  }
  csv.close();summary.close();
  std::cout<<"COMPLETED individual U/V strip near-anode numerical "
     "REFINEMENT COMPARISONS (not a claim of convergence).\n"
     <<"Maximum relative vector-gradient L2 change="<<largestGradientL2
     <<" (across requested perturbations and strip sides).\n";
  if(largestGradientL2>0.02)
    std::cout<<"NOT_CONVERGED: at least one numerical change exceeds 2% "
                "in near-anode U/V strip gradient L2; do not use fields "
                "for physical strip currents.\n";
  else
    std::cout<<"Below 2% for these selected perturbations ONLY; "
                "no independent physical strip validation.\n";
  std::cout<<"GROUNDING CAUTION: residual wire surface potential is an "
               "independent nonzero boundary error. Small gradient changes "
               "alone cannot certify strip physics.\n";
  std::cout<<"No transport or electronics. See "
      <<prefix<<"_points.csv and "<<prefix<<"_summary.csv\n";
  return 0;
 }catch(const std::exception&e){
   std::cerr<<"STEREO STRIP GRADIENT CONVERGENCE FAILED: "
            <<e.what()<<"\n";
   return 2;
 }
}

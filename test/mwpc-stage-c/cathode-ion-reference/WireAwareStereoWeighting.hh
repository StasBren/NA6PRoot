// NA60+/DiCE Stage C: 3D thin-anode-wire Dirichlet weighting-potential
// COLLOCATION MODEL for candidate two-sided stereo cathode strips.
// Base harmonic solution: analytic strip electrode between two complete plates.
// Grounded wire correction: discretized line charges along each 30um anode wire,
// with parallel-plate DIRICHLET Green function (infinite image pairs).
// All anodes are held at zero weighting potential via surface collocation;
// no ion or electron transport/electronics is computed in this class.
//
// IMPORTANT APPROXIMATIONS:
//   * Thin-wire constant line density per w segment, rather than exact curved
//     cylinder BEM; enforce mean of 4 angles on the circular wire surface.
//   * Wire sources extend over FINITE [-W,+W] along the physically long wires.
//     Output needs both segment-size and W-extent convergence tests.
//   * Infinite planar cathodes with ideal contiguous strips; no PCB dielectric,
//     interstrip gaps, guard layers, or physical finite module boundaries.
//   * Garfield's wire-aware analytic WHOLE-plane weighting potential is an
//     independent reference, not guaranteed by this approximate method.
// Coordinates and units are in mm throughout.
#pragma once
#include "StereoStripBoundaryApprox.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace NA60WireStrip {
namespace st=NA60Stereo;
constexpr double pi=st::PI;
struct Settings {
 st::Design strip{};
 double wirePitchMm=3.;
 double wireRadiusMm=0.015; // 30 micrometre diameter (NOT radius)
 int nWires=9;
 double halfSpanWmm=10.; // controlled numerical correction extent
 double segmentStepMm=1.;
 int imagePairs=12; // sum image indices -M..M
 int angularCollocationSamples=4; // average equally spaced wire-cylinder angles
};
struct Node {
 double u=0,w=0,zmin=0,zmax=0;
};
inline void Validate(const Settings&g) {
 st::Validate(g.strip);
 st::Check(g.nWires>=1&&g.nWires<=15&&g.nWires%2==1,
     "use odd 1..15 number of grounded anode wires");
 st::Check(g.wirePitchMm>g.wireRadiusMm*8. &&
     g.wireRadiusMm>0.003&&g.wireRadiusMm<=0.04,
     "wire pitch/radius outside validated numerical range");
 st::Check(g.halfSpanWmm>=4.&&g.halfSpanWmm<=30. &&
     g.segmentStepMm>=0.25&&g.segmentStepMm<=2. &&
     g.imagePairs>=6&&g.imagePairs<=64 &&
     (g.angularCollocationSamples==4 ||
      g.angularCollocationSamples==8 ||
      g.angularCollocationSamples==16),
     "invalid wire segmentation/image parameters");
 const double n=2.*g.halfSpanWmm/g.segmentStepMm;
 st::Check(std::abs(n-std::round(n))<1.e-8,
           "w length must be integer number of wire segments");
 st::Check(g.nWires*(int)std::lround(n)<=800,
           "dense wire boundary matrix exceeds 800 unknowns");
}
// Integral along [a,b] of the two-plate Dirichlet point-charge Green
// function for a unit line density at fixed (u0,v0=0), with v in [-g,g].
// The 3D Helmholtz/Laplace image series is the potential that is exactly
// zero on both infinite cathodes (when summed over all m).
inline double SegmentGreen(const Settings&g,double u,double v,double w,
                           const Node&source) {
 const double L=2.*g.strip.gapMm;
 const double du=u-source.u;
 double sum=0.;
 for(int m=-g.imagePairs;m<=g.imagePairs;++m){
   const double shift=2.*m*L;
   const double rd=std::hypot(du,v+shift);
   const double ri=std::hypot(du,v+L+shift);
   // A sample point is allowed on the cylinder surface but never exactly
   // on a charge line inside the conductor.
   st::Check(rd>1.e-11&&ri>1.e-11,"singular line-source evaluation");
   auto integral=[&](double radius){
     return std::asinh((source.zmax-w)/radius)-
            std::asinh((source.zmin-w)/radius);
   };
   sum+=integral(rd)-integral(ri);
 }
 return sum/(4.*pi);
}
enum class Electrode {SingleStrip,WholeCathode};
struct Solution {
 st::Side side=st::Side::Plus;
 Electrode electrode=Electrode::SingleStrip;
 int strip=0;
 std::vector<double> coefficient;
};
class Model {
 Settings g_;
 std::vector<Node> nodes_;
 std::vector<double> lu_;
 std::vector<int> swaps_;
 // Equally spaced angle quadrature around each 30um wire cylinder.
 // n=4 exactly reproduces the original four-angle mean equations.
 // Use 8 and 16 as explicit angular-discretization sensitivity controls.
 template<class Op>
 void ForCylinderAngles(Op&&fn) const {
   for(int a=0;a<g_.angularCollocationSamples;++a){
     const double theta=2.*pi*a/g_.angularCollocationSamples;
     fn(g_.wireRadiusMm*std::cos(theta),
        g_.wireRadiusMm*std::sin(theta));
   }
 }
 double Base(st::Side side,Electrode kind,int strip,
             double u,double v,double w) const {
   if(kind==Electrode::WholeCathode)
     return (side==st::Side::Plus?
                (v+g_.strip.gapMm)/(2.*g_.strip.gapMm):
                (g_.strip.gapMm-v)/(2.*g_.strip.gapMm));
   return st::Potential(g_.strip,side,strip,u,v,w);
 }
 void Factorize() {
   const int n=(int)nodes_.size();
   lu_.assign((std::size_t)n*n,0.);
   swaps_.assign(n,0);
   // Enforce zero mean potential on each circular wire surface at the
   // segment midpoint, sampling 4 quadrants of its cylinder.
   for(int i=0;i<n;++i)
     for(int j=0;j<n;++j){
       double sum=0.;
       ForCylinderAngles([&](double du,double dv){
         sum+=SegmentGreen(g_,nodes_[i].u+du,dv,nodes_[i].w,nodes_[j]);
       });
       lu_[(std::size_t)i*n+j]=sum/g_.angularCollocationSamples;
     }
   // Partial-pivot LU factorization; reusable for many different readout
   // strips and for complete-plane controls. One matrix for one geometry.
   for(int k=0;k<n;++k){
     int best=k;
     double maximum=std::abs(lu_[(std::size_t)k*n+k]);
     for(int i=k+1;i<n;++i){
       const double z=std::abs(lu_[(std::size_t)i*n+k]);
       if(z>maximum){maximum=z;best=i;}
     }
     st::Check(std::isfinite(maximum)&&maximum>1.e-12,
        "singular wire-collocation capacitance matrix");
     swaps_[k]=best;
     if(best!=k)for(int j=0;j<n;++j)
       std::swap(lu_[(std::size_t)k*n+j],lu_[(std::size_t)best*n+j]);
     const double pivot=lu_[(std::size_t)k*n+k];
     for(int i=k+1;i<n;++i){
       double& a=lu_[(std::size_t)i*n+k];
       a/=pivot;
       for(int j=k+1;j<n;++j)
         lu_[(std::size_t)i*n+j]-=a*lu_[(std::size_t)k*n+j];
     }
   }
 }
 std::vector<double> BackSolve(std::vector<double> rhs) const {
   const int n=(int)nodes_.size();
   st::Check((int)rhs.size()==n,"wire boundary solve RHS shape mismatch");
   for(int k=0;k<n;++k)
     if(swaps_[k]!=k)std::swap(rhs[k],rhs[swaps_[k]]);
   for(int i=0;i<n;++i)
     for(int j=0;j<i;++j)
       rhs[i]-=lu_[(std::size_t)i*n+j]*rhs[j];
   for(int i=n-1;i>=0;--i){
     for(int j=i+1;j<n;++j)
       rhs[i]-=lu_[(std::size_t)i*n+j]*rhs[j];
     rhs[i]/=lu_[(std::size_t)i*n+i];
     st::Check(std::isfinite(rhs[i]),"nonfinite wire-screening coefficient");
   }
   return rhs;
 }
public:
 explicit Model(Settings g):g_(std::move(g)){
   Validate(g_);
   const int segments=(int)std::lround(2.*g_.halfSpanWmm/g_.segmentStepMm);
   const int half=(g_.nWires-1)/2;
   nodes_.reserve((std::size_t)g_.nWires*segments);
   for(int wire=-half;wire<=half;++wire){
     const double u=wire*g_.wirePitchMm;
     for(int i=0;i<segments;++i){
       const double lo=-g_.halfSpanWmm+i*g_.segmentStepMm;
       nodes_.push_back({u,lo+0.5*g_.segmentStepMm,lo,
                        lo+g_.segmentStepMm});
     }
   }
   Factorize();
 }
 const Settings& Geometry()const{return g_;}
 int Unknowns()const{return (int)nodes_.size();}
 Solution Solve(st::Side side,Electrode kind,int index=0) const {
   std::vector<double> rhs;
   rhs.reserve(nodes_.size());
   for(const auto& node:nodes_){
     double mean=0.;
     ForCylinderAngles([&](double du,double dv){
       mean+=Base(side,kind,index,node.u+du,dv,node.w);
     });
     rhs.push_back(-mean/g_.angularCollocationSamples);
   }
   Solution out;
   out.side=side;out.electrode=kind;out.strip=index;
   out.coefficient=BackSolve(std::move(rhs));
   return out;
 }
 double Potential(const Solution&sol,double u,double v,double w) const {
   st::Check(std::isfinite(u)&&std::isfinite(v)&&std::isfinite(w)&&
       std::abs(v)<=g_.strip.gapMm,"evaluation outside cathode planes");
   st::Check(sol.coefficient.size()==nodes_.size(),"solution geometry mismatch");
   const double base=Base(sol.side,sol.electrode,sol.strip,u,v,w);
   if(std::abs(v)>=g_.strip.gapMm-1.e-12) {
     // The exact infinite-image Green function is 0 on plates; finite
     // image sums have tiny boundary residuals, so impose known Dirichlet
     // conductor values exactly at the conducting plane.
     return base;
   }
   double correction=0.;
   for(std::size_t i=0;i<nodes_.size();++i)
     correction+=sol.coefficient[i]*SegmentGreen(g_,u,v,w,nodes_[i]);
   return base+correction;
 }
 // For field differentiation choose h and observation point away from wire
 // cylinders and both conducting planes. E_w=-grad(phi_w), 1/mm.
 st::Field WeightingField(const Solution&sol,double u,double v,double w,
                          double h=0.02) const {
   st::Check(h>0. && h<0.1 &&
             std::abs(v)+h<g_.strip.gapMm,"gradient near cathode boundary");
   auto p=[&](double a,double b,double c){return Potential(sol,a,b,c);};
   return {(p(u-h,v,w)-p(u+h,v,w))/(2*h),
           (p(u,v-h,w)-p(u,v+h,w))/(2*h),
           (p(u,v,w-h)-p(u,v,w+h))/(2*h)};
 }
 // Surface averages are the enforced collocation equations, but also
 // inspect 8 azimuthal samples between the collocation directions.
 // Only check central w sites to avoid finite-z end corrections.
 double WireResidual(const Solution&sol,double maxAbsWmm=4.) const {
   double largest=0.;
   const int half=(g_.nWires-1)/2;
   for(int wire=-half;wire<=half;++wire)
     for(double w:{-maxAbsWmm,0.,maxAbsWmm}){
       for(int k=0;k<12;++k){
         const double theta=2.*pi*k/12.+0.13;
         const double u=wire*g_.wirePitchMm+
                         g_.wireRadiusMm*std::cos(theta);
         const double v=g_.wireRadiusMm*std::sin(theta);
         const double potential=Potential(sol,u,v,w);
         st::Check(std::isfinite(potential),
                   "nonfinite potential at grounded wire surface");
         largest=std::max(largest,std::abs(potential));
       }
     }
   return largest;
 }
};
} // namespace NA60WireStrip

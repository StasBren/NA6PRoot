// Candidate NA60+ simple U/V cathode stereo strip geometry, *wire-free*
// analytic parallel-plate weighting potential (LAPLACE BOUNDARY MODEL).
// WARNING: anode wires, dielectric PCB, copper gaps, finite strip ends,
// edge guards and traces are NOT present in this electrostatic approximation.
// It is NOT a validated MWPC strip weighting field and is NOT connected
// to DriftLineRKF induced-current calculations.
#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
namespace NA60Stereo {
constexpr double PI=3.14159265358979323846264338327950288;
struct Design {
  double gapMm=3.;       // distance from anode wire plane to each cathode
  double stripWidthMm=1.7; // width MEASURED along physical anode wire direction
  double stripPitchMm=1.7; // IDEAL fully tiled copper, no interstrip gap
  double tanStereo=0.1;    // strip axes at pi/2 +/- atan(0.1) to anode wires
};
enum class Side {Plus,Minus};
struct Field {
  double eu=0,ev=0,ew=0; // weighting E-field components, 1/mm
};
inline void Check(bool b,const std::string&reason) {
 if(!b)throw std::runtime_error("Stereo readout boundary model: "+reason);
}
inline void Validate(const Design&d) {
 Check(d.gapMm>=2. && d.gapMm<=5. &&
       d.stripWidthMm>0. && d.stripPitchMm>=d.stripWidthMm &&
       d.tanStereo>=0. && d.tanStereo<0.3,
       "unsupported candidate gap/width/pitch/stereo parameters");
}
// Two opposite cathodes have strips crossing the wires at pi/2 +/- alpha.
// u is coordinate across anode wires, v is electrode-normal, w is ALONG wires.
// At fixed u, the strip center of index n is at w = n*pitch +/- tan(alpha)*u.
// Hence the strip width/pitch are exactly the values seen ALONG w at fixed u.
// The first argument side distinguishes mirrored orientations, not ions.
inline double StripCoordinateMm(Side side,double uMm,double wMm,
                                const Design&d){
 const double slope=(side==Side::Plus?d.tanStereo:-d.tanStereo);
 return (wMm-slope*uMm)/std::sqrt(1.+d.tanStereo*d.tanStereo);
}
inline double StripCenterMm(int index,const Design&d){
 return index*d.stripPitchMm/std::sqrt(1.+d.tanStereo*d.tanStereo);
}
inline double HalfWidthMm(const Design&d){
 return 0.5*d.stripWidthMm/std::sqrt(1.+d.tanStereo*d.tanStereo);
}
// Harmonic measure of the selected finite-width strip at y=0 in an
// infinite parallel-plate domain 0<y<2g, with the opposite plane grounded.
// Exact for ideal 2-plate zero-thickness, unshielded straight strips.
// Complex exponential conformally maps strip domain to the upper half plane.
// Eq. Phi = (arg(z-B)-arg(z-A))/pi; A,B are mapped strip endpoints.
// IT DOES NOT INCLUDE ANODE WIRE SHIELDING at v=0.
inline double Potential(const Design&d,Side side,int index,
                        double uMm,double vMm,double wMm) {
 Validate(d);
 Check(std::isfinite(uMm)&&std::isfinite(vMm)&&std::isfinite(wMm) &&
       std::abs(vMm)<=d.gapMm,"invalid gas-side coordinate");
 const double x=StripCoordinateMm(side,uMm,wMm,d)-StripCenterMm(index,d);
 const double y=(side==Side::Plus?d.gapMm-vMm:d.gapMm+vMm);
 const double h=HalfWidthMm(d),depth=2.*d.gapMm;
 if(y<=0.)return std::abs(x)<h?1.:(std::abs(x)==h?0.5:0.);
 if(y>=depth)return 0.;
 const double k=PI/depth,xd=k*x;
 if(xd>45. || xd< -45.)return 0.; // asymptotically negligible
 const std::complex<double> z=std::exp(std::complex<double>(xd,k*y));
 const double B=std::exp(k*h),A=std::exp(-k*h);
 const double value=(std::arg(z-B)-std::arg(z-A))/PI;
 Check(std::isfinite(value)&&value>=-1.e-12&&value<=1.+1.e-12,
       "invalid harmonic potential");
 return std::max(0.,std::min(1.,value));
}
// Numerical field derivative for geometry diagnostics only.
// For physical induced currents use full wire-aware weighting fields.
inline Field WeightingField(const Design&d,Side side,int index,
                            double u,double v,double w,double hMm=0.001) {
 Validate(d);
 Check(hMm>0. && std::abs(v)+hMm<d.gapMm,
       "differentiation approaches a conductor");
 const double scale=-0.5/hMm;
 Field out;
 out.eu=scale*(Potential(d,side,index,u+hMm,v,w)-
               Potential(d,side,index,u-hMm,v,w));
 out.ev=scale*(Potential(d,side,index,u,v+hMm,w)-
               Potential(d,side,index,u,v-hMm,w));
 out.ew=scale*(Potential(d,side,index,u,v,w+hMm)-
               Potential(d,side,index,u,v,w-hMm));
 return out;
}
inline double DecodeUmm(double sPlusMm,double sMinusMm,const Design&d) {
 Check(d.tanStereo>0.,"cannot reconstruct u for zero stereo angle");
 return (sMinusMm-sPlusMm)/(2.*d.tanStereo);
}
inline double DecodeWmm(double sPlusMm,double sMinusMm) {
 return 0.5*(sPlusMm+sMinusMm);
}
} // namespace NA60Stereo

// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEPOINT__H
#define LIGHTCURVEPOINT__H

#include <iostream>
#include <iomanip>

#include "countedref.h"
#include "persistence.h"

//! lightcurve result
class LightCurvePoint: public RefCount {
public:
  CLASS_VERSION(LightCurvePoint,1);
#define LightCurvePoint__is__persistent


  LightCurvePoint();
  double julianday;
  double flux;
  double eflux;
  double mag;
  double emag_minus;
  double emag_plus;
  double zeropoint;
  void computemag(double zp) {
    zeropoint = zp;
    if(flux<1.e-12) {
      mag=99;
      emag_minus=0;
      emag_plus=0;
    }else{
      mag=-2.5*log10(flux)+zeropoint;
      if((1.-eflux/flux)<0)
	emag_minus=-99;
      else
	emag_minus=2.5*log10(1.-eflux/flux);
      emag_plus=2.5*log10(1.+eflux/flux);
    }
  }
#ifndef SWIG
  friend ostream& operator << (ostream &Stream, const LightCurvePoint &lcp) { 
    Stream << setiosflags(ios::fixed);
    Stream << setw(14) << setprecision(2) << lcp.julianday
	   << setw(15) << setprecision(3) << lcp.flux
	   << setw(15) << setprecision(3) << lcp.eflux
	   << setw(15) << setprecision(3) << lcp.mag
	   << setw(15) << setprecision(3) << lcp.emag_minus
	   << setw(15) << setprecision(3) << lcp.emag_plus
	   << setw(15) << setprecision(3) << lcp.zeropoint;
    return Stream;
  }
#endif
};


#endif // LIGHTCURVEPOINT__H

// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEPOINT__H
#define LIGHTCURVEPOINT__H

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
};


#endif // LIGHTCURVEPOINT__H

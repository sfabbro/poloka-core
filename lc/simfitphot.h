// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMFITPHOT__H
#define SIMFITPHOT__H

#include "simfit.h"
class SimFitPhot {
private:
  SimFit zeFit;

public:

  SimFitPhot(LightCurveList& Fiducials);

  void operator() (LightCurve& Lc);

};


#endif // SIMFITPHOT__H

#ifndef LSFIT__H
#define LCFIT__H



#include "reducedimage.h"
#include "vignette.h"

using namespace std;

class LcFit {

 private:
  ReducedImageRef geomRef;
  ObjectToFit objectToFit;
  VignetteList vignetteList;


 public :


  Point &PosInRef() const { return objectToFit;}

};



#endif /* LCFIT__H */



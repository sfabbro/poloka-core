// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEGURU__H
#define LIGHTCURVEGURU__H

#include <reducedimage.h>

//! 
//  \file lightcurveguru.h
//  \brief A master class to build light curves of many objects in many images.
//

class ImageSet : public ReducedImageList {


};


//! Main light curve builder that order all photometry to a 
class LightCurveGuru : public vector<ImageSet> {

public:

  LightCurveGuru() {}

  LightCurveGuru(const string& LcFileName);
 
  //! build all the light curves many different ways and take for ever
  void MonopolizeCPU();
};

/* \page lcresults Light curve results
   
Processing time depends on seeing and pixel resolution of all images. 
Output files are the following:
\arg subtracted DbImage's : <DbImage>-DbImageRef
\arg a lc-<method>.list catalog for your objects in each requested DbImage,
      where <method> is either:
      \arg "cat" for SExtractor default photometry
      \arg "sub" for subtraction aperture photometry
      \arg "sim" for simultaneous fit photometry
      \arg "apr" for N*FWHM direct aperture photometry
\arg one big ASCII file that combines all lightcurves lightcurve.list in  your current directory

*/


#endif // LIGHTCURVEGURU__H

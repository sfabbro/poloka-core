// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEBUILDER__H
#define LIGHTCURVEBUILDER__H

#include "lightcurve.h"

/*! 
  \file 
  \brief This class produces lightcurves of many objects, given images 
 */

typedef vector<LightCurve>::iterator LightCurveIterator;
typedef vector<LightCurve>::const_iterator LightCurveCIterator;

//! Lightcurve builder for one band. Execute different types of photometry
class LightCurveBuilder : public vector<LightCurve> {  
  void coincid(); // an internal routine to make sure starlists of all images are ordered the same way
  void RefStarListPhotometry(); // a routine to produce the photometry in the reference frame
  void RefAperPhotometry(const double &Nfwhm); // this sould later  helps for  absolute calibration

public:
  //! main constructor 
  LightCurveBuilder(const ReducedImageList &Images); 
  LightCurveBuilder() : Reference(NULL) {};
  //! copy constructor
  LightCurveBuilder(const LightCurveBuilder &Other);
  LightCurveBuilder& operator = (const LightCurveBuilder &Right);
  ~LightCurveBuilder();

  //! a pointer to the reference night
  Night *Reference;
  //! the list of images where each lightcurve point is refered
  NightList Nights;  
  //! a band name associated with this list
  string BandName;
  //! the number of supernovae
  unsigned nsn;
  
  //! should calibrate properly the lightcurve, but does not
  void Calibrate(const double &Nfwhm) const;
  //! init the lightcurve
  bool Init(const vector<Point> &SNe);
  //! init the lightcurve with a catalogue of standard stars
  bool Init(const vector<Point> &SNe, const string& catalogue);
  //! init all common objects
  void InitFiducials(const BaseStarList &FidList);
  //! performs a quick photometry from the starlists
  void StarListPhotometry();
  //! performs aperture photometry on the subtractions
  void SubAperPhotometry(const double &Nfwhm);
  //! build the kernels among the reference and the other images
  void BuildKernelsAndSubs();
  //! performs simultaneous photometry for each lightcurve star
  void SimFitPhotometry();
  //! write something
  void write(const string &MiddleName) const;
};

#endif // LIGHTCURVEBUILDER__H




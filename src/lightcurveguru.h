// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVEGURU__H
#define LIGHTCURVEGURU__H

#include "reducedimage.h"
#include "lightcurvebuilder.h"

class Point;

/*! \file lightcurveguru.h
  \brief A master class to run lightcurve builder in many bands
  
  The LightCurveGuru is ordering the processes to build the lightcurves
  from a set of images and points.
   
*/

struct LightCurveParams {

  double maxNeighborDist;
  unsigned maxFiducials;
  int matchOrder;
  // constructor with default values
  LightCurveParams() {
    maxNeighborDist = 5;// maximum fiducial distance to a neigbor
    matchOrder = 3;     // geometric matching order among nights
    maxFiducials = 50;  // max number of fiducials stars to produce a lightcurve simultaneously
  }
};


//! Main lightcurve builder for all bands
class LightCurveGuru {

  ReducedImage* coordRef; // the coordinate reference image
  ReducedImageList AllImages; // all images of all bands

  double jd_inf,jd_sup;  // the julian date bounds within which the SN is present
  bool check_coordinates() const;
  // an internal routine to match candidates to the coordinate reference
  vector<Point> match_candidates(const ReducedImage &GeomRef) const; 
  
  string StdCatalogue;

public:
  LightCurveGuru(): coordRef(NULL),jd_inf(0),jd_sup(1e9),fidRef(NULL) {StdCatalogue="";};
  ~LightCurveGuru();
 
  vector<Point> Supernovae;  // all candidates
  LightCurveBuilder *fidRef; // the reference lightcurve list where fiducials are expressed
  LightCurveParams lcParams; // might be useless soon
  vector<LightCurveBuilder> BandBuilders; // all lightcurve list of all bands

  //! read initial lightcurve file
  bool read(const string &FileName);
  //! check, align and stack images night by night, instrument by instrument, band by band
  bool MakeNights();
  //! select  a fiducial reference lightcurve list
  bool SelectFidRef();
  //! match fiducials objects among bands. A procedure that really needs some more thoughts.
  bool MatchFiducials();
  //! the main builder.
  bool MonopolizeCPU();
  //! Use a catalogue of standard stars for photometric calibration
  void UseStdCatalogue(const string& cata) {StdCatalogue = cata;};
  //! Select a subsample of fiducials close to supernovae
  bool SelectGoodFiducials();
};

/* \page lcresults Lightcurve results

Processing time depends on seeing and 
pixel resolution of all images. A gigantic amount of output files are
saved on disk. For each band, the following files are created:
 \arg all resampled DbImage's
 \arg co-added DbImage's
 \arg subtracted DbImage's
 \arg all vignets for your candidates (data,psf,kernel,residuals...) as FITS files
 \arg one ASCII file per band of the fiducials lightcurve results
 (simultaneous fit)    
 \arg one ASCII file per band of the fiducials lightcurve results
 (aperture on subtraction) 
 \arg one ASCII file per band for covariance matrix from simultaneous fit
 \arg all fiducial kernel and data vignets (removed when the program has
 finished running) 

Problems encountered can be numerous. Most common are:
\arg PSF DAOPHOT problem. Removed the huge core and rerun.
\arg Image registration problem. Check you have WCS in headers. Rerun.
\arg Other: email sfabbro@in2p3.fr
*/


#endif // LIGHTCURVEGURU__H

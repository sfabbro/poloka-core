// This may look like C code, but it is really -*- C++ -*-
#ifndef NIGHT__H
#define NIGHT__H

#include "reducedimage.h"

/*! \file night.h
  \brief a class to use as a ReducedImage, to avoid reopening the FitsImage 
  This is basically the same as a ReducedImage. However, if you need to use often most of the parameters of the ReducedImage, the Night class will be faster since it does not do any query on the FITS file.
  This one is used mostly to build lightcurves.
 */

//! version of the ReducedImage, quicker, but with more memory.
class Night : public ReducedImage {

public:
  //! Initialization of a night from its ReducedImage
  void init();
  Night(const string &Name);
  Night();
  ~Night(){};
  //! The saturation level
  double saturation;
  //! A photometric ratio with a given reference image
  double photomRatio;
  //! A Chi2/Dof deriving from a kernel fitting procedure
  double kernelChi2;
  //! the julian date associated with the night
  double julianDate;
  //! total exposure time for this night
  double exposure;
  //! the zero point of the night
  double zeroPoint;
  //! seeing (in pixel units) of the night
  double seeing;
  //! sigma (in pixel units) of the night PSF in x-direction
  double sigmaX;
  //! sigma (in pixel units) of the night PSF in y-direction
  double sigmaY;
  //! correlation angle of the night PSF
  double thetaXY;
  //! original sky background level of the night
  double skyLevel;
  //! current background level of the image associated with the night
  double backLevel;
  //! sky RMS of the night sky background
  double sigmaSky;
  //! gain of the image (in e-/ADU)
  double gain;
  //! an estimate of the error in flatfielding (in % of pixel intensity)
  double flatError;
  //! the readout noise in e-
  double readoutNoise;
  //! the systematic error of the chosen PSF profile
  double profileError;
  //! returns true if the night verifies a bunch of criteria
  bool IsOK() const;
  //! returns true if the night contribute as a zero point baseline for the lightcurve
  bool IsZeroRef;
  //! returns true if the night is photometric
  bool IsPhotometric;
  //! copy constructor
  ReducedImage *Clone() const;
  //! allows myNight1 < myNight2 to sort in time
  bool operator < (const Night &Right) const 
       {return (julianDate < Right.julianDate);}

#ifndef SWIG
  friend ostream& operator << (ostream & Stream, const Night& MyNight) 
       {MyNight.dump(Stream); return Stream;}
#endif

  void dump(ostream &Stream) const;
  void dump_header(ostream &Stream) const;
};

bool SameSeeing(const Night &Reference, const Night &Current, const double &Tol=0.05);

#include "imagelist.h"

typedef ImageList<Night>::iterator NightIterator;
typedef ImageList<Night>::const_iterator NightCIterator;

class NightList : public ImageList<Night> {
public:
  NightList() {};
  NightList(const string &FileName);
  NightList(const ReducedImageList &Images);
  //! the common pixel frame to all images of the list
  Frame CommonFrame();
  bool Init();
  void InitPhotomRatio(const ReducedImage &Reference);
  void FilterAllObjects(BaseStarList &RefList, const Frame &RefFrame, const double MinAssDist=1.5) const;
  void dump(ostream &Stream = cout) const;
  void write(ostream &Stream = cout) const;
  void write(const string &FileName) const;
  bool read(const string &FileName);
  bool read(istream &rd);    

#ifndef SWIG
  friend bool IncreasingNightSeeing(const Night *one, const Night *two);
  friend bool IncreasingNightJulian(const Night *one, const Night *two);
  friend bool DecreasingNightSaturation(const Night *one, const Night *two);
  friend bool DecreasingNightRatio(const Night *one, const Night *two);
  friend bool DecreasingNightArea(const Night *one, const Night *two);
#endif

};

#endif //NIGHT__H






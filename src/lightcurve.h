// This may look like C code, but it is really -*- C++ -*-
#ifndef LIGHTCURVE__H
#define LIGHTCURVE__H

#include "photstar.h"
#include "night.h"
#include "nightelement.h"

/*! \file lightcurve.h
  \brief A vector of FiducialStar of different nights
  The nights are usually geometrically aligned to use some of the
  routines of this class. 
*/

typedef NightElement<PhotStar> FiducialStar;
typedef ImageList<FiducialStar>::const_iterator FiducialStarCIterator;
typedef ImageList<FiducialStar>::iterator FiducialStarIterator;

//! The same FiducialStar monitored many nights.
class LightCurve : public ImageList<FiducialStar> {
public:
  LightCurve(const NightList &AllNights, const PhotStar &OneStar);
  LightCurve(istream &Stream);
  LightCurve(){};
  //! take the weighted mean of all positions.
  void AveragePosition();
  //! compute the relative fluxes according to the night photometric ratios.
  void RelativeFluxes();
  //! compute the absolute calibrated fluxes
  void AbsoluteFluxes();
  PhotStar refstar;
  //! returns true if the FiducialStar is a standard star
  bool IsStandard;
  //! name of the star, useful for bookeeping
  string FidName;
  void dump(ostream& Stream=cout) const;
  void write_header(ostream &Stream) const;
  void write(ostream &Stream) const;
  void write(const string &FileName) const;
  void read(istream &Stream);
  //! sorts by increasing fluxes
  friend bool ByIncreasingFlux(const FiducialStar *one, const FiducialStar *two);
  //! sorts by increasing signal to noise
  friend bool ByIncreasingSignalToNoise(const FiducialStar *one, const FiducialStar *two);
};

#endif // LIGHTCURVE__H

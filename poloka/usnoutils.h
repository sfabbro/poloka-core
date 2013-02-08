// This may look like C code, but it is really -*- C++ -*-
#ifndef USNOUTILS__H
#define USNOUTILS__H

//#include "starmatch.h"
//#include "fitsimage.h"

#include <poloka/sestar.h>
#include <poloka/gtransfo.h>
#include <poloka/vutils.h>

struct MatchCards {
  double linMatchCut; // max distance for lin match (arcsec)
  int linMatchMinCount; // min number of matches
  int distortionDegree; // degree of distortions.
  double secondMatchCut; // cut applied for collecting matches after distortion fit
  bool writeWCS;
  bool asciiWCS;
  string wcsFileName;
  string astromCatalogName;
  bool dumpMatches;
  bool ignoreSatur;
  bool ignoreBad;
  double minSigToNoise;
  
  MatchCards();
  bool ReadCards(const string &CardsFileName);

  bool cards_read;
};


#ifdef USNOUTILS__CC
MatchCards MatchPrefs;
#else
extern MatchCards MatchPrefs;
#endif


//class BaseStar;
class DbImage;
class FitsHeader;

/*! \file usnoutils.h usnoutils.cc
    \brief access the usno catalog, match image catalogs to it
    */


enum UsnoColor { RColor , BColor};



//!This routine reads the USNO catalog stars within the given bounds.
/*! The x and y coordinates of the given stars refer to RA and DEC respectively
expressed in  degrees. The location of the catalog should be set by the user
in the USNODIR environment variable. The catalog contains both
R and B magnitude. The Color argument has to be RColor or BColor. 
WARNING : The flux of the returned BaseStar's is in fact a magnitude. */

int UsnoRead(const Frame &F, UsnoColor Color, BaseStarList &ApmList);

//! same as above 
int UsnoRead(double minra, double maxra, 
             double mindec, double maxdec, 
	     UsnoColor Color, BaseStarList &ApmList);


//! replace the flux by pow(10, -flux/2.5) 
void ConvertMagToFlux(BaseStarList *List, const double Zp=0.);




StarMatchList *FindTransfoUsno(const string &FitsImageName, SEStarList &sestarlist, GtransfoLin &UsnoToPix, BaseStarList &UsnoCat);

/* uses a Guessed transfo to provide a good one. That is to say, it refines */
StarMatchList* GuessedToGood(const Gtransfo *Guess, const SEStarList *OldImag, const BaseStarList *OldCat);




int GetUsnoZeroPoint(const StarMatchList *List, UsnoColor Color,double& zeropoint, double& errzero);


bool UsnoProcess(const string &fitsFileName, const string &catalogName, 
		 DbImage *dbimage);

#endif /* USNOUTILS__H */

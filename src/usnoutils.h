// This may look like C code, but it is really -*- C++ -*-
#ifndef USNOUTILS__H
#define USNOUTILS__H

//#include "starmatch.h"
//#include "fitsimage.h"

#include "sestar.h"
#include "gtransfo.h"
#include "vutils.h"

//class BaseStar;
class DbImage;
class FitsHeader;

/*! \file usnoutils.h usnoutils.cc
    \brief access the usno catalog, match image catalogs to it
    */


typedef enum UsnoColor { RColor , BColor};



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
void ConvertMagToFlux(BaseStarList *List);




StarMatchList *FindTransfoUsno(const string &FitsImageName, SEStarList &sestarlist, GtransfoLin &UsnoToPix, BaseStarList &UsnoCat);

/* uses a Guessed transfo to provide a good one. That is to say, it refines */
StarMatchList* GuessedToGood(const Gtransfo *Guess, const SEStarList *OldImag, const BaseStarList *OldCat);




int GetUsnoZeroPoint(const StarMatchList *List, UsnoColor Color,double& zeropoint, double& errzero);


void FillMatchFile(const DbImage &Image, const Gtransfo &Pix2RaDec);

void FillMatchFile(const FitsHeader &header, const SEStarList &imageList, const Gtransfo &Pix2RaDec, 
                     const string &MatchFileName );

bool UsnoProcess(const string &fitsFileName, const string &catalogName, DbImage *dbimage, const bool write);
bool UsnoProcess(DbImage &dbimage, const bool write);

#endif /* USNOUTILS__H */

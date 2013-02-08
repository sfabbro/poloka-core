// This may look like C code, but it is really -*- C++ -*-
#include <iostream>
#include <string>

#include <poloka/datacards.h>
#include <poloka/sestar.h>

#ifndef SEEINGBOX_SEEN
#define SEEINGBOX_SEEN





// Structure pour le calcul du seeing:
// - selection des etoiles
// - calcul de la moyenne de la fwhm

struct DatSeeing {
  int prlevel  ;
  double saturation  ; 
  float prctage  ;
  int Nmax  ;
  DatSeeing();
  DatSeeing(const double Saturlev);
  DatSeeing(const string &FileName, const double Saturlev);
  void Print(ostream & s=cout) const ;
  void Default(const double saturlevel)  ;
  void LitDataCard(DataCards  & data );

};


// structure de sortie du calul du seeing
struct SortieSeeing{
  // nbre etoiles gardees pour les histos
  int Ngardees ; 
  // nbre etoiles gardees pour le calcul
  int Ncalcul ;
  float seeing ;
  float seeing_histo ;


  SortieSeeing() ;
  void Print(ostream & s=cout) const ;
};


// procedures de calcul de seeing



/*DOCF Seeing computation on SExtractor parameters 
  
  It computes a 2D histogram : fwhm of the object on one axis, shape
  of the object on the other.  
  The histogram is filled with cutted sextractor list, all the saturated
  and sextractor flagged are removed.
  The shape is defined to be -2.5*log10(flux/fluxmax) (fluxmax is the
  higher pixel of the object).
  This histogram exhibits a maximum which correspond to the best stars
  and then the fwhm of this object give a good and robust estimation
  of the seeing.

*/

int  
CalculeSeeingSE(DatSeeing const & dat, SortieSeeing & sortie, 
		SEStarList const & stlse, 
		SEStarList & seestar);


#endif

// This may look like C code, but it is really -*- C++ -*-
#ifndef ALIGNBOX_SEEN
#define ALIGNBOX_SEEN

#include <cmath>
#include <iostream>
#include <string>
#include <cstdio>

#include "image.h"
#include "datacards.h"
#include "frame.h"
#include "gtransfo.h"
#include "starmatch.h"

#ifdef STORAGE

#define NLIM_ALIGN 5


enum {ALIGN_SIMPLE = 0, ALIGN_USE_TRUNC };

struct DatAlign{
  //int prlevel ;
  //double marge ;
  //  int flag ;
  double saturation   ;
  //double  b_min ;
  //int degre ;
  //int type_align_phot ;
  double ph_fluxmax_prct ;
  int ph_fluxmax_N  ;
  //double n_sigma ;
  //double lim_n1n2 ;


  DatAlign() ;
  void Saturation(const double saturlevel) ;
  void Print(ostream & s=cout) const ;
  void  Default(const double saturlevel) ;
  void LitDataCard(DataCards  & data );
};


/*DOCF selection for geometric alignement */

int
GoodForGeomAlign(DatAlign const & dat, SEStarList const & stli, 
		   SEStarList & stlo);





/*DOCF selection for photometric alignement */
int
GoodForPhotomAlign(StarMatchList & liste, double saturation, 
		   double pcrtage, int N_max);

int
GoodForPhotomAlign(DatAlign const & dat, 
		   StarMatchList & liste);
/*DOCF second selection for photometric alignement after first fit*/

#ifdef IS_IT_USEFUL
int
CutForPhotomAlign(Poly & pl, double n_sigma, 
		   StarMatchList & liste);
#endif

/*DOCF routine for photom alignement (still in progress).
If an error occur, the returned polynomial is zero. Its degree
obtained by the method {\tt Degre()} (see poly.h) can be
tested to be different from zero. */


void SearchPhotomAlign(double & A, double & B, const SEStarList & stl1,
		       const SEStarList & stl2, const Gtransfo * tf,
		       const double saturation, const double Nmax);


#ifdef IS_IT_USEFUL
Poly
SearchPhotomAlign(DatAlign & datal, SortieAlign & sortie,  
		  SEStarList & stl1,
		  SEStarList & stl2, Gtransfo * tf);

int MeanSigAl(double *mag,double *err,int *ndata,
	      double nsigma,double *mean,double *sigma);
/*DOCF fit of flux2 = alpha * flux1 + 0.  on  StarMatchList */
Poly 
PhotomAlign2_(StarMatchList & liste);
/*DOCF fit of flux2 = Pl(flux1) on  StarMatchList */
Poly 
PhotomAlign_(StarMatchList & liste, int degre, double & chi2);


/*DOCF The polynome Pl is applied on every pixel */
void 
Poly_Apply(Poly const & pl, Image & img);

/*DOCF The polynome Pl is applied on all computed fluxes */
void 
Poly_Apply(Poly const & pl, SEStarList & stl);


void 
Permutation_XY(SEStarList & stl);

StarMatchList *
SimpleAlignment( BaseStarList *lr1, BaseStarList *lr2, double distance);
#endif /* IS_IT_USEFUL */


void 
SEStarListCutEdges(SEStarList *L, const Frame &F, const double MinDist);

#endif



#endif

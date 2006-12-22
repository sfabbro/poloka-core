#ifndef SEXTRACTOR_BOX__H
#define SEXTRACTOR_BOX__H

#include "sestar.h"

#include "image.h"


/* DOC The two routines that follows enable to invoke sextractor
and get the resulting catalog in a SEStarList. */

/* DOCF Structures to get every parameters right before calling for SExtractor */
class ForSExtractor{
public:
string FitsFileName ;
string FitsMaskName ; 
string FitsWeightName  ;
string FitsBackName  ;
 string FitsSegmentationName;
string FitsMiniBackName ;
 string TempDir;
 string UniqueName;

double saturation  ;
bool back_type_manual ; // if true, a constant value for the background is taken.
double sigma_back ;
double backmean ;
public:
 ForSExtractor(){back_type_manual = false ;sigma_back=-1;backmean=0.;};
 void DecompressIfNeeded();
 void Print();

 ~ForSExtractor();
 private :
  string ToRemove;
};

class AllForSExtractor: public ForSExtractor{

public:
string SexConfigFileName  ;
string SexParamName  ; 
string SexNNWName ;
string SexFilterName ;
public:
AllForSExtractor(ForSExtractor const & data0) 
  : ForSExtractor(data0){}
 void FillFromEnvironnement();
 void Print();
};




/*! make a SEStarList on an image. 
Background sigma and value are recovered
in the  Fond and SigmaFond variables. 
Whether FitsMaskName is given or not, or is empyt or not,
a mask is used to flag the objects.
default values are taken for specific SExtractor datacards (in TOADSCARDS) */
int
SEStarListMake(const ForSExtractor & shortdata,
	       SEStarList &List, double & Fond, 
	       double &SigmaFond,
	       Image * pmask_sat);

// To get the background map from the mini background map

Image *BackFromMiniBack(Image const & minib, int Nx, int Ny, 
			int back_meshx, int back_meshy);





#endif /* SEXTRACTOR_BOX__H */

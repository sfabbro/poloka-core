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

//--------------------- In Case of 2 images

string FitsFileName_0 ;
string FitsWeightName_0  ;
string FitsFileName_1 ;
string FitsWeightName_1  ;
string UniqueName_0 ;
string UniqueName_1 ;
string TempDir_0;
string TempDir_1;

//---------------------





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
 bool FillFromEnvironnement();
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
int
SEStarListMake_2(const ForSExtractor & shortdata,
	       SEStarList &List, double & Fond_0, 
	       double &SigmaFond_0, double & Fond_1, 
	       double &SigmaFond_1,bool weight_from_measurement_image);

// To get the background map from the mini background map

Image *BackFromMiniBack(Image const & minib, int Nx, int Ny, 
			int back_meshx, int back_meshy);

bool SubtractMiniBack(Image& Im, const string& minibackfile);



#endif /* SEXTRACTOR_BOX__H */

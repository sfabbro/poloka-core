#ifndef CALIBRATEDSTAR__H
#define CALIBRATEDSTAR__H


#include "basestar.h"

class DicStar;

struct CalibratedStar : public BaseStar {
 public:
 //  CalibratedStar() {};
 // CalibratedStar(BaseStar& toto) : BaseStar(toto) {};
 CalibratedStar(const DicStar &D);
  double ra,dec;
  double u,g,r,i,z,x,y;
  double ue,ge,re,ie,ze;
  double neighborDist,neighborFlux;
  double neighborFluxContamination;
  double neighborNsigma;
  int id;
};


class Frame;
class Gtransfo;

#include "reducedimage.h"

class CalibratedStarList : public StarList<CalibratedStar>
{
 public :
   // check ra dec catalog is in frame and return x y list
  CalibratedStarList(const string &CatalogName, const Gtransfo *WCS, 
		     const Frame& ImageFrame);

 // check ra dec catalog is in refimage frame, compares with refimage apersecatalog to get contamination, and return x y list
  CalibratedStarList(const string &CatalogName,const  ReducedImageRef & refimage);

  CalibratedStarList() {};

};

typedef CalibratedStarList::iterator CalibratedStarIterator;
typedef CalibratedStarList::const_iterator CalibratedStarCIterator;



#endif /* CALIBRATEDSTAR__H */

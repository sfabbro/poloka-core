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
  int id;
};


class Frame;
class Gtransfo;

class CalibratedStarList : public StarList<CalibratedStar>
{
 public :
  CalibratedStarList(const string &CatalogName, const Gtransfo *WCS, 
		     const Frame& ImageFrame);

  CalibratedStarList() {};

};

typedef CalibratedStarList::iterator CalibratedStarIterator;
typedef CalibratedStarList::const_iterator CalibratedStarCIterator;



#endif /* CALIBRATEDSTAR__H */

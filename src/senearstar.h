#ifndef SENEARSTAR__H
#define SENEARSTAR__H

#include <string>

#include "sestar.h"
#include "starlist.h"

class Image;
class FitsImage;
class Gtransfo;

class SENearStar : public SEStar 
{

 private:
  
  int xShift, yShift;
  double photFactor;
  double galFlux;
  
 public:
  
  SENearStar(const SEStar &AStar, const int Dx = 0, const int Dy = 0, const double PhotFactor=0, const double GalFlux=0);
  SENearStar();
  void AddToFitsImage(Image &image, Image & dest, 
		      const Gtransfo *Transfo) const;
  void AddToFitsImage(Image &image, const Gtransfo *Transfo) const{
      AddToFitsImage(image,image,Transfo);};
  SENearStar *ActualFake() const;
  void writen(ostream & pr = cout) const;
  std::string WriteHeader_(ostream &pr = cout, const char* i = NULL) const;
  static SENearStar* read(istream& r, const char* Format);
  

};

class SENearStarList : public StarList<SENearStar> 
{
  public :
  void AddToFitsImage(FitsImage &image, const Gtransfo *Transfo) const;
  void ActualFakes(SENearStarList &Result) const;
 
};

#ifndef SWIG
BaseStarList* SENear2Base(SENearStarList * This);
const BaseStarList* SENear2Base(const SENearStarList * This);
#endif

typedef SENearStarList::const_iterator SENearStarCIterator;
typedef SENearStarList::iterator SENearStarIterator;


#endif /* SENEARSTAR__H */

#ifndef FIDSTAR__H
#define FIDSTAR__H

#include "persistence.h"


#include "photstar.h"

//using namespace std;
class SEStar;

class FidStar : public PhotStar {
  CLASS_VERSION(FidStar,1);
  #define FidStar__is__persistent
public:

  FidStar();
  FidStar(const BaseStar &BStar);
  FidStar(const SEStar &SStar);
  
  int nd;
  double jd;
  static FidStar* read(istream & Stream, const char *Format);
  friend ostream& operator << (ostream &Stream, const FidStar &PStar)
    {PStar.dump(Stream); return Stream;}
  virtual void dumpn(ostream  &Stream = cout) const;
  virtual void writen(ostream &Stream = cout) const;
  string WriteHeader_(ostream & Stream, const char* Num=NULL) const;
  static const char* TypeName() {return "FidStar";}

};

#include "starlist.h"

typedef StarList<FidStar> FidStarList;
typedef StarList<class FidStar>::const_iterator FidStarCIterator;
typedef StarList<class FidStar>::iterator FidStarIterator;

#endif // FIDSTAR__H

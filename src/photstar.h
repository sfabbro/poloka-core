// This may look like C code, but it is really -*- C++ -*-
#ifndef PHOTSTAR__H
#define PHOTSTAR__H

#include "persistence.h"

#include "basestar.h"
class SEStar;

class PhotStar : public BaseStar {
  CLASS_VERSION(PhotStar,1);
  #define PhotStar__is__persistent

public:
  PhotStar();
  PhotStar(const BaseStar &BStar);
  PhotStar(const SEStar &SStar);
  double sky,varx,vary,covxy,varflux,varsky;
  static PhotStar* read(istream & Stream, const char *Format);
  friend ostream& operator << (ostream &Stream, const PhotStar &PStar)
  {PStar.dump(Stream); return Stream;}
  virtual void dumpn(ostream  &Stream = cout) const;
  virtual void writen(ostream &Stream = cout) const;
  string WriteHeader_(ostream & Stream, const char* Num=NULL) const;
  static const char* TypeName() {return "PhotStar";}
};

#include "starlist.h"

typedef StarList<PhotStar> PhotStarList;
typedef list<class PhotStar*>::const_iterator PhotStarCIterator;
typedef list<class PhotStar*>::iterator PhotStarIterator;

#endif // PHOTSTAR__H


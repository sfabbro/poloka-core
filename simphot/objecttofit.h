#ifndef OBJECTTOFIT__H
#define OBJECTTOFIT__H


#include <string>

#include "basestar.h"

using namespace std;

class ObjectToFit : public BaseStar
{
 private :
  double jdmin,jdmax; // begin and end of "on"
  string name; // for IO purposes
  string band;
  int type; // fit type;


 public:
  ObjectToFit(const string& Line);

  int FitType() const { return type;}
  const string Name() const { return name;}
  const string &Band() const { return band;}
  const double JdMin() const { return jdmin;}
  const double JdMax() const {return jdmax;}

};


#include <starlist.h>

typedef StarList<ObjectToFit> ObjectToFitList;
typedef ObjectToFitList::const_iterator ObjectToFitCIterator;
typedef ObjectToFitList::iterator ObjectToFitIterator;


#endif /* OBJECTTOFIT__H */

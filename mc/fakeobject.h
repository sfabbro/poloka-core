#ifndef FAKEOBJECT__H
#define FAKEOBJECT__H

#include <iostream>
//#include "basestar.h"
//#include "sestar.h"
#include "point.h"
#include "countedref.h"
#include <cstdio>

using namespace std;


class FakeObject : public Point , public RefCount 

#ifdef USE_ROOT
	, public TObject
#endif
{
protected:
double mag;

 public:
 
 double Ra() const {return x;}
 double& Ra()  {return x;}
 double Dec() const {return y;}
 double& Dec()  {return y;}
 double Mag() const {return mag;}
 double& Mag()  {return mag;}

FakeObject(double Ra, double Dec, double fluxflux);
FakeObject();

  /* DOC \noindent {\bf For read & write}: */
  
  //! for dump with NO end-of-line
  virtual void    dumpn(ostream& s = cout) const;

  //! for dump
  virtual void    dump(ostream& s = cout) const ;
  
  virtual void write(ostream& s = cout) const ;
  
  //! for write with NO end-of-line
  virtual void    writen(ostream& s = cout) const ;

  //! to read once the object is created 
  void   read_it(istream& r, const char *Format); 

   //! to read and create the object  

  static FakeObject* read(istream& r, const char *Format); 

  /* DOCF  to write the FakeList header with the string ${}^{\star}i$
     appended to every ntuple variable (with no end)  */

  std::string WriteHeader_(ostream & pr = cout, const char* i = NULL) const;
  
  virtual void WriteHeader (ostream & stream =cout) const; 

  static const char *TypeName() { return "FakeObject";}

};


#include "fakelist.h"
#include <list>

#ifdef USE_ROOT

typedef FakeListWithRoot<FakeObject> FakeObjectList;

#else

typedef FakeList<FakeObject> FakeObjectList;

#endif /* USE_ROOT */



typedef FakeObjectList::const_iterator FakeObjectCIterator;
typedef FakeObjectList::iterator FakeObjectIterator;
typedef CountedRef<FakeObject> FakeObjectRef;

#endif

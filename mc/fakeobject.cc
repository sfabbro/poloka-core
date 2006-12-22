#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "fakeobject.h"
#include "basestar.h" 
#include "fakelist.h"


#ifdef STORAGE
//! object aperture data is incomplete or corrupted
#define APER_CORRUPT(m)    ((bool) ((m) & APER_COR_FLAG) )
// object isophotal data are incomplete or corrupted
#define ISO_CORRUPT(m)    ((bool) ((m) & ISO_COR_FLAG ) )
// memory overflow during deblending
#define MEM_OVFLOW1(m)    ((bool) ((m) & MEM_OV1_FLAG) )
// memory overflow during extraction
#define MEM_OVFLOW2(m)    ((bool) ((m) & MEM_OV2_FLAG) )
#endif

FakeObject::FakeObject()
  : Point(0.0,0.0)
{mag=0.;}

FakeObject::FakeObject(double RaRa, double DecDec, double fluxflux)
  : Point(RaRa,DecDec) {mag=fluxflux;}




void
FakeObject::dumpn(ostream& s) const
{
  s << " Ra : " << x ;
  s << " Dec : " << y ;
s << " Mag : " << mag ;

}
 


void
FakeObject::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}

void
FakeObject::write(ostream& s)  const
{
	writen(s);
	s << endl;
}

void
FakeObject::writen(ostream& s)  const
{

  s << Ra() << " ";
  s << Dec() << " ";
  s << Mag() << " ";

 

}


void
FakeObject::read_it(istream& s, const char * Format)
{
  int format = DecodeFormat(Format, "FakeObject");
  if (format >=1)
    { 
    s >> Ra();
    s >> Dec();
    s >> Mag();

    }
 
  return ;
}

FakeObject*  FakeObject::read(istream& r, const char *Format)
{
  FakeObject *pstar = new FakeObject();  
  pstar->read_it(r, Format);
  return(pstar);
}


void FakeObject::WriteHeader(ostream & stream) const
{
	string format =WriteHeader_(stream);
	stream <<"# format " << format << endl;
	stream << "# end" << endl;
	
}

std::string FakeObject::WriteHeader_(ostream & pr, const char *i) const
{
  if (i== NULL) i= "";
  
  pr    << "# Ra"<< i <<" : " << endl 
  	<< "# Dec"<< i <<" : " << endl 
	<< "# Mag"<<i <<" : " << endl;

    
/* 1 is the current format id for FakeObjects (when being written) it must correspond
to the right behaviour of the read routine ( and match what write does ! ) */
return " FakeObject 1 ";
}






#ifdef USE_ROOT
ClassImp(FakeObject)

/* To Generate the sestardict.cc file :
LINKDEF_CONTENT : #pragma link C++ class FakeObject+;
*/
#endif /* USE_ROOT */


//********************   FINDEFINITION FakeObject   *********************



#include "fakelist.cc" /* since starlist is a template class */


#ifdef USE_ROOT
template class FakeListWithRoot<FakeObject>;
ClassImpT(FakeListWithRoot,FakeObject);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class CountedRef<FakeObject>-;
LINKDEF_CONTENT : #pragma link C++ class list<CountedRef<FakeObject> >;
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeObject> >::unique();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeObject> >::sort();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<FakeObject> >::merge(list <CountedRef<FakeObject> >)&;
LINKDEF_CONTENT : #pragma link C++ class StarList<FakeObject>-;
LINKDEF_CONTENT : ostream& operator << (ostream&, const StarList<FakeObject>&);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<FakeObject>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<FakeObject>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<FakeObject>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef FakeObjectIterator;
*/
#include "root_dict/sestardict.cc"
#endif /* USE_ROOT */

template class FakeList<FakeObject>; // because StarListWithRoot<> derives from StarList<>



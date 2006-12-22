#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
 
#include "host.h"
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

Host::Host()
  : Point(0.0,0.0)
{
  z=0.0;
  err_z=0.0;
  gal_type=0.0;
  a=0.0;
  b=0.0;
  theta=0.0;
  i_mag=0.0;
  }

Host::Host(double RaRa, double DecDec)
  : Point(RaRa,DecDec) 
  {
  z=0.0;
  err_z=0.0;
  gal_type=0.0;
  a=0.0;
  b=0.0;
  theta=0.0;
  i_mag=0.0;
  }
  





void
Host::dumpn(ostream& s) const
{
  s << " Ra : " << x ;
  s << " Dec : " << y ;
  s << " Z : " << z ;
  s << " err_Z : " << err_z ;
  s << " gal_type : " << gal_type ;
  s << " A : " << a ;
  s << " B : " << b ;
  s << " Theta : " << theta ;
  s << " I_Mag : " << i_mag ;
}
 


void
Host::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}

void
Host::write(ostream& s)  const
{
	writen(s);
	s << endl;
}

void
Host::writen(ostream& s)  const
{

  s << Ra() << " ";
  s << Dec() << " ";
  s << Z() << " ";
  s << Err_Z() << " ";
  s << Gal_Type() << " ";
  s << A() << " ";
  s << B() << " ";
  s << Theta() << " ";
  s << I_Mag() << " ";

 

}


void
Host::read_it(istream& s, const char * Format)
{
  int format = DecodeFormat(Format, "Host");
  if (format >=1)
    { 
    s >> Ra();
    s >> Dec();
    s >> Z();
    s >> Err_Z();
    s >> Gal_Type();    
    s >> A();
    s >> B();
    s >> Theta();
    s >> I_Mag(); 
    }
 
  return ;
}

Host*  Host::read(istream& r, const char *Format)
{
  Host *pstar = new Host();  
  pstar->read_it(r, Format);
  return(pstar);
}


void Host::WriteHeader(ostream & stream) const
{
	string format =WriteHeader_(stream);
	stream <<"# format " << format << endl;
	stream << "# end" << endl;
	
}

std::string Host::WriteHeader_(ostream & pr, const char *i) const
{
  if (i== NULL) i= "";
  
  pr    << "# Ra"<< i <<" : " << endl 
  	<< "# Dec"<< i <<" : " << endl 
	<< "# Z"<<i <<" : " << endl
  	<< "# Err_Z"<< i <<" : " << endl 
	<< "# Gal_Type"<<i <<" : " << endl
	<< "# A"<< i <<" : " << endl 
	<< "# B"<<i <<" : " << endl
  	<< "# Theta"<< i <<" : " << endl 
	<< "# I_Mag"<<i <<" : " << endl;
    
/* 1 is the current format id for Hosts (when being written) it must correspond
to the right behaviour of the read routine ( and match what write does ! ) */
return " Host 1 ";
}






#ifdef USE_ROOT
ClassImp(Host)

/* To Generate the sestardict.cc file :
LINKDEF_CONTENT : #pragma link C++ class Host+;
*/
#endif /* USE_ROOT */


//********************   FINDEFINITION Host   *********************



#include "fakelist.cc" /* since starlist is a template class */


#ifdef USE_ROOT
template class FakeListWithRoot<Host>;
ClassImpT(FakeListWithRoot,Host);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class CountedRef<Host>-;
LINKDEF_CONTENT : #pragma link C++ class list<CountedRef<Host> >;
LINKDEF_CONTENT : #pragma link off function list<CountedRef<Host> >::unique();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<Host> >::sort();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<Host> >::merge(list <CountedRef<Host> >)&;
LINKDEF_CONTENT : #pragma link C++ class StarList<Host>-;
LINKDEF_CONTENT : ostream& operator << (ostream&, const StarList<Host>&);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<Host>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<Host>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<Host>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef HostIterator;
*/
#include "root_dict/sestardict.cc"
#endif /* USE_ROOT */

template class FakeList<Host>; // because StarListWithRoot<> derives from StarList<>



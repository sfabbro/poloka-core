#include <cstdio>
#include <cstdlib>

#include <iomanip>
#include "basestar.h"


using namespace std;

/* comment on coordinates value writen to disk : to follow IRAF-FITS convention
we decided to write coodinates values to disk which are expressed in pixels 
starting at (1,1), with cebter of the leftmots bottom pixel = (1.,1.). 
Internally, Toads will keep on with the natural C behaviour, where the leftmost
pixle is (0,0) and its center (0.,0.).
So when we read, we subtract 1, when we write, we add 1. In order to
be able to locate all occurences of this shift, we #define MEMPIX2DISK 1*/


void BaseStar::read_it(istream & rd, const char *format)
{
 int formatValue = 0;
 if (format) 
   formatValue = DecodeFormat(format,"BaseStar");
 rd >> x >> y >> flux;
 if (formatValue >= 1) // only shift back if shifted when written
   {
     x -= MEMPIX2DISK;
     y -= MEMPIX2DISK;
   }
}  


BaseStar* BaseStar::read(istream & rd, const char *format)
{
  BaseStar *p = new BaseStar();
  p->read_it(rd, format);
  return p;
}


string BaseStar::WriteHeader_(ostream & stream, const char*i) const
{
  if (i==NULL) i = "";
  stream << "# x"<< i <<" : x position (pixels)" << endl 
	 << "# y"<< i <<" : y position (pixels)" << endl 
	 << "# flux"<< i <<" : flux en unites du pixel" << endl ;
  return " BaseStar 1 "; 
}

void BaseStar::WriteHeader(ostream & stream) const
{ 
  string format = WriteHeader_(stream); 
  stream << "# format " << format << endl;
  stream << "# end " << endl ;
};


void BaseStar::writen(ostream &s) const 
{
  s << x + MEMPIX2DISK << " " << y + MEMPIX2DISK << " " << flux << " " ;
}


void BaseStar::write(ostream &s) const 
{
  writen(s);
 s << endl ;
}

bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2)
{
return (S1->flux > S2->flux);
}

bool IncreasingMag(const BaseStar *S1, const BaseStar *S2)
{
return (S1->flux < S2->flux);
}

#include "rootstuff.h"

ClassImp(BaseStar)

/* To Generate the %dict.cc file :
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class Point;
LINKDEF_CONTENT : #pragma link C++ function operator<<(ostream &, const Point&);
LINKDEF_CONTENT : #pragma link C++ class BaseStar-;
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream &, const BaseStar &);

use a custom Streamer (-) in the previous statement to avoid
making "Point" known to Root
*/

#ifdef USE_ROOT

void BaseStar::Streamer(TBuffer &R__b)
{
   // Stream an object of class BaseStar. hand made
 
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> x;
      R__b >> y;
      R__b >> flux;
      R__b.CheckByteCount(R__s, R__c, BaseStar::IsA());
   } else {
      R__c = R__b.WriteVersion(BaseStar::IsA(), kTRUE);
      R__b << x;
      R__b << y;
      R__b << flux;
      R__b.SetByteCount(R__c, kTRUE);
   }
}
#endif /*USE_ROOT */


/**************** BaseStarList ******************/


#include "starlist.h"

int DecodeFormat(const char *FormatLine, const char *StarName)
{
if (!FormatLine || !StarName) return 0;
const char *p= strstr(FormatLine, StarName);
 if (!p) return  0;
return atoi( p + strlen(StarName));
}

#include "starlist.cc" /* since starlist is a template class */

#ifdef USE_ROOT
template class StarListWithRoot<BaseStar>; /* to force instanciation */
ClassImpT(StarListWithRoot,BaseStar)

  /* all ancestors of StarListWithRoot<T>  ave to be bound explicitely to cint,
     because I did not find what triggers that by default (because semetimes, it 
     just happens by itself) */

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT

LINKDEF_CONTENT : #pragma link C++ class CountedRef<BaseStar>;
LINKDEF_CONTENT : #pragma link C++ class list<CountedRef<BaseStar> >;
LINKDEF_CONTENT : #pragma link off function list<CountedRef<BaseStar> >::unique();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<BaseStar> >::sort();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<BaseStar> >::merge(list <CountedRef<BaseStar> >)&;
LINKDEF_CONTENT : #pragma link C++ class StarList<BaseStar>;
LINKDEF_CONTENT : #pragma link C++ class StarList<BaseStar>::iterator;
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream &, const StarList<BaseStar> &);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<BaseStar>-;
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream &, const BaseStar &);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream &, const BaseStarList &);
*/
#include "root_dict/basestardict.cc"
#endif /* ifdef USE_ROOT */

// whatever happens (root or not), StarList<BaseStar> has to be instanciated
template class StarList<BaseStar>; /* to force instanciation */


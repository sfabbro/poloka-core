#include <assert.h>

//#include <iomanip>
#include "basestar.h"
#include "starlistexception.h"


using namespace std;

/* comment on coordinates value writen to disk : to follow IRAF-FITS convention
we decided to write coodinates values to disk which are expressed in pixels 
starting at (1,1), with cebter of the leftmots bottom pixel = (1.,1.). 
Internally, Toads will keep on with the natural C behaviour, where the leftmost
pixle is (0,0) and its center (0.,0.).
So when we read, we subtract 1, when we write, we add 1. In order to
be able to locate all occurences of this shift, we #define MEMPIX2DISK 1*/


/* comment to the comment :
this proves to be a bad idea : BaseStar's may contain coordinates
which do NOT represent pixels. It could be degrees for example.
In this case subtracting 1 is a very bad idea. So it is no longer the case.
Old files (BaseStar format = 1) are correctly interpreted, but we no
longer use this trick of having different coordinate origins on disk and
in memory.
*/

static double sq(const double &x) { return x*x;}


#include "fastifstream.h"

void BaseStar::read_it(fastifstream & rd, const char *format)
{
 int formatValue = 0;
 if (format) 
   formatValue = DecodeFormat(format,"BaseStar");
 if (formatValue == 3)
   {
     rd >> x >> y >> vx >> vy >> vxy >> flux;
     /* write sig(x) sig(y), rho ... */
     vxy *= (vx*vy);
     vx *= vx;
     vy *= vy;
   }
 else if (formatValue == 2 || formatValue == 0)
   rd >> x >> y >> flux;
 else if (formatValue == 1) // only shift back if shifted when written
   {
     x -= MEMPIX2DISK;
     y -= MEMPIX2DISK;
   }
 else throw(StarListException(" Unknown format value for BaseStar "));
}  

/* this routine is used by DicStar. It *HAS* to be consistent with BaseStar::read_it, above */
unsigned NValsBaseStar(const char *Format)
{
  int formatValue = DecodeFormat(Format, "BaseStar");
  if (formatValue <= 2) return 3;
  if (formatValue == 3) return 6;
  return 0;
}



BaseStar* BaseStar::read(fastifstream & rd, const char *format)
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
	 << "# sx"<< i <<" : x position r.m.s " << endl
	 << "# sy"<< i <<" : y position r.m.s " << endl
	 << "# rhoxy"<< i <<" : xy correlation " << endl
	 << "# flux"<< i <<" : flux in image ADUs" << endl ;
  return " BaseStar 3 "; 
}

void BaseStar::WriteHeader(ostream & stream) const
{ 
  string format = WriteHeader_(stream); 
  stream << "# format " << format << endl;
  stream << "# end " << endl ;
};


void BaseStar::writen(ostream &s) const 
{
  assert(vx>0 && vy>0 && sq(vxy)<vx*vy);
  /* write (sigx,sigy,rho) rather than (vx,vy,vxy). 
     This limits shortcomings of truncation, and is more useful 
     when reading or plotting.
  */
  double sx = sqrt(vx);
  double sy = sqrt(vy);
  s << x << ' ' << y << ' ' 
    << sx << ' ' << sy << ' ' << vxy/sx/sy << ' ' 
    << flux << ' ' ;
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



// whatever happens (root or not), StarList<BaseStar> has to be instanciated
template class StarList<BaseStar>; /* to force instanciation */


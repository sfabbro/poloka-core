#include <iomanip>
#include "photstar.h"
#include "sestar.h"

PhotStar::PhotStar()
  : sky(0), varx(0), vary(0), covxy(0), varflux(0),varsky(0)
{}

PhotStar::PhotStar(const BaseStar &BStar)
  : BaseStar(BStar),sky(0),varx(0),vary(0),covxy(0),varflux(0),varsky(0)
{}

PhotStar::PhotStar(const SEStar &SStar)
  : BaseStar(SStar),sky(SStar.Fond()),varx(0),vary(0),covxy(0),
    varflux(SStar.EFlux()*SStar.EFlux()),varsky(0)
{}

void PhotStar::writen(ostream &Stream) const
{
  int old = Stream.precision();
  BaseStar::writen(Stream); 
  Stream << ' ' << sky << ' '
	 << varx << ' '
	 << vary << ' '
	 << covxy << ' '
	 << varflux << ' '
	 << varsky << ' ';
  Stream << setprecision(old);
}

PhotStar* PhotStar::read(istream & Stream, const char *Format)
{
  PhotStar *pstar = new PhotStar();
  pstar->BaseStar::read(Stream, Format);
  if (Stream >> pstar->sky >> pstar->varx >> pstar->vary >> pstar->covxy >> pstar->varflux >> pstar->varsky) 
    return pstar;
  return NULL;    
}

void PhotStar::dumpn(ostream &Stream) const
{
  BaseStar::dumpn(Stream);
  Stream << " sky " << sky << " varx " << varx << " vary " << vary 
	 << " covxy " << covxy << " varflux " << varflux;    
}

string PhotStar::WriteHeader_(ostream & Stream, const char* Num) const
{
  if (Num==NULL) Num = "";
  const string baseStarFormat = BaseStar::WriteHeader_(Stream,Num);
  Stream << "# sky"<< Num <<" : mean sky value per pixel" << endl 
	 << "# varx"<< Num <<" : variance in x position (pixels)" << endl 
	 << "# vary"<< Num <<" : variance in y position (pixels)" << endl 
	 << "# covxy"<< Num <<" : covariance in x and y position (pixels)" << endl 
	 << "# varflux"<< Num <<" : variance of measured flux" << endl
	 << "# varsky" << Num <<" : variance of measured sky" << endl;  
  return baseStarFormat + "PhotStar 0 ";  
}

#include "starlist.cc" /* since starlist is a template class */
template class StarList<PhotStar>; /* to force instanciation */


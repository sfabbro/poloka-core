#include "fidstar.h"
#include <iomanip>


FidStar::FidStar()
  : nd(0), jd(0)
{}

FidStar::FidStar(const BaseStar &BStar)
  : PhotStar(BStar),nd(0),jd(0)
{}

FidStar::FidStar(const SEStar &SStar)
  : PhotStar(SStar),nd(0),jd(0)
{}


void FidStar::writen(ostream &Stream) const
{
  Stream << resetiosflags(ios::scientific) ;
  Stream << setiosflags(ios::fixed) ;
  int old = Stream.precision();
  Stream << setprecision(12) ;
  Stream << ' ' << nd << ' ' << jd << ' ';
  PhotStar::writen(Stream);
  Stream << ' 0';
  Stream << setprecision(old) ;
}

FidStar* FidStar::read(istream & Stream, const char *Format)
{
  FidStar *pstar = new FidStar();
  Stream >> pstar->nd >> pstar->jd;
  int toto;
  if (Stream >> pstar->x>> pstar->y >> pstar->flux >> pstar->sky >> pstar->varx >> pstar->vary >> pstar->covxy >> pstar->varflux >> pstar->varsky >> toto) 
    return pstar;
  return NULL;
}


void FidStar::dumpn(ostream &Stream) const
{
  Stream << resetiosflags(ios::scientific) ;
  Stream << setiosflags(ios::fixed) ;
  int old = Stream.precision();
  Stream << setprecision(12) ;
  Stream << " nd " << nd << " jd " << jd;
  Stream << setprecision(old) ;
    
  PhotStar::dump(Stream);
}

string FidStar::WriteHeader_(ostream & Stream, const char* Num) const
{
  if (Num==NULL) Num = "";
  //const string baseStarFormat = BaseStar::WriteHeader_(Stream,Num);
  Stream << "# n"<< Num <<" : fiducial tag number" << endl 
	 << "# jd"<<Num <<" : reduced julian date " << endl 
	 << "# x"<<Num <<" : x position (pixels)" << endl 
    	 << "# y"<<Num <<" : y position (pixels)" << endl 
    	 << "# flux"<<Num <<" : flux en unites du pixel" << endl 
	 << "# sky"<< Num <<" : mean sky value per pixel" << endl 
	 << "# varx"<< Num <<" : variance in x position (pixels)" << endl 
	 << "# vary"<< Num <<" : variance in y position (pixels)" << endl 
	 << "# covxy"<< Num <<" : covariance in x and y position (pixels)" << endl 
	 << "# varflux"<< Num <<" : variance of measured flux" << endl
	 << "# varsky" << Num <<" : variance of measured sky" << endl
	 << "# format" << Num <<" : format " << endl; 
  return "PhotStar 0 ";  
}



#include "starlist.cc"
template class StarList<FidStar>; // To force instanciation




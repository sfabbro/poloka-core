#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "datdetec.h"
#include "sestar.h"
#include "image.h"
#include "yquemstar.h"

//********************   DEFINITION  YquemStar   *********************
// Converter :
BaseStarList* Yquem2Base(YquemStarList * This)
{ return (BaseStarList*) This;} 

const BaseStarList* Yquem2Base(const YquemStarList * This)
{ return (BaseStarList*) This;} 


YquemStar::YquemStar(CandidateStar const & star) 
  : CandidateStar(star)
{
  dass_ref=0;
}

void
YquemStar::dumpn(ostream& s) const
{
  CandidateStar::dumpn(s);
  s << " objet associe sur ref : "  ;
  StarRef.dumpn(s);
  s << "dist. assoc. " << dass_ref ;
}


void
YquemStar::dump(ostream& s) const
{
  dumpn(s);
  s << endl ;
}


void
YquemStar::writen(ostream& s)  const
{
  CandidateStar::writen(s); 
  s  << resetiosflags(ios::scientific) ;
  s  << setiosflags(ios::fixed) ;
  s  << setprecision(4) ;
  StarRef.writen(s);
  s << " " << dass_ref << " " ;

}

void
YquemStar::writen_scan(ostream& s)  const
{
  s  << resetiosflags(ios::scientific) ;
  s  << setiosflags(ios::fixed) ;
  s  << setprecision(4) ;
  s << x << " "; 
  s << y << " "; 
  s << Fluxmax() << " ";
  s << flux << " ";
  s << EFlux() << " ";
  s << Noise()   << " ";
  s << SigToNoise()   << " "; 
  s << Faper_ref()   << " "; 
  s << PrctIncrease()   << " "; 
  s << StarRef.flux << " "; 
  s << StarRef.x << " " ;
  s << StarRef.y << " "; 
  s << " " << dass_ref << " " ;
}





// le format specifique de la YquemStar est celui de la SEStar StarRef pour l'instant



void
YquemStar::read_it(istream& r, const char * Format)
{
  CandidateStar::read_it(r,Format);
  StarRef.read_it(r,Format);
  r >> dass_ref ;
  return ;
}

YquemStar*  YquemStar::read(istream& r, const char * Format)
{
  YquemStar *pstar = new YquemStar();  
  pstar->read_it(r, Format);
  return(pstar);
}


string YquemStar::WriteHeader_(ostream & pr, const char*i) const
{
  if (i==NULL) i = "";
  string cStarFormat =  CandidateStar::WriteHeader_(pr,i);
  string s = i ;
  string u = "ref" + s ;
  string formatref = SEStar::WriteHeader_(pr, u.c_str());
  pr << "# dass_r"<< i << " : dist. with nearest object on ref. " << endl;
  return (cStarFormat + formatref + " YquemStar 1 ");
}

const char *YquemStar::WriteHeader_scan(ostream & pr, const char*i) const
{
  if (i==NULL) i = "";
  pr << "# x"<< i <<" : x position (pixels)" << endl 
     << "# y"<< i <<" : y position (pixels)" << endl 
     << "# amp"<< i <<" : Peak pixel value above background on sub" << endl 
     << "# flux"<< i <<" : flux in pixel units on sub " << endl 
     << "# eflux"<< i <<" : error on the flux" << endl;
  pr << "# noise"<< i <<" : noise (for sub) " << endl;  
  pr << "# stnoise"<< i <<" : S/N for obj. detecetd on sub " << endl;
  pr << "# fap_r"<< i <<" : aperture flux computed on ref at candidate position " << endl; 
   pr << "# prcti"<< i <<" : percentage increase ie flux/fap_r " << endl; 
  pr << "# hostflux"<< i <<" : Flux of the nearest object on the ref " << endl; 
  pr << "# xref"<< i <<" : X position of the nearest object on the ref " << endl; 
  pr << "# yref"<< i <<" : Y position of the nearest object on the ref " << endl; 
  pr << "# dass_r"<< i << " : dist. with nearest object on ref. " << endl;
  static char format[256];
  sprintf(format,"YquemStar %d",1 );
  return format; 
}



int YquemStar::write_nice(ostream & pr) const
{ 
  pr << resetiosflags(ios::scientific) ;
  pr << setiosflags(ios::fixed) ;
  pr  << setprecision(3) ;

  pr << endl << "Candidat " << Numero() << endl ;
  pr << "--------------------------" << endl ;
  pr << "x,y,flux, noise, S/N : " << x << " " << y << " " << flux
     << " " << Noise() << " " << SigToNoise() << endl ;


  double u1 =0., u3=0 ;
  if ( fabs(StarRef.Fluxmax()) > 1.e-10)
    u1 =  Fluxmax()/StarRef.Fluxmax() ;

  if ( fabs(StarRef.flux ) > 1.e-10)
    u3 = flux/StarRef.flux ;
  

  pr << endl <<  "dass_r, fluxmax, fluxmaxr, fluxmax/fluxmaxr : " 
     << dass_ref << " | " 
     << StarRef.Fluxmax() << " "  << StarRef.Fluxmax() << " " 
     << u1 << endl ;
  pr <<  "flux, flux_ref ouverture , flux/flux_ref ouverture: " 
     << flux << " "  << Faper_ref()  << " " 
     << PrctIncrease()  << endl ;
  pr <<  "flux, flux_ref  , flux/flux_ref: " 
     << flux << " "  <<  StarRef.flux << " " 
     << u3  << endl ;


 

  pr << endl << "candidat : " << endl;
  dumpn(pr);
  pr << endl ;

return 1;
}






#include "starlist.cc" /* since starlist is a template class */

template class StarList<CandidateStar>;
template class StarList<YquemStar>;



void YquemStar::write_scan(ostream &pr) const
{
  pr << x << " " 
     << y << " " 
     << Fluxmax() << " "
     << flux << " "
     << EFlux() << " " 
     << Noise()   << ' '
     << SigToNoise()   << ' ' 
     << Faper_ref()   << ' ' 
     << PrctIncrease()   << ' ' 
     << StarRef.flux << " " 
     << StarRef.x << " " 
     << StarRef.y << " " 
     << dass_ref << endl;
  
}


void YquemStar::WriteHeader_Scan(ostream &pr)
{
  pr << "# x : x position (pixels)" << endl 
     << "# y : y position (pixels)" << endl 
     << "# amp : Peak pixel value above background on sub" << endl 
     << "# flux : flux in pixel units on sub " << endl 
     << "# eflux : error on the flux" << endl;
  pr << "# noise : noise (for sub) " << endl;  
  pr << "# stnoise : S/N for obj. detecetd on sub " << endl;
  pr << "# fap_r : aperture flux computed on ref at candidate position " << endl; 
  pr << "# prcti : percentage increase ie flux/fap_r " << endl; 
  pr << "# hostflux : Flux of the nearest object on the ref " << endl; 
  pr << "# xref : X position of the nearest object on the ref " << endl; 
  pr << "# yref : Y position of the nearest object on the ref " << endl; 
  pr << "# dass_r : dist. with nearest object on ref. " << endl;
  pr << "# end" << endl;
}
//********************   FIN DEFINITION YquemStar   *********************


void YquemStarList::write_scan(const string &FileName) 
{
  ofstream pr(FileName.c_str());
  
  YquemStar::WriteHeader_Scan(pr);
  for(YquemStarCIterator it= this->begin(); it!= this->end(); it++)
    {
      (*it)->write_scan(pr);
    }
  pr.close();
}  


void
YquemStarList::ComputeFAperRef(Image & imgref, double rayonref, double fd_ref)
{
  for (YquemStarIterator it= this->begin(); it!= this->end(); ++it )
    (*it)->ComputeFaperRef(imgref,rayonref,fd_ref  );
  
}


#include "fastfinder.h"

void
YquemStarList::Construct_List(CandidateStarList & stl, SEStarList & stlref)
{
  // ttes les etoiles sont dans le systeme Candidate ie SE
  int n = 0 ;

  // Add the fastfinder... Tested. It works...
  FastFinder finder(*((BaseStarList*) &stlref));
  
  for (CandidateStarIterator it= stl.begin(); it!= stl.end(); it++ )
    { 
      BaseStar *pstar =  *it ;
      
      // la plus proche sur ref
      const BaseStar *pr = finder.FindClosest(*pstar,100);
      
      YquemStar *ystar = new YquemStar(**it);
      if (pr != NULL)
	{
	  ystar->StarRef  = **it ;
	  ystar->dass_ref = pr->Distance(*pstar);
	}
      else
	ystar->dass_ref = -1;
      ystar->Numero() = n;
      this->push_back(ystar); n++;
    }
    
}

// utilise lors de la detection sur 1 sub
void
YquemStarList::Construct(CandidateStarList & stl, SEStarList & stlref, 
			 Image & imgref, double rayon_ref, double fd_ref)
{
  this->Construct_List(stl, stlref);
  this->ComputeFAperRef(imgref, rayon_ref, fd_ref);
}



// APPLY CUTS ON THE LISTS
void
YquemStarList::Cut(YquemStarList & stlcut, DatDetec & dat)
{  
  int n = 0 ;
  for (YquemStarIterator it= this->begin(); it!= this->end(); ++it )
    {
      
      YquemStar *star =  *it ;
      if ( star->IsGood(dat.sigtonoise, dat.frac_faper_ref) )
	{
	  star->Numero() = n;
	  stlcut.push_back(new YquemStar(*star)); n++;
	}
    }
}

YquemStarList::YquemStarList(const string FileName)
{
  this->read(FileName);
}

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "sestar.h"

//********************   DEFINITION  SEStar   *********************


// Converter :
BaseStarList* SE2Base(SEStarList * This)
{ return (BaseStarList*) This;} 

const BaseStarList* SE2Base(const SEStarList * This)
{ return (BaseStarList*) This;} 


// CONSTRUCTORS

#define BAD_FLAG 1
#define COSMIC_FLAG 2

#define NEIGBHOUR_FLAG 1
#define BLEND_FLAG 2
#define SATUR_FLAG 4
#define TRUNC_FLAG 8
#define APER_COR_FLAG 16
#define ISO_COR_FLAG 32
#define MEM_OV1_FLAG 64
#define MEM_OV2_FLAG 128
#define DAOFLAG 256


#define NOTSATUR_FLAG ~SATUR_FLAG


//! The object has neighbours, bright and close enough to bias aperture photometry
bool SEStar::HasNeighbours() const
{
return (Flag() & NEIGBHOUR_FLAG); 
}

//! The object is blended with another one 
bool SEStar::IsBlended() const
{
return (Flag() & BLEND_FLAG);
}

//! at least one pixel is saturated

bool SEStar::IsSaturated(const double saturation) const
{
return (Fluxmax() + Fond() >= saturation );
}

bool SEStar::IsSaturated() const
{
return (Flag() & SATUR_FLAG);
}

bool SEStar::IsCosmic() const
{
return (FlagBad() & COSMIC_FLAG);
}

bool SEStar::IsBad() const
{
return (FlagBad() > 0);
}


//! the object is truncated
bool SEStar::IsTruncated() const
{
return (Flag() & TRUNC_FLAG);
}

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


//! To un-flag a star because it is not saturated
void SEStar::FlagAsNotSaturated()
{
  Flag() = Flag() & (NOTSATUR_FLAG);
}

//! To flag a star as saturated
void SEStar::FlagAsSaturated()
{
  Flag() = Flag() | SATUR_FLAG;
}

//! To flag a star as cosmic
void SEStar::FlagAsCosmic()
{
  FlagBad() = FlagBad() | COSMIC_FLAG;
}



SEStar::SEStar()
: BaseStar(0.,0.,0.)
{
  Set_to_Zero();
}

SEStar::SEStar(double xx, double yy, double ff)
: BaseStar(xx,yy,ff)
{  
 Set_to_Zero();
}

void
SEStar::Set_to_Zero()
{
  
  // x_image = x  de BaseStar
  // y_image = y  de BaseStar
  // flux_best = flux de BaseStar  

  xpeak = 0 ;  
  ypeak = 0 ;
  fluxmax = 0 ; 
  e_flux  = 0 ;  
  fond = 0. ;

  flux_aper  = 0 ;
  e_flux_aper  = 0 ;
  flux_fixaper  = 0 ;
  e_flux_fixaper  = 0 ;
  flux_iso = 0 ; 
  e_flux_iso = 0 ;
  flux_isocor = 0 ; 
  e_flux_isocor = 0 ;

  fwhm = 0 ;
  kronradius = 0 ;
  isoarea = 0 ;
 
  mxx = 0 ;
  myy = 0 ;
  mxy = 0 ;

  a = 0 ; 
  b = 0 ; 
  gyr_angle= 0;
  flag= 0;
  flagbad= 0;
  cstar= 0;
  xtrunc = 0;
  ytrunc = 0;
  num = 0 ;
  iter = 0;
  chi = 0;
  sharp = 0;

}





void
SEStar::dumpn(ostream& s) const
{
  s << " x : " << x ;
  s << " y : " << y ;
  s << " flux : " << flux ;
  s << " Fluxmax : " << Fluxmax() ;
  s << " EFlux : " << EFlux() ;
  s << " Fond : " << Fond() ;
  s << " Flux_aper : " << Flux_aper() ;
  s << " Eflux_aper : " << Eflux_aper() ;
  s << " Flux_iso : " << Flux_iso() ;
  s << " Eflux_iso : " << Eflux_iso() ;
  s << " Flux_isocor : " << Flux_isocor() ;
  s << " Eflux_isocor : " << Eflux_isocor() ;
  s << " Kronradius : " << Kronradius() ;
  s << " Isoarea : " << Isoarea() ;
  s << " Fwhm : " << Fwhm() ;
  s << " Mxx : " << Mxx() ;
  s << " Myy : " << Myy() ;
  s << " Mxy : " << Mxy() ;
  s << " A : " << A() ;
  s << " B : " << B() ;
  s << " Gyr_Angle : " << Gyr_Angle() ;
  s << " Flag : " << Flag() ;
  s << " FlagBad : " << FlagBad() ;
  s << " Cstar : " << Cstar() ;
  s << " xpeak : " << xpeak ;
  s << " ypeak : " << ypeak ;
  s << " Flux_fixaper : " << Flux_fixaper() ;
  s << " Eflux_fixaper : " << Eflux_fixaper() ;
  s << " Number : " << N(); 
  s << " Iterations: " << Iter();
  s << " Chi: " << Chi();
  s << " Sharp: " << Sharp();
}


void
SEStar::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}


void
SEStar::writen(ostream& s)  const
{

  BaseStar::writen(s);
  s << Fluxmax()  << " " ;
  s << EFlux()   << " " ; 
  s << Fond()  << " " ;
  s << Flux_aper ()  << " " ;
  s << Eflux_aper()   << " " ;
  s << Flux_iso()  << " " ; 
  s << Eflux_iso()  << " " ;
  s << Flux_isocor()  << " " ; 
  s << Eflux_isocor()  << " " ;
  s << Kronradius()  << " " ;
  s << Isoarea () << " " ;
  s << Fwhm()  << " " ; 
  s << Mxx()  << " " ;
  s << Myy()  << " " ;
  s << Mxy() << " " ;
  s << A()  << " " ; 
  s << B()  << " " ; 
  s << Gyr_Angle() << " " ;
  s << Flag() << " " ; 
  s << FlagBad() << " " ; 
  s << Cstar()  << " " ;
  s << Xtrunc()  << " " ;
  s << Ytrunc()  << " " ; 
  s << X_Peak() << " " ;
  s << Y_Peak() << " " ; 
  s << Flux_fixaper ()  << " " ;
  s << Eflux_fixaper()   << " " ; 
  s << N()   << " " ; 
  s << Iter() << " ";
  s << Chi() << " ";
  s << Sharp() << " ";
}


void
SEStar::read_it(istream& r, const char * Format)
{
  BaseStar::read_it(r, Format);
  r >> Fluxmax()  ;
  r >> EFlux()   ; 
  r >> Fond()  ;
  r >> Flux_aper ()  ;
  r >> Eflux_aper()   ;
  r >> Flux_iso()  ; 
  r >> Eflux_iso()  ;
  r >> Flux_isocor()  ; 
  r >> Eflux_isocor()  ;
  r >> Kronradius()  ;
  r >> Isoarea () ;
  r >> Fwhm()  ; 
  r >> Mxx()  ;
  r >> Myy()  ;
  r >> Mxy() ;
  r >> A()  ; 
  r >> B()  ; 
  r >> Gyr_Angle() ;
  r >> Flag() ; 
  r >> FlagBad() ; 
  r >> Cstar()  ;
  int format = DecodeFormat(Format, "SEStar");
  if (format >=1)
    {
    r >> Xtrunc()  ; 
    r >> Ytrunc()  ; 
    r >> X_Peak() ;
    r >> Y_Peak() ;
    }
  if (format >=2)
    {
      r >> Flux_fixaper ()  ;
      r >> Eflux_fixaper()   ;
    }
  if (format < 3)
    {
      x -= 1.0;
      y -= 1.0;
    }
  if (format >=4)
    {
      r >> N()  ;
    }
  if (format >=5)
    {
      r  >> Iter();
      r  >> Chi();
      r  >> Sharp();
    }
  return ;
}

SEStar*  SEStar::read(istream& r, const char *Format)
{
  SEStar *pstar = new SEStar();  
  pstar->read_it(r, Format);
  return(pstar);
}


std::string SEStar::WriteHeader_(ostream & pr, const char *i) const
{
  if (i== NULL) i= "";
  string baseStarFormat =  BaseStar::WriteHeader_(pr, i);
  pr    << "# fluxmax"<< i <<" : Peak pixel value above background " << endl 
	<< "# with SExtractor flux = flux_best = fisoc if the object is not  crowded, = faper otherwise. Very close to total flux (5%)" << endl 
	<< "# eflux"<< i <<" : RMS error for  flux " << endl 
	<< "# fond"<< i <<" :  background " << endl 
        << "# faper"<< i <<" : Flux within a Kron-like elliptical aperture. The radius is 2.5 * Kron-Radius (see krad param.)" << endl 
	<< "# efaper"<< i <<" : RMS error for  faper " << endl 
	<< "# fiso"<< i <<" : isophotal flux " << endl
	<< "# efiso"<< i <<" : RMS error for  isophotal flux " << endl 
	<< "# fisoc"<< i <<" : isophotal corrected flux = total flux of object " << endl 
	<< "# efisoc"<< i <<" : RMS error for  isophotal corrected flux" << endl 
        << "# krad"<< i <<" : extension of the  aperture , in units of A or B" 
	<< endl 
	<< "# isoar"<< i <<" : area of lowest isophote en pixels " << endl 
	<< "# fwhm"<< i <<" : " << endl 
        << "# mxx"<< i <<" : x second order moment = <x^2> - <x>^2" << endl 
	<< "# myy"<< i <<" : y second order moment = <y^2> - <y>^2" << endl 
	<< "# mxy"<< i <<" : <xy> - <y><x>" << endl 
	<< "# a"<< i <<" : Profile RMS along major axis " << endl 
	<< "# b"<< i <<" : Profile RMS along minor axis" << endl 
	<< "# angle"<< i <<" :gyration angle of the major axis 0.0 = axis en degres "
	<< endl
	<< "# flag"<< i <<" : Extraction flags (flag=0 = object is unsaturated and uncrowded)" << endl 
        << "# flagbad"<< i <<" : number of bad pixels in iso area" << endl 
	<< "# cstar"<< i <<" : star/galaxy classification in [0,1], 0=gal and 1=star"  << endl  
	<< "# xtrunc"<< i <<" : x [peak/2,peak]"  << endl 
	<< "# ytrunc"<< i <<" : y [peak/2,peak]"  << endl
	<< "# xpeak"<< i <<" : x-coordinate of the brightest pixel" << endl 
	<< "# ypeak"<< i <<" : y-coordinate of the brightest pixel" << endl 
        << "# ffixa"<< i <<" : Flux within a fixed aperture" << endl 
	<< "# effixa"<< i <<" : RMS error for  ffixaper " << endl 
	<< "# num"<< i <<" : star number " << endl 
	<< "# iter" << i << "  : number of iterations" << endl 
	<< "# chi" << i << "  : robustified chi per star" << endl 
	<< "# sharp" << i << "  : sharp index" << endl; 


/* 1 is the current format id for SEStars (when being written) it must correspond
to the right behaviour of the read routine ( and match what write does ! ) */
return baseStarFormat + " SEStar 5 ";
}




// test si etoile pas sat
bool SEStar::IsOK(const double &saturation) const 
{
  return ((FlagBad() == 0 ) &&
	  ( Flag() == 0 ) &&
	  ( fabs(Fluxmax()) > 1e-10)  &&
	  ( fabs(flux) > 1e-10)  &&
	  ( Fluxmax() + Fond() < saturation ));
}




bool DecFluxMax(const SEStar *S1, const SEStar *S2)
{
return (S1->Fluxmax() > S2->Fluxmax());
}


#ifdef USE_ROOT
ClassImp(SEStar)

/* To Generate the sestardict.cc file :
LINKDEF_CONTENT : #pragma link C++ class SEStar+;
*/
#endif /* USE_ROOT */


//********************   FINDEFINITION SEStar   *********************


#include "rootstuff.h"
#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */


#ifdef USE_ROOT
template class StarListWithRoot<SEStar>;
ClassImpT(StarListWithRoot,SEStar);

/* comments to drive the Makefile part that runs rootcint
RUN_ROOTCINT
LINKDEF_CONTENT : #pragma link C++ class CountedRef<SEStar>-;
LINKDEF_CONTENT : #pragma link C++ class list<CountedRef<SEStar> >;
LINKDEF_CONTENT : #pragma link off function list<CountedRef<SEStar> >::unique();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<SEStar> >::sort();
LINKDEF_CONTENT : #pragma link off function list<CountedRef<SEStar> >::merge(list <CountedRef<SEStar> >)&;
LINKDEF_CONTENT : #pragma link C++ class StarList<SEStar>-;
LINKDEF_CONTENT : ostream& operator << (ostream&, const StarList<SEStar>&);
LINKDEF_CONTENT : #pragma link C++ function operator << (ostream&, const StarList<SEStar>&);
LINKDEF_CONTENT : #pragma link C++ class StarListWithRoot<SEStar>-;
LINKDEF_CONTENT : #pragma link C++ class StarList<SEStar>::iterator;
LINKDEF_CONTENT : #pragma link C++ typedef SEStarIterator;
*/
#include "root_dict/sestardict.cc"
#endif /* USE_ROOT */

template class StarList<SEStar>; // because StarListWithRoot<> derives from StarList<>

//! Flag objects near pixels  > 0
//! OK for images with very few pixels > 0
int 
FlagCosmicFromImage(SEStarList & stl, 
		    Image const & image, const double dist)
{
  int nx = image.Nx();
  int ny = image.Ny();
  Pixel *a = image.begin();
  int nobj = 0;
  int n=0;
  for (int j=0; j<ny; ++j) for (int i=0; i<nx; ++i, ++a) 
      if (*a > 0) 
	{
	  SEStar *cosmic = stl.FindClosest(double(i),double(j));
	  if (cosmic->Distance(Point(i,j)) < dist) 
	    {
	      if (! (cosmic->IsCosmic()) ) 
		{
		  cosmic->FlagAsCosmic();
		  if (cosmic->IsBad()) n++;
		  nobj++ ;
		}
	    }
	}
  cout << nobj << "cosmics detected " << endl ;
  cout << n << " obj already flagged as bad by SExtractor" << endl ;
  return(nobj);
}

int 
FlagSaturatedStars(SEStarList & stl, const double saturation)
{
  int n = 0 ;
  for (SEStarIterator it= stl.begin(); it!=stl.end(); it++)
    {
      if ( !((*it)->IsSaturated()) && (*it)->IsSaturated(saturation) )
	{ 
	  (*it)->FlagAsSaturated() ; 
	  n++ ;
	}
    }
  return n;
}

int 
FlagSatFromImage(SEStarList & stl, Image const & image)
{
  int n=0;
  for (SEStarIterator it= stl.begin(); it!=stl.end(); it++)
    {
      if ( image(int((*it)->x),int((*it)->y)) > 0 )
	{
	  (*it)->FlagAsSaturated();n++;
	}
    }
  return n;
}

int KeepIt(double saturation ,SEStarList const & stli,SEStarList & temp)
{
  int Nok = 0;
  for (SEStarCIterator it= stli.begin(); it!=stli.end(); it++)
    {
      
      if ( (*it)->IsOK(saturation) && ((*it)->Fwhm() > 1.e-5) 
	   && ((*it)->flux > 1.e-5) )
	{
	  SEStar * starse = new SEStar(**it) ;
	  temp.push_back(starse); Nok++;
	}

    } 
  return Nok;
}

void 
SetStarsBackground(SEStarList & stl, double const background)
{
  SEStarIterator it ;
  int nflag=0;
  for(it = stl.begin() ; it != stl.end() ; it++)
    {
      SEStar *star = *it;
      star->Fond() =  background ;
      if (star->FlagBad()>0) nflag++;
    }
  cout << nflag << " etoiles flagees. " << endl ;
}

#include "histo2d.h"
int StarFinder(SEStarList & stlse, SEStarList & PrettyStars, 
	       double saturation)
{
  Histo2d histo(100,0,10,50,-10,0); //HC
  
  SEStarList cutlist;
  for (SEStarCIterator it=stlse.begin(); it!=stlse.end(); ++it) 
    { 
      if (((*it)->Mxx() > 1.e-10) &&
	  ((*it)->Myy() > 1.e-10) && ( (*it)->IsOK(saturation)) )
	{
	  histo.Fill((*it)->Fwhm() ,-2.5*log10((*it)->flux/(*it)->Fluxmax()), 1 ); 
	  cutlist.push_back(new SEStar(**it));
	}
    }

  // recherche du mod de l'histo
  double X, Y;
  histo.MaxBin(X, Y);
  double BinX, BinY;
  histo.BinWidth(BinX, BinY);
  int nsigma = 10;
  double moyenne = X + 0.5 * BinX;
  double sigma = 0.5 * BinX;
  for (SEStarIterator it= cutlist.begin(); it!=cutlist.end(); ++it)
    {
      SEStar *pstar = *it ;
      double flux = pstar->flux;
      double fluxmax = pstar->Fluxmax();
      double fwhm = pstar->Fwhm();
      double shape = -2.5*log10(flux/fluxmax);
      if ( (fabs(fwhm-moyenne) < nsigma * sigma) &&  
	   (Y-BinY < shape) && (shape < Y + BinY))
	{
	  PrettyStars.push_back(new SEStar(*pstar));
	}
    }
  
  return 1;
}

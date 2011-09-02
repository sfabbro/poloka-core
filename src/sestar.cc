#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>

#include "sestar.h"
#include "fastifstream.h"

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
#define SAT_ON_MASK 512 // le centroid tombe dans une zone flaggee comme saturee


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

  local_seeing = 0 ;
   aper_flux  = 0 ;
   err_aper_flux = 0  ;
   aper_flux_other = 0  ;

  xpeak = 0 ;  
  ypeak = 0 ;
  fluxmax = 0 ; 
  e_flux  = 0 ;  
  fond = 0. ;

  flux_auto  = 0 ;
  e_flux_auto  = 0 ;
  flux_petro = 0 ;
  e_flux_petro  = 0 ;
  flux_circ_aper  = 0 ;
  e_flux_circ_aper  = 0 ;
  flux_iso = 0 ; 
  e_flux_iso = 0 ;
  flux_isocor = 0 ; 
  e_flux_isocor = 0 ;

  fwhm = 0 ;
  kronradius = 0 ;
  petroradius = 0 ;
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
  s << " Flux_auto : " << Flux_auto() ;
  s << " Eflux_auto : " << Eflux_auto() ;
  s << " Flux_petro : " << Flux_petro() ;
  s << " Eflux_petro : " << Eflux_petro() ;
  s << " Flux_iso : " << Flux_iso() ;
  s << " Eflux_iso : " << Eflux_iso() ;
  s << " Flux_isocor : " << Flux_isocor() ;
  s << " Eflux_isocor : " << Eflux_isocor() ;
  s << " Kronradius : " << Kronradius() ;
  s << " Petroradius : " << Petroradius() ;
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
  s << " Flux_fixaper : " << Flux_circ_aper() ;
  s << " Eflux_fixaper : " << Eflux_circ_aper() ;
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
  s << Flux_auto ()  << " " ;
  s << Eflux_auto()   << " " ;
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
  s << Flux_circ_aper ()  << " " ;
  s << Eflux_circ_aper()   << " " ; 
  s << N()   << " " ; 
  s << Iter() << " ";
  s << Chi() << " ";
  s << Sharp() << " ";
  s << Flux_petro ()  << " " ;
  s << Eflux_petro()   << " " ;
  s << Petroradius()   << " " ;
}


void
SEStar::read_it(fastifstream& r, const char * Format)
{
  BaseStar::read_it(r, Format);
  r >> Fluxmax()  ;
  r >> EFlux()   ; 
  r >> Fond()  ;
  r >> Flux_auto ()  ;
  r >> Eflux_auto()   ;
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
      r >> Flux_circ_aper ()  ;
      r >> Eflux_circ_aper()   ;
    }
  /* I don't know why this is there, but in the absence of format, it just 
     messes up the whole thing :
  if (format < 3)
    {
      x -= 1.0;
      y -= 1.0;
    }
  */

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
  if (format >=6)
    { 
      r >> Flux_petro()  ;
      r >> Eflux_petro() ;
      r >> Petroradius() ;
    }
  return ;
}

SEStar*  SEStar::read(fastifstream& r, const char *Format)
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
	<< "# eflux"<< i <<" : RMS error for  flux " << endl 
	<< "# fond"<< i <<" :  background " << endl 
        << "# faper"<< i <<" : Flux within a Kron-like elliptical aperture. The radius is Kron-Radius (see krad param.)" << endl 
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
	<< "# sharp" << i << "  : sharp index" << endl
	<< "# fpet" << i << "  : f petro" << endl 
	<< "# efpet" << i << "  : err f petro" << endl  
	<< "# rpet" << i << "  : petro radius" << endl; 


/* 1 is the current format id for SEStars (when being written) it must correspond
to the right behaviour of the read routine ( and match what write does ! ) */
return baseStarFormat + " SEStar 6 ";
}



// test si etoile ok
bool SEStar::IsOK() const 
{
  return ((FlagBad() == 0 ) &&
	  ( Flag() == 0 ) &&
	  ( Fluxmax() > 1e-10)  &&
	  ( flux > 1e-10) && 
	  ( Fwhm() > 1.e-5 ) );
}

// test si etoile ok et pas sat
bool SEStar::IsOK(double saturation) const 
{
  return ( IsOK() &&
	  ( Fluxmax() + Fond() < saturation ));
}


static int NBad(float x, float y , Image const & image, float demi_cote) 
{
  int xstart = int(x - demi_cote + 0.5);
  int ystart = int(y - demi_cote + 0.5);
  int xend = min(int(xstart + 2.*demi_cote + 0.5),image.Nx());
  int yend = min(int(ystart + 2.*demi_cote+0.5),image.Ny());
  xstart = max(xstart,0);
  ystart = max(ystart,0);
  float seuil = 1.e-20 ;
  int nbad = 0 ;
  for (int j = ystart; j < yend; ++j)
    {
      for (int i = xstart; i < xend; ++i)
	{
	  if (image(i,j) <= seuil )
	    nbad++;
	}
    }
  return(nbad);
  
}
// to be able to eliminate a star for which the stamp is contaminated by a dead column
bool SEStar::IsOK(Image const & image) const 
{
  if ( ! IsOK() )
    return false ;
  int nbad = NBad(x,y,image,StampSize());
  if ( nbad > StampSize() )
    return false ;
  else
    return true ;
}

// to be able to eliminate a star with too big a neighbour

bool
SEStar::HasBigCloseNeighbor(BaseStarList const & stl, double dist, double ratio) const
{
  return (stl.HasCloseNeighbor(x,y, dist, 0.5, ratio*flux) );
}



bool DecFluxMax(const SEStar *S1, const SEStar *S2)
{
return (S1->Fluxmax() > S2->Fluxmax());
}

bool IncNumber(const SEStar *S1, const SEStar *S2)
{
  return (S1->N() < S2->N());
}





//********************   FINDEFINITION SEStar   *********************


#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */


template class StarList<SEStar>; // because StarListWithRoot<> derives from StarList<>

//! Flag objects near pixels  > 0
//! OK for images with very few pixels > 0
//UNUSED ?? (DH)
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
      if ( image(int((*it)->x+0.5),int((*it)->y+0.5)) > 0 )
	{
	  (*it)->FlagAsSaturated();n++;
	}
    }
  return n;
}

int KeepOK(double saturation ,SEStarList const & stli,SEStarList & temp)
{
  int Nok = 0;
  for (SEStarCIterator it= stli.begin(); it!=stli.end(); it++)
    {
      
      if ( (*it)->IsOK(saturation) )
	{
	  SEStar * starse = new SEStar(**it) ;
	  temp.push_back(starse); Nok++;
	}

    } 
  return Nok;
}
void RemoveNonOKObjects(SEStarList &List)
{
  cout << "Nombre d'objects avant nettoyage " << List.size()<<endl;
  // On degage les objets a probleme
  for (SEStarIterator it= List.begin(); it != List.end();) 
    {
      SEStar *star = *it;
      if ( !(star->IsOK()) ) 
	it = List.erase(it);
      else ++it;
    }
  cout << "Nombre d'objects apres nettoyage " << List.size()<<endl;
}

void 
SetStarsBackground(SEStarList & stl, double const background)
{
  SEStarIterator it ;
  for(it = stl.begin() ; it != stl.end() ; it++)
    {
      SEStar *star = *it;
      star->Fond() =  background ;
    }
}

/*************** Star Finder, Galaxy Finder  *********/
/* Finder needs a clean list */
#include "histo2d.h"
#include "histo1d.h"

static
void HistoFwhmStarFinder(SEStarList const & stlse,  
	       double & mfwhm, double  & bin_fwhm)
{
  
  double minfwhm = 1. ;
  double  maxfwhm = 10. ;
  int NBin = (int)((maxfwhm -minfwhm)/(1.*bin_fwhm) + 0.5)  ;
  cerr << "N bin " << NBin << endl ;
  Histo1d histo(NBin,minfwhm,maxfwhm); //HardCoded
  
  for (SEStarCIterator it=stlse.begin(); it!=stlse.end(); ++it) 
    { 
      histo.Fill((*it)->Fwhm() , 1. );
    }

  // recherche du mode de l'histo
  histo.MaxBin(mfwhm);
  cout << "Mode en fwhm: " << mfwhm << " (Bin: " 
       << histo.BinWidth() << " = " << bin_fwhm << " )"
       << endl ; 
  bin_fwhm = histo.BinWidth();


  double minfwhm2 = mfwhm - 5*bin_fwhm ;
  double  maxfwhm2 = mfwhm + 5*bin_fwhm ;
  double S =0., S2=0. ;
  int n = 0 ;
  for (SEStarCIterator it=stlse.begin(); it!=stlse.end(); ++it) 
    { 
      double fwhm = (*it)->Fwhm() ;
      if ( (fwhm > minfwhm2 ) && (fwhm < maxfwhm2 ))
	{
	  S += fwhm ;
	  S2 += fwhm*fwhm ;
	  n++;
	}
    }
  double mean = S/(n * 1.);
  double sigma = S2/(n * 1.) - mean*mean ;
  sigma = sqrt(sigma);
  cerr << "Fwhm mean, sigma : " << mean << " " << sigma << endl ;


  return ;
}

static
void HistoShapeStarFinder(SEStarList const & stlse,  
	       double & mshape, double & bin_shape)
{
  double  minshape = -6. ;
  double  maxshape = -2. ;
  int NBin = (int)((maxshape -minshape)/bin_shape + 0.5 )  ;
  cerr << "N bin " << NBin << endl ;
  Histo1d histo(NBin,minshape ,maxshape); //HardCoded
  
  for (SEStarCIterator it=stlse.begin(); it!=stlse.end(); ++it) 
    { 
      double flux = (*it)->flux ;
      double fluxmax = (*it)->Fluxmax() ;
      if ( ( flux > 1.e-5 ) && ( fluxmax > 1.e-5 ) )
	{
	  double shape = -2.5*log10(flux/fluxmax);
	  histo.Fill(shape , 1. );
	}
    }

  // recherche du mode de l'histo
  histo.MaxBin(mshape);
  cout << "Mode en shape: " << mshape << " (Bin: " 
       << histo.BinWidth() << " = " << bin_shape << " )"
       << endl ; 
  bin_shape = histo.BinWidth();
  return ;
}

 


// The StarFinder needs a clean list. cuts are adaptated to SNLS ref
void StarFinder(SEStarList const & stlse, SEStarList & PrettyStars)
{
  SEStarList stl1;
  stlse.CopyTo(stl1);
  stl1.FluxSort();
  stl1.CutTail(200);
  double mfwhm1, bin_fwhm1 = 0.2 ;
  HistoFwhmStarFinder(stl1, mfwhm1, bin_fwhm1);
  double mshape1, bin_shape1 = 0.1;
  HistoShapeStarFinder(stl1, mshape1, bin_shape1);

  // pour info
  double Mfwhm, Bin_fwhm, Mshape, Bin_shape;
  HistoStarFinder(stlse, Mfwhm, Bin_fwhm, 
		     Mshape, Bin_shape)  ;

  double min_fwhm = mfwhm1 - 2. * bin_fwhm1;
  double max_fwhm = mfwhm1 + 3. * bin_fwhm1;

  // galaxies: shape tres neg, fwhm grande
  // stars : shape petit (pas tres neg), fwhm petit
  double min_shape = mshape1  - 2. * bin_shape1;
  double max_shape = mshape1  + 3.  * bin_shape1 ;

  cerr << "fwhm min, max: " << min_fwhm << " " << max_fwhm << endl ;
  cerr << "shape min, max: " << min_shape << " " << max_shape << endl ;

  for (SEStarCIterator it= stlse.begin(); it!= stlse.end(); ++it)
    {
      const SEStar *pstar = *it ;
      double flux = pstar->flux;
      double fluxmax = pstar->Fluxmax();
      double fwhm = pstar->Fwhm();
      double shape = -2.5*log10(flux/fluxmax);
      if ( (fwhm < max_fwhm ) && ( fwhm > min_fwhm ) &&
	   (shape < max_shape ) && ( shape > min_shape))
	{
	  PrettyStars.push_back(new SEStar(*pstar));
	}
    }

}
void GalaxyFinder(SEStarList const & stlse, 
		  SEStarList & PrettyGal, 
		  float maglim, float zerop)
{
  SEStarList stl1;
  stlse.CopyTo(stl1);
  stl1.FluxSort();
  stl1.CutTail(200);
  double mfwhm1, bin_fwhm1 = 0.2 ;
  HistoFwhmStarFinder(stl1, mfwhm1, bin_fwhm1);
  double mshape1, bin_shape1 = 0.1;
  HistoShapeStarFinder(stl1, mshape1, bin_shape1);

  // pour info
  double Mfwhm, Bin_fwhm, Mshape, Bin_shape;
  HistoStarFinder(stlse, Mfwhm, Bin_fwhm, 
		     Mshape, Bin_shape)  ;

  double shape_thresh = mshape1  - 2. * bin_shape1 ;
  double fwhm_thresh = mfwhm1 + 3. * bin_fwhm1;

  cerr << " fwhm_thresh, shape_thresh: " << fwhm_thresh << " " << shape_thresh << endl ;

  for (SEStarCIterator it= stlse.begin(); it!= stlse.end(); ++it)
    {
      const SEStar *pstar = *it ;
      double flux = pstar->flux;
      double fluxmax = pstar->Fluxmax();
      double fwhm = pstar->Fwhm();
      double shape = -2.5*log10(flux/fluxmax);
      bool okmag = false ;
      // objects fainter than maglim are galaxies
      // anyway
      if (( maglim > 0) && (zerop > 0))
	{
	  double mag =  -2.5*log10(flux) + zerop ;
	  okmag = (mag > maglim);
	}
      if (  (( fwhm  > fwhm_thresh) &&  
	     ( shape < shape_thresh) ) || okmag )
	{
	  PrettyGal.push_back(new SEStar(*pstar));
	}
    }

}

void HistoStarFinder(SEStarList const & stlse,  
	       double & mfwhm, double & bin_fwhm, 
	       double & mshape, double & bin_shape)
{
  double nstarperbin = 10. ;
  int NBin = max(min(100,int(stlse.size()/nstarperbin)),1);
  Histo2d histo(NBin,0,10,max(NBin/2,1),-6,-2); //HC
  
  for (SEStarCIterator it=stlse.begin(); it!=stlse.end(); ++it) 
    { 
      histo.Fill((*it)->Fwhm() ,-2.5*log10((*it)->flux/(*it)->Fluxmax()), (*it)->Fluxmax() );
    }

  // recherche du mode de l'histo
  histo.MaxBin(mfwhm, mshape);
  histo.BinWidth(bin_fwhm, bin_shape);
  cout << "Mode en fwhm: " << mfwhm << " (Bin: " << bin_fwhm
       << "), Mode en mag-mu: " << mshape <<  " (Bin: " 
       << bin_shape <<  ")" << endl ;
  return ;
}



#include <cmath>
#include "matvect.h"

static double sq(const double &x) { return x*x;}



/* one more star finder:
1) locate a peak (the maximum) in fwhm%log(flux/fluxmax) histogram

2) try to figure out the span of the peak by looking at the bin contents
   around the peak.

3) "fit" a 2D gaussian on the data itself, using an iterative approach:

   using a guess of the peak extension (vxx, vyy, vxy),
   you may compute weighted second moments (sum(x**2*w)/sum(w)),
   with a gaussian w. You don't get the right answer, but if w
   has the same width as the data, the outcome is a variance
   which is half of the true one. There are several inversions
   of 2D symetric matrices done on the flight to switch from covariance 
   to weight matrix. 


   Once we have located the peak and its extension
   in the plane, we can finally select objects whithin an elliptic
   aperture around the star peak, and call them stars.
   If you don't need the stars, you may however consider the Fwhm
   that comes out from this routine: it is somehow better that the
   clipped average done in seeing_box.cc

*/


#include "histopeakfinder.h"

void RoughStarFinder(const SEStarList &In, double &Fwhm,
		     const double SigCut, SEStarList *Out)
{

  Histo2d histo(100,0,10.,30,2.,6.);
  StarScoresList scores;
  for (SEStarCIterator i = In.begin(); i != In.end(); ++i)
    {
      const SEStar &s = **i;
      if (s.Flag() !=0 || s.FlagBad() != 0) continue;
      if (s.flux < 10*s.EFlux()) continue; // no resolution on shape
      double xx = s.Fwhm();
      double yy = log(s.flux/s.Fluxmax());
      histo.Fill(xx,yy, s.Fluxmax());
      scores.push_back(new StarScores(xx,yy, &s));
    }
  double fwhmMax, shapeMax;
  histo.MaxBin(fwhmMax, shapeMax);
  Ellipse ellipse;
  if (!HistoPeakFinder(scores, histo, fwhmMax, shapeMax, ellipse))
    {
      cout << " WARNING : HistoPeakFinder did not work properly and " 
	   << " RoughStarFinder does not handle it " << endl;
    }

  double meanShape;
  ellipse.GetCenter(Fwhm, meanShape);

  if (Out)
    {
      Out->clear();
      for (StarScoresCIterator i = scores.begin(); i != scores.end(); ++i)
	{
	   const StarScores &star = **i;
	   if (star.nSig<SigCut) // it lies within the star cluster
	     {
	       const SEStar &s = dynamic_cast< const SEStar &>(*star.star);
	       Out->push_back(&s);
	     }
	}
    }
}

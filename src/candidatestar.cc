#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "sestar.h"
#include "image.h"
#include "candidatestar.h"
#include "dodetection.h"


#ifdef STORAGE
double
Noise_Aperture(Image  const & img, double x, double y, double radius, 
	       double Fond, double SigFond)
{
  /*ddouble r1 = 2.4* radius;
  double r2 = 3.2* radius;
  double fraclow = 0.15 ;
  double frachigh = 0.15 ;
  double fd_loc = FondLocal( img, x, y, r1, r2, fraclow, frachigh);
  cout << "FondLoc: " << fd_loc << " " << Fond << endl ; */
  int tx = img.Nx();
  int ty = img.Ny();
  double S = 0 ;
  double Sig2 = SigFond* SigFond ;
  int dtaille = (int) ( radius + 1. );
  int ix = (int) (x+0.5);
  int iy = (int) (y+0.5);
  int ntot = 0 ;
  for(int i = -dtaille; i <=dtaille; i++)
    {

      for(int j = -dtaille; j <=dtaille; j++)
	{
	  int xx = ix + i ;
	  int yy = iy + j ;
	  double d = sqrt( i*i + j*j );
	  if ( ( xx > 0 ) && ( xx < tx ) &&
	       ( yy > 0 ) && ( yy < ty )&&
	       (d < radius) )
	    {

	      double v = img(xx,yy)- Fond;
	      if (  v < 1.e-10 )
		v = 0. ;
	      S += v ;
              ntot++;
	    }

	}

    } 

  S += ntot *  Sig2 ;
  if ( S > 1.e-10 )
    S = sqrt(S);
  else
    S = 0. ;
  /*cerr << "Surface " << 3.14*radius*radius << " " << ntot << endl ;
  cerr << "Sigma " << SigFond  << endl ;
  cerr << "S = " << SigFond*sqrt(3.14)*radius<< " " << S<< endl << endl ;*/

  return(S);
}

#endif
//********************   DEFINITION  CandidateStar   *********************

// Converter :
BaseStarList* Candidate2Base(CandidateStarList * This)
{ return (BaseStarList*) This;} 

const BaseStarList* Candidate2Base(const CandidateStarList * This)
{ return (BaseStarList*) This;} 




CandidateStar::CandidateStar()
{
  Set_to_Zero();
}

CandidateStar::CandidateStar(SEStar const & sestar) 
  : SEStar(sestar)
{
  Set_to_Zero();
}


void
CandidateStar::Set_to_Zero()
{
  numero = 0;
  fluxcv = 0 ;
  noise = 0 ;
  sigtonoise = 0;
  faper_ref = 0;
  prctincrease = 0;
  sigx=0;
  sigy=0;
}




void
CandidateStar::dumpn(ostream& s) const
{
  SEStar::dumpn(s);
  s << " numero : " <<   numero ;
  s << " flux cv : " << FluxCv()  ;
  s << " noise : " << noise ;
  s << " sigtonoise   : " <<  sigtonoise ;
  s << " faper_ref    : " <<  faper_ref ;
  s << " prctincrease : " <<  prctincrease ;     
  s << " Sigx : " << sigx;  
  s << " Sigy : " << sigy;
  s << " Sigxy : "<< sigxy;
  s << " Chi2 : "<< chi2;
}


void
CandidateStar::dump(ostream& s) const
{
 dumpn(s);
 s << endl ;
}


void
CandidateStar::writen(ostream& s)  const
{
  SEStar::writen(s); 
  s  << resetiosflags(ios::scientific) ;
  s  << setiosflags(ios::fixed) ;
  s  << setprecision(8) ;
  s  << numero << " " ;
  s  << FluxCv() << " " ;
  s  << noise << " " ;
  s  << sigtonoise << " " ;
  s  << Faper_ref() << " " ;
  s  << prctincrease << " "  ;  
  s  << sigx << " " ;
  s  << sigy << " " ;
  s  << sigxy << " ";
  s  << chi2 << " ";
}


void
CandidateStar::read_it(istream& r, const char * Format)
{
  SEStar::read_it(r, Format); // lu avec format SEStar, y compris alignement
  r >> numero ;
  r >> FluxCv();
  r >>  noise ;
  r >>  sigtonoise ;
  r >>  faper_ref ;
  r >>  prctincrease ;
  r >> sigx;
  r >> sigy;
  r >> sigxy;
  r >> chi2;
  //int format = DecodeFormat(Format, "CandidateStar");
  return ;
}

CandidateStar*  CandidateStar::read(istream& r, const char *Format)
{
  CandidateStar *pstar = new CandidateStar();  
  pstar->read_it(r, Format);
  return(pstar);
}


string CandidateStar::WriteHeader_(ostream & pr, const char*i) const
{
  if (i==NULL) i = "";
  string  seStarFormat =  SEStar::WriteHeader_(pr, i);
  pr    << "#numero"<< i <<" : candidate number " << endl 
	<< "#fluxcv"<< i <<" : flux estime par convolution " << endl     
	<< "#noise"<< i <<" : noise" << endl    
	<< "#stnoise"<< i <<" : signal to noise " << endl 
	<< "#fap_r"<< i <<": aperture flux on ref" << endl
	<< "#prcti"<< i <<" : perctage increase" << endl
	<< "#sigx"<< i <<": x error " << endl
	<< "#sigy"<< i <<": y error " << endl 
	<< "#sigxy"<< i <<": correlation between x and y " << endl 
	<< "#chi2"<< i <<": chi2 for flux estimation " << endl ;
  
  
  string format = format +  " CandidateStar 1 ";
  /* 1 is the current format id for CandidateStars (when being written) it must correspond
     to the right behaviour of the read routine ( and match what write does ! ) */
  return format; 
}



void CandidateStar::write_scan(ostream &pr) const
{
  pr << x << " " 
     << y << " " 
     << Fluxmax() << " "
     << flux << " "
    //<< FluxCv()   << ' ' pas dans le header ! 
     << Noise()   << ' '
     << SigToNoise()   << ' ' 
     << Faper_ref()   << ' ' 
     << PrctIncrease()   << endl;

}



void CandidateStar::WriteHeader_Scan(ostream &pr)
{
  pr << "# x : x position (pixels)" << endl 
     << "# y : y position (pixels)" << endl 
     << "# amp : Peak pixel value above background on sub" << endl 
     << "# flux : flux in pixel units on sub " << endl ;
  pr << "# noise : noise (for sub) " << endl;  
  pr << "# stnoise : S/N for obj. detecetd on sub " << endl;
  pr << "# fap_r : aperture flux computed on ref at candidate position " << endl; 
  pr << "# prcti : percentage increase ie flux/fap_r " << endl; 
  pr << "# end" << endl;
}




bool  
CandidateStar::IsGood(double cut_sigtonoise)
{
  //cut signal-to-noise
  if (SigToNoise()< cut_sigtonoise)
    return false ;
  return(true);
}

bool  
CandidateStar::IsGood(double cut_sigtonoise, double frac_faper_ref)
{
  
  // Cut on signal to noise
  if (!IsGood(cut_sigtonoise))
    return false;

  // Cut on the second momentum
  if ( Mxx() < 1.e-10 || Myy() < 1.e-10 )
    return false ;

  // Test on the flag on the object (saturated, bad pixels,...)
  if ( FlagBad() != 0 )
    return false ;  

  //cut prctage increase . first check that the flux on ref is > 0.  
  if (Faper_ref() >0 &&  flux < frac_faper_ref * Faper_ref())
    return false ;

  return(true);
}


void  
CandidateStar::ComputeFaperRef(Image & imgref, double rayon,double  fd_ref)
{

  // on prend le fond moyen
  // a 0 de ttes facons normalement
  double X = x - DECALAGE_IJ_SE ;
  double Y = y - DECALAGE_IJ_SE ;
  int npix;
  Faper_ref() = Flux_Aperture(imgref, X, Y, rayon, fd_ref,npix  );

  if (fabs(Faper_ref()) > 1.e-10)
    PrctIncrease() = Flux_fixaper()/Faper_ref();
  
  return;

}

void
CandidateStar::ComputeNoise(Image & imgref, Image & imgnew, double rayon_ref, double rayon_new,
			    double sigma_ref, double sigma_new, double fd_ref, double fd_new)
{
  // passer des coorrdonnees SExtracor aux coordonnees images
  double X = x - DECALAGE_IJ_SE ; 
  double Y = y - DECALAGE_IJ_SE ;
  double Nr = Noise_Aperture(imgref, X, Y, rayon_ref, 
			     fd_ref, sigma_ref);
  double Nn = Noise_Aperture(imgnew, X,Y, rayon_new, 
			     fd_new, sigma_new);
  Noise() = sqrt(Nr*Nr + Nn*Nn);
  
  if (Noise() > 1.e-10)
    SigToNoise() = flux/Noise();
  
}


//********************   FIN DEFINITION CandidateStar   *********************

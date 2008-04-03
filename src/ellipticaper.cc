#include <iostream>
#include <cmath> // asin, sqrt


#include "ellipticaper.h"
#include "sestar.h"
#include "apersestar.h"
#include "image.h"






void Elliptic_Aperture::WriteHeader_(ostream & pr, const char* i) const 
{
	if (i== NULL) i= "";
	pr << "# xc" << i << " : " << endl 
	   << "# yc" << i << " : " << endl 
	   << "# a" << i << " : " << endl 
	   << "# b"<<i<< " :  " << endl 
	   << "# t" << i << " : " << endl 
	   << "# cxx" << i << " : " << endl 
	   << "# cyy" << i << " : " << endl 
	   << "# cxy" << i << " : " << endl 
	   << "# krf" << i << " : " << endl 
	   << "# fixr" << i << " : " << endl 
	   << "# bck" << i << " : " << endl 
	    << "# flux" << i << " : " << endl 
	   << "# eflux" << i << " : " << endl 
	   << "# nbad" << i << " : " << endl
	   << "# isc" << i << " : " << endl
	   << "# isg" << i << " : " << endl;}

 
void Elliptic_Aperture::writen(ostream & pr) const {pr << " " << xc << " " << yc << " " << a << " " << b << " " << angle <<  " " << cxx << " " << cyy << " " << cxy << " " << kron_factor << " " << fixradius << " " << background << " " << flux << " " << eflux << " " << nbad << " " << is_circle << " "  << is_good << " " ;}


void Elliptic_Aperture::read(fastifstream& rd ) 
{rd >>xc >>yc >>a >>b >>angle >> cxx >>cyy >>cxy >>kron_factor >>fixradius >>background >>flux >>eflux >>nbad >> is_circle  >> is_good ;}





void Elliptic_Aperture::SetParameters(const SEStar & star, 
				      const double dilatation, 
				      const double RadMin,
				      const double Radius, bool fromabtheta)
{
  xc = star.X();
  yc = star.Y() ;
  a = star.A() ;
  b = star.B() ;
  angle = star.Gyr_Angle() ;
  if ( star.Gyr_Angle() > 0 )
    angle = star.Gyr_Angle();
  else
    angle =star.Gyr_Angle() + 180. ;

  // parametres de l'ellipse d'ouverture: repris de scan.c 
  // resume : les moments mxx, myy, mxy sont calcules d'apres donnees.
  // on definit x' et y' dans le ref rotate de theta. quel theta maximise < x'^2> --> donne theta, A, et B
  // l'ellispe de definition de l'objet est alors, avec X = x-xc et Y = y - yc :
  // cxx . X^2 + cyy . Y^2 + cxy XY = 1, ellipse de 1/2 grd axe A, de 1/2 petit axe B,  et d'inclinaison theta
  // cela donne  :
  if ( !fromabtheta)
    {
      double temp2 =  star.Mxx()* star.Myy()- star.Mxy()* star.Mxy(); // normalement protege dans SExtractor (scan.c) contre valeur trop petite
      cxx =  star.Myy() /temp2;
      cyy =  star.Mxx()/temp2;
      cxy = -2.* star.Mxy()/temp2;
    }
  else
  // si on a que a, b, et theta pour calculer (e.g. catalogue terapix)
    {
      double thet = star.Gyr_Angle()*M_PI/180.;
      cxx = cos(thet)*cos(thet)/(star.A()*star.A()) + sin(thet)*sin(thet)/(star.B()*star.B()) ;
      cyy = cos(thet)*cos(thet)/(star.B()*star.B()) + sin(thet)*sin(thet)/(star.A()*star.A()) ;
      cxy = 2*cos(thet)*sin(thet)*(1./(star.A()*star.A()) - 1./(star.B()*star.B())) ;
    }

  kron_factor =  star.Kronradius()* dilatation; 
  fixradius = -1 ;
  is_circle=-1;
  double radmin = 0 ;
  if (RadMin > radmin)
    radmin = RadMin;
  radmin *= dilatation;
  if (kron_factor * sqrt( star.A()* star.B()) <=  radmin)
    {
      cxx = 1. ;
      cyy = 1. ;
      cxy = 0. ;
      fixradius = Radius*dilatation ;
      is_circle=1;
    }
}
// repris de  photom.c (SExtractor)
  // cxx . X^2 + cyy . Y^2 + cxy XY = 1, ellipse de 1/2 grd axe A, de 1/2 petit axe B,  et d'inclinaison theta
  // on va integrer sur cette ellipse dilatee, ie blabla = R^2 avec R >1.
  // definition de R = ce que E.B. appelle le kronfactor, qui est 
// "le kron radius calcule en unite de A ou B"  MULTIPLIE  PAR 
// le facteur PHOT_AUTOPARAM[0] definit dans la datacard (2.5 ds notre exemple et dans le default)
// PHOT_AUTOPARAMS 2.5 3.5
//NB : PHOT_AUTOPARAM[1] est une valeur minimale utilisee lors du calcul 
// du kron_factor par SExtractor. Si kronfactor < PHOT_AUTOPARAM[1], alors
// kronfactor = PHOT_AUTOPARAM[1].

// kron_factor = 2.5 * kron_radius, si < 3.5, alors = 3.5

// le kron_factor est bien ce qu'on recupere de SExtractor mais qu'on appelle malencontreusement Kronradius() dans SEStar.

// Apres le calcul du kron_factor, il y a eu dans SE la verification suivante de faite :
// if (kron_factor * sqrt(A*B) <=  PHOT_AUTOAPER_[1] * 0.5 ) alors kron_factor =0.0 
// specifies dans datacard PHOT_AUTOAPERS : mag_auto minimum circ. aperture diameter: estimation disk, measurement disk
// (non precise dans notre datacard donc valeur par defaut : 0., 0.)
// Il est alors prevu dans ce cas d'integrer sur  un cercle de rayon Radius = PHOT_AUTOAPER[1] * 0.5 en pixel.
// ie avec une equation X^2 + Y^2 = Radius^2, attention Radius est alors en pixels !
// ici on teste donc sur kron_factor > RadMin, mis a 0 par defaut 

//NB: PHOT_AUTOAPER_[0] sert sur le meme genre de test lors de l'estimation de la taille du dsique qui servira a calculer le kronfactor.

// dans SE, le calcul est donc fait avec dilatation =1., et si PHOT_AUTOAPERS
// n'est pas precise, alors Radius =0. 


void Elliptic_Aperture::computeflux(const Image& I, const Image& W, 
				    const Image *pC, const Image *pS, 
				    const double Gain, int segmentation_number, double scale_fact)
{
  
  double dxlim =0., dylim =0., klim2 = 0. ;
  int  nncos=0, nnbad =0; // a modifier si oversampling en fraction de pixel
  double ffcos=0;
  double var = 0., fflux = 0. , ffother =0.;
  double area =0. ; // en preparation de 1) back !=0, 2) oversampling en fraction de pixel


  int thisStarNumber = segmentation_number; // the S (Segmentation) image contains the "num" value of pixels that Sextractor attributed to objects

  if (! IsCircle())
    //if (kron_factor * sqrt(A()*B()) >  RadMin)
    //if (kron_factor > 0. )
    {
      // ellipse inscrite dans un rectangle (-dxlim, +dxlim, -dylim, +dylim)
      dxlim = cxx - cxy * cxy /(4. * cyy) ;
      if (dxlim > 0)
	  dxlim = kron_factor*scale_fact / sqrt(dxlim);
      else
	dxlim = 0. ;
      dylim = cyy - cxy * cxy /(4. * cxx) ;
      if (dylim > 0)
	  dylim = kron_factor*scale_fact / sqrt(dylim);
      else
	dylim = 0. ;
      klim2 = kron_factor*kron_factor*scale_fact*scale_fact ;
    }
  else // il etait prevu d'integrer sur un cercle de rayon  PHOT_AUTOAPER_1 * 0.5, ici fixradius en pixels
    {
      // on doit avoir cxx = cyy = 1., cxy = 0.0 ; 
      dxlim = dylim = fixradius*scale_fact;
      klim2 = dxlim * dylim ;
    }
  int xmin = (int) (xc - dxlim) ;
  int xmax = (int) (xc + dxlim +1) ;
  int Nx = I.Nx() ;
  int ymin = (int) (yc - dylim) ;
  int ymax = (int) (yc + dylim +1) ;
  int Ny = I.Ny() ;
  
   for (int y1 = ymin; y1 < ymax ; y1++)
    {
      for (int x1 = xmin; x1 < xmax ; x1++)
	{
	  
	  double dx = x1 - xc ;
	  double dy = y1 - yc ;
	  if ( ( cxx * dx * dx + cyy * dy *dy + cxy * dx *dy ) <= klim2 ) // on est dans l'ellipse
	    {
	      // calcul de la fraction de pixel par oversampling ici si on veut
	      if ( ( x1 < 0 ) || (y1 < 0 ) || ( x1 >= Nx ) || (y1 >= Ny ) )
		{
		  nnbad += 1;
		  continue;
		  }
	      double w = W(x1,y1);
	      double pixVal =  I(x1,y1) ;
	      
	      if (w > 0 )
		{
		  int seg =  0 ;
		  if (pS != NULL) seg = int( (*pS)(x1,y1) ) ;
		  if ( (seg ==0) || (seg == thisStarNumber) )
		    {
		      
		      var += 1. / w ;
		      fflux += pixVal ;
		      area += 1. ;
		    }
		  else
		    ffother+= pixVal;
		}
	      else
		{
		  nnbad += 1;
		  if( (pC != NULL) && ((*pC)(x1,y1)>0 ) ) 
		    {
		      nncos += 1;
		      ffcos += pixVal ;
		    }  
	      
		}
	    }
	}
    }
  fflux -= area * background ; 
  ffother -= area * background ; 
  if (Gain > 0 && fflux > 0) // j'ajoute la condition fflux >0 par rapport a SExtractor
    var += fflux/Gain ;
  if (var > 0. )
    eflux = sqrt(var) ;
  else
   eflux  = 0. ;
  flux = fflux ;
  //fother = ffother ;
  nbad = nnbad ;
  ncos = nncos ;
  fcos = ffcos ;
  return ;
}

#define N_SIGMAS 5.  
#define Q_AMIN 4.
#define MAX_BINS 4096
#define MIN_PIX_FRAC 0.1
#include "imageback.h"
#include "histo1d.h"

static void WClippedMean(double *Itab, double *Wtab, double &mean, 
			double &sigma, int & ntab)
{
  double sum = 0;
  double sum2 = 0;
  double sumw = 0;
  double value,w;
  int ntot =0 ;
  double mmean =0., sig=0. ; ;
  for (int i=0; i < ntab; ++i)
    {
      w = Wtab[i];
      value = Itab[i];
      sum += value*w;
      sum2 += value*value*w;
      sumw += w;
    }
  if (sumw)
    {
      mmean  = sum/sumw;
      sig = (sum2/sumw - mmean*mmean); if (sig>0) sig = sqrt(sig); else sig = 0;
    }
  else { sigma = 1; mean = 0; ntab = 0; return;}
  int loop;
  for (loop=2; loop>0; loop--)
    {
      double lcut =  -2.0*sig;
      double hcut =  2.0*sig;
      sum = sum2 = 0; sumw = 0; 
      for (int i=0; i < ntab; ++i)
	 {
	   value = Itab[i]-mmean;
	   w = Wtab[i];
	   if (w == 0) continue;
	   if (value > hcut || value < lcut) continue;
	   sum += value*w;
	   sum2 += value*value*w;
	   sumw += w;
	   ntot++;
	 }
     if (sumw)
       {
	 sum =  sum/sumw;
	 sum2 = sum2/sumw - sum*sum;
	 mmean += sum;
	 sig = (sum2>0) ? sqrt(sum2) : 0.0;
       }
     else break; // keep "old" values
   }
  mean = mmean;
  sigma = sig;
  ntab = ntot ;
}

static void ClippedMean(double *Itab,double &mean, 
			double &sigma, int & ntab)
{
  double sum = 0;
  double sum2 = 0;
  double value;
  double mmean =0., sig=0. ; ;
  for (int i=0; i < ntab; ++i)
    {
      value = Itab[i];
      sum += value;
      sum2 += value*value;
    }
  if (ntab)
    {
      mmean  = sum/ntab;
      sig = (sum2/ntab - mmean*mmean); 
      if (sig>0) sig = sqrt(sig); else sig = 0;
    }
  else { sigma = 1; mean = 0; ntab = 0; return;}
  int loop;  
  int ntot =0 ;
  for (loop=2; loop>0; loop--)
    {
      ntot =0 ;
      double lcut =  -2.0*sig;
      double hcut =  2.0*sig;
      sum = sum2 = 0; 
      for (int i=0; i < ntab; ++i)
	 {
	   value = Itab[i]-mmean;
	   if (value > hcut || value < lcut) continue;
	   sum += value;
	   sum2 += value*value;
	   ntot++;
	 }
     if (ntot)
       {
	 sum =  sum/ntot;
	 sum2 = sum2/ntot - sum*sum;
	 mmean += sum;
	 sig = (sum2>0) ? sqrt(sum2) : 0.0;
       }
     else break; // keep "old" values
   }
  mean = mmean;
  sigma = sig;
  ntab = ntot ;
}



double  Elliptic_Aperture::SqEllipticDistance(double x1, double y1)
{	  
  return( SqEllipticDistance(xc,yc,x1,y1));
}

double Elliptic_Aperture::SqEllipticDistance(double x_or, double y_or, double x1, double y1)
{
  double dx = x1 - x_or ;
  double dy = y1 - y_or ;
  return( cxx * dx * dx + cyy * dy *dy + cxy * dx *dy);
}



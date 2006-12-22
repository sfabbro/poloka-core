#include "agaussian.h"
#include "image.h"
#include "basestar.h"

#include <cmath> 

 
AGaussian::AGaussian()
{
  xc=0.; 
  yc=0.; 
  sigma_x=0.;
  sigma_y=0.; 
  rho=0.;
  flux =0.;
  fond=0.;
}

double 
AGaussian::ExpValue(double xin, double yin) const
{ 
  double ux = (xin-xc) / sigma_x;
  double uy = (yin-yc) / sigma_y ; 
  double res = exp( -0.5 * (ux * ux + uy * uy - 2. * rho * ux * uy) ) ;
  return(res);
}


//! integrated value over pixel(i,j)
// ie. int_{x=i-0.5, x=i+0.5}_{y=j-0.5, x=j+0.5}
#define NPAS 10
double 
AGaussian::ExpIntegValue(int i, int j)  const
{ 
  double I = 0 ;
  double dl = 1./NPAS ;
  // on divise en NPAS sous-intervalles [i-0.5,i+0.5], de largeur dl
  // on prend pour chaque sous-intervalle la valeur au milieu x dl.
  double xbeg = i-0.5 + 0.5 * dl ;
  double xx = xbeg ; 
  double yy = j-0.5 + 0.5 * dl ;
  for(int ii=0; ii< NPAS; ii++) 
    {
      for(int jj=0; jj< NPAS; jj++) 
	{
	  I += ExpValue(xx,yy);
	  xx += dl ; 
	}
      yy+=dl ;
      xx = xbeg ;
    }
  I *= dl*dl;
  return(I);
}

double 
AGaussian::window_size() const
{
  double prec = 0.01 ;
  double sigma = max(sigma_x,sigma_y); // on ajoute un carre.
  // true formula is : sqrt(2)*sigma*sqrt(log(flux/2*pi*sigma^2*prec))
  double u = flux/(6.3 * sigma * sigma * prec);
  // le flux est tellement petit que meme la valeur au max est 
  // plus petite que prec.
  if ( u < 1) 
    return(0.);
  else
    {
      double d = 1.4 * sigma * sqrt(log(u)) ;
      return(d);
    }
}

void AGaussian::AddToImage(Image &image, 
			  Image * psat,
			  double saturation) const
{
  // pour ne pas ecrire le message de saturation 36 fois
 int flagsaturated = 0 ;
 double w = window_size();
 //cerr << " window_size " << w << endl ;
  int xstart = int(xc - w + 0.5);
  int ystart = int(yc - w + 0.5);

  int taille =  int(2*w + 1. + 0.5);

  int xend = min(xstart + taille,image.Nx()) ;
  int yend = min(ystart + taille,image.Ny()) ;

  xstart = max(xstart,0);
  ystart = max(ystart,0);

  //   gauss = new AGaussian(XF,YF,Fwhm_i()/(sqrt(2.*log(2.))*2.),Flux_i()*photFactor);
  double normalisat = Norma();

  for (int j = ystart; j < yend; ++j)
    {     
      for (int i = xstart; i < xend; ++i)
	{
	  double val = ExpIntegValue(i,j) * normalisat + fond ;

	  if ( (saturation > 0 ) && 
	       (image(i,j)+ val > saturation ) )
	    {
	      image(i,j) = saturation ;
	      if (psat != NULL)
		{
		  if ((flagsaturated == 0 ) && ((*psat)(i,j) == 0 ))
		    // not saturated before, saturated because of SNe
		      {
			cerr << "SN is saturated (i= " 
			     << i << ", j= " 
			     << j << "), saturation = " 
			     << saturation << ", image = " 
			     << image(i,j) << ", new image val = " 
			     << image(i,j)+ val  << endl ;
			flagsaturated = 1;
		      }
	
		    
		  (*psat)(i,j) = 1 ;
		}
	      else
		if (flagsaturated == 0 )
		  {
		    cerr << "SN or SN zone  is saturated: " << endl ;
		    flagsaturated = 1;
		  }
	      
	    }
	  else
	    image(i,j) += val ;
	}
    }
}



void AddListWGaussianToImage(double sigmax, double sigmay, double rho,
			     BaseStarList *list,
			     Image & dest,   
			     Image * psat, 
			     double satlevel) 
{
  AGaussian gauss;
  gauss.sigma_x = sigmax ;
  gauss.sigma_y = sigmay ;
  gauss.rho = rho  ;
  for (BaseStarCIterator i = list->begin(); i != list->end(); ++i)
    {
      gauss.xc = (*i)->x ;
      gauss.yc =  (*i)->y ;
      gauss.flux = (*i)->flux ;
      gauss.AddToImage(dest, psat, satlevel);
    }
}

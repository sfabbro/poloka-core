
#include "sestar.h"
#include "imageback.h"
#include "fitsimage.h"
#include "fileutils.h"
#include "convolution.h"
#include "candstar.h"
#include "frame.h" 
#include "datacards.h"
#include "dodetection.h"

static double sqr(double x){return x*x;}

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif
// ----------------------------------- 
// ------ routines generalistes ------
// ----------------------------------- 
#ifdef STORAGE
void
Sig_Position(Image const & img ,double const &rad_flux, 
	     double & xbar, double & ybar, const double & fond, 
	     double & Sigx, double & Sigy, double & Sigxy)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = (int) ((rad_flux + 1.));
  double sx=0;
  double sy=0;
  double S=0;
  double SumSigx = 0, SumSigy = 0, SumSigxy = 0, 
    Sx = 0, Sy = 0, Sxy;
  double rad_flux2 = rad_flux * rad_flux ;
  int ix = (int) (xbar + 0.5 ) ;
  int iy = (int) (ybar + 0.5 ) ;
  for (int l=-d; l <= d; l++)
    {
      for (int m=-d; m <=d ;  m++)
	{
	  int xu = ix + l ; 
	  int yu = iy + m ;
	  double dist2 = l*l+m*m;
	  if ( ( xu < tx ) &&  (xu >= 0 )&&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist2 <rad_flux2) )
	    {
	      double val = img(xu,yu) ;
	      sx += l*val;
	      sy += m*val;
	      SumSigx += sqr(l) * val;
	      SumSigy += sqr(m) * val;
	      SumSigxy += l * m * val;
	      S +=  val;
	    }
	}
    }
      Sx = SumSigx / S - sqr(sx/S);
      Sy = SumSigy / S - sqr(sy/S);
      Sxy = SumSigxy / S - sx*sy/sqr(S);
  if (Sy>0 && Sy>0)
    {
      Sigx = Sx;
      Sigy = Sy;
      Sigxy = Sxy;
    }
  return;	  


}
#endif
void Barycentre(Image const & img , double x, double y , double rad_flux,
		double fond, double & xbar, double & ybar, double & flux )
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = int(rad_flux + 1);
  double sx=0;
  double sy=0;
  double s=0;
  for (int l=-d; l < d+1; l++)
    {
      for (int m=-d; m < d + 1;  m++)
	{
	  int xu , yu ;
	  xu = (int) ( x + l ); // sert pas a grd chose mais plus clair
	  yu = (int) ( y + m );
	  double dist = sqrt(double(l*l+m*m));
	  if ( ( xu < tx ) &&  (xu >= 0 )&&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist <rad_flux) )
	    {
	      double val =(double) (img(xu,yu) - fond );
	      
	      s += val ;
	      sx += xu*val;
	      sy += yu*val;
	    }
	}
    }
  flux = s ;
  if (s>0) 
    {
      xbar   = sx/s;
      ybar   = sy/s;
    }  
  return;
}




bool
IsNotBadPix( Image const & img , Image const & mask, int x , int y , 
	     float  rad_bad)
{

  int i , j ;

  double rad_bad2 = rad_bad * rad_bad ;
  int irad_bad = (int) ( rad_bad + 1. );

  
  int tx = img.Nx() ;
  int ty = img.Ny() ;


  for (i=-irad_bad; i<=irad_bad; i++)
    for (j=-irad_bad; j<=irad_bad; j++)
      {
	double d2 = i*i + j*j ;
	int xx = x + i ;
	int yy = y + j ;
	if ( ( xx > 0 ) && ( xx < tx ) &&
	     ( yy > 0 ) && ( yy < ty ) &&
	     ( d2 < rad_bad2 ) )
	  if (mask(xx,yy) > 0 ) 
	      return  false;
      }

  return true;
}


// test si (x,y) ( qui est entier, on boucle sur des pixels)
// si c'est un maximum local ds un rayon rad_locmax, 

static bool IsLocalMax( Image const & img , int x , int y , int half_square_size)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  float val;
  if ( ( x > 0 ) && ( x < tx ) &&
       ( y > 0 ) && ( y < ty ) )
    val = img(x,y);
  int i , j ;
  
  for (i=-half_square_size; i<=half_square_size; i++)
    for (j=-half_square_size; j<=half_square_size; j++)
      {
	if ((i==0) && (j==0)) continue;
	int xx = x + i ;
	int yy = y + j ;
	if ( ( xx > 0 ) && ( xx < tx ) &&
	     ( yy > 0 ) && ( yy < ty ) )
	  if (img(xx,yy)>val) 
	    return false;
      }
  return true;
}

double Flux_Aperture(const Image   & img, const double x, const double y, const double radius, const double Fond, int &Npix)
{
  int tx = img.Nx();
  int ty = img.Ny();
  double S = 0 ;
  int dtaille =(int) ( radius + 1. );
  int ix = (int) (x+0.5);
  int iy = (int) (y+0.5);
  Npix = 0;
  for(int i = -dtaille; i <=dtaille; i++)
    {
 
      for(int j = -dtaille; j <=dtaille; j++)
        {
          int xx = ix + i ;
          int yy = iy + j ;
          double d = sqrt( double(i*i + j*j) );
          if ( ( xx > 0 ) && ( xx < tx ) &&
               ( yy > 0 ) && ( yy < ty )&&
               (d < radius) )
            {
 
              double v = img(xx,yy)- Fond;
              S += v ;
              ++Npix;
            }
 
        }
 
    }
  return(S);
}              

// attention: ne prd que les pixels au essus du fond
// hyp: Fond + (Signal >0) sur l'image
// + pas de gain !

double
Noise_Aperture(Image  const & img, double x, double y, double radius, 
	       double Fond, double SigFond)
{
  /*
   double r1 = 2.4* radius;
   double r2 = 3.2* radius;
   double fraclow = 0.15 ;
   double frachigh = 0.15 ;
   double fd_loc = FondLocal( img, x, y, r1, r2, fraclow, frachigh);
   cout << "FondLoc: " << fd_loc << " " << Fond << endl ; */
  int tx = img.Nx();
  int ty = img.Ny();
  double S = 0 ;
  double Sig2 = SigFond* SigFond ;
  int dtaille =(int) ( radius + 1. );
  for(int i = -dtaille; i <=dtaille; i++)
    {

      for(int j = -dtaille; j <=dtaille; j++)
	{
	  int xx = (int) (x + i) ;
	  int yy = (int) (y + j) ;
	  double d = sqrt( double(i*i + j*j) );
	  if ( ( xx > 0 ) && ( xx < tx ) &&
	       ( yy > 0 ) && ( yy < ty )&&
	       (d < radius) )
	    {

	      double v = img(xx,yy)- Fond;
	      if (  v < 1.e-10 )
		v = 0. ;
	      S += Sig2 + v ;
	    }

	}

    }  
  if ( S > 1.e-10 )
    S = sqrt(S);
  else
    S = 0. ;

  return(S);
}

#ifdef STORAGE
double
FondLocal(Image const & image ,float xcentre , float ycentre,
	  float rayon1 , float rayon2 , float fraclow,
	  float frachigh , int lp)
{  

  int dtaille = (int) ( rayon2 + 1.) ;
  int Ntot = (2*dtaille + 1) * (2*dtaille + 1);
  double *tab = new double[Ntot] ;
  int l = 0 ;
  for (int i=-dtaille; i<=dtaille; i++)
    {
      
      for (int j=-dtaille; j<=dtaille; j++)
	{
      
	  int x = (int) ( xcentre + i ); 
	  int y = (int) (ycentre + j );
	  double d = sqrt(double(i*i + j*j));
	  if ( (x <  image.Nx() ) && (x  >= 0 )&&
	       ( y < image.Ny() )&& (y >= 0)  &&
	       (d > rayon1) && (d < rayon2) )
	    {
	      tab[l]= image(x,y) ;
	      l++ ;
	    }
	}
    }
  SimpleSort T(l) ;
  for (int k=0; k<l; k++)
    T.Set(k , tab[k]) ;
  T.Sort();
  
  if ( lp > 9 )
    {
      cout << " l = " << l << endl ;
      for (int k=0; k<l; k++)
        cout << k << " " << T.Value(k) << endl ;
    }
  
  int badlow = (int) ( fraclow * l ) ;
  int badhigh = (int) ( frachigh * l ) ;
  badhigh = l - badhigh ;
  
  if ( lp > 9 )
    {
      cout << " badlow  = " << badlow << endl ;
      cout << " badhigh = " << badhigh << endl;    }
  
  
  double F = 0 ; 
  int count = 0 ;
  
  for (int k=badlow; k<badhigh; k++)
    {
      F += T.Value(k) ;
      count++ ;
    }
  
  if ( count > 0 )
    {
      if ( lp > 9 )
        {
          cout << "count = " << count << " , F = " 
               << F << " , F/count = " << F / count << endl ;
        }
      F = F / count ;
      if (lp > 4 )
        cout << "Fond = " << F << endl ;
    }
  else
    {
      if ( lp > 9 )
        {
          cout << "count = " << count << " , F = " 
               << F << endl ;
        }
      if ( lp > 4 )
        {
          cout << "Pb dans FondLocaL pour l'etoile ( " << xcentre 
               << " , " << ycentre << " ) : count = 0 " << endl ;
        }
    }
  
  delete [] tab;
  return(F) ;
  
}
#endif

//Calcul le flux pese par une PSF de sig seeing d'un objet.

static double 
Weighted_Flux_Aperture(const Image & img, const double & xbar, 
		       const double & ybar, const double & radius, 
		       int &Npix, const double & seeing, 
		       double const & sigweigth, 
		       const double & Fond)
{
  int tx = img.Nx();
  int ty = img.Ny();
  double S = 0 ;
  int dtaille =(int) ((radius + 0.5));
  int ix = (int) (xbar+0.5);
  int iy = (int) (ybar+0.5);
  Npix = 0; 
  double SumDiv=0, Wi=0, PSFi=0, SumPSFi=0;
  for(int i = -dtaille; i <=dtaille; i++)
    {
      for(int j = -dtaille; j <=dtaille; j++)
	{
	  int xx = ix + i ;
	  int yy = iy + j ;
	  double d = sqrt( double(i*i + j*j) );
	  if ( ( xx > 0 ) && ( xx < tx ) &&
	       ( yy > 0 ) && ( yy < ty )&&
	       (d < radius) )
	    {
	      double v = img(xx,yy)-Fond;
	      double r2 = sqr(double(xx) - xbar) + sqr(double(yy)-ybar);
	      Wi = exp(-(r2)/(2*sqr(sigweigth)));
	      PSFi = exp(-(r2)/(2*sqr(seeing)));
	      S += v * Wi ;

	      SumPSFi += PSFi;
	      
	      SumDiv += PSFi * Wi;
              ++Npix;
	    }

	}

    }
  
  double flux=S*SumPSFi/SumDiv;

  return(flux);
}


static double 
Chi2(const Image & img, const double & xbar, 
     const double & ybar, const double & radius, 
     int &Npix, const double & seeing, const double SkySigma,
     const double & flux, const double & Fond)
{
  int tx = img.Nx();
  int ty = img.Ny();
  double Chi2 = 0 ;
  int dtaille =(int) ((radius + 0.5));
  int ix = (int) (xbar+0.5);
  int iy = (int) (ybar+0.5);
  Npix = 0; 
  double   PSFi=0;
  double norm = 1/(2 * M_PI * sqr(seeing));
  for(int i = -dtaille; i <=dtaille; i++)
    {
      for(int j = -dtaille; j <=dtaille; j++)
	{
	  int xx = ix + i ;
	  int yy = iy + j ;
	  double d = sqrt( double(i*i + j*j) );
	  if ( ( xx > 0 ) && ( xx < tx ) &&
	       ( yy > 0 ) && ( yy < ty )&&
	       (d < radius) )
	    {
	      double r2 = sqr(double(xx) - xbar) + sqr(double(yy)-ybar);
	      double v = img(xx,yy)-Fond;

	      PSFi = norm * exp(-(r2)/(2*sqr(seeing)));
	      
	      Chi2 +=  sqr((PSFi*flux- v)/(SkySigma));
	      
	    }

	}

    }
  return(Chi2);

  //return(Chi2/SumPsf);
}




// calcule le barycentre pese par une PSF de sig seeing dans une ouverture 

static void 
Weighted_Barycentre(Image const & img, double rad_flux, 
		    double & xbar, double & ybar,double & flux, 
		    const double & seeing, const double &fond)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = (int) (floor(rad_flux));
  double sx=0;
  double sy=0;
  double s=0;
  double SumWi = 0, Wi = 0;
  double rad_flux2 = rad_flux * rad_flux ;
  int ix = (int) (floor(xbar)) ;
  int iy = (int) (floor(ybar)) ;
  for (int l=-d; l < d+1; l++)
    {
      for (int m=-d; m < d +1;  m++)
	{
	  int xu = ix + l ; 
	  int yu = iy + m ;
	  double dist2 = l*l+m*m;
	  if ( ( xu < tx ) &&  (xu >= 0 ) &&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist2 <rad_flux2) )
	    {
	      double r2=sqr(double(xu) - xbar) + sqr(double(yu)-ybar);
	      double val = img(xu,yu)- fond ;
	      Wi = exp(-(r2)/(2*seeing*seeing));
	      s += val*Wi ;
	      sx += l*val*Wi;
	      sy += m*val*Wi;
	      SumWi += Wi;
	    }
	}
    }
  flux = s/SumWi ;
  if (s>0) 
    {
      sx  = sx/s;
      sy  = sy/s;
      xbar = sx + ix;
      ybar = sy + iy;
     }
  return; 
}

//#ifdef STORAGE   
static void
Var_Position(Image const & img , double const &rad_flux, 
	     double & xbar, double & ybar, const double & fond, 
	     double & varx, double & vary, double & varxy, double & seeing)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = (int) ((rad_flux + 1.));
  double S=0;
  double Sigma = 0;
  double rad_flux2 = rad_flux * rad_flux ;
  int ix = (int) (xbar + 0.5 ) ;
  int iy = (int) (ybar + 0.5 ) ;
  for (int l=-d; l < d + 1; l++)
    {
      for (int m=-d; m < d +1;  m++)
	{
	  int xu = ix + l ; 
	  int yu = iy + m ;
	  double dist2 = l*l+m*m;
	  if ( ( xu < tx ) &&  (xu >= 0 )&&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist2 <rad_flux2) )
	    {
	      double val = img(xu,yu);
	      double r2 = sqr(double(xu) - xbar) + sqr(double(yu)-ybar);
	      double wi = exp(-(r2)/(2*sqr(seeing)));
	      Sigma = fond + val;
	      double wi2 = sqr(wi);
	      varx += sqr(xu - xbar) * Sigma * wi2;
	      vary += sqr(yu - ybar) * Sigma * wi2;
	      varxy+= (xu - xbar) * (yu - ybar) * Sigma * wi2;
	      S += val * wi;
	    }
	}
    }
  S=sqr(S);
  varx = varx / S;
  vary = vary / S;
  varxy= varxy/ S;

  return;	  

}

//#endif
void
Flux_Coord(Image const & img , double x, double y , double rad_flux,
	   float fond, double & xbar, double & ybar, double & flux )
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = int( rad_flux + 1);
  double sx=0;
  double sy=0;
  double s=0;
  for (int l=-d; l <= d; l++)
    {
      for (int m=-d; m <=d ;  m++)
	{
	  int xu , yu ;
	  xu = (int) ( x + l ); // sert pas a grd chose mais plus clair
	  yu = (int) ( y + m );
	  double dist = sqrt(double(l*l+m*m));
	  if ( ( xu < tx ) &&  (xu >= 0 )&&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist <rad_flux) )
	    {
	      double val =(float) (img(xu,yu) - fond );
	      s += val ;
	      sx += xu*val;
	      sy += yu*val;
	    }
	}
    }
  flux = s ;
  if (s>0) 
    {
      xbar   = sx/s;
      ybar   = sy/s;
    }


  return;
}

//#ifdef STORAGE
static void
Sig_Weighted_Position(Image const & img ,double const &rad_flux, 
		      double & xbar, double & ybar, const double & seeing, 
		      const double & fond, 
		      double & Sigx, double & Sigy, double & Sigxy)
{
  int tx = img.Nx() ;
  int ty = img.Ny() ;
  int d = (int) ((rad_flux + 1.));
  double sx=0;
  double sy=0;
  double S=0;
  double SumWi = 0, Wi = 0, SumSigx = 0, SumSigy = 0, SumSigxy = 0, 
    Sx = 0, Sy = 0, Sxy;
  double rad_flux2 = rad_flux * rad_flux ;
  int ix = (int) (xbar + 0.5 ) ;
  int iy = (int) (ybar + 0.5 ) ;
  for (int l=-d; l <= d; l++)
    {
      for (int m=-d; m <=d ;  m++)
	{
	  int xu = ix + l ; 
	  int yu = iy + m ;
	  double dist2 = l*l+m*m;
	  if ( ( xu < tx ) &&  (xu >= 0 )&&
	       ( yu < ty ) &&  (yu >= 0 ) && (dist2 <rad_flux2) )
	    {
	      double r2=sqr(double(xu) - xbar) + sqr(double(yu)-ybar);
	      double val = img(xu,yu)- fond ;
	      Wi = exp(-(r2)/(2*seeing*seeing));
	      sx += l*val*Wi;
	      sy += m*val*Wi;
	      SumSigx += sqr(l) * val * Wi;
	      SumSigy += sqr(m) * val * Wi;
	      SumSigxy += l * m * val* Wi;
	      SumWi += Wi; 
	      S += Wi * val;
	    }
	}
    }
      Sx = SumSigx / S - sqr(sx/S);
      Sy = SumSigy / S - sqr(sy/S);
      Sxy = SumSigxy / S - sx*sy/sqr(S);
  if (Sy>0 && Sy>0)
    {
      Sigx = Sx;
      Sigy = Sy;
      Sigxy = Sxy;
    }
  return;	  

}
//#endif


static void
Seuils_CvDetection(Image const & img, Image & mask, Image & imgcv, 
		   double seeing, double & SigWeigth, double &mean_amp, 
		   double & sig_amp, double & mean_fluxcv, 
		   double  &sig_fluxcv, double & Norma)
{
  Image *imgcv1 = new Image(img.Nx(),img.Ny()); 
  double precision = 1.e-3 ;
  int l = Demi_Fenetre(-1 ,SigWeigth, precision);
  Convole1DX(img, SigWeigth, l , *imgcv1 ) ;
  Convole1DY(*imgcv1, SigWeigth, l , imgcv ) ;
  delete imgcv1 ;

  //Image noyau = Convole_MakeNoyau(seeing, l) ;
  Image poids = Convole_MakeNoyau(SigWeigth, l);

  double S2 = poids.SumSquaredPixels()/sqr(poids.SumPixels()) ;
  //double S2 = ComputeS2(noyau, poids); 
  S2 = 1. / S2 ;
  cout << "S2 " << S2 << " " << 4.*M_PI*seeing*seeing << endl ;

  imgcv *= S2 ;
  
  Frame frame(2*l,2*l,img.Nx()-2*l,img.Ny()-2*l );
  mask.Masking(frame);
  double mean = 0., sigma = 0., meancv=0., sigmacv = 0. ;

  img.ClippedMeanSigmaValue(mean,sigma,&mask);
  imgcv.ClippedMeanSigmaValue(meancv,sigmacv,&mask);
  
  cout << " Mean, sigma image:  " <<  mean  << " " << sigma 
       << " Mean, sigma imagecv: " <<  meancv << " " << sigmacv << endl;

  sig_amp = sigma ;
  sig_fluxcv =  sigmacv ;
  mean_amp = mean ;
  mean_fluxcv =  meancv ;
  Norma = S2 ;
}
// enleve les multiples detections

static bool IgnoreIt(const BaseStar *star )
{
  return (star->flux == 0);
}

#include "fastfinder.h"
/* The new cleaner routine uses the fastfinder.
It looks for the nearest star from the current star, if the neighbour is in 
the cut radius and have a flux greater than the currrent star, then the current
star is set to zero and then deleted from the list*/
static void
CleanCandidateList(CandidateStarList & stl, double delta)
{
  //cerr << "ATTENTION L'ANCIEN DEDOUBLONNEUR EST ENCORE EN MARCHE " << endl;
  cout << "Size of the list before cleaning: " << stl.size() << endl ;
  cout << "Cut radius for double detectiom " << delta << endl;
  // sort and copy of the list for old cleaner
  stl.FluxSort() ;
  CandidateStarList Newstl;
  stl.CopyTo(Newstl);
  
  FastFinder finder(*((BaseStarList*) &stl));
  
  for (CandidateStarIterator it = stl.begin(); it != stl.end(); it++)
    {
      CandidateStar * star = *it ;
      // To avoid having neighbour = current star using FindClosest
      // The flux is reaffected to its initial value after.
      double KeepFlux = star->flux;
      star->flux = 0;
      bool keep_it = true;
      const BaseStar *neighbour = finder.FindClosest(*star, delta, &IgnoreIt); 
      if ( neighbour != NULL)
	{
	  if ( star->flux < neighbour->flux)
	    keep_it = false;
	}
      if (keep_it) star->flux = KeepFlux;
    }
  
  // to erase the star with 0 flux.
  for (CandidateStarIterator it = stl.begin(); it != stl.end(); it++)
    {
      CandidateStar * star = *it ; 
      if ( star->flux == 0) // on erase
	it = stl.erase(it);
    }
  
  
  
  
  // ANCIEN DEDOUBLONEUR EN PARRALLELE POUR CHECK
  
#ifdef STORAGE
  for (it = Newstl.begin(); it != Newstl.end(); it++)
    {
      CandidateStar * star = *it ; 
      CandidateStarIterator itplus1 = it ;
      itplus1++;
      for (CandidateStarIterator iti = itplus1  ; iti != Newstl.end(); )
	{
	  CandidateStar * stari = *iti ;
	  double d =  star->Distance(*stari);
	
	  if ( d < delta ) // on erase
	    iti = Newstl.erase(iti);
	  else
	    iti++;
	}
      
    }



  if ( stl.size() != Newstl.size() )
    {
      cout << "PROBLEME AVEC LE DEDOUBLONEUR" << endl;
      cout << "Taille liste sortie fast: " << stl.size() << endl ;
      cout << "Taille liste sortie slow: " << Newstl.size() << endl ;
      stl.write("dedoub_fast.list");
      Newstl.write("dedoub_slow.list");
    }
#endif
  cout << "Size of the list after cleaning: " << stl.size() << endl ;
  return ;
  
}



int
NewCvDetection( Image & img , Image & mask, CandidateStarList & stl,  
		DatDetec & datdet, double & seeing, double & sigweigth)
{
  cerr << "Running CvDetection " << endl ; 
  cout << "MEASURED seeing : " << seeing << endl ;
  // On va chercher les coupures dans la datacard
  //seeing=1.6;  
  datdet.ComputeRadius(seeing);
  //  datdet.Print();
  
  double rad_flux = datdet.rayon;
  double rad_bad = datdet.rad_bad;
  int rad_locmax = datdet.rad_locmax;
  float nsig_flux = datdet.n_seuil_flux;
  
  cout << "Valeurs des coupures " << endl;
  cout << "rad_flux : " << rad_flux << endl;
  cout << "rad_bad : " << rad_bad << endl;
  cout << "rad_locmax : " << rad_locmax << endl;
  cout << "nsig_flux : " <<nsig_flux << endl; 
 
   
  const Image *pimage = &img ;

  // On retire un fond global et on ne s'en occupe plus
  bool GlobBack = true;
  if (GlobBack) { img.Surface(40,img);}
  
  //  FitsImage("back.fits",img);
  //declaration de l'image convoluee
  Image imgcv(img.Nx(),img.Ny());
  
  // Declaration et calcul des coupures pour la detection
  double mean=0., meancv=0., sigma = 0. , sigmacv = 0., S2=0. ;
  
  Seuils_CvDetection(*pimage,mask, imgcv, seeing, sigweigth, mean, 
		     sigma , meancv, sigmacv, S2 );
  
  double fluxmin = nsig_flux * sigmacv + meancv;
  cout << "Seuils flux image vonvoluee :  " << " " << fluxmin << endl ;

  // Declaration des limites pour les calculs de flux et iterateur pour la boucle de detection
  int irad_locmax = rad_locmax;
  double deuxrad_flux = 1.6*rad_flux;
  int tx = pimage->Nx();
  int ty = pimage->Ny();
  int N = 0 ;
  Pixel * p = pimage->begin() ;
  Pixel * pcv = imgcv.begin() ;
  Pixel * pm = mask.begin() ;
  p += irad_locmax + tx * irad_locmax ;
  pcv += irad_locmax + tx * irad_locmax ;
  pm += irad_locmax + tx * irad_locmax ;
  int delta = 2 * irad_locmax ;

  for( int y = irad_locmax ; y < ty - irad_locmax ; y++)
    {
      for( int x = irad_locmax ; x < tx - irad_locmax ; x++)
	{
	  // test sur valeur du pixel central, 
	  // puis si pas de pixel bad pres (a rad_bad), 
	  // si c'est un max. loc. ds un rayon rad_locmaxsur l'image convolue, 
	  //	  Pixel valpix = *p;
	  Pixel valpixcv = *pcv  ;
	  double valpix2=0., valpixcv2 = 0. ;
	  
	  // MAX LOCAL sur l'image convolue
	  if ( ( valpixcv  > fluxmin )  && 
	       IsLocalMax(imgcv,x,y,rad_locmax) &&
	       IsNotBadPix(img, mask, x,y,rad_bad) )
	    {
	      // Xbar et ybar coord du barycentre de l'objet
	      double xbar=x, ybar=y, flux=0.,  fd_loc = 0 ;
	      
	      // Calcul de la position du barycentre
	      double xloop=0, yloop=0, deltax = 1, deltay = 1;
	      int compteur=0;
	      
	      while (deltax>0.005 || deltay>0.005)
		{
		  flux=0;
		  Weighted_Barycentre(*pimage, deuxrad_flux, 
				      xbar, ybar, flux, seeing, fd_loc);
		  deltax = fabs(xloop-xbar);
		  deltay = fabs(yloop-ybar);
		  xloop=xbar;
		  yloop=ybar;
		  compteur++;
		  if (compteur>10) break;
		}
	      CandidateStar * star = new CandidateStar;
	      
	      double Sigx = 0, Sigy = 0, Sigxy = 0;
	      
	      Sig_Weighted_Position(*pimage, deuxrad_flux, xbar, ybar, 
				    seeing, 0, Sigx, Sigy, Sigxy);
	      
	      int  area;
	      // Calcul de la variance sur la position du barycentre
	      
	      double varx = 0, vary = 0, varxy = 0;
	      
	      Var_Position(*pimage, deuxrad_flux, 
			   xbar, ybar, sigma*sigma, 
			   varx, vary, varxy, seeing);
	      

#ifdef STORAGE	      	      
	      // calcul fond local un peu ole
	      double rad1 = 6.*rad_flux/2.5 ; 
	      double rad2 = 8.*rad_flux/2.5 ; 
	      double fond = FondLocal(img, xbar, ybar, rad1, rad2, 0.15, 0.15); 
#endif     	      

	      // Calcul du flux pese
	      //Attention tres sensible au seeing a flux plus grand
	      flux = Weighted_Flux_Aperture(img, xbar, ybar,deuxrad_flux,
					    area, seeing, sigweigth, 
					    /*fond*/ 0);
	      double chi2 = Chi2(img, xbar, ybar, deuxrad_flux, area, 
				 seeing, sigma, flux, /*fond*/0);
	      
	      if (xbar >0 && xbar < tx && ybar > 0  && ybar < ty )
		{ 
		  valpixcv2 = imgcv((int)(xbar+0.5), (int)(ybar+0.5))  ;
		  valpix2 = (*pimage)((int)(xbar+0.5), (int)(ybar+0.5))  ; 
		}
	      
	      //Position du barycentre
	      star->x = xbar + DECALAGE_IJ_SE ;
	      star->y = ybar + DECALAGE_IJ_SE ;
	      
	      //Position du pixel allume
	      star->X_Peak() = x;
	      star->Y_Peak() = y;
	      int npix=0;
	      
	      double aper_flux = Flux_Aperture(img, x, y, rad_flux,0,npix);
	      
	      
	      if (flux>10*aper_flux || flux<0.1*aper_flux) flux=aper_flux;
	      
	      
	      
	      star->Mxx() = Sigx ;
	      star->Myy() = Sigy ;
	      star->Mxy() = Sigxy;
	      

	      star->flux = flux;
	      star->Flux_fixaper() = aper_flux ;
              star->Isoarea() = area;
	      
	      star->FluxCv() = (valpixcv2-meancv) ;
	      star->Noise() = sigmacv;
	      star->EFlux() = sigmacv;
	      star->Fond() = sigmacv ;

	      //Variance sur la poistion du barycentre

	      star->Fwhm() = sqrt(Sigx*Sigx + Sigy*Sigy);
	      star->Sigx() = varx ;
	      star->Sigy() = vary ;
	      star->Sigxy() = varxy;
	      star->Chi2() = chi2;

	      if ( star->Noise() > 1.e-10 )
		star->SigToNoise() = star->flux / star->Noise() ;
	      
	      star->Fluxmax() = valpix2 /*0- ffd */; 
	      //star->Fond() = ffd ;
	      
	      N++ ;
	      
	      stl.push_back(star);
	    }
	  p++ ;
	  pcv++ ;
	  pm++ ;
	}
      p += delta ;
      pcv += delta ;
      pm += delta ;
    }
  
  // on cleane des doubles detections
  double distmin = 0.5*rad_flux;
  CleanCandidateList(stl, distmin);
  
  cerr << N << " " << stl.size()
       << " etoiles detectees par CvDetection " << endl ;
  
  return(N) ;
}


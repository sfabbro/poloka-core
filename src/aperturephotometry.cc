#include <iostream>

#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif

//#define DEBUG
//#define TEST

#include "aperturephotometry.h"
#include "starlist.h"
#include "starlist.cc" /* since starlist is a template class */

// premier probleme, quelle est la fraction de pixel couverte par un disque de rayon R ?
// x_pix, y_pix coord du centre du pixel (de taille 1x1 par convention)

static double fraction(double x_centre, double y_centre, double R, double x_pix, double y_pix)
{
  double x_min, x_max;
	// Je m'arrange pour que le pixel soit dans le premier quadrant par rapport au centre du cercle
	
	if (x_pix < x_centre) {x_pix = 2.*x_centre - x_pix;};
	if (y_pix < y_centre) {y_pix = 2.*y_centre - y_pix;};
       
	// Je m'arrange même pour qu'il soit dans la zone y>x
	
	if ((x_pix-x_centre)>(y_pix-y_centre))
	{
		double prov = x_pix; x_pix = y_pix; y_pix = prov;
		prov = x_centre; x_centre = y_centre; y_centre = prov;
	}
	
#ifdef DEBUG
	cout << x_centre << " " << y_centre << " " << R << " " << x_pix << " " << y_pix << endl;
#endif

	// je definis les coordonnees par rapport au coin inf gauche du pixel

	x_centre -= x_pix - 0.5;
	y_centre -= y_pix - 0.5;

#ifdef DEBUG
	cout << "coord centre utilisees : " << x_centre << " " << y_centre << endl;
#endif
	// si la distance entre centre et bord inf gauch du pix est > R, le pixel est dehors (faux pour R petit)

	if (pow(x_centre,2.)+pow(y_centre,2.) > pow(R,2.)) return 0.;

	// si la distance entre centre et bord sup droit du pix est < R, le pixel est dedans (faux pour R petit)

	if (pow(x_centre-1.,2.)+pow(y_centre-1.,2.) < pow(R,2.)) return 1.;

	// L'equation du cercle est y = y_centre + sqrt(R^2-(x-x_centre)^2)
	
	// 
	// On calcule les intersections de ce cercle avec les bords du pixel.
	//
	
		
	// en quels points du bord inferieur du pixel passe le cercle.

	x_min = x_centre - sqrt(pow(R,2.) - pow(y_centre,2.));
	x_max = x_centre + sqrt(pow(R,2.) - pow(y_centre,2.));

	if (x_min > x_max) { double prov = x_min; x_min=x_max; x_max=prov;}
	
#ifdef DEBUG
	cout << "intersections bord inf : " << x_min << " " << x_max << endl;
#endif

	if (x_min < 0.) x_min = 0.;
	if (x_max > 1.) x_max = 1.;


	// Si le point UL du pixel n'est pas dans le cercle, il suffit de calculer l'integrale

	double complement;
	if (pow(x_centre,2.)+pow(y_centre+1.,2.) < pow(R,2.))
	  {
#ifdef DEBUG
	    cout << "cas 1 : complement = 0" << endl;
#endif
	    complement = 0.;
	  }
	else
	// Sinon il faut recalculer x_min, le point ou le cercle coupe le bord sup du pixel
	  {
	    double x_min1 = x_centre - sqrt(pow(R,2.) - pow(y_centre-1.,2.));
	    double x_min2 = x_centre + sqrt(pow(R,2.) - pow(y_centre-1.,2.));

#ifdef DEBUG
	    cout << "intersections bord sup : " << x_min1 << " " << x_min2 << endl;
#endif

	    x_min = x_min1;
	    if (x_min<x_min2) x_min=x_min2;
	    complement = x_min;
#ifdef DEBUG
	    cout << "cas 2 :complement = " << x_min << endl;
#endif
	  }

	//int(sqrt(R^2-x^2),x) = 1/2 x sqrt(R^2-x^2) + 1/2 R^2 arcsin (x/R)

	double integrale = (x_max-x_centre) * sqrt(pow(R,2.) - pow(x_max-x_centre,2.))/2. 
	  + pow(R,2.)* asin ((x_max-x_centre)/R)/2. + y_centre*(x_max-x_centre);
	integrale -= (x_min-x_centre) * sqrt(pow(R,2.) - pow(x_min-x_centre,2.))/2. 
	  + pow(R,2.)* asin ((x_min-x_centre)/R)/2. + y_centre*(x_min-x_centre);
	integrale += complement;
	  
#ifdef DEBUG
	cout << integrale << endl;
#endif
	
	return x_min + integrale;	
}



// ==============================================
//
// FillFlux : calcule les flux d'ouverture et le barycentre de l'objet (mais ne remet pas a jour les valeurs)
//
// ==============================================

bool AperturePhotomSEStar::FillFlux(FitsImage& fitsimage, FitsImage& fitsweight, int star_id=0, int image_id=0)
{
  bool valid=true;
  for (int i=0; i<nb_rayons; ++i) aperture_flux[i] = 0.;

  double rayon_max = rayon[nb_rayons-1];
  
  int xmin = -1 - (int) rayon_max;
  int xmax =  1 + (int) rayon_max;
  int ymin = -1 - (int) rayon_max;
  int ymax = 1 + (int) rayon_max;

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin<=1.)
      ||( (x+delta_x+xmax+1.) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin <= 1.)
      ||(y+delta_y+ymax+1. >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  float provmax = 0.;
  float provmax2 = 0.;

  //
  // Je remplis les flux et je repere les 2 pixels de flux maximal.
  //

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x), (int) (yy+y+delta_y));
	if  (fitsweight((int) (xx+x+delta_x), (int) (yy+y+delta_y))==0.)  valid=false;
	if ((xx*xx+yy*yy<20.)&& (val_pixel > provmax)) 
	  {
	    provmax2 = provmax;
	    provmax = val_pixel;	    
	  }
	for (int i=0; i<nb_rayons; ++i) aperture_flux[i] += fraction(0., 0., rayon[i], xx, yy)*val_pixel;
	
      }
  
  fluxmax  = provmax;
  fluxmax2 = provmax2;
  
  return valid;
}


bool AperturePhotomSEStar::FillFlux_noweight(FitsImage& fitsimage, double mean, int star_id, int image_id)
{
  bool valid=true;
  for (int i=0; i<nb_rayons; ++i) aperture_flux[i] = 0.;

  int dec_x = 32;
  int dec_y = 1;

  double rayon_max = rayon[nb_rayons-1];
  
  int xmin = -1 - (int) rayon_max;
  int xmax =  1 + (int) rayon_max;
  int ymin = -1 - (int) rayon_max;
  int ymax = 1 + (int) rayon_max;

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin+dec_x<=1.)
      ||( (x+delta_x+xmax+1.+dec_x) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin+dec_y <= 1.)
      ||(y+delta_y+ymax+1.+dec_y >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  float provmax = 0.;
  float provmax2 = 0.;

  //
  // Je remplis les flux et je repere les 2 pixels de flux maximal.
  //

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x+dec_x), (int) (yy+y+delta_y+dec_y)) - mean;

	if ((xx*xx+yy*yy<20.)&& (val_pixel > provmax)) 
	  {
	    provmax2 = provmax;
	    provmax = val_pixel;	    
	  }
	
	for (int i=0; i<nb_rayons; ++i) 
	  {
	    aperture_flux[i] += fraction(0., 0., rayon[i], xx, yy)*val_pixel;
	  }
      }
  
  fluxmax  = provmax;
  fluxmax2 = provmax2;
  
  return valid;
}

//
// calcule les flux d'ouverture et recalcule le barycentre de l'objet.
//

bool AperturePhotomSEStar::ComputeBarycenter(FitsImage& fitsimage, FitsImage& fitsweight)
{
  bool valid=true;
  float total = 0.;
  float x_moyen = 0., y_moyen=0.;

  int dec_x = 32;
  int dec_y = 1;

  int xmin = -1 - (int) rayon[5];
  int xmax =  1 + (int) rayon[5];
  int ymin = -1 - (int) rayon[5];
  int ymax =  1 + (int) rayon[5];

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin+dec_x<=1.)
      ||( (x+delta_x+xmax+1.+dec_x) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin +dec_y<= 1.)
      ||(y+delta_y+ymax+1.+dec_y >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x+dec_x), (int) (yy+y+delta_y+dec_y));
	if  (fitsweight((int) (xx+x+delta_x+dec_x), (int) (yy+y+delta_y+dec_y))==0.)  valid=false;
	float frac = fraction(0., 0., rayon[5], xx, yy);
	total += frac*val_pixel;
	x_moyen += frac*val_pixel*xx;
	y_moyen += frac*val_pixel*yy;
      }

  delta_x = x_moyen/total;
  delta_y = y_moyen/total;

  if (delta_x*delta_x+delta_y*delta_y<25.) return valid;
  else return false;

}
bool AperturePhotomSEStar::ComputeBarycenter_noweight(FitsImage& fitsimage, int dec_x = 0, int dec_y=0)
{
  bool valid=true;
  float total = 0.;
  float x_moyen = 0., y_moyen=0.;

  int xmin = -1 - (int) rayon[5];
  int xmax =  1 + (int) rayon[5];
  int ymin = -1 - (int) rayon[5];
  int ymax =  1 + (int) rayon[5];

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin+dec_x<=1.)
      ||( (x+delta_x+xmax+1.+dec_x) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin +dec_y<= 1.)
      ||(y+delta_y+ymax+1.+dec_y >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x+dec_x), (int) (yy+y+delta_y+dec_y));

	float frac = fraction(0., 0., rayon[5], xx, yy);
	total += frac*val_pixel;
	x_moyen += frac*val_pixel*xx;
	y_moyen += frac*val_pixel*yy;
      }

  delta_x = x_moyen/total;
  delta_y = y_moyen/total;

  if (delta_x*delta_x+delta_y*delta_y<25.) return valid;
  else return false;

}

//
// calcule la racine de l'ecart-type des mesures entre rayon[8] et 
// rayon[nb_rayons], normalisee a la valeur moyenne.
//

float AperturePhotomSEStar::compute_platitude()
{
  float mean = 0., var = 0.;
  
  int i_debut = 8;
  for (int i_rayon=i_debut; i_rayon<nb_rayons; ++i_rayon)
    {
      mean += aperture_flux[i_rayon];
      var += pow(aperture_flux[i_rayon],(float) 2.);
    }
  mean /= (nb_rayons-i_debut);
  var /= (nb_rayons-i_debut);
  if (mean>0.)  return sqrt(var-mean*mean)/mean;
  else return -1;
}


void AperturePhotomSEStar::write(ostream& s, const bool write_header)
{
  if (write_header)
    {
      WriteHeader_(s);
      s << "# ra : "   << endl
	<< "# dec : "  << endl
	<< "# r3 : "   << endl
	<< "# r6 : "   << endl
	<< "# r10 : "  << endl
	<< "# r19 : "  << endl
	<< "# f3 : "   << endl
	<< "# f6 : "   << endl
	<< "# f10 : "  << endl
	<< "# f19 : "  << endl
	<< "# plat : " << endl
	<< "# end" << endl;
    }
    else
      {
	writen(s);
	s << ra << " " 
	  << dec << " "
	  << rayon[3] << " " 
	  << rayon[6] << " " 
	  << rayon[10] << " " 
	  << rayon[19] << " " 
	  << aperture_flux[3] << " " 
	  << aperture_flux[6] << " " 
	  << aperture_flux[10] << " " 
	  << aperture_flux[19] << " " 
	  << compute_platitude() << " "
	  << endl;
      }
  return;
}

void AperturePhotomSEStar::write_short(ostream& s, const bool write_header)
{
  if (write_header)
      s << "# x : "     << endl 
	<< "# y : "     << endl
	<< "# ra : "    << endl
	<< "# dec : "   << endl
	<< "# flux : "  << endl
	<< "# eFlux : " << endl
	<< "# fmax : "  << endl
	<< "# fmax2 : " << endl
	<< "# fwhm : "  << endl
	<< "# r3 : "    << endl
	<< "# r6 : "    << endl
	<< "# r10 : "   << endl
	<< "# r19 : "   << endl
	<< "# f3 : "    << endl
	<< "# f6 : "    << endl
	<< "# f10 : "   << endl
	<< "# f19 : "   << endl
	<< "# plat : "  << endl
	<< "# end"      << endl;
  else
    s << x << " "
      << y << " "
      << ra << " "
      << dec << " "
      << flux << " "
      << EFlux() << " "
      << fluxmax << " "
      << fluxmax2 << " "
      << Fwhm() << " "
      << rayon[3] << " " 
      << rayon[6] << " " 
      << rayon[10] << " " 
      << rayon[19] << " " 
      << aperture_flux[3] << " " 
      << aperture_flux[6] << " " 
      << aperture_flux[10] << " " 
      << aperture_flux[19] << " " 
      << compute_platitude() << " " 
      << endl;
  return;
}

void AperturePhotomSEStar::write_profile(ostream& s, const bool write_header)
{
  if (write_header)
    s << "# x :"     << endl 
      << "# y :"     << endl 
      << "# r :"     << endl 
      << "# f :"     << endl 
      << "# fmax :"  << endl 
      << "# fmax2 :" << endl 
      << "# fwhm :"  << endl 
      << "# f3 :"    << endl 
      << "# f6 :"    << endl 
      << "# f10 :"   << endl
      << "# f19 :"   << endl 
      << "# plat :"  << endl 
      << "# end"     << endl;
  else
    for (int i=0; i<nb_rayons; ++i)
      s << x << " "
	<< y << " "
	<< rayon[i] << " "
	<< aperture_flux[i] << " "
	<< fluxmax << " "
	<< fluxmax2 << " "
	<< Fwhm() << " "
	<< aperture_flux[3] << " " 
	<< aperture_flux[6] << " " 
	<< aperture_flux[10] << " " 
	<< aperture_flux[19] << " " 
	<< compute_platitude() << " " 
	<< endl;
  return;
}



// ==============================================
//
// FillFlux : calcule les flux d'ouverture et le barycentre de l'objet (mais ne remet pas a jour les valeurs)
//
// ==============================================

bool AperturePhotomBaseStar::FillFlux(FitsImage& fitsimage, FitsImage& fitsweight, int star_id=0, int image_id=0)
{
  bool valid=true;
  float total[nb_rayons];
  for (int i=0; i<nb_rayons; ++i) total[i] = 0.;

  double x_moyen = 0., y_moyen=0.;

  double rayon_max = rayon[nb_rayons-1];
  
  int xmin = -1 - (int) rayon_max;
  int xmax =  1 + (int) rayon_max;
  int ymin = -1 - (int) rayon_max;
  int ymax = 1 + (int) rayon_max;

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin<=1.)
      ||( (x+delta_x+xmax+1.) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin <= 1.)
      ||(y+delta_y+ymax+1. >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  float provmax = 0.;
  float provmax2 = 0.;

  //
  // Je remplis les flux et je repere les 2 pixels de flux maximal.
  //

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x), (int) (yy+y+delta_y));
	if  (fitsweight((int) (xx+x+delta_x), (int) (yy+y+delta_y))==0.)  valid=false;
	if ((xx*xx+yy*yy<20.)&& (val_pixel > provmax)) 
	  {
	    provmax2 = provmax;
	    provmax = val_pixel;	    
	  }
	for (int i=0; i<nb_rayons; ++i) 
	 {
	   total[i] += fraction(0., 0., rayon[i], xx, yy)*val_pixel;
	   cout << setprecision(1) << fraction(0., 0., rayon[i], xx, yy) << " ";
	 }
	 cout << setprecision(12) << endl;
	
	x_moyen += fraction(0., 0., rayon[14], xx, yy)*val_pixel*xx;
	y_moyen += fraction(0., 0., rayon[14], xx, yy)*val_pixel*yy;
      }
  
  fluxmax  = provmax;
  fluxmax2 = provmax2;
  
  for (int i=0; i<nb_rayons; ++i) aperture_flux[i] = total[i];
  
  //on ne remets pas les valeurs de delta_x et delta_y a jour
  //delta_x = x_moyen/total15;
  //delta_y = y_moyen/total15;
  
  return valid;
}


bool AperturePhotomBaseStar::FillFlux_noweight(FitsImage& fitsimage, double mean, int star_id, int image_id)
{
  bool valid=true;
  for (int i=0; i<nb_rayons; ++i) aperture_flux[i] = 0.;

  double rayon_max = rayon[nb_rayons-1];
  
  int xmin = -1 - (int) rayon_max;
  int xmax =  1 + (int) rayon_max;
  int ymin = -1 - (int) rayon_max;
  int ymax = 1 + (int) rayon_max;

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin<=1.)
      ||( (x+delta_x+xmax+1.) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin <= 1.)
      ||(y+delta_y+ymax+1. >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  float provmax = 0.;
  float provmax2 = 0.;

  //
  // Je remplis les flux et je repere les 2 pixels de flux maximal.
  //

  for (int yy = ymin; yy <= ymax; ++yy)
    {
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x), (int) (yy+y+delta_y));
	
	val_pixel -= mean;

	if ((xx*xx+yy*yy<20.)&& (val_pixel > provmax)) 
	  {
	    provmax2 = provmax;
	    provmax = val_pixel;	    
	  }
	for (int i=0; i<nb_rayons; ++i) 
	 {
	   aperture_flux[i] += fraction(0., 0., rayon[i], xx, yy)*val_pixel;
	 }
      }
    }
  
  fluxmax  = provmax;
  fluxmax2 = provmax2;
  
  //  for (int i=0; i<nb_rayons; ++i) aperture_flux[i] = total[i];
  
  return valid;
}

//
// calcule les flux d'ouverture et recalcule le barycentre de l'objet.
//

bool AperturePhotomBaseStar::ComputeBarycenter(FitsImage& fitsimage, FitsImage& fitsweight)
{
  bool valid=true;
  float total = 0.;
  float x_moyen = 0., y_moyen=0.;

  cout << "(" << x << "," << y << ") ";

  int xmin = -1 - (int) rayon[5];
  int xmax =  1 + (int) rayon[5];
  int ymin = -1 - (int) rayon[5];
  int ymax =  1 + (int) rayon[5];

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin<=1.)
      ||( (x+delta_x+xmax+1.) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin <= 1.)
      ||(y+delta_y+ymax+1. >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x), (int) (yy+y+delta_y));
	if  (fitsweight((int) (xx+x+delta_x), (int) (yy+y+delta_y))==0.)  valid=false;
	float frac = fraction(0., 0., rayon[5], xx, yy);
	total += frac*val_pixel;
	x_moyen += frac*val_pixel*xx;
	y_moyen += frac*val_pixel*yy;
      }

  delta_x = x_moyen/total;
  delta_y = y_moyen/total;

  if (delta_x*delta_x+delta_y*delta_y<25.) return valid;
  else return false;

}
bool AperturePhotomBaseStar::ComputeBarycenter_noweight(FitsImage& fitsimage)
{
  bool valid=true;
  float total = 0.;
  float x_moyen = 0., y_moyen=0.;

  int xmin = -1 - (int) rayon[5];
  int xmax =  1 + (int) rayon[5];
  int ymin = -1 - (int) rayon[5];
  int ymax =  1 + (int) rayon[5];

  //  cout << "test : (" << x << "," << delta_x << "," << xmin << "," << xmax << ") ";
  //cout << "test : (" << y << "," << delta_y << "," << ymin << "," << ymax << ") ";

  // Si on est trop près du bord on arrete
  if ((x+delta_x+xmin<=1.)
      ||( (x+delta_x+xmax+1.) >= (double) fitsimage.KeyVal("NAXIS1"))
      ||(y+delta_y+ymin <= 1.)
      ||(y+delta_y+ymax+1. >= (double) fitsimage.KeyVal("NAXIS2")))
    return false;

  for (int yy = ymin; yy <= ymax; ++yy)
    for (int xx = xmin; xx <= xmax ; ++xx)
      {
	float val_pixel =  fitsimage((int) (xx+x+delta_x), (int) (yy+y+delta_y));

	float frac = fraction(0., 0., rayon[5], xx, yy);
	total += frac*val_pixel;
	x_moyen += frac*val_pixel*xx;
	y_moyen += frac*val_pixel*yy;
      }

  delta_x = x_moyen/total;
  delta_y = y_moyen/total;

  if (delta_x*delta_x+delta_y*delta_y<25.) return valid;
  else return false;

}

// ==============================
//
// calcule la racine de l'ecart-type des mesures entre rayon[8] et 
// rayon[nb_rayons], normalisee a la valeur moyenne.
//
// ==============================

float AperturePhotomBaseStar::compute_platitude(double &mean, double &var, int i_debut)
{
  mean = 0.; var = 0.;
  
  for (int i_rayon=i_debut; i_rayon<nb_rayons; ++i_rayon)
    {
      mean += aperture_flux[i_rayon];
      var += pow(aperture_flux[i_rayon],(float) 2.);
    }
  mean /= (nb_rayons-i_debut);
  var /= (nb_rayons-i_debut);
  if (mean>0.)  return sqrt(var-mean*mean)/mean;
  else return -1;
}

float AperturePhotomBaseStar::compute_platitude(int i_debut)
{
  double mean = 0., var = 0.;
  return compute_platitude(mean, var, i_debut) ;
}


float AperturePhotomBaseStar::compute_indep_platitude(double &mean, double &var, int i_debut)
{
  mean = 0.; var = 0.;
  double flux_diff[nb_rayons];

  if (i_debut < 2) return -1.;

  //
  // flux_diff is the mean flux per square pixel in a ring between two consectutive radii
  //

  int nb=0;

  //cout << "\n" << endl;

  for (int i_rayon=i_debut; i_rayon<nb_rayons; ++i_rayon)
    {
      nb++;
      flux_diff[i_rayon] = (aperture_flux[i_rayon] - aperture_flux[i_rayon-1])
	/(i_rayon*i_rayon*M_PI - (i_rayon-1.)*(i_rayon-1.)*M_PI) ;
      mean += flux_diff[i_rayon];
      var += pow(flux_diff[i_rayon], 2.);
      //cout << i_rayon << " " <<  flux_diff[i_rayon] << endl;
    }
  
  mean /= nb;
  var /= nb;
  if (mean>0.)  return sqrt(var-mean*mean)/mean;
  else return -1;
}

float AperturePhotomBaseStar::compute_indep_platitude(int i_debut)
{
  double mean = 0., var = 0.;
  return compute_indep_platitude(mean, var, i_debut);
}


void AperturePhotomBaseStar::write(ostream& s, const bool write_header)
{
  if (write_header)
    {
      WriteHeader_(s);
      s << "# ra : "   << endl
	<< "# dec : "  << endl
	<< "# r3 : "   << endl
	<< "# r6 : "   << endl
	<< "# r10 : "  << endl
	<< "# r19 : "  << endl
	<< "# f3 : "   << endl
	<< "# f6 : "   << endl
	<< "# f10 : "  << endl
	<< "# f19 : "  << endl
	<< "# plat : " << endl
	<< "# end" << endl;
    }
    else
      {
	writen(s);
	s << ra << " " 
	  << dec << " "
	  << rayon[3] << " " 
	  << rayon[6] << " " 
	  << rayon[10] << " " 
	  << rayon[19] << " " 
	  << aperture_flux[3] << " " 
	  << aperture_flux[6] << " " 
	  << aperture_flux[10] << " " 
	  << aperture_flux[19] << " " 
	  << compute_platitude() << " "
	  << endl;
      }
  return;
}

void AperturePhotomBaseStar::write_short(ostream& s, const bool write_header)
{
  if (write_header)
      s << "# x : "     << endl 
	<< "# y : "     << endl
	<< "# ra : "    << endl
	<< "# dec : "   << endl
	<< "# flux : "  << endl
	<< "# fmax : "  << endl
	<< "# fmax2 : " << endl
	<< "# r3 : "    << endl
	<< "# r6 : "    << endl
	<< "# r10 : "   << endl
	<< "# r19 : "   << endl
	<< "# f3 : "    << endl
	<< "# f6 : "    << endl
	<< "# f10 : "   << endl
	<< "# f19 : "   << endl
	<< "# plat : "  << endl
	<< "# end"      << endl;
  else
    s << x << " "
      << y << " "
      << ra << " "
      << dec << " "
      << flux << " "
      << fluxmax << " "
      << fluxmax2 << " "
      << rayon[3] << " " 
      << rayon[6] << " " 
      << rayon[10] << " " 
      << rayon[19] << " " 
      << aperture_flux[3] << " " 
      << aperture_flux[6] << " " 
      << aperture_flux[10] << " " 
      << aperture_flux[19] << " " 
      << compute_platitude() << " " 
      << endl;
  return;
}

void AperturePhotomBaseStar::write_profile(ostream& s, const bool write_header)
{
  if (write_header)
    s << "# x :"     << endl 
      << "# y :"     << endl 
      << "# r :"     << endl 
      << "# f :"     << endl 
      << "# fmax :"  << endl 
      << "# fmax2 :" << endl 
      << "# f3 :"    << endl 
      << "# f6 :"    << endl 
      << "# f10 :"   << endl
      << "# f19 :"   << endl 
      << "# plat :"  << endl 
      << "# end"     << endl;
  else
    for (int i=0; i<nb_rayons; ++i)
      s << x << " "
	<< y << " "
	<< rayon[i] << " "
	<< aperture_flux[i] << " "
	<< fluxmax << " "
	<< fluxmax2 << " "
	<< aperture_flux[3] << " " 
	<< aperture_flux[6] << " " 
	<< aperture_flux[10] << " " 
	<< aperture_flux[19] << " " 
	<< compute_platitude() << " " 
	<< endl;
  return;
}



void MultiMeasuredStar::AddInfo(float rayon_mesure[], float aperture_flux[], float eflux, bool val, float fluxmax, float fluxmax2)
{
  if (nb_measurements>99) {cerr << "too many measurements" << endl; abort();}
  
  for (int i=0; i<nb_rayons; ++i)
    {
      mesures[nb_measurements].measured_aperture_flux[i] = aperture_flux[i];
      mesures[nb_measurements].measurement_radius[i] = rayon_mesure[i];
    }
  mesures[nb_measurements].valid = val;
  mesures[nb_measurements].fluxmax = fluxmax;
  mesures[nb_measurements].fluxmax2 = fluxmax2;
  mesures[nb_measurements].measured_errors = eflux;
  ++nb_measurements;
  
  return;
}

void MultiMeasuredStar::AddInfo(float rayon_mesure[], bool must_be_false)
{
  
  
  for (int i=0; i<nb_rayons; ++i)
    {
      mesures[nb_measurements].measured_aperture_flux[i] = -1.;
      mesures[nb_measurements].measurement_radius[i] = rayon_mesure[i];
    }
  mesures[nb_measurements].valid = false;
  mesures[nb_measurements].measured_errors = -1.;
  ++nb_measurements;
  
  return;
}

void MultiMeasuredStar::AddInfo(bool must_be_false)
{
  mesures[nb_measurements].valid = false;
  ++nb_measurements;
  
  return;
}

void MultiMeasuredStar::AddInfo(AperturePhotomSEStar *aps)
{
  if (nb_measurements>99) {cerr << "too many measurements" << endl; abort();}
  
  for (int i=0; i<nb_rayons; ++i)
    {
      mesures[nb_measurements].measured_aperture_flux[i] = aps->aperture_flux[i];
      mesures[nb_measurements].measurement_radius[i] = aps->rayon[i];
    }
  mesures[nb_measurements].valid = 1;
  mesures[nb_measurements].fluxmax = aps->fluxmax;
  mesures[nb_measurements].fluxmax2 = aps->fluxmax2;
  ++nb_measurements;
  
  return;
}

// calcule la moyenne des mesures pour un rayon donne

float MultiMeasuredStar::compute_mean(int i_rayon)
{
  float mean = 0.;
  int nb=0;
  for (int measurement=0; measurement<nb_measurements; ++measurement)
    {
      if (mesures[measurement].valid) 
	{
	  mean += mesures[measurement].measured_aperture_flux[i_rayon];
	  ++nb;
	}
    }
  if (nb>0) return mean/nb;
  else return -1.;
}

// calcule la moyenne des mesures au carre (pour l'ecart-type)

float MultiMeasuredStar::compute_mean2(int i_rayon)
{
  float mean2 = 0.;
  int nb=0;
  for (int measurement=0; measurement<nb_measurements; ++measurement)
    {
      if (mesures[measurement].valid) 
	{
	  mean2 += pow(mesures[measurement].measured_aperture_flux[i_rayon],(float) 2.);
	  ++nb;
	}
    }
  if (nb>0) return mean2/nb;
  else return -1.;
}

float MultiMeasuredStar::compute_platitude(int measurement)
{
  float mean = 0., var = 0.;
  
  if (! mesures[measurement].valid) return -1.;
  
  int i_debut = 8;
  for (int i_rayon=i_debut; i_rayon<nb_rayons; ++i_rayon)
    {
      mean += mesures[measurement].measured_aperture_flux[i_rayon];
      var += pow(mesures[measurement].measured_aperture_flux[i_rayon],(float) 2.);
    }
  mean /= (nb_rayons-i_debut);
  var /= (nb_rayons-i_debut);
  
  if (mean>0.)  return sqrt(var-mean*mean)/mean;
  else return -1;
}

//
// J'applique un critere sur la platitude de la courbe f(ouverture)
// pour selectionner les objets interessants
//

void MultiMeasuredStar::finalize(int measurement)
{
  nb_valid = 0;
  if (! mesures[measurement].valid)  {mesures[measurement].is_final=false; return;}

  if (compute_platitude(measurement)<0.01) 
    { 
      mesures[measurement].is_final=true;
      mesures[measurement].final_flux = compute_mean(10);
      ++nb_valid;
    }
  return;
}

//
// Ecrit les infos sur les mesures valides
//

void MultiMeasuredStar::write_valid(ostream& s, int star_id, bool write_header=false)
{
  if (write_header)
    s << "# x : "       << endl
      << "# y : "       << endl
      << "# star_id : " << endl
      << "# i : "       << endl
      << "# j : "       << endl
      << "# valid : "   << endl
      << "# fluxmax : " << endl
      << "# fluxmax2 : " << endl
      << "# plat : "    << endl
      << "# f  : "      << endl
      << "# r  : "      << endl
      << "# f2  : "     << endl
      << "# f4  : "     << endl
      << "# f5  : "     << endl
      << "# f6  : "     << endl
      << "# f7  : "     << endl
      << "# f8  : "     << endl
      << "# f10  : "    << endl
      << "# f12  : "    << endl
      << "# f15  : "    << endl
      << "# f19  : "    << endl
      << "# end"        << endl;
  else
    {
      int total_valid=0;
      for(int i=0; i< nb_measurements; ++i) if (mesures[i].valid) ++total_valid;
      
      for(int i=0; i< nb_measurements; ++i)
	{
	  if (mesures[i].valid)
	    for (int j=0; j<nb_rayons; ++j)
	      s << " " << setprecision(12) << x << " " << y << setprecision(8)
		<< " " << star_id
		<< " " << i
		<< " " << j
		<< " " << total_valid
		<< " " << mesures[i].fluxmax
		<< " " << mesures[i].fluxmax2
		<< " " << compute_platitude(i)
		<< " " << mesures[i].measured_aperture_flux[j] 
		<< " " << mesures[i].measurement_radius[j] 
		<< " " << mesures[i].measured_aperture_flux[2] 
		<< " " << mesures[i].measured_aperture_flux[4] 
		<< " " << mesures[i].measured_aperture_flux[5] 
		<< " " << mesures[i].measured_aperture_flux[6] 
		<< " " << mesures[i].measured_aperture_flux[7] 
		<< " " << mesures[i].measured_aperture_flux[8] 
		<< " " << mesures[i].measured_aperture_flux[10] 
		<< " " << mesures[i].measured_aperture_flux[12] 
		<< " " << mesures[i].measured_aperture_flux[15] 
		<< " " << mesures[i].measured_aperture_flux[19] 
		<< endl;
	}
    }
  return;
}

//
// Ecrit les infos sur toutes les mesures
//

void MultiMeasuredStar::write_all(ostream& s, int star_id, bool write_header=false)
{
  if (write_header)
    s << "# x : "       << endl
      << "# y : "       << endl
      << "# star_id : " << endl
      << "# i : "       << endl
      << "# j : "       << endl
      << "# valid : "   << endl
      << "# fluxmax : " << endl
      << "# plat : "    << endl
      << "# f  : "      << endl
      << "# r  : "      << endl
      << "# f2  : "     << endl
      << "# f4  : "     << endl
      << "# f5  : "     << endl
      << "# f6  : "     << endl
      << "# f7  : "     << endl
      << "# f8  : "     << endl
      << "# f10  : "    << endl
      << "# f12  : "    << endl
      << "# f15  : "    << endl
      << "# f19  : "    << endl
      << "# end"        << endl;
  else
    {
      int total_valid=0;
      for(int i=0; i< nb_measurements; ++i) if (mesures[i].valid) ++total_valid;
      
      for(int i=0; i< nb_measurements; ++i)
	{
	  for (int j=0; j<nb_rayons; ++j)
	    s << " " << setprecision(12) << x << " " << y << setprecision(8)
	      << " " << star_id
	      << " " << i
	      << " " << j
	      << " " << total_valid
	      << " " << mesures[i].fluxmax
	      << " " << compute_platitude(i)
	      << " " << mesures[i].measured_aperture_flux[j] 
	      << " " << mesures[i].measurement_radius[j] 
	      << " " << mesures[i].measured_aperture_flux[2] 
	      << " " << mesures[i].measured_aperture_flux[4] 
	      << " " << mesures[i].measured_aperture_flux[5] 
	      << " " << mesures[i].measured_aperture_flux[6] 
	      << " " << mesures[i].measured_aperture_flux[7] 
	      << " " << mesures[i].measured_aperture_flux[8] 
	      << " " << mesures[i].measured_aperture_flux[10] 
	      << " " << mesures[i].measured_aperture_flux[12] 
	      << " " << mesures[i].measured_aperture_flux[15] 
	      << " " << mesures[i].measured_aperture_flux[19] 
	      << endl;
	}
    }
  
  return;
  
}

//
// cherche l'etoile de la liste laliste qui correspond a la multimesured
// star en cours
//

bool MultiMeasuredStar::match_aps(const SEStarList& laliste,  Gtransfo * pix2radec, AperturePhotomSEStar* aps)
{

#define  ASSOCIATION_LANDOLT_DISTANCE 0.0001

  const double ra1 = x;
  const double dec1 = y;

  // je cherche l'etoile du catalogue correspondant a celle que je cherche
  // Je suis oblige car le se.list a les coord en pixels alors que ma liste
  // de depart les a en ra-dec.

  for (SEStarCIterator si = laliste.begin(); si != laliste.end(); ++si)
    {
      const SEStar *se = *si;

      const double x0 = se->x;
      const double y0 = se->y;
      double ra,dec;
      
      pix2radec->apply(x0, y0, ra, dec);

      double distance = sqrt( pow(ra1 - ra,2.)+pow(dec1 - dec,2.));
      
      if (distance < ASSOCIATION_LANDOLT_DISTANCE) 
	{
	  aps->x = se->x;
	  aps->y = se->y;
	  aps->flux = se->flux;
	  return true;
	}
    }
  return false;
}



typedef StarList<MultiMeasuredStar> MultiMeasuredStarList;
typedef MultiMeasuredStarList::const_iterator MultiMeasuredStarCIterator;
typedef MultiMeasuredStarList::iterator MultiMeasuredStarIterator;
typedef CountedRef<MultiMeasuredStar> MultiMeasuredStarRef;

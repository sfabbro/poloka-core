/* Calcule le flux en ADU par photometrie d'ouverture  */

#define CHECKBOUNDS

//#define DEBUG
//#define TEST
//#define WRITE_ALL
//#define WRITE_FINAL_LIST
#define WRITE_COMPLETE_LIST
#define WRITE_SE2_LIST

#include <iostream>
#include <fstream>
#include "fileutils.h"
#include "dbimage.h"
#include <math.h>
#include "wcsutils.h"  // pour passer des coordonnees pixel aux RA-DEC
#include "aperturephotometry.h"


#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif


int burn_object(FitsImage& fitsimage, FitsImage& masque_objets, int x, int y, BaseStar *star, double seuil);
bool validate_neighbour(FitsImage& fitsimage, FitsImage& masque_objets, int x_neighbour, int y_neighbour, double seuil);
double fraction(double, double, double, double, double);
void extract_seeing_exptime(FitsImage& limage, float &seeing, float &exptime);
void extract_airmass(FitsImage& fitsimage, float &phot_c, float &phot_k, float &airmass);
bool extract_noweight(FitsImage& fitsimage,  float seeing, double mean, AperturePhotomBaseStar * se1);
bool test_no_other_same(BaseStarList& apstarlist,const BaseStar *star);
void dist_corners(Gtransfo *pix2radec, double ra, double dec, double x0, double y0, double &dist_seuil, double &x, double &y);
void wcs_to_pixel(Gtransfo *pix2radec, double ra, double dec, double &x, double &y);

// =================================================================
// =================================================================
//
// MAIN
//
// =================================================================
// =================================================================

int main(int argc, char ** argv)
{

  if (argc < 2)
    {
      std::cerr << "usage : myphotometry -l list <dbimage>" << std::endl;
      std::cerr << "  process coordinates listed in ra-dec in list" << std::endl;
      std::cerr << "  extracts a rudimentary list if none is provided" << std::endl;      exit(1);
    }

  // Je decortique les arguments, et je cherche l'image fits

  BaseStarList starlist;
  bool list_is_provided=false;
  char coordinate_system = 'n';

  for (int i=1;  i<argc; ++i) // argv[0] is the executable name
    {      
      char *arg = argv[i];

      if (arg[0] == '-')
	{
	  switch (arg[1])
	    {
	    case 'l' : ++i; list_is_provided = true; 
	      starlist = BaseStarList(argv[i]);
	      break;
	      
	    case 'r' : coordinate_system='r'; break;
	    case 'p' : coordinate_system='p'; break;

	    default: std::cerr << "don't understand " << arg << std::endl; 
	    }
	}
      else
	{

	  if ((coordinate_system!='r')&&(coordinate_system!='p')&&list_is_provided)
	    {
	      cerr << "you must specify the coordinate format for lists\nuse -r (radec) or -p (pixels) option" << endl;
	    }

	  if (!list_is_provided) coordinate_system='p';

	  char *fitsname = argv[i];
	  
	  FitsImage fitsimage(fitsname);
	  FitsHeader fitsheader(fitsname);
	  
	  Gtransfo *pix2radec;
	  WCSFromHeader(fitsname, pix2radec);
	  
	  Pixel mean, sigma;
	  fitsimage.SkyLevel(&mean, &sigma);
	  
	  //	  cout << "skylev : " << mean << " " << sigma << endl;
	  
#define THRESHOLD 5
	  
	  //
	  // Si pas de liste fournie, on va chercher les objets
	  //
	  
	  if (!list_is_provided)
	    {
	      FitsImage masque_objets("mask.fits", fitsheader );
	      for (int y=0; y<fitsimage.Ny(); ++y)
		for (int x=0; x<fitsimage.Nx(); ++x)
		  {
		    if ((masque_objets(x,y) < 0.5 )&&(fitsimage(x,y) > mean + THRESHOLD*sigma))
		      {
			BaseStar *star = new BaseStar();

			int npixels = burn_object(fitsimage,masque_objets,x,y,star,mean + THRESHOLD*sigma*0.5 );

			if (npixels>10)
			  {
			    star->flux -= npixels*mean;
			    starlist.push_back(star);
			  }
		      }
		  }
	      BaseStarList starlist_radec;
	      for (BaseStarCIterator bs = starlist.begin(); bs != starlist.end(); ++bs)
		{
		  const BaseStar *tagada = *bs;
		  BaseStar *prov = new BaseStar(*tagada);
		  double ra, dec;
		  pix2radec->apply(prov->x, prov->y, ra, dec);
		  prov->x = ra;
		  prov->y = dec;
		  starlist_radec.push_back(prov);
		}
	      starlist_radec.write("mylist_radec.list");
	      starlist.write("mylist.list");
	    }
	  
	  //
	  // Je vais lire le seeing et le temps d'exposition de l'image
	  //
	  

	  float seeing, exptime, phot_k, phot_c, airmass;
	  extract_seeing_exptime(fitsimage, seeing, exptime);
	  extract_airmass(fitsimage, phot_c, phot_k, airmass);

	  seeing = 1.;

	  ofstream output_file("test.list");
	  ofstream output_profile("test_profil.list");
	  AperturePhotomBaseStar apstarlist;
	  apstarlist.write(output_file, true);
	  apstarlist.write_profile(output_profile, true);
	  
	  //	  cout << "# x :\n# y :\n# ra :\n# dec :\n# mag5 :\n# mag10 :\n# mag15 :\n# mag18 :\n# flux10 :\n# fmax :\n# fmax2 :\n# plat :\n# end\n" << endl; 

	  //cout << "# x :\n# y :\n# ra :\n# dec :\n# flux5 :\n# flux10 :\n# flux15 :\n# flux1 :\n# fmax :\n# fmax2 :\n# plat :\n# end\n" << endl; 



	  for (BaseStarCIterator bs = starlist.begin(); bs != starlist.end(); ++bs)
	    {
	      const BaseStar *prov = *bs;
	      AperturePhotomBaseStar pabs(*prov);

	      //
	      // Si les coordonnees de la liste sont en radec, je les transforme en pixels
	      // Commenter la ligne suivante si la liste est en pixels
	      //
	      double xx=pabs.x,yy=pabs.y;
	      if (coordinate_system=='r') wcs_to_pixel(pix2radec, prov->x, prov->y, xx, yy);
	      pabs.x = xx;
	      pabs.y = yy;

	      //	      cout << "pos in pix : " << pabs.x << " " << pabs.y << endl;

	      pabs.SetSeeing(seeing);
	      
	      //
	      // Compute the aperture fluxes for this star
	      //

	      bool valid = extract_noweight(fitsimage, seeing, mean, &pabs);

	      // 
	      // teste si des objets de la liste sont en double... 
	      // obsolete avec l'agorithme burn-out
	      //
	      
	      // if (!list_is_provided && valid) valid = test_no_other_same(starlist,prov);

	      if (valid||list_is_provided) 
		{
		  double ra, dec;
		  pix2radec->apply(pabs.x, pabs.y, ra, dec);
		  pabs.ra = ra;
		  pabs.dec = dec;
		  pabs.write(output_file);
		  pabs.write_profile(output_profile);

		  double local_mean, local_var;
		  //		  float plat = pabs.compute_indep_platitude(local_mean, local_var);
		  float plat = pabs.compute_platitude(local_mean, local_var);
		  float mag_5 = phot_c + phot_k*(airmass-1.) - 2.5*log10(pabs.aperture_flux[5])+2.5*log10(exptime);		  
		  float ZP = phot_c + phot_k*(airmass-1.) + 2.5*log10(exptime);
		  float mag_10 = ZP - 2.5*log10(pabs.aperture_flux[10]);
		  float mag_15 = ZP - 2.5*log10(pabs.aperture_flux[15]);		  
		  float mag_18 = ZP - 2.5*log10(pabs.aperture_flux[18]);		  
		  /*		  
		  cout << setprecision(6) 
		       << pabs.x  << " "
                       << pabs.y << " "
		       << setprecision(12) 
		       << pabs.ra  << " " 
		       << pabs.dec << " " 
		       << setprecision(6) 
		       << mag_5   << " " 
		       << mag_10  << " " 
		       << mag_15  << " " 
		       << mag_18  << " " 
		       << pabs.aperture_flux[10]  << " " 		    
		       << pabs.fluxmax  << " " 		    
		       << pabs.fluxmax2  << " " 
		       << plat<< " "
		       << endl;
		  */

		  /*
		    cout << setprecision(6) 
		       << pabs.x  << " "
                       << pabs.y << " "
		       << setprecision(12) 
		       << pabs.ra  << " " 
		       << pabs.dec << " " 
		       << setprecision(6) 
		       << pabs.aperture_flux[5]   << " " 
		       << pabs.aperture_flux[10]  << " " 
		       << pabs.aperture_flux[15]  << " " 
		       << pabs.aperture_flux[19]  << " " 
		       << pabs.fluxmax  << " " 		    
		       << pabs.fluxmax2  << " " 
		       << plat<< " "
		       << endl;
		  */
		}
	    }
	  
	  output_file.close();
	  output_profile.close();

	}
	  
    }

  
  return 0;
}


/* ================================================================== */
/* ================================================================== */
/* ================================================================== */

int burn_object(FitsImage& fitsimage, FitsImage& masque_objets, int x, int y, BaseStar *star, double seuil)
{
#define Nmax 100000
  int x_to_burn[Nmax]; 
  int y_to_burn[Nmax];
  int last_item = 0;
  int current_item = 0;
  int npixels = 0;
  x_to_burn[0] = x;
  y_to_burn[0] = y;

  int x_neighbour, y_neighbour;

  double x_moyen=0., y_moyen=0., total=0.;

  while ( (current_item <= last_item)&&(current_item <Nmax))
    {
      if (masque_objets(x_to_burn[current_item],y_to_burn[current_item]) < 0.5)
	{
	  double val = fitsimage(x_to_burn[current_item],y_to_burn[current_item]);
	  npixels++;
	  total += val;
	  x_moyen += x_to_burn[current_item]*val;
	  y_moyen += y_to_burn[current_item]*val;
	  masque_objets(x_to_burn[current_item],y_to_burn[current_item]) = 1.;
	  
	  x_neighbour = x_to_burn[current_item] + 1;
	  y_neighbour = y_to_burn[current_item] + 0;
	  	  
	  if (validate_neighbour(fitsimage, masque_objets, x_neighbour, y_neighbour, seuil))
	    {
	      last_item++;
	      x_to_burn[last_item] = x_neighbour; 
	      y_to_burn[last_item] = y_neighbour;	  
	    }

	  x_neighbour = x_to_burn[current_item] - 1;
	  y_neighbour = y_to_burn[current_item] + 0;
	  	  
	  if (validate_neighbour(fitsimage, masque_objets, x_neighbour, y_neighbour, seuil))
	    {
	      last_item++;
	      x_to_burn[last_item] = x_neighbour; 
	      y_to_burn[last_item] = y_neighbour;	  
	    }

	  x_neighbour = x_to_burn[current_item] + 0;
	  y_neighbour = y_to_burn[current_item] - 1;
	  	  
	  if (validate_neighbour(fitsimage, masque_objets, x_neighbour, y_neighbour, seuil))
	    {
	      last_item++;
	      x_to_burn[last_item] = x_neighbour; 
	      y_to_burn[last_item] = y_neighbour;	  
	    }

	  x_neighbour = x_to_burn[current_item] + 0;
	  y_neighbour = y_to_burn[current_item] + 1;
	  	  
	  if (validate_neighbour(fitsimage, masque_objets, x_neighbour, y_neighbour, seuil))
	    {
	      last_item++;
	      x_to_burn[last_item] = x_neighbour; 
	      y_to_burn[last_item] = y_neighbour;	  
	    }

	}
      else current_item++;
    }

  star->x = x_moyen/total;
  star->y = y_moyen/total;
  star->flux = total;
  
  return npixels;
}

/* ================================================================== */

bool validate_neighbour(FitsImage& fitsimage, FitsImage& masque_objets, int x_neighbour, int y_neighbour, double seuil)
{

  if ((x_neighbour < 1) || (x_neighbour >=  fitsimage.Nx())
      || (y_neighbour < 1) || (y_neighbour >=  fitsimage.Ny())) return false; // out of image
  
  if (masque_objets(x_neighbour,y_neighbour) > 0.5) return false; // already burnt
  if (fitsimage(x_neighbour,y_neighbour) < seuil) return false; // too low to burn
	  
  return true;
}

/* ================================================================== */

bool extract_noweight(FitsImage& fitsimage, float seeing, double mean, AperturePhotomBaseStar * se1)
{
  bool valid = se1->ComputeBarycenter_noweight(fitsimage);

  //    cout << "1> " << valid;
  
  if (valid) valid = se1->FillFlux_noweight(fitsimage, mean);
  
  // cout << " 2> " << valid;

  for (int i=0; i<nb_rayons; ++i) 
    {
      if (se1->aperture_flux[i]< 0.) valid = false; 
    }
  
  // cout << " 3> " << valid << endl;
  return valid;	    
}

/* ================================================================== */

/* je veux calculer le flux contenu dans un cercle de rayon R, sur une image composee
   de nx*ny pixels */

//#include <math>
#include <iostream>
#include <fstream>

using namespace std; 

// premier probleme, quelle est la fraction de pixel couverte par un disque de rayon R ?
// x_pix, y_pix coord du centre du pixel (de taille 1x1 par convention)

double fraction(double x_centre, double y_centre, double R, double x_pix, double y_pix)
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

/* ================================================================== */

void extract_seeing_exptime(FitsImage& fitsimage, float &seeing, float &exptime)
{

  if (fitsimage.HasKey("SESEEING"))
    seeing = fitsimage.KeyVal("SESEEING");
  else
    {
      cerr << "Image has no SESEEING key" << endl;
      seeing = 1.5;
    }

  if (fitsimage.HasKey("EXPTIME"))
    exptime = fitsimage.KeyVal("EXPTIME");
  else
    {
      cerr << "Image has no EXPTIME key" << endl;
      exptime = 1.;
    }

  return;
}
/* ================================================================== */

void extract_airmass(FitsImage& fitsimage, float &phot_c, float &phot_k, float &airmass)
{

  if (fitsimage.HasKey("PHOT_K"))
    phot_k = fitsimage.KeyVal("PHOT_K");
  else
    {
      cerr << "Image has no PHOT_K key" << endl;
      phot_k = 0.;
    }
  if (fitsimage.HasKey("PHOT_C"))
    phot_c = fitsimage.KeyVal("PHOT_C");
  else
    {
      cerr << "Image has no PHOT_C key" << endl;
      phot_c = 0.;
    }
  if (fitsimage.HasKey("AIRMASS"))
    airmass = fitsimage.KeyVal("AIRMASS");
  else
    {
      cerr << "Image has no AIRMASS key" << endl;
      airmass = 1.;
    }

  return;
}

/* ================================================================== */

bool test_no_other_same(BaseStarList& starlist,const BaseStar *star)
{
  bool result = true;

  for (BaseStarCIterator bs = starlist.begin(); bs != starlist.end(); ++bs)
    {
      const BaseStar *prov = *bs;

      double distance = sqrt( (prov->x - star->x)*(prov->x - star->x)
			      - (prov->y - star->y)*(prov->y - star->y));
 
      if ((distance < 8.)&&(distance > 0.0001)) result = false;
    }

  //if (result == false) cout << "rejected" << endl;
  return result;
}

/* ================================================================== */

void wcs_to_pixel(Gtransfo *pix2radec, double ra, double dec, double &x, double &y)
{
#define NPAS 30

  double taille_x=2112.;
  double taille_y=4644.;

  double x_nearest=taille_x/2.;
  double y_nearest=taille_y/2.;

  for (int niter=1; niter<=NPAS; ++niter)
    {
      double x_11 = x_nearest;
      double y_11 = y_nearest;
      double d_nearest=10000.;

      dist_corners(pix2radec,ra,dec,x_11-taille_x/2., y_11-taille_y/2., d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11-taille_x/2., y_11            , d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11-taille_x/2., y_11+taille_y/2., d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11            , y_11-taille_y/2., d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11            , y_11            , d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11            , y_11+taille_y/2., d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11+taille_x/2., y_11-taille_y/2., d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11+taille_x/2., y_11            , d_nearest, x_nearest, y_nearest);
      dist_corners(pix2radec,ra,dec,x_11+taille_x/2., y_11+taille_y/2., d_nearest, x_nearest, y_nearest);

      //      cout << "dist : " << d_nearest << " " << x_nearest << " " << y_nearest << endl;

      taille_x *= 0.75;
      taille_y *= 0.75;

    }

  x = x_nearest;
  y = y_nearest;

  return;
}

void dist_corners(Gtransfo *pix2radec, double ra, double dec, double x0, double y0, double &dist_seuil, double &x, double &y)
{
  double ra_0,dec_0;
  pix2radec->apply(x0, y0, ra_0, dec_0);
  double distance = sqrt( (ra_0-ra)*(ra_0-ra) + (dec_0-dec)*(dec_0-dec));
  
  if (distance < dist_seuil) 
    {
      x = x0;
      y = y0;
      dist_seuil = distance;
    }
  
  return;
}

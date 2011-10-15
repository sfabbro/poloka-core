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
#include <cmath>
#include "wcsutils.h"  // pour passer des coordonnees pixel aux RA-DEC
#include "aperturephotometry.h"
#include "frame.h"


#ifndef M_PI
#define  M_PI           3.14159265358979323846  /* pi */
#endif


int burn_object(FitsImage& fitsimage, FitsImage& masque_objets, int x, int y, BaseStar *star, double seuil);
bool validate_neighbour(FitsImage& fitsimage, FitsImage& masque_objets, int x_neighbour, int y_neighbour, double seuil);
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

static void usage(const char *prog)
{
  std::cerr << "usage :" << std::endl
	    << prog << " [options] <fitsimage>" << std::endl
	    << "  process coordinates listed in ra-dec in list" << std::endl
	    << "  extracts a rudimentary list if none is provided" << std::endl
	    << " options: " << std::endl
	    << "  -l <list>: do photometry of objects in this external catalogue" << std::endl
	    << "  -p       : coordinates of external catalog in pixels" << std::endl
	    << "  -r       : coordinates of external catalog in RA/Dec" << std::endl;
    exit(1);
}  


int main(int argc, char ** argv)
{

  if (argc < 2) usage(argv[0]);

  // Je decortique les arguments, et je cherche l'image fits

  BaseStarList starlist;
  bool list_is_provided=false;
  char coordinate_system = 'n';
  string fitsname;

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
	    case 'h' : usage(argv[0]); break;

	    default: std::cerr << "don't understand " << arg << std::endl;
	      usage(argv[0]);
	    }
	}
      else
	{
	  fitsname = argv[i];
	}
    }
  if (fitsname == "")
    {
      cerr << " I need a fitsimage !! " << endl;
      usage(argv[0]);
    }


  if ((coordinate_system!='r')&&(coordinate_system!='p')&&list_is_provided)
    {
      cerr << "you must specify the coordinate format for lists\nuse -r (radec) or -p (pixels) option" << endl;
    }
  
  if (!list_is_provided) coordinate_system='p';

	  
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
		else delete star;
	      }
	  }
      
      BaseStarList starlist_radec;
      starlist.CopyTo(starlist_radec);
      starlist_radec.ApplyTransfo(*pix2radec);
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
  apstarlist.WriteHeader(output_file);
  apstarlist.write_profile(output_profile, true);
  
  
  Gtransfo *radec2pix = pix2radec->InverseTransfo(0.01,Frame(fitsheader));
  for (BaseStarCIterator bs = starlist.begin(); bs != starlist.end(); ++bs)
    {
      const BaseStar *prov = *bs;
      cout << *prov << endl;
      AperturePhotomBaseStar pabs(*prov);
      
      //
      // Si les coordonnees de la liste sont en radec, je les transforme en pixels
      // Commenter la ligne suivante si la liste est en pixels
      //
      
      if (coordinate_system=='r') pabs.Apply(*radec2pix);
      
      //	      cout << "pos in pix : " << pabs.x << " " << pabs.y << endl;
      
      pabs.SetSESeeing(seeing);
      
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



#include "daophotpsf.h"
#include "dimage.h"
#include "fileutils.h"

static void usage(char *progName)
{
  cerr << progName << " <PSF file> <x> <y> <flux> \n" 
       << "  produces a FITS image of a star from a PSF, given a PSF file, daophot style \n"
       << "    <x> <y> : if given, position of the star (default is PSF image center) \n"
       << "    <flux>  : if given the flux of the star  (default is 1) \n"
       << " The FITS image of the star is hardcoded of psf_star.fits \n"
       << " PS: do not change the order of arguments. I am lasy. \n";
  exit(-1);
}

int main(int nargs, char **args)
{
  if ((nargs < 2) || (nargs > 5)) {usage(args[0]);}

  string name = args[1];
  DaoPsf psf(name);
  
  double x = psf.Xpsf();
  if (nargs >2) x = atof(args[2]);

  double y = psf.Ypsf();
  if (nargs >3) y = atof(args[3]);

  double flux = 1.;
  if (nargs >4) flux = atof(args[4]);

  Kernel star_image;
 
  BaseStar star(x,y,flux);
  psf.MakeStar(star_image, star);
  cout << " Pixel sum = " << star_image.sum() << endl;

  star_image.writeFits("psf_star.fits");

  return EXIT_SUCCESS;
}

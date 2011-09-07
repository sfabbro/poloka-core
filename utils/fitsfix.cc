#include "fitsimage.h"
#include "fitstoad.h" 
#include "imageutils.h"

static bool fitsfix_process(const char* fitsname) {
  
  FitsImage image(fitsname, RW);
  if (!image.IsValid()) {
    cerr << " invalid FITS file : "  << fitsname << endl;
    return false;
  }

  // trim
  Frame Illu = TotalIlluRegion(image);
  image.Trim(Illu);
  
  // set saturation level  
  image.AddOrModKey("SATURLEV", ComputeSaturation(image), "Current saturation level");

  // compute stats
  Pixel sky,sig;
  image.SkyLevel(&sky, &sig);
  double rdnoise = image.KeyVal("TOADRDON");
  double gain = image.KeyVal("TOADGAIN");
  double sigth = sqrt(sky*gain+rdnoise*rdnoise)/gain;
  image.AddOrModKey("SKYLEV",sky," Original sky level");
  image.AddOrModKey("SKYSIGEX",sig," Original sky sigma");
  image.AddOrModKey("SKYSIGTH",sigth," app. Theoretical sigma (sqrt(sky*gain+rdnoise^2)/gain");
  
  // fix the masked 0 pixels with sky (assume non subtracted sky)
  Pixel *p = image.begin();
  Pixel *pend = image.end();
  for ( ; p < pend ; ++p) if (fabs(*p) < 0.1)  *p = sky;
  return true;
}

static void usage(const char *progname) {
  cout << progname << " <FITS image>...<FITS image>\n"
       << " Prepare an image:\n" 
       << "    - trim the overscan regions\n"
       << "    - compute the saturation and sky levels\n"
       << "    - fix the masked pixels if 0\n"
       << endl;
  exit(-1);
}

int main(int argc,char **argv) {

  if(argc <=1) usage(argv[0]);
  
  bool ok = true;

  for (int i=1; i< argc; ++i)
    ok &= fitsfix_process(argv[i]);

  return (ok ? EXIT_SUCCESS: EXIT_FAILURE);
}

